require 'nokogiri'	
require 'mspire/mzml'
require 'pry'
#require_relative 'binary_search'

GLOBALJOIN = ";"

PpmThreshold = 100
RTThreshold = 400 # amount forward and backward to search
IntensityThreshold = 500
Matches = Struct.new(:retention_times, :mzs, :intensities)
# determine_stats returns [min, max, mean, sample_variance, standard_deviation, sum, size]
PeakIDMap = Struct.new(:aaseq, :proteins, :mods, :charge, :mh, :ion_score, :ppm_error, :spectrum_file, :match_file, :rt, :mz, :int, :rt_array, :mz_array, :int_array) # Data = [rt_array, mz_array, int_array]
HEADERLINE = %w{Sequence Proteins RTmin RTmax RTmean RTvariance RTstdev RTsum RTsize Modifications Charge MZmin MZmax MZmean MZvariance MZstdev MZsum MZsize MH+ deltaM(ppm) Ion_Score SpectrumFile MatchFile Imin Imax Imean Ivariance Istdev Isum Isize RT MZ Intensities RT_array MZ_array Int_array}.join(GLOBALJOIN)
PepID = Struct.new(:precursor_neutral_mass, :aaseq, :mods, :charge, :scan_num, :retention_time, :proteins, :match_file, :ion_score, :ppm_error, :mh)

def ppm(m1,m2)
  (m2-m1)/m1*1e6
end
def ppm_range(mass, ppm = MatchThreshold)
  range = (ppm*mass)/1e6
  (mass-range..mass+range)
end

module Enumerable
  def sum
    self.inject(:+)
  end
  def mean
    self.sum/self.length.to_f
  end
  def sample_variance
    m = self.mean
    sum = self.inject(0) {|prod, i| prod + (i-m)**2 }
    sum/(self.length - 1).to_f
  end
  def standard_deviation
    Math.sqrt(self.sample_variance)
  end
end

class PeakIDMapper
  StandardWindow = [-60, 200]
  def self.determine_stats(array)
    min, max = array.minmax
    [min, max, array.mean, array.sample_variance, array.standard_deviation, array.sum, array.size]
  end
  def self.peakIDs_to_csv(peakids, file = nil)
    file ||= peakids.map{|a| [a.spectrum_file, a.match_file].map{|f| File.basename(f).gsub(File.extname(f),"")}}.flatten.uniq.join("_") + '.csv'
    File.open(file, 'w+') do |outstream|
      outstream.puts HEADERLINE
      peakids.each do |pk|
        outstream.puts [pk.aaseq, pk.proteins, pk.rt, pk.mods, pk.charge, pk.mz, pk.mh, pk.ppm_error, pk.ion_score, pk.spectrum_file, pk.match_file, pk.int, pk.rt.join(","), pk.mz.join(","), pk.int.join(","), pk.rt_array.join(","), pk.mz_array.join(","), pk.int_array.join(",")].flatten.join(GLOBALJOIN)
      end
    end
  end
  class PepxmlParser
    # This is from the run_compare code I wrote
    def self.parse(file)
      doc = Nokogiri.XML(File.open(file))
      pepids = []
      search_sum = doc.xpath("//xmlns:search_summary") # search for aminoacid_modification and terminal_modification
      doc.xpath('//xmlns:search_hit[@hit_rank=$value]', nil, {value: '1'}).each do |search_hit|
        spec_query = search_hit.parent.parent
        output = PepID.new
        output.scan_num = spec_query.attributes["start_scan"].value.to_i 
        output.charge = spec_query.attributes["assumed_charge"].value.to_i
        mod_info = search_hit.css "modification_info"
        unless mod_info.empty?
          output.mods = mod_info.map do |mod| 
            tmp = mod.attributes.entries
            tmp.map {|f,l| [f, mod.attributes[f].value.to_f]}
          end
        end
        output.precursor_neutral_mass = spec_query.attributes["precursor_neutral_mass"].value.to_f
        output.retention_time = spec_query.attributes["start_scan"].value.to_i
        output.proteins = search_hit.attributes["protein"].value
        output.aaseq = search_hit.attributes["peptide"].value
        output.match_file = file
        output.ion_score = search_hit.css("search_score").to_s[/"ionscore" value="(.*?)"\/>/,1].to_f
        output.mh = search_hit.attributes["calc_neutral_pep_mass"].value.to_f
        output.ppm_error = search_hit.attributes["massdiff"].value.to_f
        pepids << output
      end
      pepids#.uniq {|a| a.aaseq} ## PROBABLY a bad idea
    end
  end

  class MzmlParser
    SpectralObject = Struct.new(:peaks, :retention_time)
    def self.join_pepxml_with_mzml_file(pepdata, file)
      output_map = []
      # This is the file which you must cache.  
      # Basically, pre-read the file into memory and then play with the structures you've cached inside the loop
      # The problem otherwise is that the random access abilities of Mspire reek havoc on the CPU cycles.
      
      # I would do something like this:
      spectra = []
      Mspire::Mzml.open(file) do |mzml|
        mzml.each_spectrum do |spectrum|
          spectra << SpectralObject.new(spectrum.peaks, spectrum.retention_time)
        end
      end
      # Now,do the following loops using the previous data instead...
      
      pepdata.each do |pepid|
        rt_array, mz_array, int_array = [],[],[]
        mass = pepid.precursor_neutral_mass
        mass_range = ppm_range(mass, PpmThreshold)
        time_range = (pepid.retention_time-RTThreshold)..(pepid.retention_time+RTThreshold)
        spectra.each do |spectralobject| # |peaks, retention_time|
          next if spectralobject.peaks.empty?
          # ... continued.
            next unless time_range.include?(spectralobject.retention_time)
            list = spectralobject.peaks
            rrange = list.transpose.first.bsearch_range do |a| 
              mass_range.include?(a) ? 0 : a <=> mass
            end
            resp = spectralobject.peaks[rrange].flatten #select {|a| puts "mz: #{a.first}"; puts "ppm: #{ppm(a.first,mass)}"; ppm(a.first,mass) < PpmThreshold}
            next unless mass_range.include?(resp.first)
            next if resp.first == nil
            next if resp.last < IntensityThreshold
            rt_array << spectralobject.retention_time
            mz_array << resp.first
            int_array << resp.last
          end
          next if mz_array.empty?
          output_map << PeakIDMap.new(pepid.aaseq, pepid.proteins, pepid.mods, pepid.charge, pepid.mh, pepid.ion_score, pepid.ppm_error, file, pepid.match_file, PeakIDMapper.determine_stats(rt_array), PeakIDMapper.determine_stats(mz_array), PeakIDMapper.determine_stats(int_array), rt_array, mz_array, int_array)
        end
      output_map
    end
  end
  module CommandLine
    def self.run(pepxml, mzml, opts = {})
      pepdata = PepxmlParser.parse(pepxml)
      peakIDs = MzmlParser.join_pepxml_with_mzml_file(pepdata, mzml)
    end
  end
end


if __FILE__ == $0
  args = ARGV
  if args.size < 2
    puts "Error."
    puts "USAGE: pepxml_file, mzml_file, [future_opts]"
    exit
  end
  peak_ids = PeakIDMapper::CommandLine.run(*args)
  PeakIDMapper.peakIDs_to_csv(peak_ids)
end
