require 'nokogiri'	
require 'mspire/mzml'
require 'pry'
#require_relative 'binary_search'

def ppm(m1,m2)
  (m2-m1)/m1*1e6
end
def ppm_range(mass, ppm = MatchThreshold)
  range = (ppm*mass)/1e6
  (mass-range..mass+range)
end


PepID = Struct.new(:precursor_neutral_mass, :aaseq, :mods, :charge, :scan_num, :retention_time, :protein)
module PeakIDMapper
  StandardWindow = [-60, 200]
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
        output.precursor_neutral_mass = spec_query.attributes["precursor_neutral_mass"].value.to_f
        output.retention_time = nil
        output.protein = search_hit.attributes["protein"].value
        output.aaseq = search_hit.attributes["peptide"].value
        pepids << output
      end
      pepids.uniq {|a| a.aaseq} ## PROBABLY a bad idea
    end
  end

  class MzmlParser
    Matches = Struct.new(:retention_times, :mzs, :intensities)
    PpmThreshold = 100
    def self.join_pepxml_with_mzml_file(pepdata, file)
      output_map = []
      # This is the file which you must cache.  
      # Basically, pre-read the file into memory and then play with the structures you've cached inside the loop
      # The problem otherwise is that the random access abilities of Mspire reek havoc on the CPU cycles.
      
      # I would do something like this:
      SpectralObject = Struct.new(:peaks, :retention_time)
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
        spectra.each do |spectralobject| # |peaks, retention_time|
          next if spectralobject.peaks.empty?
          # ... continued.
      ## Previous code starts here here
      Mspire::Mzml.open(file) do |mzml|
        pepdata.each do |pepid|
          rt_array, mz_array, int_array = [],[],[]
          mass = pepid.precursor_neutral_mass
          mzml.each_spectrum do |spectrum|
            next if spectrum.peaks.empty?
            list = spectrum.peaks
            rrange = list.transpose.first.bsearch_range do |a| 
              ppm_range(a, PpmThreshold).include?(mass) ? 0 : a <=> mass
            end
            resp = spectrum.peaks[rrange].flatten
            next if resp.first == nil 
            rt_array << spectrum.retention_time
            mz_array << resp.first
            int_array << resp.last
          end
          next if mz_array.empty?
          output_map << PeakIDMap.new(pepid.aaseq, pepid.charge, nil, mz_array.inject(:+)/mz_array.size.to_f, int_array.inject(:+), [rt_array, mz_array, int_array])
        end
      end
      output_map
    end
  end

  PeakIDMap = Struct.new(:aaseq, :charge, :mods, :mz_mean, :total_intensity, :data) # Data = [rt_array, mz_array, int_array]

  module CommandLine

    def self.run(pepxml, mzml, opts = {})
      pepdata = PepxmlParser.parse(pepxml)
      peakIDs = MzmlParser.join_pepxml_with_mzml_file(pepdata, mzml)
      p peakIDs
    end
  end
end


if __FILE__ == $0
  args = ARGV
  PeakIDMapper::CommandLine.run(*args)
end
