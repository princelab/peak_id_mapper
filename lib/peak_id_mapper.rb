require 'nokogiri'	
require 'mspire/mzml'
require 'mspire/mass/all'
require 'pry'
require 'gnuplot'

#require_relative 'binary_search'

$PIMVERBOSE = true

def putsv(thing)
  puts thing if $PIMVERBOSE
end

GLOBALJOIN = ";"

PpmThreshold = 100
RTThreshold = 400 # amount forward and backward to search
IntensityThreshold = 500
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
def ppm_mz_range(mz, ppm = MatchThreshold)
  range = (ppm*mz)/1e6
  (mz-range..mz+range)
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

  def self.plot_peakIDs(peak_ids)
    peak_ids.each_with_index do |peak_id,i|
      Gnuplot.open do |gp|
        gp << "set term svg enhanced\n"
        gp << %Q{set output "spectrum#{i}.svg"\n}
        Gnuplot::Plot.new( gp ) do |plot|
          title_string = [:aaseq, :charge, :mh, :ion_score, :ppm_error, :mz, :rt].map do |key| 
            val = peak_id.send(key)
            [key, val.is_a?(Float) ? val.round(4) : val].join(": ") 
          end.join(", ")
          File.write("spectrum#{i}.txt", title_string)

          plot.title title_string
          plot.xlabel "time (s)"
          plot.ylabel "intensity"
          plot.y2label "m/z"

          plot.data << Gnuplot::DataSet.new( [peak_id.rt_array, peak_id.mz_array] ) do |ds|
            ds.axes = "x1y1"
            ds.with = "lines"
            ds.title = "ion intensity"
          end
          plot.data << Gnuplot::DataSet.new( [peak_id.rt_array, peak_id.int_array] ) do |ds|
            ds.axes = "x1y2"
            ds.title = "m/z"
          end
        end
      end
    end
  end

  def self.peakIDs_to_csv(peakids, file = nil)
    putsv "Loading the peakIDS to csv"
    file ||= peakids.map{|a| [a.spectrum_file, a.match_file].map{|f| File.basename(f).gsub(File.extname(f),"")}}.flatten.uniq.join("_") + '.csv'
    putsv "Writing file: #{file}"
    File.open(file, 'w') do |outstream|
      outstream.puts HEADERLINE
      peakids.each do |pk|
        outstream.puts [pk.aaseq, pk.proteins, pk.rt, pk.mods, pk.charge, pk.mz, pk.mh, pk.ppm_error, pk.ion_score, pk.spectrum_file, pk.match_file, pk.int, pk.rt.join(","), pk.mz.join(","), pk.int.join(","), pk.rt_array.join(","), pk.mz_array.join(","), pk.int_array.join(",")].flatten.join(GLOBALJOIN)
      end
    end
    putsv "Failure" unless File.exists?(file)
  end
  class PepxmlParser
    # returns an array of pep identifications (PepID's).
    def self.parse(file)
      putsv "Parsing the pepxml #{file}"
      doc = Nokogiri.XML(File.open(file))
      search_sum = doc.xpath("//xmlns:search_summary") # search for aminoacid_modification and terminal_modification
      doc.xpath('//xmlns:search_hit[@hit_rank=$value]', nil, {value: '1'}).map do |search_hit|
        spec_query = search_hit.parent.parent
        pep_id = PepID.new
        pep_id.scan_num = spec_query.attributes["start_scan"].value.to_i 
        pep_id.charge = spec_query.attributes["assumed_charge"].value.to_i
        mod_info = search_hit.css "modification_info"
        unless mod_info.empty?
          pep_id.mods = mod_info.map do |mod| 
            tmp = mod.attributes.entries
            tmp.map {|f,l| [f, mod.attributes[f].value.to_f]}
          end
        end
        pep_id.precursor_neutral_mass = spec_query.attributes["precursor_neutral_mass"].value.to_f
        pep_id.retention_time = spec_query.attributes["start_scan"].value.to_i
        pep_id.proteins = search_hit.attributes["protein"].value
        pep_id.aaseq = search_hit.attributes["peptide"].value
        pep_id.match_file = file
        pep_id.ion_score = search_hit.css("search_score").to_s[/"ionscore" value="(.*?)"\/>/,1].to_f
        pep_id.mh = search_hit.attributes["calc_neutral_pep_mass"].value.to_f
        pep_id.ppm_error = search_hit.attributes["massdiff"].value.to_f
        pep_id
      end
    end
  end

  class MzmlParser
    # spectrum_array is [mzs, intensities]
    SpectralObject = Struct.new(:spectrum_array, :retention_time) do
      def mzs
        spectrum_array[0]
      end

      def intensities
        spectrum_array[1]
      end

      def peaks(&block)
        spectrum_array[0].zip(spectrum_array[1], &block) 
      end
    end

    def self.join_pepxml_with_mzml_file(pep_ids, file)
      putsv "Joining the pepids with the mzml"
      output_map = []
      # This is the file which you must cache.  
      # Basically, pre-read the file into memory and then play with the structures you've cached inside the loop
      # The problem otherwise is that the random access abilities of Mspire reek havoc on the CPU cycles.

      # I would do something like this:
      spectra = []
      Mspire::Mzml.foreach(file) do |spectrum|
        if spectrum.ms_level == 1
          spectra << SpectralObject.new(spectrum.mzs_and_intensities, spectrum.retention_time)
        end
      end
      putsv "Finished reading the file, and have now closed the MZML"
      # Now,do the following loops using the previous data instead...

      cnt = 0
      pep_ids.each do |pep_id|
        rt_array, mz_array, int_array = [],[],[]

        mass = pep_id.precursor_neutral_mass
        mass_range = ppm_range(mass, PpmThreshold)

        neutral_mass = pep_id.precursor_neutral_mass - Mspire::Mass::Element[:H] # this is (M+H) !!!
        observed_mz = (neutral_mass + (Mspire::Mass::PROTON * pep_id.charge)) / pep_id.charge

        mz_range = ppm_mz_range(observed_mz, PpmThreshold)
        time_range = (pep_id.retention_time-RTThreshold)..(pep_id.retention_time+RTThreshold)

        spectra.each do |spectralobject|
          next if spectralobject.first.empty?
          next unless time_range.include?(spectralobject.retention_time)

          (mzs, intensities) = spectralobject.spectrum_array
          
          rrange = mzs.bsearch_range do |a| 
            mz_range.include?(a) ? 0 : a <=> observed_mz
          end

          indices = rrange.to_a
          best_index = 
            case indices.size
            when 0
              next
            when 1
              indices.first
            else
              # grab the closest peak if there are multiple peaks within the range
              # in case of tie, choose the largest intensity
              indices.sort_by {|index| [(mzs[index] - observed_mz).abs, -intensities[index], mzs[index]] }.first
            end

          next if intensities[best_index] < IntensityThreshold

          rt_array << spectralobject.retention_time
          mz_array << mzs[best_index]
          int_array << intensities[best_index]
        end
        next if mz_array.empty?
        output_map << PeakIDMap.new(pep_id.aaseq, pep_id.proteins, pep_id.mods, pep_id.charge, pep_id.mh, pep_id.ion_score, pep_id.ppm_error, file, pep_id.match_file, PeakIDMapper.determine_stats(rt_array), PeakIDMapper.determine_stats(mz_array), PeakIDMapper.determine_stats(int_array), rt_array, mz_array, int_array)
      end
      output_map
    end
  end
  module CommandLine
    def self.run(pepxml, mzml, opts = {})
      pep_ids = PepxmlParser.parse(pepxml)
      peakIDs = MzmlParser.join_pepxml_with_mzml_file(pep_ids, mzml)
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
  #PeakIDMapper.peakIDs_to_csv(peak_ids)
  PeakIDMapper.plot_peakIDs(peak_ids)
end
