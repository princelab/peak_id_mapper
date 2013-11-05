require 'nokogiri'	
require 'pry'

PepID = Struct.new(:precursor_neutral_mass, :aaseq, :charge, :scan_num, :retention_time, :protein)
module PeakIDMapper
  class PepxmlParser
    # This is from the run_compare code I wrote
    def self.parse(file)
      doc = Nokogiri.XML(File.open(file))
      pepids = []
      binding.pry
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
      pepids
    end

  end
  class MzmlParser
    def self.parse(file)
    end
  end
  PeakIDMap = Struct.new(:aaseq, :charge, :mods, :mz_mean, :total_intensity, :data) # Data = [rt_array, mz_array, int_array]
  module CommandLine
    def self.run(pepxml, mzml, opts = {})
      pepdata = PepxmlParser.parse(pepxml)
      mzdata = MzmlParser.parse(mzml)
      join_pepxml_and_mzml_data(pepdata, mzdata)
    end
  end
end


if __FILE__ == $0
  args = ARGV
  PeakIDMapper::CommandLine.run(*args)
end
