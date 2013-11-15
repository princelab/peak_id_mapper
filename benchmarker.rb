require 'benchmark'

Pfile = "ptest.pepxml"
Mfile = "test.mzML"


puts RUBY_DESCRIPTION
puts Benchmark.measure {5.times { system "ruby lib/peak_id_mapper.rb #{Pfile} #{Mfile}" } }
