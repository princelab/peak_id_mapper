require 'benchmark'

Pfile = "ptest.pepxml"
Mfile = "test.mzML"

REPEATS = 10

puts RUBY_DESCRIPTION
puts Benchmark.measure {REPEATS.times { system "ruby lib/peak_id_mapper.rb #{Pfile} #{Mfile}" } }
