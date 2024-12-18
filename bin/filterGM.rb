#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# filterGM
FILTERGMVER = "0.3.3"
# Michael G. Campana and Ellie E. Armstrong, 2020-2023
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2023. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. Bioinformatics. 39: btac784.
# 10.1093/bioinformatics/btac784.
#----------------------------------------------------------------------------------------

# Script to filter GenMap bed output for inclusion/exclusion of regions using VCFtools
# ARGV[0] is input GenMap bed, ARGV[1] is cutoff float, ARGV[2] is optional, but 'exclude' converts to exclude BED (from default include)
# This script has been modified from the RatesTools version by moving gz_file_open and format_splash from denovolib.rb here

#-----------------------------------------------------------------------------------------------
# From BaitsTools 1.6.8: Campana 2018
def gz_file_open(file) # Determine whether input file is gzipped or not and set method to open it
	if file[-3..-1] == ".gz"
		yield Zlib::GzipReader.open(file)
	else
		yield File.open(file)
	end
end
#-----------------------------------------------------------------------------------------
# From RatesTools 1.2.1: Armstrong & Campana 2023
def format_splash(cmd, version, cmdline) # Format output for basic script help screens
	puts "\033[1m#{cmd} #{version}\033[0m"
	puts "\nUsage: ruby #{cmd}.rb #{cmdline}"
end
#-----------------------------------------------------------------------------------------
if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('filterGM', FILTERGMVER, '<in_GenMap.bed[.gz]> <cutoff> [exclude] > <out.bed>')
else
	gz_file_open(ARGV[0]) do |f1|
		while line = f1.gets
			line_arr = line.split
			if ARGV[2] == "exclude"
				puts line if line_arr[4].to_f < ARGV[1].to_f
			else
				puts line if line_arr[4].to_f >= ARGV[1].to_f
			end
		end
	end
end
