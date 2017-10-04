#! /usr/bin/env python

#This script creates a hi-c matrix out of a mix of input files.
#It will not only get the narmalized signal values, but will also track the features at each pixel

import argparse
import numpy
import scipy
import matplotlib
import math


#fill in the arguements
parser = argparse.ArgumentParser(description='Create a data matrix out of multiple HI-C inputs')

parser.add_argument('-n', '--norm_file', type=str, required=True,
					help='Normalization file that describes if a bin can be used or not')
parser.add_argument('-c', '--count_file', type=str, required=True,
					help='File containing the normalized reads between bin')
parser.add_argument('-o', '--out_dir', type=str, required=True,
					help='The directory where you would like the output to be placed')
parser.add_argument('-l', '--loop_file', type=str, required=True,
					help='File that defines every loop')
parser.add_argument('-chr', '--chromosome',type=int, required=True,
					help='The chromosome you want to look at. File name will include this too')
parser.add_argument('-s', '--bin_size', type=int, required=True,
					help='The size of the bins you are using')

#create an arguement object after reading in the arguments
args = parser.parse_args()

#Add a dictionary that contains all the chromome lengths for humans
#Data obtained from genome hg19 on the ucsc genome browser
human_chrom_sizes_dict = {
	'1':249250621,
	'2':243199373,
	'3':198022430,
	'4':191154276,
	'5':180915260,
	'6':180915260,
	'7':159138663,
	'8':146364022,
	'9':141213431,
	'10':135534747,
	'11':135006516,
	'12':133851895,
	'13':115169878,
	'14':107349540,
	'15':102531392,
	'16':90354753,
	'17':81195210,
	'18':78077248,
	'19':59128983,
	'20':63025520,
	'21':48129895,
	'22':51304566,
	'X':155270560,
	'Y':59373566,
}

### FUNCTIONS ###

def create_matrix( chrom, bin_size, human_chrom_sizes_dict ):
	num_bins = int(math.ceil(human_chrom_sizes_dict[str(chrom)]/bin_size))
	main_dict = {}
	#create a list within a dictionary within a dictionary. The list will be formatted as (signal, #loops, #domains, #f0, #f1, #f2)
	for i in range(num_bins+1) :
		main_dict[i*bin_size] = {}
		for j in range(201) :
			spot = (j*bin_size)+(i*bin_size)
			if spot > human_chrom_sizes_dict[str(chrom)]:
				break
			#this creates the link between one bin and the other "end" bin
			main_dict[i*bin_size][spot] = [0,0,0,0,0,0]
		
	
	return main_dict

def mark_unusable_rows(norm_f, bin_size, chrom, matrix_dict):
	counter = 0;
	norm_file = open(norm_f)
	
	#loop through file
	line = norm_file.readline()
	while line:
		if 'Inf' in line:
			spot = counter*bin_size
			#Gets everything associated with it.
			for ends in matrix_dict[spot].keys():
				matrix_dict[spot][ends][0] = "NA"
			#Now have to go back and look at the bins before it
			for starts in reversed(range(counter+1)):
				starts = starts*bin_size
				if (abs(starts - spot) <= 200*bin_size):
					matrix_dict[starts][spot][0] = "NA"
				
		counter += 1
		line = norm_file.readline()
	
	#A sanity check
	if len(matrix_dict) != counter:
		print 'length of matrix is ', len(matrix_dict), '\nlength of the norm file is ', counter, '\nFilling the remaining bins with NAs'
		for i in range( abs(len(matrix_dict) - counter) ):
			for ends in matrix_dict[(counter+i)*bin_size].keys():
				matrix_dict[(counter+i)*bin_size][ends][0] = "NA"
		

	return(matrix_dict)

def add_signal_values( bin_size, signal_file, matrix_dict):
	signal_file = open(signal_file)
	 
	 #loop through file and add signals
	line = signal_file.readline()
	while line:
		line_vals = line.split("\t")
		start = int(line_vals[0])*bin_size
		end = int(line_vals[1])*bin_size
		if abs(end - start) <= 200*bin_size:
			matrix_dict[start][end][0] = float(line_vals[2])
		line = signal_file.readline()
	
	return(matrix_dict)

def add_features( chrom, loops_f, matrix_dict, bin_size ):
	#hold the two points that define the location of the loop
	loop_list_start = []
	loop_list_end = []
	#get the chrom string to test if the chrom is there
	chrom = 'chr' + str(chrom)
	
	#go through the loop file making sure to only grab loops from the correct chromosome
	loops_file = open(loops_f)
	line = loops_file.readline()
	while line:
		line_vals = line.split("\t")
		if( chrom in line_vals[0] and len(chrom) == len(line_vals[0]) and
		   chrom in line_vals[3] and len(chrom) == len(line_vals[3]) ):
			loop_list_start.append(int(line_vals[1]))
			loop_list_end.append(int(line_vals[4]))
	
		line = loops_file.readline()
		
	#go through and add loops to the matrix
	for i in range(len(loop_list_start)):
		if abs(loop_list_end[i] - loop_list_start[i])/bin_size <= 200:
			matrix_dict[loop_list_start[i]][loop_list_end[i]][1] += 1
	
	#get and set domains
	dom_starts, dom_ends = get_domain_locations(loop_list_start, loop_list_end, bin_size)
	for i in range(len(dom_starts)):
		matrix_dict[dom_starts[i]][dom_ends[i]][2] += 1 
	
	#get and set flares ( the return from get flare location 0&1 is two away 2&3 is one away and 4&5 is 0 away)
	lists_of_lists_starts_stops = get_flare_locations(loop_list_start, loop_list_end, bin_size)
	count = 5;
	for i in range(0, len(lists_of_lists_starts_stops), 2):
		for j in range(len(lists_of_lists_starts_stops[i])):
			matrix_dict[lists_of_lists_starts_stops[i][j]][lists_of_lists_starts_stops[i+1][j]][count] += 1
		count -= 1
	
	return(matrix_dict)

def get_domain_locations( starts, ends, bin_size ):
	domain_position_starts = []
	domain_position_ends = []
	
	#Get domain structure for each loop
	for i in range(len(starts)):
		#I add bin size because range does not include the end value
		for j in range(starts[i], ends[i]+bin_size, bin_size):
			for k in reversed(range(j,ends[i]+bin_size, bin_size)):
				domain_position_starts.append(j)
				domain_position_ends.append(k)
				
	
	return( domain_position_starts, domain_position_ends)

def get_flare_locations( starts, ends, bin_size ):
	two_away_flares_s = []
	two_away_flares_e = []
	one_away_flares_s = []
	one_away_flares_e = []
	line_flares_s = []
	line_flares_e = []
	
	
	for i in range(len(starts)):
		#get the two away flares
		two_out_s = starts[i] - 2*bin_size
		two_out_e = ends[i] + 2*bin_size
		two_in_s = starts[i] + 2*bin_size
		two_in_e = ends[i] - 2*bin_size
		
		two_s, two_e = get_flare_vals_in_between(two_out_s, two_out_e, bin_size)
		in_s, in_e = get_flare_vals_in_between(two_in_s, two_in_e, bin_size)
		two_away_flares_s.extend(two_s.extend(in_s))
		two_away_flares_e.extend(two_e.extend(in_e))
		
		#get the one away flares
		one_out_s = starts[i] - 1*bin_size
		one_out_e = ends[i] + 1*bin_size
		one_in_s = starts[i] + 1*bin_size
		one_in_e = ends[i] - 1*bin_size
		
		one_s, one_e = get_flare_vals_in_between(one_out_s, one_out_e, bin_size)
		in_s, in_e = get_flare_vals_in_between(one_in_s, one_in_e, bin_size)
		one_away_flares_s.extend(one_s.extend(in_s))
		one_away_flares_e.extend(one_e.extend(in_e))
		
		#get the flare at 0 away
		flare_s, flare_e = get_flare_vals_in_between(starts[i], ends[i], bin_size)
		line_flares_s.extend(flare_s)
		line_flares_e.extend(flare_e)
		
		
	return([two_away_flares_s, two_away_flares_e, one_away_flares_s, one_away_flares_e, line_flares_s, line_flares_e])

def get_flare_vals_in_between( start, end, bs ):
	start_tot = []
	end_tot = []
	for i in range(start, end+bs, bs):
		start_tot.append(start)
		end_tot.append(i)
	
	for i in reversed(range(start,end, bs)):
		start_tot.append(i)
		end_tot.append(end)
	
	return(start_tot, end_tot)

def print_out_matrix( out_dir, matrix_dict, chrom ):
	out_file = out_dir + "/" + "count_matrix_chr" + str(chrom) + ".txt"
	out = open(out_file, 'w')
	
	out.write("Chr\tBin1\tBin2\tSignal\tLoops\tDomains\tFlare_0\tFlare_1\tFlare_2\n")
	
	for start in matrix_dict.keys():
		for end in matrix_dict[start].keys():
			out.write(chrom, "\t", start, "\t", end, "\t", matrix_dict[start][stop][0], "\t", matrix_dict[start][stop][1],
			"\t", matrix_dict[start][stop][2], "\t", matrix_dict[start][stop][3], "\t", matrix_dict[start][stop][4], "\t",
			matrix_dict[start][stop][5], "\n")
	
	return(1)

### MAIN ###
#Create a dictionary the size of the chromosome passed to a max distance off 200 bins
mat_dict = create_matrix(args.chromosome, args.bin_size, human_chrom_sizes_dict)

#put NA's in the locations where the reads were not enough
mat_dict = mark_unusable_rows(args.norm_file, args.bin_size, args.chromosome, mat_dict)

#update the signal at each pixel
mat_dict = add_signal_values(args.bin_size, args.count_file, mat_dict)

#add the features at each pixel
mat_dict = add_features(args.chromosome, args.loop_file, mat_dict, args.bin_size)

#print out the matrix and features ( 1 row per pixel )
print_out_matrix(args.out_dir, mat_dict)
