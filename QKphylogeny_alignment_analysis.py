#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] -a phylip_file
Assesses an alignment and generates a non-redundant curated alignment based on user parameters
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Imports an alignment in Phylip format
	2. Assesses for redundancy in the alignment, reports to user
	3. Reduces the alignment to only polymorphic sites (removing positions where only missing data exists (N or X))
	4. Exports file

Tested with Python 3.9.7 on MacOS Monterey 12.6
"""

## modules
import optparse
from optparse import OptionParser 

import subprocess


## global variables
nucleotides = ['A', 'C', 'G', 'T']
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


# functions
def sequence_selection(sequence, position_index):
	selected_sequence = ''

	for index in range(len(sequence)):
		if index in position_index:
			selected_sequence += sequence[index]

	return selected_sequence


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-a", "--alignment", action="store", type="string", dest="alignment", default='', help="Alignment file in Phylip format")
parser.add_option("-b", "--breadth", action="store", type="float", dest="breadth", default=-1.0, help="Breadth coverage required for inclusion of sequence, coverage above 0.0 but equal to or lower than 1.0. Breadth refers to the start and end of a sequence relative to the entire alignment.")
parser.add_option("-d", "--depth", action="store", type="float", dest="depth", default=-1.0, help="Depth coverage required for inclusion of sequence, coverage above 0.0 but equal to or lower than 1.0. Depth coverage refers to the number of characters with data (i.e. no missing data or gaps).")
parser.add_option("-m", "--missing", action="store", type="float", dest="missing", default=-1.0, help="Alternative to breadth coverage, exclude sequences based on degree of missing data with a range from 0.0 to 1.0.")
parser.add_option("-n", "--nonredundant", action="store_true", dest="nr", default=False, help="Remove redundant sequences from alignment")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='', help="Output file for tree with updated identifiers")
parser.add_option("-r", "--reduce", action="store_true", dest="reduce", default=False, help="Reduce alignment to polymorphic sites")
parser.add_option("-t", "--type", action="store", type="string", dest="type", default='', help="Specify nucleotide or protein")
parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="Don't print status messages to stdout")
(options, args) = parser.parse_args()


## QC input parameters
if options.type not in ['nucleotide', 'protein']:
	print('Type parameter is incorrect, select either nucleotide or protein')
	exit()

## Import phylip file
# check for identical sequence and notify user
if options.verbose:
	print('Import Phylip file')

phylip_input = open(options.alignment, 'r')

truth = False
sequence_ID = {}
ID_sequence = {}
IDs = []

for line in phylip_input.readlines():
	sline = line.split()

	if truth:	
		ID_sequence[sline[0]] = sline[1]

		if sline[1] not in sequence_ID.keys():
			sequence_ID[sline[1]] = [sline[0]]
		else:
			sequence_ID[sline[1]].append(sline[0])

			if options.verbose:
				print('\t' + 'Identical sequence for: ' + sequence_ID[sline[1]][0] + ' ' + sline[0])

		IDs.append(sline[0])

	truth = True

IDs.sort()

phylip_input.close()

if options.verbose:
	print('Unique sequences    : ' + str(len(sequence_ID)))
	print('Number of genes     : ' + str(len(set(IDs))))
	print('Length of alignment : ' + str(len(list(sequence_ID.keys())[0])))
	print()

## Assess alignment breadth coverage, requires user-specified breadth coverage of the alignment
if options.verbose:
	print('Assess alignment breath coverage')

alignment_length = len(list(sequence_ID.keys())[0])

truth = False
position_sequence = []

threshold_IDs = []
nr_threshold_IDs = []

present_threshold_IDs = []
nr_present_threshold_IDs = []

for sequence in sequence_ID.keys():
	first_sequence = -1
	last_sequence = -1

	if len(position_sequence) == 0:
		for position_index in range(len(sequence)):
			position_sequence.append([])

	for position_index in range(len(sequence)):
		position_sequence[position_index].append(sequence[position_index])

		if first_sequence < 0:
			if options.type == 'nucleotide':
				if sequence[position_index] in nucleotides:
					first_sequence = position_index
			elif options.type == 'protein':
				if sequence[position_index] in amino_acids:
					first_sequence = position_index

		if options.type == 'nucleotide':
			if sequence[position_index] in nucleotides:
				last_sequence = position_index
		elif options.type == 'protein':
			if sequence[position_index] in amino_acids:
				last_sequence = position_index

	if ((float(last_sequence) - float(first_sequence)) / float(alignment_length)) > options.breadth:
		nr_threshold_IDs.append(sequence_ID[sequence][0])

		for gene in sequence_ID[sequence]:
			threshold_IDs.append(gene)

	if (sequence.count('-') / float(alignment_length)) < (1.0 - options.missing):
		nr_present_threshold_IDs.append(sequence_ID[sequence][0])

		for gene in sequence_ID[sequence]:
			present_threshold_IDs.append(gene)

# if coverage used, limit to genes that meet threshold
if options.breadth > 0:
	selected_genes = threshold_IDs

	if options.reduce:
		selected_genes = nr_threshold_IDs

	if options.verbose:
		print('\t' + 'Number of sequences meeting breadth threshold : ' + str(len(threshold_IDs)))
		print('\t' + 'Number of non-redundant sequences meeting breadth threshold : ' + str(len(nr_threshold_IDs)))
		print()
elif options.missing > 0:
	selected_genes = present_threshold_IDs

	if options.reduce:
		selected_genes = nr_present_threshold_IDs

	if options.verbose:
		print('\t' + 'Number of sequences meeting missing threshold : ' + str(len(present_threshold_IDs)))
		print('\t' + 'Number of non-redundant sequences meeting missing threshold : ' + str(len(nr_present_threshold_IDs)))
		print()
else:
	selected_genes = []

	for sequence in sequence_ID.keys():
		for ID in sequence_ID[sequence]:
			selected_genes.append(ID)

# Export alignment coverage 
if options.verbose:
	print('Export alignment coverage')

alignment_coverage_file = open(options.alignment + '_alignment_coverage.txt', 'w')
alignment_coverage_file.write('position' + '\t' + 'coverage' + '\n')

position_selection = []

for position_index in range(alignment_length):
	positive_sites = 0

	for singlesequence in position_sequence[position_index]:
		if options.type == 'nucleotide':
			if singlesequence in nucleotides:
				positive_sites += 1

		elif options.type == 'protein':
			if singlesequence in amino_acids:
				positive_sites += 1
	
	alignment_coverage_file.write(str(position_index + 1) + '\t' + str(float(positive_sites) / float(len(position_sequence[position_index]))) + '\n')

	if (float(positive_sites) / float(len(position_sequence[position_index]))) > options.depth:
		position_selection.append(position_index)

if options.verbose:
	print('Number of positions meeting depth threshold of ' + str(options.depth) + ' is ' + str(len(position_selection)))

# Export present data coverage
present_coverage_file = open(options.alignment + '_present_coverage.txt', 'w')
present_coverage_file.write('ID' + '\t' + 'coverage' + '\n')

for sequence in sequence_ID.keys():
	present_coverage_file.write(str(list(sequence_ID.keys()).index(sequence)) + '\t' + str(alignment_length - sequence.count('-')) + '\n')

present_coverage_file.close()

# Visualization of alignment coverage
if options.verbose:
	print('Visualization of alignment coverage')

data_visualization_file = open(options.alignment + '_alignment_coverage.R', 'w')
data_visualization_file.write('library(ggplot2)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('alignment = read.table(file="' + options.alignment + '_alignment_coverage.txt", header=T)' + '\n')
data_visualization_file.write('alignment = data.frame(alignment)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('png(file="' + options.alignment + '_alignment_coverage.png", height=600, width=600)' + '\n')
data_visualization_file.write('ggplot(alignment, aes(coverage)) + geom_histogram(binwidth = 0.05) + scale_x_continuous(limits = c(-0.1, 1.1))' + '\n')
data_visualization_file.write('dev.off()' + '\n')
data_visualization_file.write('' + '\n')
data_visualization_file.close()

data_visualization_file = open(options.alignment + '_present_coverage.R', 'w')
data_visualization_file.write('library(ggplot2)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('present = read.table(file="' + options.alignment + '_present_coverage.txt", header=T)' + '\n')
data_visualization_file.write('present = data.frame(present)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('png(file="' + options.alignment + '_present_coverage.png", height=600, width=600)' + '\n')
data_visualization_file.write('ggplot(present, aes(coverage)) + geom_histogram() + scale_x_continuous(limits = c(-5, ' + str(alignment_length * 1.05) + '))' + '\n')
data_visualization_file.write('dev.off()' + '\n')
data_visualization_file.write('' + '\n')
data_visualization_file.close()

subprocess.getstatusoutput('R --vanilla < ' + options.alignment + '_alignment_coverage.R')
subprocess.getstatusoutput('R --vanilla < ' + options.alignment + '_present_coverage.R')


## Reduce evaluated sites to polymorphic sites
# check site frequency, evaluate if variation exists
if options.reduce:
	if options.verbose:
		print('Reduce evaluated sites')

	ID_sequence_polymorphic_sites = {}
	
	polymorphic_sites = 0
	monomorphic_sites = 0

	for gene in selected_genes:
		ID_sequence_polymorphic_sites[gene] = ''
	
	for position_index in range(alignment_length):
		evaluated_site_sequence = []
	
		for gene in selected_genes:
			evaluated_site_sequence.append(ID_sequence[gene][position_index].upper())
		
		if len(set(evaluated_site_sequence) - set(['X', 'N', '-', '?', 'O'])) > 1:
			polymorphic_sites += 1
	
			for gene in selected_genes:
				ID_sequence_polymorphic_sites[gene] += ID_sequence[gene][position_index].upper()
		else:
			monomorphic_sites += 1
	
	sequence_ID_polymorphic_sites = {}
	
	for ID in ID_sequence_polymorphic_sites.keys():
		if ID_sequence_polymorphic_sites[ID] not in sequence_ID_polymorphic_sites.keys():
			sequence_ID_polymorphic_sites[ID_sequence_polymorphic_sites[ID]] = []
	
		sequence_ID_polymorphic_sites[ID_sequence_polymorphic_sites[ID]].append(ID)
	
	if options.verbose:
		print('Number of sites that were:')
		print('\t' + 'Polymorphic: ' + '\t' + polymorphic_sites)
		print('\t' + 'Monomorphic: ' + '\t' + monomorphic_sites)
	
	selected_genes = []
	
	for sequence in sequence_ID_polymorphic_sites.keys():
		identical_sequence_IDs = sequence_ID_polymorphic_sites[sequence]
		identical_sequence_IDs.sort()
	
		selected_genes.append(identical_sequence_IDs[0])
	
	if options.verbose:
		print('Unique sequences : ' + str(len(sequence_ID_polymorphic_sites.keys())))
		print('Number of genes  : ' + str(len(selected_genes)))
		print()
else:
	ID_sequence_polymorphic_sites = ID_sequence
	sequence_ID_polymorphic_sites = sequence_ID


## output files
# digitial: non-redundant correction, optional coverage requirement
# crosslist: redundant identifiers
if options.verbose:
	print('output files')

phylip_digital = open(options.output, 'w')

if options.verbose:
	print('Genotypes in alignment: ' + str(len(sequence_ID_polymorphic_sites.keys())))
	print('Length of alignment:    ' + str(len(sequence_selection(list(sequence_ID_polymorphic_sites.keys())[0], position_selection))))

# export identifier and sequence, final check of gene identifier length
if options.nr:
	# assess redundancy in final alignment based on breadth and depth calculation
	sequence_ID_output = {}

	for gene in selected_genes:
		sequence = sequence_selection(ID_sequence_polymorphic_sites[gene], position_selection)
		
		if sequence not in sequence_ID_output.keys():
			sequence_ID_output[sequence] = []

		sequence_ID_output[sequence].append(gene)

	phylip_digital.write(' ' + str(len(selected_genes)) + ' ' + str(len(sequence_ID_output.keys()[0])) + '\n')
	phylip_digital_crosslist = open(options.output + '_crosslist.txt', 'w')
	
	for sequence in sequence_ID_output.keys():
		phylip_digital.write(sequence_ID_output[sequence][0] + '   ' + sequence + '\n')

		phylip_digital_crosslist.write(sequence_ID_output[sequence][0])

		if len(sequence_ID_output[sequence]) > 1:
			for ID in sequence_ID_output[sequence][1:]:
				phylip_digital_crosslist.write('\t' + ID)
	
		phylip_digital_crosslist.write('\n')
	
	phylip_digital_crosslist.close()
else:
	phylip_digital.write(' ' + str(len(selected_genes)) + ' ' + str(len(sequence_selection(ID_sequence_polymorphic_sites[selected_genes[0]], position_selection))) + '\n')

	for gene in selected_genes:
		phylip_digital.write(gene + '   ' + sequence_selection(ID_sequence_polymorphic_sites[gene], position_selection) + '\n')


phylip_digital.close()
