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
"""

## modules
import commands
import optparse
from optparse import OptionParser 
import sets
import string


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
parser.add_option("-n", "--nonredundant", action="store_true", dest="nr", default=False, help="Remove redundant sequences from alignment")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='', help="Output file for tree with updated identifiers")
parser.add_option("-r", "--reduce", action="store_true", dest="reduce", default=False, help="Reduce alignment to polymorphic sites")
parser.add_option("-t", "--type", action="store", type="string", dest="type", default='nucleotide', help="Specify nucleotide or protein")
(options, args) = parser.parse_args()


## Import phylip file
# check for identical sequence and notify user
print 'Import Phylip file'
phylip_input = open(options.alignment, 'r')

truth = False
sequence_ID = {}
ID_sequence = {}
IDs = []

for line in phylip_input.readlines():
	sline = string.split(line)

	if truth:	
		ID_sequence[sline[0]] = sline[1]

		if sline[1] not in sequence_ID.keys():
			sequence_ID[sline[1]] = [sline[0]]
		else:
			sequence_ID[sline[1]].append(sline[0])
			print '\t' + 'Identical sequence for:', sequence_ID[sline[1]][0], sline[0] 

		IDs.append(sline[0])

	truth = True

IDs.sort()

phylip_input.close()

print 'Unique sequences    :', len(sequence_ID)
print 'Number of genes     :', len(sets.Set(IDs))
print 'Length of alignment :', len(sequence_ID.keys()[0])
print

## Assess alignment breadth coverage, requires user-specified breadth coverage of the alignment
print 'Assess alignment breath coverage'
alignment_length = len(sequence_ID.keys()[0])

truth = False
position_sequence = []

threshold_IDs = []
nr_threshold_IDs = []
	
for sequence in sequence_ID.keys():
	first_sequence = -1
	last_sequence = -1

	alignment = sequence

	if len(position_sequence) == 0:
		for position_index in range(len(alignment)):
			position_sequence.append([])

	for position_index in range(len(alignment)):
		position_sequence[position_index].append(alignment[position_index])

		if first_sequence < 0:
			if options.type == 'nucleotide':
				if alignment[position_index] in nucleotides:
					first_sequence = position_index
			elif options.type == 'protein':
				if alignment[position_index] in amino_acids:
					first_sequence = position_index

		if options.type == 'nucleotide':
			if alignment[position_index] in nucleotides:
				last_sequence = position_index
		elif options.type == 'protein':
			if alignment[position_index] in amino_acids:
				last_sequence = position_index

	if ((float(last_sequence) - float(first_sequence)) / float(alignment_length)) > options.breadth:
		nr_threshold_IDs.append(sequence_ID[sequence][0])

		for gene in sequence_ID[sequence]:
			threshold_IDs.append(gene)


# if coverage used, limit to genes that meet threshold
if options.breadth > 0:
	selected_genes = threshold_IDs

	print '\t' + 'Number of sequences meeting breadth threshold :', len(threshold_IDs)
	print '\t' + 'Number of non-redundant sequences meeting breadth threshold :', len(nr_threshold_IDs)
	print
else:
	selected_genes = []

	for sequence in sequence_ID.keys():
		for ID in sequence_ID[sequence]:
			selected_genes.append(ID)

# Export alignment coverage 
print 'Export alignment coverage'
alignment_coverage_file = open(options.alignment + '_alignment_coverage.txt', 'w')
	
alignment_coverage_file.write('position' + '\t' + 'coverage' + '\n')

position_selection = []

for position_index in range(len(alignment)):
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


print 'Number of positions meeting depth threshold of', options.depth, 'is', len(position_selection)


# Visualization of alignment coverage
print 'Visualization of alignment coverage'

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

commands.getstatusoutput('R --vanilla < ' + options.alignment + '_alignment_coverage.R')


## Reduce evaluated sites to polymorphic sites
# check site frequency, evaluate if variation exists
if options.reduce:
	print 'Reduce evaluated sites'

	ID_sequence_polymorphic_sites = {}
	
	polymorphic_sites = 0
	monomorphic_sites = 0

	for gene in selected_genes:
		ID_sequence_polymorphic_sites[gene] = ''
	
	for position_index in range(alignment_length):
		evaluated_site_sequence = []
	
		for gene in selected_genes:
			evaluated_site_sequence.append(ID_sequence[gene][position_index].upper())
		
		if len(sets.Set(evaluated_site_sequence) - sets.Set(['X', 'N', '-', '?', 'O'])) > 1:
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
	
	print 'Number of sites that were:'
	print '\t' + 'Polymorphic:', polymorphic_sites
	print '\t' + 'Monomorphic:', monomorphic_sites
	
	selected_genes = []
	
	for sequence in sequence_ID_polymorphic_sites.keys():
		identical_sequence_IDs = sequence_ID_polymorphic_sites[sequence]
		identical_sequence_IDs.sort()
	
		selected_genes.append(identical_sequence_IDs[0])
	
	print 'Unique sequences :', len(sequence_ID_polymorphic_sites.keys())
	print 'Number of genes  :', len(selected_genes)
	print
else:
	ID_sequence_polymorphic_sites = ID_sequence
	sequence_ID_polymorphic_sites = sequence_ID


## output files
# digitial: non-redundant correction, optional coverage requirement
# crosslist: redundant identifiers
print 'output files'

phylip_digital = open(options.output, 'w')

print 'Genotypes in alignment:', len(sequence_ID_polymorphic_sites.keys())
print 'Length of alignment:   ', len(sequence_selection(sequence_ID_polymorphic_sites.keys()[0], position_selection))

# export header for phylip file

# export identifier and sequence, final check of gene identifier length
if options.nr:
	phylip_digital.write(' ' + str(len(selected_genes)) + ' ' + str(len(sequence_selection(ID_sequence[selected_genes[0]], position_selection))) + '\n')
	phylip_digital_crosslist = open(options.output + '_crosslist.txt', 'w')
	
	for sequence in sequence_ID_polymorphic_sites.keys():
		phylip_digital.write(sequence_ID_polymorphic_sites[sequence][0] + '   ' + sequence_selection(sequence, position_selection) + '\n')

		phylip_digital_crosslist.write(sequence_ID_polymorphic_sites[sequence][0])

		if len(sequence_ID_polymorphic_sites[sequence]) > 1:
			for ID in sequence_ID_polymorphic_sites[sequence][1:]:
				phylip_digital_crosslist.write('\t' + ID)
	
		phylip_digital_crosslist.write('\n')
	
	phylip_digital_crosslist.close()
else:
	phylip_digital.write(' ' + str(len(selected_genes)) + ' ' + str(len(sequence_selection(ID_sequence_polymorphic_sites[selected_genes[0]], position_selection))) + '\n')

	for gene in selected_genes:
		phylip_digital.write(gene + '   ' + sequence_selection(ID_sequence_polymorphic_sites[gene], position_selection) + '\n')


phylip_digital.close()
