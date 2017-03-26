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
import optparse
from optparse import OptionParser 
import sets
import string


## global variables
nucleotides = ['A', 'C', 'G', 'T']
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-a", "--alignment", action="store", type="string", dest="alignment", default='', help="Alignment file in Phylip format")
parser.add_option("-c", "--coverage", action="store", type="float", dest="coverage", default=-1.0, help="Coverage required for inclusion of sequence, coverage above 0.0 but equal to or lower than 1.0")
parser.add_option("-m", "--missing", action="store", type="float", dest="missing", default=1.0, help="Missing data allowed (0.0 (none) to 1.0 (all))")
parser.add_option("-n", "--nonredundant", action="store_true", dest="nr", default=False, help="Remove redundant sequences from alignment")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='', help="Output file for tree with updated identifiers")
parser.add_option("-r", "--reduce", action="store_true", dest="reduce", default=False, help="Reduce alignment to polymorphic sites")
parser.add_option("-t", "--type", action="store", type="string", dest="type", default='nucleotide', help="Specify nucleotide or protein")
(options, args) = parser.parse_args()


## import phylip file
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

print 'Unique sequences :', len(sequence_ID)
print 'Number of genes  :', len(sets.Set(IDs))
print

## assess alignment coverage, requires user-specified coverage of the alignment
print 'assess alignment coverage'
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

	if ((float(last_sequence) - float(first_sequence)) / float(alignment_length)) > options.coverage:
		nr_threshold_IDs.append(sequence_ID[sequence][0])

		for gene in sequence_ID[sequence]:
			threshold_IDs.append(gene)


# if coverage used, limit to genes that meet threshold
if options.coverage > 0:
	selected_genes = threshold_IDs

	print '\t' + 'Number of sequences meeting threshold :', len(threshold_IDs)
	print '\t' + 'Number of non-redundant sequences meeting threshold :', len(nr_threshold_IDs)
	print
else:
	selected_genes = []

	for sequence in sequence_ID.keys():
		for ID in sequence_ID[sequence]:
			selected_genes.append(ID)


# export alignment coverage 
print 'export alignment coverage'
alignment_coverage_file = open(options.alignment + '_alignment_coverage.txt', 'w')
	
alignment_coverage_file.write('position' + '\t' + 'coverage' + '\n')

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


## reduce evaluated sites
# check site frequency, evaluate if variation exists
if options.reduce:
	print 'reduce evaluated sites'

	ID_sequence_polymorphic_sites = {}
	
	polymorphic_sites = 0
	monomorphic_sites = 0

	missing_data_threshold_pass = 0
	missing_data_threshold_fail = 0
	
	for gene in selected_genes:
		ID_sequence_polymorphic_sites[gene] = ''
	
	for position_index in range(alignment_length):
		evaluated_site_sequence = []
	
		for gene in selected_genes:
			evaluated_site_sequence.append(ID_sequence[gene][position_index].upper())
		
		if len(sets.Set(evaluated_site_sequence) - sets.Set(['X', 'N'])) > 1:
			polymorphic_sites += 1
	
			missing_data = 0.0

			if options.type == 'nucleotide':
				for gene in selected_genes:
					if ID_sequence[gene][position_index].upper() in ['N', 'O', 'X', '?', '-']:	
						missing_data += 1
			elif options.type == 'protein':
				for gene in selected_genes:
					if ID_sequence[gene][position_index].upper() in ['X', '?', '­']:
						missing_data += 1

			if ((float(len(evaluated_site_sequence)) - missing_data) / float(len(evaluated_site_sequence))) >= (1.0 - options.missing):
				missing_data_threshold_pass += 1

				for gene in selected_genes:
					ID_sequence_polymorphic_sites[gene] += ID_sequence[gene][position_index].upper()
			else:
				missing_data_threshold_fail += 1
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
	
	print 'Sites that passed missing data threshold of ' + str(options.missing * 100.0) + '%:'
	print '\t' + 'Pass:', missing_data_threshold_pass
	print '\t' + 'Fail:', missing_data_threshold_fail

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
print 'Length of alignment:   ', len(sequence_ID_polymorphic_sites.keys()[0])

# export header for phylip file

# export identifier and sequence, final check of gene identifier length
if options.nr:
	phylip_digital.write(' ' + str(len(selected_genes)) + ' ' + str(len(ID_sequence[selected_genes[0]])) + '\n')
	phylip_digital_crosslist = open(options.output + '_crosslist.txt', 'w')
	
	for sequence in sequence_ID_polymorphic_sites.keys():
		phylip_digital.write(sequence_ID_polymorphic_sites[sequence][0] + '   ' + sequence + '\n')

		phylip_digital_crosslist.write(sequence_ID_polymorphic_sites[sequence][0])

		if len(sequence_ID_polymorphic_sites[sequence]) > 1:
			for ID in sequence_ID_polymorphic_sites[sequence][1:]:
				phylip_digital_crosslist.write('\t' + ID)
	
		phylip_digital_crosslist.write('\n')
	
	phylip_digital_crosslist.close()
else:
	phylip_digital.write(' ' + str(len(selected_genes)) + ' ' + str(len(ID_sequence_polymorphic_sites[selected_genes[0]])) + '\n')

	for gene in selected_genes:
		phylip_digital.write(gene + '   ' + ID_sequence_polymorphic_sites[gene] + '\n')


phylip_digital.close()
