#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Renames identifiers in a phylogenetic tree based on a user-defined translation table
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Reads phylogenetic tree in Newick format
	2. Reads label translation table
	3. Replaces labels within the tree
	4. Exports tree
"""

## modules
import optparse
from optparse import OptionParser 
import string


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-l", "--labels", action="store", type="string", dest="labels", default='', help="Tab-delimited file containing replacement identifier")
parser.add_option("-t", "--tree", action="store", type="string", dest="tree", default='', help="Tree file in Newick format")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='', help="Output file for tree with updated identifiers")
(options, args) = parser.parse_args()


# read identifer translation file
node_label_file = open(options.labels, 'r')

gene_identifier = {}

for line in node_label_file.readlines():
	sline = string.split(line)

	if len(sline) > 1:
		gene_identifier[sline[0]] = sline[1]

node_label_file.close()


# read phylogenetic tree
tree_file = open(options.tree, 'r')

tree = ''

for line in tree_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line)

	for element in sline:
		tree += element

tree_file.close()


# replace identifiers
for gene in gene_identifier.keys():
	tree = string.replace(tree, gene, gene_identifier[gene])


# export tree
tree_out = open(options.output, 'w')

tree_out.write(tree)

tree_out.close()
