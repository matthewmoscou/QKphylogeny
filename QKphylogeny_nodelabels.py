#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Extracts node identifiers from a phylogenetic tree in Newick format
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Reads phylogenetic tree in Newick format
	2. Export node labels

"""

## modules
import optparse
from optparse import OptionParser 
import string


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--tree", action="store", type="string", dest="tree", default='', help="Tree file in Newick format")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='', help="Output file for node labels")
(options, args) = parser.parse_args()


tree_file = open(options.tree, 'r')

tree = ''

for line in tree_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line)

	for element in sline:
		tree += element

tree_file.close()

tree = string.replace(tree, '(', '')
nodes = string.split(tree, ',')

node_label_file = open(options.output, 'w')

for node in nodes:
	ID = string.split(node, ':')[0]
	
	node_label_file.write(ID + '\n')

node_label_file.close()
