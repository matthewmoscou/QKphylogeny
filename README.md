# QKphylogeny
A set of scripts for phylogenetic tree assessment and editing.

## Scripts
<i>QKphylogeny_nodelabels.py</i> exports the genes/protein identifiers in a phylogenetic tree (Newick format).

<i>QKphylogeny_rename_nodes.py</i> will convert genes/protein identifiers in a phylogenetic tree based on a tab-delimited translation table.

<i>QKphylogeny_alignment_analysis.py</i> assesses the quality of a multiple sequence alignment. Parameters exist for the removal of redundant sequences, reduction of the alignment based on breadth and depth coverage, and reduction of alignment to polymorphic sites.

## Examples
Export the proteins in a phylogenetic tree in Newick format. In this case, using the neighbor joining tree for thioredoxins in <i>Arabidopsis thaliana</i>.
```bash
python QKphylogeny_nodelabels.py -t examples/AT_TRX_phylogenetic_tree_NJ.newick -o AT_TRX_phylogenetic_tree_NJ_proteins.txt
```

Convert the protein identifiers from one format to another. In this case, convert long gene format (i.e. AT1G03680.1 to AtTRXm1).
```bash
python QKphylogeny_rename_nodes.py -t examples/AT_TRX_phylogenetic_tree_NJ.newick -l examples/AT_TRX_abbreviations.txt -o AT_TRX_phylogenetic_tree_NJ_renamed.newick
```

Restrict the nucleotide alignment to evaluated sites that have more than 80% of data present
```bash
python QKphylogeny_alignment_analysis.py -a test.phy -d 0.2 -t nucleotide -o test_missing20.phy
```
