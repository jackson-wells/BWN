# Bio-Sequence Search using the Burrows-Wheeler Transform

## Current Version:
 - Index and search protein or nucleotide sequences
 - Index a single or multi-fasta file, or a string input
 - Search index using string, single or multi-fasta

## Additions To Come:
 - Induced sorting of SA for calculation in "linear" time
 - Support for multiple scoring matrices

### To Compile:
 1. Clone this repository to your local machine.
 2. Run "make" in the protein or nucl directory. The decision between nucl and protein searches should depend upon the format of your sequence files.
 3. Excutables will now available in the directory that "make" was run in
    - Use the "-h" flag during execution of programs for additional help
