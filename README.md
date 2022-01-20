# FastTree reproduction
A reproduction of the Fast Tree algorithm by Price et al. (2009), which generates a phylogenetic tree from a list of nucleotide sequences.

## Installation
To install the software, navigate to the fasttree directory and run `pip install .`. This registers

## How to run the program

### Using a Python compatible IDE
* Update the constant INPUT_FILE_PATH in main.py. This should be the relative path of the .aln input file containing the sequences based on which a topology should be generated. Default is '/data.test-small.aln'.
* Update the constant OUTPUT_FILE_PATH in main.py. This is the file location where the resulting tree is stored in Newick format.
* Run main.py

### Using command line interface
Our Fast Tree reproduction can be run as a command line program. To run using the command line, 

Note that the constants in main.py (INPUT_FILE PATH and OUTPUT_FILE_PATH) are ignored when running via the command line.

* Required parameters:
	* file: path to the input file (*.aln) containing the input sequences.
* Optional parameters:
	* output_file: path where the resulting tree is stored in Newick format.
	* print_tree: flag whether or not to print the output to the terminal, defaults to False

## Output
The program writes the generated topology to a string in [Newick format](https://en.wikipedia.org/wiki/Newick_format), which is written by default to OUTPUT_FILE_PATH.
