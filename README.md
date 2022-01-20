# FastTree reproduction
A reproduction of the Fast Tree algorithm by [Price et al. (2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2693737/#bib8), which generates a phylogenetic tree from a list of nucleotide sequences.

## Installation
To install the software, navigate to the `fasttree_reproduction` directory and run `pip install .`. This registers the command line interface utility.

## How to run the program
The software can be executed as a command line interface (CLI) or inside an IDE.
### Using CLI (preferred)
Our FastTree reproduction can be run as a command line program. To run using the command line, first ensure that the CLI is correctly installed (see Installation).
The only required argument is the input filename, which should be a file of *.aln type or any text file with similar structure.
Note that the constants in main.py (INPUT_FILE PATH and OUTPUT_FILE_PATH) are ignored when running via the command line.

* Required parameters:
	* file: path to the input file (*.aln) containing the input nucleotide sequences.
* Optional parameters:
	* output_file: path where the resulting tree is stored in Newick format.
	* print_tree: flag whether or not to print the output to the terminal, defaults to False

### Using a Python compatible IDE
The `main.py` module contains constants which can be manually updated to run the program.
* Update the constant `INPUT_FILE_PATH` in the main module. This should be the relative path of the .aln input file containing the sequences based on which a topology should be generated. Default is '/data.test-small.aln'.
* Update the constant `OUTPUT_FILE_PATH` in the main module. This is the file location where the resulting tree is stored in Newick format.
* Run `main.py` using the python version provided in `venv/bin/python`.


## Output
The program writes the generated topology to a string in [Newick format](https://en.wikipedia.org/wiki/Newick_format), which is written by default to OUTPUT_FILE_PATH.
It is also possible to print a text view of the tree using the `--print_tree` flag.
