import click

from main import run

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', '--output_file', default="output_tree.txt",
              help='filename of the output to be generated, defaults to output_tree.txt')
@click.option('-p', '--print_tree', default=False, is_flag=True,
              help='flag whether or not to print the output to the terminal, defaults to False')
def run_fasttree(output_file, print_tree, input_file):
    """
    Runs the FastTree algorithm to generate a phylogenetic tree from a list of nucleotide sequences. Requires a *.aln file as input.
    See Price MN, Dehal PS, Arkin AP. FastTree: computing large minimum evolution trees with profiles instead of a distance matrix. Mol Biol Evol. 2009;26(7):1641-1650. doi:10.1093/molbev/msp077
    for more information about the FastTree algorithm.
    Reproduction of the FastTree algorithm, by Eljo Dorrestijn, Jurrian de Boer and Frank te Nijenhuis
    """
    run(input_file, output_file, print_tree)


if __name__ == '__main__':
    run_fasttree()
