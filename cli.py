import click
import newick

from main import run, to_newick, parse_input


@click.command()
@click.argument('file', type=click.File())
def run_fasttree(filename):
    sequences, sequence_length = parse_input(filename)
    tree, final_tree = run(sequences, sequence_length)
    newick_string = to_newick(tree, named_parent_nodes=False)
    newick_string_final = to_newick(final_tree, named_parent_nodes=False)
    n_tree = newick.loads(newick_string)[0]
    f_tree = newick.loads(newick_string_final)[0]
    print(n_tree.ascii_art())
    print(f_tree.ascii_art())
