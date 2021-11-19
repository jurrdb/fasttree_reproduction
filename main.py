
def parse_input():
    with open('data/test-small.aln') as f:
        lines = f.readlines()
        sequences = {}
        for i in range(0, len(lines), 2):
            sequences[lines[i].strip()] =  lines[i + 1].strip()
        return sequences


def run(sequences):
    pass


if __name__ == '__main__':
    sequences = parse_input()
    run(sequences)
