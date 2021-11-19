
def parse_input():
    with open('data/test-small.aln') as f:
        lines = f.readlines()
        sequences = []
        for i in range(0, len(lines)):
            if i % 2 == 1:
                # Uneven lines are the sequences
                sequences.append(lines[i].strip())
        return sequences


def run(sequences):
    pass


if __name__ == '__main__':
    sequences = parse_input()
    run(sequences)
