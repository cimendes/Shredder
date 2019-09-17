import random


BASES = ['A', 'T', 'C', 'G', 'N']
COMPLEMENT = str.maketrans("ACTGN", "TGACN")


def reverse_complement(seq):
    seq = seq.upper()  # make sure it's uppercase
    return seq.translate(COMPLEMENT)[::-1]


def simulate_errors(error_rate, sequence):
    sequence_with_errors = ''
    num_errors = 0

    if error_rate > 0:
        threshold = 1000 * error_rate
        for i in range(0, len(sequence)):
            if random.randint(0, 1000) < threshold:
                sequence_with_errors += random.choice(BASES)
                num_errors += 1
            else:
                sequence_with_errors += sequence[i]
    else:
        sequence_with_errors = sequence

    return sequence_with_errors, num_errors

