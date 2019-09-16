#!/usr/bin/env python3
import sys
import os
import glob
from itertools import groupby
from Shredder.generator import PairedEndReads
from Shredder.assembler import Assembly


def main():

    try:
        # TODO add a read argparse...
        # TODO check dependencies
        reference_dir = sys.argv[1]
        read_length = int(sys.argv[2])
        coverage = [int(i) for i in sys.argv[3:]]
        if len(coverage) == 1:
            coverage = [coverage[0]] * len(os.listdir(reference_dir))
        else:
            if len(coverage) != len(os.listdir(reference_dir)):
                print("number of coverage values and reference sequences don't match.")
                raise IndexError

    except IndexError:
        print("Shred your genomes! Usage:\n")
        print("\t shredder.py <reference> <read_length> <target coverage>")
        sys.exit(0)

    print("Saving shredded sample reads in: shredded_R1.fastq, shredded_R2.fastq")
    if os.path.isfile('shredded_R1.fastq') or os.path.isfile('shredded_R2.fastq'):
        print("Output files already exist... Overwriting")
        os.remove('shredded_1.fastq')
        os.remove('shredded_2.fastq')

    count = 0

    for file in os.listdir(reference_dir):
        target_coverage = coverage[os.listdir(reference_dir).index(file)]
        print("Generating reads for {} with {}x coverage:".format(file, target_coverage))
        fh = open(os.path.join(reference_dir, file), "r")
        entry = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        header_str = ""
        seq = ""
        for header in entry:
            header_str = header.__next__()[1:].strip().split(' ')[0]
            seq += "".join(s.strip() for s in entry.__next__())

        gen = PairedEndReads(coverage=target_coverage, fasta_name=header_str,
                             fasta_seq=seq, short_read_length=read_length)
        f, r = gen.generate()

        assembly = Assembly(filename=file.split('.')[0], forward=f, reverse=r, read_length=read_length)
        assembly_file = assembly.do_SPAdes_assembly()

        # assemble each read file


if __name__ == '__main__':
    main()

