#!/usr/bin/env python3
import sys
import os
import argparse
from multiprocessing import Pool
from itertools import groupby
from Shredder.generator import PairedEndReads
from Shredder.assembler import Assembly

# TODO - SPAdes not working
# TODO - Add spades thread option, with verification to ensure that threads * spades_threads < total_cpu

def run_assembly(forward, reverse, filename, read_length):

    assembly = Assembly(filename=filename.split('.')[0], forward=forward, reverse=reverse, read_length=read_length)
    assembly_file = assembly.do_SPAdes_assembly()

    return assembly_file


def gen_reads(config):
    print("Generating reads for {} with {}x coverage:".format(config["reference"], config["coverage"]))

    fh = open(config["reference"], "r")
    entry = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    header_str = ""
    seq = ""
    for header in entry:
        header_str = header.__next__()[1:].strip().split(' ')[0]
        seq += "".join(s.strip() for s in entry.__next__())

    gen = PairedEndReads(coverage=config["coverage"], fasta_name=header_str,
                         fasta_seq=seq, short_read_length=config["read-length"],
                         insert_size=config["insert-size"], insert_size_std=config["insert-size-std"])
    f, r = gen.generate()
    return f, r


def run_shredder(config):
    forward, reverse = gen_reads(config)
    assembly = run_assembly(forward, reverse, config["reference"], config["read-length"])

    return forward, reverse, assembly


def main():

    parser = argparse.ArgumentParser(prog='shredder.py', description="Genome goes in, bits and pieces come out.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v', '--version', help="Print version and exit.", action='version', version='alpha-0.1')

    parser.add_argument('-r', '--reference', help='Path of a file or a folder of files', nargs='+')
    parser.add_argument('-l', '--read_length', help='Target read length', type=int, required=True)
    parser.add_argument('-c', '--coverage', help='Coverage, or list of coverages for each reference file (alphabetical)',
                        nargs='+', required=True, type=int)
    parser.add_argument('-t', '--threads', help="Number of threads to use", required=False, default=2)

    # Optional arguments
    parser.add_argument('--insert_size', help="Paired-end read insert size.", type=int, default=140, required=False)
    parser.add_argument('-insert_size_std', help="Standard deviation for paired-end read insert size.", type=int,
                        default=60, required=False)

    args = parser.parse_args()

    if len(args.coverage) == 1:
        coverage = [args.coverage[0]] * len(args.reference)
    else:
        if len(args.coverage) != len(args.reference):
            print("number of coverage values and reference sequences don't match.")
            sys.exit(1)
        else:
            coverage = args.coverage

    print("Saving shredded sample reads in: shredded_R1.fastq, shredded_R2.fastq")
    if os.path.isfile('shredded_R1.fastq') or os.path.isfile('shredded_R2.fastq'):
        print("Output files already exist... Overwriting")
        os.remove('shredded_1.fastq')
        os.remove('shredded_2.fastq')

    to_multiprocess = []
    for file in args.reference:
        to_multiprocess.append({"reference": file,
                                "coverage": coverage[args.reference.index(file)],
                                "read-length": args.read_length,
                                "insert-size": args.insert_size,
                                "insert-size-std": args.insert_size_std})

    reads = []
    p = Pool(processes=args.threads)
    r = p.map_async(run_shredder, to_multiprocess, callback=reads.extend)
    r.wait()

    print(reads)


if __name__ == '__main__':
    main()

