#!/usr/bin/env python3
import sys
import shutil
import os
import argparse
import tempfile
import subprocess
from multiprocessing import Pool, cpu_count
from subprocess import PIPE
from itertools import groupby
from Shredder.generator import PairedEndReads
from Shredder.assembler import Assembly


def check_dependencies():

    try:
        cli = ["shuffle.sh", "-h"] #BBTools?
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        print("shuffle.sh found.")

    except OSError:
        print("shuffle.sh not found. Exiting...")
        sys.exit(1)


def shuffle_reads(reads_1, reads_2):

    cli = ["shuffle.sh",
           "in={}".format(reads_1),
           "in2={}".format(reads_2),
           "out={}".format(os.path.join(os.path.dirname(reads_1), "shuffled_" + os.path.basename(reads_1))),
           "out2={}".format(os.path.join(os.path.dirname(reads_2), "shuffled_" + os.path.basename(reads_2)))]

    print("\t Running shuffle subprocess with command: \n\t{}".format(cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        print("\terror while running shuffle.sh")
        print(stderr)
        return None, None
    else:
        os.remove(reads_1)
        os.remove(reads_2)
        return os.path.join(os.path.dirname(reads_1), "shuffled_" + os.path.basename(reads_1)), \
               os.path.join(os.path.dirname(reads_2), "shuffled_" + os.path.basename(reads_2))


def run_assembly(forward, reverse, filename, read_length, cpus, outdir, memory):

    dirpath = tempfile.mkdtemp(dir=os.path.join(outdir, os.getcwd()))

    assembly = Assembly(filename=filename.split('.')[0], forward=forward, reverse=reverse,
                        read_length=read_length, output_dir=dirpath, cpus=cpus, current_dir=outdir,
                        memory=memory)
    assembly_file = assembly.do_spades_assembly()

    shutil.rmtree(dirpath)

    return assembly_file


def gen_reads(config):
    print("\tGenerating reads for {} with {}x coverage.".format(config["reference"], config["coverage"]))

    fh = open(config["reference"], "r")
    entry = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    header_str = ""
    seq = ""
    for header in entry:
        header_str = header.__next__()[1:].strip().split(' ')[0]
        seq += "".join(s.strip() for s in entry.__next__())

    gen = PairedEndReads(coverage=config["coverage"], fasta_name=header_str,
                         fasta_seq=seq, short_read_length=config["read-length"],
                         insert_size=config["insert-size"], insert_size_std=config["insert-size-std"],
                         outdir =config["outdir"])
    f, r = gen.generate()
    return f, r


def main():

    parser = argparse.ArgumentParser(prog='shredder.py', description="Genome goes in, bits and pieces come out.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v', '--version', help="Print version and exit.", action='version', version='alpha-0.1')

    parser.add_argument('-r', '--reference', help='Path of a file or a folder of files', nargs='+')
    parser.add_argument('-l', '--read_length', help='Target read length', type=int, required=True)
    parser.add_argument('-c', '--coverage', help='Coverage, or list of coverages for each reference file (alphabetical)',
                        nargs='+', required=True)
    parser.add_argument('-t', '--threads', help="Number of threads to use", required=False, default=cpu_count()-2,
                        type=int)
    parser.add_argument('-m', '--memory', help="Memory limit in GB", required=False, default=8)

    # Optional arguments
    parser.add_argument('--insert_size', help="Paired-end read insert size.", type=int, default=140, required=False)
    parser.add_argument('-insert_size_std', help="Standard deviation for paired-end read insert size.", type=int,
                        default=60, required=False)
    parser.add_argument('-o', '--outdir', help="Output directory", default=os.getcwd(), required=False)

    args = parser.parse_args()

    check_dependencies()

    try:
        args.coverage = args.coverage[0].split(',')
    except AttributeError:
        args.coverage = args.coverage[0]

    print(args.coverage)
    print(args.reference)

    if len(args.coverage) == 1:
        coverage = [args.coverage[0]] * len(args.reference)
    else:
        if len(args.coverage) != len(args.reference):
            print("Number of coverage values and reference sequences don't match.")
            sys.exit(1)
        else:
            coverage = args.coverage

    if args.threads > cpu_count():
        print("Number of cpus requested ({}) larger than the available cpus ({}).".format(args.threads, cpu_count()))
        sys.exit(1)

    reads_r1_filename = os.path.join(args.outdir, 'shredded_R1.fastq')
    reads_r2_filename = os.path.join(args.outdir, 'shredded_R2.fastq')

    print("Saving shredded sample reads in: {}, {}".
          format(reads_r1_filename, reads_r2_filename))

    if os.path.isfile(reads_r1_filename) or os.path.isfile(reads_r2_filename):
        print("Output files already exist... Overwriting")
        os.remove(reads_r1_filename)
        os.remove(reads_r2_filename)

    read_file_1 = open(reads_r1_filename, "w")
    read_file_2 = open(reads_r2_filename, "w")

    print("Generating reads..")
    to_multiprocess = []
    for ref_file in args.reference:
        to_multiprocess.append({"reference": ref_file,
                                "coverage": coverage[args.reference.index(ref_file)],
                                "read-length": args.read_length,
                                "insert-size": args.insert_size,
                                "insert-size-std": args.insert_size_std,
                                "outdir": args.outdir
                                })

    reads = []
    p = Pool(processes=args.threads)
    r = p.map_async(gen_reads, to_multiprocess, callback=reads.extend)
    r.wait()

    read_file_1.close()
    read_file_2.close()

    print("Shuffling reads...")
    shuffle_reads(reads_r1_filename, reads_r1_filename)
    print("Finished!")


if __name__ == '__main__':
    main()
