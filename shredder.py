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
        cli = ["spades.py", "--version"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        print("SPades version {} found.".format(stdout.strip().split()[-1][1:].decode("utf8")))

    except OSError:
        print("SPAdes not found. Exiting...")
        sys.exit(1)

    try:
        cli = ["shuffle.sh", "-h"]
        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, _ = p.communicate()

        print("shuffle.sh found.")

    except OSError:
        print("shuffle.sh not found. Exiting...")
        sys.exit(1)


def shuffle_reads(reads_1, reads_2, outdir):

    cli = ["shuffle.sh",
           "in={}".format(reads_1),
           "in2={}".format(reads_2),
           "out={}".format(os.path.join(outdir, "shuffled_" + os.path.basename(reads_1))),
           "out2={}".format(os.path.join(outdir, "shuffled_" + os.path.basename(reads_2)))]
    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        print("\terror while running SPAdes")
        print(stderr)
        return None, None
    else:
        os.remove(reads_1)
        os.remove(reads_2)
        return os.path.join(outdir, "shuffled_" + os.path.basename(reads_1)), \
               os.path.join(outdir, "shuffled_" + os.path.basename(reads_2))


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
                        nargs='+', required=True, type=int)
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

    print("Saving shredded sample reads in: {}/shredded_R1.fastq, {}/shredded_R2.fastq".
          format(args.outdir, args.outdir))
    print("Saving shredded assembly file in: {}/shredded_assembly.fasta".format(args.outdir))
    if os.path.isfile(os.path.join(args.outdir, 'shredded_R1.fastq')) or \
            os.path.isfile(os.path.join(args.outdir, 'shredded_R2.fastq')) or \
            os.path.isfile(os.path.join(args.outdir, 'shredded_assembly.fasta')):
        print("Output files already exist... Overwriting")
        os.remove(os.path.join(args.outdir, 'shredded_R1.fastq'))
        os.remove(os.path.join(args.outdir, 'shredded_R2.fastq'))
        os.remove(os.path.join(args.outdir, 'shredded_assembly.fasta'))

    read_file_1 = open(os.path.join(args.outdir, 'shredded_R1.fastq'), "w")
    read_file_2 = open(os.path.join(args.outdir, 'shredded_R2.fastq'), "w")
    assembly_file = open(os.path.join(args.outdir, "shredded_assembly.fasta"), "w")

    print("Saving bins in: {}/Bins/".format(args.outdir))
    try:
        os.mkdir(os.path.join(args.outdir, "Bins"))
    except OSError:
        print("Bin folder already exists... Overwriting")
        shutil.rmtree(os.path.join(args.outdir, "Bins"))
        os.mkdir(os.path.join(args.outdir, "Bins"))

    print("Generating reads..")
    to_multiprocess = []
    for file in args.reference:
        to_multiprocess.append({"reference": file,
                                "coverage": coverage[args.reference.index(file)],
                                "read-length": args.read_length,
                                "insert-size": args.insert_size,
                                "insert-size-std": args.insert_size_std,
                                "outdir": args.outdir
                                })

    reads = []
    p = Pool(processes=args.threads)
    r = p.map_async(gen_reads, to_multiprocess, callback=reads.extend)
    r.wait()

    assemblies = []

    print("Generating assemblies...")
    for read_pair in reads:
        assemblies.append(run_assembly(forward=read_pair[0], reverse=read_pair[1],
                                       filename=read_pair[0].split('_R1.fastq')[0],
                                       read_length=args.read_length, cpus=args.threads,
                                       outdir=args.outdir, memory=args.memory))
        with open(read_pair[0], 'r') as infile1:
            shutil.copyfileobj(infile1, read_file_1)
        with open(read_pair[1], 'r') as infile2:
            shutil.copyfileobj(infile2, read_file_2)

    for file in assemblies:
        if file is not None:
            with open(file, "r") as infile_assembly:
                for line in infile_assembly:
                    if line.startswith(">"):
                        line = line.strip() + "_{}\n".format(os.path.basename(file).split(".")[0])
                        assembly_file.write(line)
                    else:
                        assembly_file.write(line)

            shutil.move(file, os.path.join(args.outdir, "Bins",
                                           "Bin_" + str(assemblies.index(file)) + "_" + os.path.basename(file)))
    read_file_1.close()
    read_file_2.close()

    print("Shuffling reads...")
    shuffle_reads(read_file_1, read_file_2, args.outdir)

    print("Finished!")


if __name__ == '__main__':
    main()
