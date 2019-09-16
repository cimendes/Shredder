import os
import subprocess
from subprocess import PIPE


class Assembly:
    def __init__(self, forward, reverse, filename, read_length):
        self.forward_reads = forward
        self.reverse_reads = reverse
        self.ouput = filename
        if read_length >= 175:
            self.kmers = [55, 77, 99, 113, 127]
        else:
            self.kmers = [21, 33, 55, 67, 77]

    def do_SPAdes_assembly(self):

        # TODO add check for spades!
        # TODO add multithreading
        cli = [
            "spades.py",
            "--careful",
            "--only-assembler",
            "--threads",
            "4"
            "-o",
            "."
        ]

        cli += ["-k {}".format(",".join([str(x) for x in self.kmers]))]

        cli += [
            "-1",
            self.forward_reads,
            "-2",
            self.reverse_reads
        ]

        print("\t Running SPAdes subprocess for with command: \n\t{}".format(cli))

        p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
        # string
        try:
            stderr = stderr.decode("utf8")
            stdout = stdout.decode("utf8")
        except (UnicodeDecodeError, AttributeError):
            stderr = str(stderr)
            stdout = str(stdout)

        if p.returncode != 0:
            print("\terror while running SPAdes")

        assembly_file = "{}_spades.fasta".format(self.ouput)
        os.rename("contigs.fasta", assembly_file)

        return assembly_file

