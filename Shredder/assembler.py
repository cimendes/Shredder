import os
import subprocess
import shutil
from subprocess import PIPE


class Assembly:
    def __init__(self, forward, reverse, filename, read_length, output_dir, cpus, current_dir, memory):
        self.forward_reads = forward
        self.reverse_reads = reverse
        self.ouput_name = filename
        self.output_dir = output_dir
        self.threads = cpus
        self.current_dir = current_dir
        self.mem_limit = memory
        if read_length >= 175:
            self.kmers = [55, 77, 99, 113, 127]
        else:
            self.kmers = [21, 33, 55, 67, 77]

    def do_spades_assembly(self):

        cli = [
            "spades.py",
            "--careful",
            "--only-assembler",
            "--threads",
            str(self.threads),
            "-o",
            self.output_dir,
            "-m",
            str(self.mem_limit)
        ]

        cli += ["-k {}".format(",".join([str(x) for x in self.kmers]))]

        cli += [
            "-1",
            self.forward_reads,
            "-2",
            self.reverse_reads
        ]

        print("\t Running SPAdes subprocess with command: \n\t{}".format(cli))

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
            print(stderr)
            print(stdout)
            return None
        else:
            assembly_file = os.path.join(self.current_dir, "{}_spades.fasta".format(self.ouput_name))
            shutil.move(os.path.join(self.output_dir, "contigs.fasta"), assembly_file)

            return assembly_file

