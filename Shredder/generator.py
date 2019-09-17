import random
from Shredder.utils import simulate_errors, reverse_complement


class PairedEndReads:
    def __init__(self, short_read_length, fasta_seq, coverage, fasta_name, insert_size, insert_size_std):
        self.short_read_length = short_read_length
        self.insert_size = insert_size
        self.insert_size_std = insert_size_std
        self.output_sample_name = fasta_name
        self.fasta_seq = fasta_seq
        self.coverage = coverage
        self.Q_str = ''
        self.error_rate = 0.05

    def generate(self):

        output_r1 = open('{}_R1.fastq'.format(self.output_sample_name), 'w+')
        output_r2 = open('{}_R2.fastq'.format(self.output_sample_name), 'w+')

        x = self.short_read_length
        self.Q_str = ''.join(['A'] * x)

        c = int(self.coverage)

        total_r1_errors = 0
        total_r2_errors = 0
        total_num_reads = 0

        L = len(self.fasta_seq)

        av_num_short_reads_needed = int(L / x * c)
        total_num_reads += av_num_short_reads_needed

        av_num_pairs_needed = int(av_num_short_reads_needed / 2)

        for index_pair in range(0, av_num_pairs_needed):
            I = int(round(random.gauss(self.insert_size, self.insert_size_std)))
            if L - ((x * 2) + I) > 0:
                start_pos = random.randint(0, L - ((x * 2) + I))
            else:
                start_pos = random.randint(0, L - (x * 2))

            read_1_start = start_pos
            read_1_stop = read_1_start + x

            read_2_start = read_1_stop + I
            read_2_stop = read_2_start + x

            read_1, num_errors_r1 = simulate_errors(self.error_rate, self.fasta_seq[read_1_start:read_1_stop])
            read_2, num_errors_r2 = simulate_errors(self.error_rate, self.fasta_seq[read_2_start:read_2_stop])

            total_r1_errors += num_errors_r1
            total_r2_errors += num_errors_r2

            c1, c2 = random.randint(1, 10000), random.randint(1, 10000)
            output_r1.write('@%s:23:B02CBACXX:8:2315:%d:%d 1:N:0:GATCAG\n' % (self.output_sample_name, c1, c2))
            output_r1.write(read_1 + '\n')
            output_r1.write('+source:%s; start:%d; stop:%d; insert_size:%d\n' % (self.output_sample_name, read_1_start,
                                                                                 read_1_stop, I))
            output_r1.write('%s\n' % self.Q_str)

            output_r2.write('@%s:23:B02CBACXX:8:2315:%d:%d 2:N:0:GATCAG\n' % (self.output_sample_name, c1, c2))
            output_r2.write(reverse_complement(read_2) + '\n')
            output_r2.write('+source:%s; start:%d; stop:%d; insert_size:%d\n' % (self.output_sample_name, read_2_start,
                                                                                 read_2_stop, I))
            output_r2.write('%s\n' % self.Q_str)

        total_num_errors = total_r1_errors + total_r2_errors

        # print("\tWrote {} pairs into {} and {} with {} errors.".format(av_num_pairs_needed,
        #                                                               self.output_sample_name + '_R1.fastq',
        #                                                               self.output_sample_name + '_R2.fastq',
        #                                                               total_num_errors))

        output_r1.close()
        output_r2.close()
        return '{}_R1.fastq'.format(self.output_sample_name), '{}_R2.fastq'.format(self.output_sample_name)
