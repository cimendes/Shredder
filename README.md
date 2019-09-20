# Shredder 
:octopus: A (tiny) script to randomly shred complete genomes into metagenomic-like data :octopus:


To properly benchmark metagenomic software a set of golden standards are necessary for fair comparison. Shredder 
is a very simple python tool that aims to recreate the best possible outcome out of a metagenomic analysis process: 
reads, assembly and bins. 


When provided with complete reference sequences, shredder will generate artificial sequencing reads, in Illumina 
paired-end form, with a given target coverage. The reads generates will cover the reference genome uniformly. 
Each read dataset, prior to being (TODO: shuffled and) concatenated in a single metagenomic sequencing dataset, is 
assembled with SPAdes. Each resulting assembly is saved as a Bin file, and all assemblies are concatenated (TODO: 
ordered by length) into a single metagenomic assembled contigs file. 


## Dependencies & Installation
Right now the only dependency is [SPAdes](http://cab.spbu.ru/software/spades/). I hope to remove this in the future. 
The easiest way to install this is through conda with the command `conda install -c bioconda spades`. 


You can download Shredder by cloning the GitHub repository and adding it to your PATH. 
`git clone https://github.com/cimendes/Shredder.git`
`export "PATH=$PWD/Shredder:$PATH"`


## Usage
To use shredder you need to provide at least one reference file (`-r`), a read length (`-l`) and a target 
coverage (`-c`). Several references can be passed, using the wildcard `*`. A coverage value must be provided
for each of the sequences in `-r` or, alternatively, a single value that will be used for all files. 


    shredder.py --help
    usage: shredder.py [-h] [-v] [-r REFERENCE [REFERENCE ...]] -l READ_LENGTH -c
                       COVERAGE [COVERAGE ...] [-t THREADS] [-m MEMORY]
                       [--insert_size INSERT_SIZE]
                       [-insert_size_std INSERT_SIZE_STD] [-o OUTDIR]
    
    Genome goes in, bits and pieces come out.
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         Print version and exit.
      -r REFERENCE [REFERENCE ...], --reference REFERENCE [REFERENCE ...]
                            Path of a file or a folder of files (default: None)
      -l READ_LENGTH, --read_length READ_LENGTH
                            Target read length (default: None)
      -c COVERAGE [COVERAGE ...], --coverage COVERAGE [COVERAGE ...]
                            Coverage, or list of coverages for each reference file
                            (alphabetical) (default: None)
      -t THREADS, --threads THREADS
                            Number of threads to use (default: 6)
      -m MEMORY, --memory MEMORY
                            Memory limit in GB (default: 8)
      --insert_size INSERT_SIZE
                            Paired-end read insert size. (default: 140)
      -insert_size_std INSERT_SIZE_STD
                            Standard deviation for paired-end read insert size.
                            (default: 60)
      -o OUTDIR, --outdir OUTDIR
                            Output directory (default: current directory)


## TODO
* Shuffle the resulting *in silico* metagenomic reads (possibly with [seqkit](https://github.com/shenwei356/seqkit)? 
Full python alternative preferred).
* Order resulting MAGs by sequence length. 