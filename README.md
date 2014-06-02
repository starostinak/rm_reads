Description
----------------------

A tool for filtering reads by given kmers set. Initially created to filter reads with primers and adapters from a large datasets.

Requirements
----------------------

- make
- gcc 4.7 at least

Usage
----------------------

./rm_reads [-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq] -o output_dir --polyG 13 --length 50 --adapters adapters.dat --dust_cutoff cutoff --dust_k k -e 1 -N

    -i              input file
    -1              first input file for paired reads
    -2              second input file for paired reads
    -o              output directory
    --polyG, -p     length of polyG/polyC tails
    --adapters, -a  file with adapter kmers
    --dust_k, -k    window size for dust filter
    --dust_cutoff, -c   cutoff by dust score
    --errors, -e    maximum error count in match, possible values - 1, 2
    --filterN, -N   allow filter by N's in reads

Input files
--------------------

Tools takes files with reads in fastq format as input. You can also use paired end reads. In case if you are using paired end reads, please, make sure that all reads from first file have correct pairs in second file.

Output files
--------------------

Tool creates following files in output directory:
input_prefix.ok.fastq       file with correct reads
input_prefix.fitered.fastq  file with reads, containing adapter kmers, N's, polyG/polyC tails or filtered by dust filter. Reason why read was filtered is given in the read id.
input_prefix.se.fastq       for paired reads only. File with correct reads which have incorrect pair.

Project page
--------------------

Last version of project can be found at https://github.com/allivi/rm_reads
