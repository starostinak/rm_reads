#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <algorithm>
#include <locale>
#include <getopt.h>
#include <stdlib.h>
#include <unordered_map>

#include "search.h"
#include "stats.h"
#include "seq.h"
#include "rm_reads.h"

#define LENGTH_CUTOFF 50
#define DUST_K 4
#define POLYG 13

struct FilterCmd {
    FilterCmd()
        : reads1_fp(nullptr), reads2_fp(nullptr),
          ok1_fp(nullptr), ok2_fp(nullptr), bad1_fp(nullptr), bad2_fp(nullptr),
          se1_fp(nullptr), se2_fp(nullptr),
          stats1(), stats2(), root(nullptr),
          length(LENGTH_CUTOFF), dust_k(DUST_K), dust_cutoff(0),
          errors(0) {}

private:

    ReadType check_read(std::string const & read, std::vector <std::pair<std::string, Node::Type> > const & patterns)
    {
        if (length && read.size() < length) {
            return ReadType::length;
        }
        if (dust_cutoff && get_dust_score(read, dust_k) > dust_cutoff) {
            return ReadType::dust;
        }

        if (errors) {
            return (ReadType)search_inexact(read, root, patterns, errors);
        } else {
            return (ReadType)search_any(read, root);
        }
    }

    void filter_single_reads(std::vector <std::pair<std::string, Node::Type> > const & patterns)
    {
        Seq read;

        std::ifstream & reads_f = *reads1_fp;
        std::ofstream & ok_f = *ok1_fp;
        std::ofstream & bad_f = *bad1_fp;

        while (read.read_seq(reads_f)) {
            ReadType type = check_read(read.get_seq(), patterns);
            stats1.update(type);
            if (type == ReadType::ok) {
                read.write_seq(ok_f);
            } else {
                read.update_id(type);
                read.write_seq(bad_f);
            }
        }
    }

    void filter_paired_reads(std::vector <std::pair<std::string, Node::Type> > const & patterns)
    {
        Seq read1;
        Seq read2;

        std::ifstream & reads1_f = *reads1_fp;
        std::ifstream & reads2_f = *reads2_fp;
        std::ofstream & ok1_f = *ok1_fp;
        std::ofstream & ok2_f = *ok2_fp;
        std::ofstream & se1_f = *se1_fp;
        std::ofstream & se2_f = *se2_fp;
        std::ofstream & bad1_f = *bad1_fp;
        std::ofstream & bad2_f = *bad2_fp;

        while (read1.read_seq(reads1_f) && read2.read_seq(reads2_f)) {
            ReadType type1 = check_read(read1.get_seq(), patterns);
            ReadType type2 = check_read(read2.get_seq(), patterns);
            if (type1 == ReadType::ok && type2 == ReadType::ok) {
                read1.write_seq(ok1_f);
                read2.write_seq(ok2_f);
                stats1.update(type1, true);
                stats2.update(type2, true);
            } else {
                stats1.update(type1, false);
                stats2.update(type2, false);
                if (type1 == ReadType::ok) {
                    read1.write_seq(se1_f);
                    read2.update_id(type2);
                    read2.write_seq(bad2_f);
                } else if (type2 == ReadType::ok) {
                    read1.update_id(type1);
                    read1.write_seq(bad1_f);
                    read2.write_seq(se2_f);
                } else {
                    read1.update_id(type1);
                    read2.update_id(type2);
                    read1.write_seq(bad1_f);
                    read2.write_seq(bad2_f);
                }
            }
        }
    }

public:

    void filter_reads(std::vector <std::pair<std::string, Node::Type> > const & patterns) {
        if (reads2_fp == nullptr) {
            filter_single_reads(patterns);
        } else {
            filter_paired_reads(patterns);
        }
    }

    std::ifstream * reads1_fp;
    std::ifstream * reads2_fp;
    std::ofstream * ok1_fp;
    std::ofstream * ok2_fp;
    std::ofstream * bad1_fp;
    std::ofstream * bad2_fp;
    std::ofstream * se1_fp;
    std::ofstream * se2_fp;
    Stats stats1;
    Stats stats2;
    Node * root;
    size_t length;
    int dust_k;
    int dust_cutoff;
    int errors;
};

void build_patterns(std::ifstream & kmers_f, std::vector <std::pair <std::string, Node::Type> > & patterns, int polyG, bool filterN)
{
    std::string tmp;
    while (!kmers_f.eof()) {
        std::getline(kmers_f, tmp);
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        if (!tmp.empty()) {
            size_t tab = tmp.find('\t');
            if (tab == std::string::npos) {
                patterns.push_back(std::make_pair(tmp, Node::Type::adapter));
            } else {
                patterns.push_back(std::make_pair(tmp.substr(0, tab), Node::Type::adapter));
            }
        }
    }
    kmers_f.close();

    if (filterN) {
        patterns.push_back(std::make_pair("N", Node::Type::n));
    }
    if (polyG) {
        patterns.push_back(std::make_pair(std::string(polyG, 'G'), Node::Type::polyG));
        patterns.push_back(std::make_pair(std::string(polyG, 'C'), Node::Type::polyC));
    }
}

double get_dust_score(std::string const & read, int k)
{
    std::unordered_map <int, int> counts;
    static std::unordered_map <char, int> hashes = {{'N', 1},
                                          {'A', 2},
                                          {'C', 3},
                                          {'G', 4},
                                          {'T', 5}};
    unsigned int hash = 0;
    unsigned int max_pow = pow(10, k - 1);
    for (auto it = read.begin(); it != read.end(); ++it) {
        char c = std::toupper(*it);
        hash = hash * 10 + hashes[c];
        if (it - read.begin() >= k - 1) {
            ++counts[hash];
            hash = hash - (hash / max_pow) * max_pow;
        }
    }
    double score = 0;
    double total = 0;
    for (auto it = counts.begin(); it != counts.end(); ++it) {
        score += it->second * (it->second - 1) / 2;
        total += score;
    }
//    std::cout << (total / (read.size() - k + 1)) << std::endl;
    return (total / (read.size() - k + 1));
}

std::string basename(std::string const & path)
{
    std::string res(path);
    size_t pos = res.find_last_of('/');
    if (pos != std::string::npos) {
        res.erase(0, pos);
    }
    pos = res.find('.');
    if (pos != std::string::npos) {
        res.erase(pos, path.size());
    }
    return res;
}

void print_help() 
{
    std::cerr << "Usage:\n./rm_reads <-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq> --adapters adapters.dat [-o output_dir --polyG POLYG --length LENGTH_CUTOFF --dust_cutoff cutoff --dust_k k -errors 0 -filterN]\n"
        << "\nOptions:\n"
        << "\t-i\t\tinput file \n"
        << "\t-1\t\tfirst input file for paired reads\n"
        << "\t-2\t\tsecond input file for paired reads\n"
        << "\t-o\t\toutput directory (current directory by default)\n"
        << "\t--polyG, -p\tlength of polyG/polyC tails (13 by default)\n"
        << "\t--length, -l\tminimum length cutoff (50 by default)\n"
        << "\t--adapters, -a\tfile with adapter kmers\n"
        << "\t--dust_k, -k\twindow size for dust filter (not used by default)\n"
        << "\t--dust_cutoff, -c\tcutoff by dust score (not used by default)\n"
        << "\t--errors, -e\tmaximum error count in match, possible values - 0, 1, 2 (by default 0)\n"
        << "\t--filterN, -N\tallow filter by N's in reads" << std::endl;
}

int main(int argc, char ** argv)
{
    Node root('0');
    std::vector <std::pair<std::string, Node::Type> > patterns;

    std::string kmers, reads, out_dir;
    std::string reads1, reads2;
    char rez = 0;
    int polyG = POLYG;
    bool filterN = false;
    FilterCmd cmd;

    const struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"length",required_argument,NULL,'l'},
        {"polyG",required_argument,NULL,'p'},
        {"adapters",required_argument,NULL,'a'},
        {"dust_k",required_argument,NULL,'k'},
        {"dust_cutoff",required_argument,NULL,'c'},
        {"errors", required_argument, NULL, 'e'},
        {"filterN", no_argument, NULL, 'N'},
        {NULL,0,NULL,0}
    };

    while ((rez = getopt_long(argc, argv, "hN1:2:l:p:a:i:o:e:", long_options, NULL)) != -1) {
        switch (rez) {
        case 'l':
            cmd.length = std::atoi(optarg);
            break;
        case 'p':
            polyG = std::atoi(optarg);
            // polyG = boost::lexical_cast<int>(optarg);
            break;
        case 'a':
            kmers = optarg;
            break;
        case 'i':
            reads = optarg;
            break;
        case '1':
            reads1 = optarg;
            break;
        case '2':
            reads2 = optarg;
            break;
        case 'o':
            out_dir = optarg;
            break;
        case 'c':
            cmd.dust_cutoff = std::atoi(optarg);
            break;
        case 'k':
            cmd.dust_k = std::atoi(optarg);
            break;
        case 'e':
            cmd.errors = std::atoi(optarg);
            break;
        case 'N':
            filterN = true;
            break;
        case '?':
        case 'h':
            print_help();
            return -1;
        }
    }

    if (cmd.errors < 0 || cmd.errors > 2) {
        std::cerr << "Possible errors count are 0, 1, 2" << std::endl;
        return -1;
    }

    if (out_dir.empty()) {
        out_dir = ".";
    }

    if (kmers.empty() || (reads.empty() &&
            (reads1.empty() || reads2.empty()))) {
        std::cerr << "Please, specify reads and kmers files" << std::endl;
        print_help();
        return -1;
    }

    std::ifstream kmers_f (kmers.c_str());
    if (!kmers_f.good()) {
        std::cerr << "Cannot open kmers file" << std::endl;
        print_help();
        return -1;
    }

    init_type_names(length, polyG, cmd.dust_k, cmd.dust_cutoff);

    build_patterns(kmers_f, patterns, polyG, filterN);

    if (patterns.empty()) {
        std::cerr << "patterns are empty" << std::endl;
        return -1;
    }

    build_trie(root, patterns, cmd.errors);
	add_failures(root);

    cmd.root = &root;

    if (!reads.empty()) {
        std::string reads_base = basename(reads);
        std::ifstream reads_f (reads.c_str());
        std::ofstream ok_f((out_dir + "/" + reads_base + ".ok.fastq").c_str(), std::ofstream::out);
        std::ofstream bad_f((out_dir + "/" + reads_base + ".filtered.fastq").c_str(), std::ofstream::out);

        if (!reads_f.good()) {
            std::cerr << "Cannot open reads file, please, make sure that it exists" << std::endl;
            print_help();
            return -1;
        }

        if (!ok_f.good() || !bad_f.good()) {
            std::cerr << "Cannot open output files, please, make sure that output directory exists and you can write there" << std::endl;
            print_help();
            return -1;
        }

        cmd.stats1 = Stats(reads);
        cmd.reads1_fp = &reads_f;
        cmd.ok1_fp = &ok_f;
        cmd.bad1_fp = &bad_f;

        cmd.filter_reads(patterns);

        std::cout << cmd.stats1;

        ok_f.close();
        bad_f.close();
        reads_f.close();
    } else {
        std::string reads1_base = basename(reads1);
        std::string reads2_base = basename(reads2);
        std::ifstream reads1_f(reads1.c_str());
        std::ifstream reads2_f(reads2.c_str());
        std::ofstream ok1_f((out_dir + "/" + reads1_base + ".ok.fastq").c_str(),
                            std::ofstream::out);
        std::ofstream ok2_f((out_dir + "/" + reads2_base + ".ok.fastq").c_str(),
                            std::ofstream::out);
        std::ofstream se1_f((out_dir + "/" + reads1_base + ".se.fastq").c_str(),
                            std::ofstream::out);
        std::ofstream se2_f((out_dir + "/" + reads2_base + ".se.fastq").c_str(),
                            std::ofstream::out);
        std::ofstream bad1_f((out_dir + "/" + reads1_base + ".filtered.fastq").c_str(),
                            std::ofstream::out);
        std::ofstream bad2_f((out_dir + "/" + reads2_base + ".filtered.fastq").c_str(),
                            std::ofstream::out);

        if (!reads1_f.good() || !reads2_f.good()) {
            std::cerr << "Cannot open reads file, please, make sure that it exists" << std::endl;
            print_help();
            return -1;
        }

        if (!ok1_f.good() || !ok2_f.good() || !bad1_f.good() || !bad2_f.good() ||
                !se1_f.good() || !se2_f.good()) {
            std::cerr << "Cannot open output files, please, make sure that output directory exists and you can write there" << std::endl;
            print_help();
            return -1;
        }

        cmd.stats1 = Stats(reads1);
        cmd.stats2 = Stats(reads2);
        cmd.reads1_fp = &reads1_f;
        cmd.reads2_fp = &reads2_f;
        cmd.ok1_fp = &ok1_f;
        cmd.ok2_fp = &ok2_f;
        cmd.se1_fp = &se1_f;
        cmd.se2_fp = &se2_f;
        cmd.bad1_fp = &bad1_f;
        cmd.bad2_fp = &bad2_f;

        cmd.filter_reads(patterns);

        std::cout << cmd.stats1;
        std::cout << cmd.stats2;

        ok1_f.close();
        ok2_f.close();
        bad1_f.close();
        bad2_f.close();
        reads1_f.close();
        reads2_f.close();
    }
    return 0;
}
