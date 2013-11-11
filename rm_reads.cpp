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
#include <boost/lexical_cast.hpp>

#include "search.h"
#include "stats.h"
#include "seq.h"

void build_patterns(std::ifstream & kmers_f, int polyG, std::vector <std::pair <std::string, Node::Type> > & patterns)
{
    std::string tmp;
    while (!kmers_f.eof()) {
        std::getline(kmers_f, tmp);
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        if (!tmp.empty()) {
            size_t tab = tmp.find('\t');
            if (tab == std::string::npos) {
                patterns.push_back(std::make_pair(tmp, Node::Type::adapter));
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
                patterns.push_back(std::make_pair(tmp, Node::Type::adapter));
            } else {
                patterns.push_back(std::make_pair(tmp.substr(0, tab), Node::Type::adapter));
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
                patterns.push_back(std::make_pair(tmp.substr(0, tab), Node::Type::adapter));
            }
        }
    }
    kmers_f.close();

    patterns.push_back(std::make_pair("N", Node::Type::n));
    patterns.push_back(std::make_pair("n", Node::Type::n));
    if (polyG) {
        patterns.push_back(std::make_pair(std::string(polyG, 'G'), Node::Type::polyG));
        patterns.push_back(std::make_pair(std::string(polyG, 'C'), Node::Type::polyC));
        patterns.push_back(std::make_pair(std::string(polyG, 'g'), Node::Type::polyG));
        patterns.push_back(std::make_pair(std::string(polyG, 'c'), Node::Type::polyC));
    }
}

ReadType check_read(std::string const & read, Node * root, unsigned int length)
{
    if (read.size() < length) {
        return ReadType::length;
    } 
    return (ReadType)search_any(read, root); 
}

void filter_single_reads(std::ifstream & reads_f, std::ofstream & ok_f, std::ofstream & bad_f, 
                         Stats & stats, Node * root, int length)
{
    Seq read;

    while (read.read_seq(reads_f)) {
        ReadType type = check_read(read.seq, root, length);
        stats.update(type);
        if (type == ReadType::ok) {
            read.write_seq(ok_f);
        } else {
            read.update_id(type);
            read.write_seq(bad_f);
        }
    }
}

void filter_paired_reads(std::ifstream & reads1_f, std::ifstream & reads2_f,
                         std::ofstream & ok1_f, std::ofstream & ok2_f,
                         std::ofstream & bad1_f, std::ofstream & bad2_f,
                         std::ofstream & se1_f, std::ofstream & se2_f,
                         Stats & stats1, Stats & stats2,
                         Node * root, int length)
{
    Seq read1;
    Seq read2;

    while (read1.read_seq(reads1_f) && read2.read_seq(reads2_f)) {
        ReadType type1 = check_read(read1.seq, root, length);
        ReadType type2 = check_read(read2.seq, root, length);
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
    std::cout << "./rm_reads [-i raw_data.fastq | -1 raw_data1.fastq -2 raw_data2.fastq] -o output_dir --polyG 13 --length 50 --adapters adapters.dat" << std::endl;
}

int main(int argc, char ** argv)
{
	Node root('0');
	std::vector <std::pair<std::string, Node::Type> > patterns;

    std::string kmers, reads, out_dir;
    std::string reads1, reads2;
    char rez = 0;
    int length = 0;
    int polyG = 0;

    const struct option long_options[] = {
        {"length",required_argument,NULL,'l'},
        {"polyG",required_argument,NULL,'p'},
        {"adapters",required_argument,NULL,'a'},
        {NULL,0,NULL,0}
    };

    while ((rez = getopt_long(argc, argv, "1:2:l:p:a:i:o:", long_options, NULL)) != -1) {
        switch (rez) {
        case 'l':
            length = boost::lexical_cast<int>(optarg);
            break;
        case 'p':
            polyG = boost::lexical_cast<int>(optarg);
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
        case '?':
            print_help();
            return -1;
        }
    }

    if (kmers.empty() || out_dir.empty() || (
            reads.empty() &&
            (reads1.empty() || reads2.empty()))) {
        print_help();
        return -1;
    }

    std::ifstream kmers_f (kmers.c_str());
    if (!kmers_f.good()) {
        std::cout << "kmers file is bad" << std::endl;
        print_help();
        return -1;
    }

    init_type_names(length, polyG);

    build_patterns(kmers_f, polyG, patterns);

/*
    for (std::vector <std::string> ::iterator it = patterns.begin(); it != patterns.end(); ++it) {
        std::cout << *it << std::endl;
    }
    */

    if (patterns.empty()) {
        std::cout << "patterns are empty" << std::endl;
        return -1;
    }

	build_trie(root, patterns);
	add_failures(root);

    if (!reads.empty()) {
        std::string reads_base = basename(reads);
        std::ifstream reads_f (reads.c_str());
        std::ofstream ok_f((out_dir + "/" + reads_base + ".ok.fastq").c_str(), std::ofstream::out);
        std::ofstream bad_f((out_dir + "/" + reads_base + ".bad.fastq").c_str(), std::ofstream::out);

        if (!reads_f.good()) {
            std::cout << "reads file is bad" << std::endl;
            print_help();
            return -1;
        }

        if (!ok_f.good() || !bad_f.good()) {
            std::cout << "out file is bad" << std::endl;
            print_help();
            return -1;
        }

        Stats stats(reads);

        filter_single_reads(reads_f, ok_f, bad_f, stats, &root, length);

        std::cout << stats;

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
        std::ofstream bad1_f((out_dir + "/" + reads1_base + ".bad.fastq").c_str(),
                            std::ofstream::out);
        std::ofstream bad2_f((out_dir + "/" + reads2_base + ".bad.fastq").c_str(),
                            std::ofstream::out);

        if (!reads1_f.good() || !reads2_f.good()) {
            std::cout << "reads file is bad" << std::endl;
            print_help();
            return -1;
        }

        if (!ok1_f.good() || !ok2_f.good() || !bad1_f.good() || !bad2_f.good() ||
                !se1_f.good() || !se2_f.good()) {
            std::cout << "out file is bad" << std::endl;
            print_help();
            return -1;
        }

        Stats stats1(reads1);
        Stats stats2(reads2);

        filter_paired_reads(reads1_f, reads2_f, ok1_f, ok2_f,
                            bad1_f, bad2_f, se1_f, se2_f,
                            stats1, stats2,
                            &root, length);

        std::cout << stats1;
        std::cout << stats2;

        ok1_f.close();
        ok2_f.close();
        bad1_f.close();
        bad2_f.close();
        reads1_f.close();
        reads2_f.close();
    }
}
