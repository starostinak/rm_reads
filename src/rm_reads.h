#ifndef RM_READS_H
#define RM_READS_H

void build_patterns(std::ifstream & kmers_f, std::vector <std::pair <std::string, Node::Type> > & patterns, int polyG, bool filterN);
double get_dust_score(std::string const & read, int k);
std::string basename(std::string const & path);
void print_help();

#endif // RM_READS_H
