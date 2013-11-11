#include "stats.h"

void Stats::update(ReadType type, bool paired)
{
    auto it = reads.find(type);
    if (it == reads.end()) {
        reads[type] = 1;
    } else {
        ++it->second;
    }
    ++complete;
    if (type == ReadType::ok) {
        if (paired) {
            ++pe;
        } else {
            ++se;
        }
    }
}

std::ostream & operator << (std::ostream & out, const Stats & stats)
{
    out << stats.filename << std::endl;
    unsigned int bad = 0;
    for (auto it = stats.reads.begin(); it != stats.reads.end(); ++it) {
        out << "\t" << get_type_name(it->first) << "\t" << it->second << std::endl;
        if (it->first != ReadType::ok) {
            bad += it->second;
        }
    }
    out << "\t" << "fraction" << (double)(stats.complete - bad)/stats.complete << std::endl;
    if (stats.pe) {
        out << "\t" << "se\t" << stats.se << std::endl;
        out << "\t" << "pe\t" << stats.pe << std::endl;
    }
    return out;
}
