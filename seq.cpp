#include <map>
#include "seq.h"

std::map <ReadType, std::string> type_names;
void init_type_names(int length, int polyG)
{
    type_names[ReadType::ok] = "ok";
    type_names[ReadType::adapter] = "adapter";
    type_names[ReadType::n] = "n";
    type_names[ReadType::polyG] = "polyG" + std::to_string(polyG);
    type_names[ReadType::polyC] = "polyC" + std::to_string(polyG);
    type_names[ReadType::length] = "length" + std::to_string(length);
}

const std::string & get_type_name (ReadType type) {
    return type_names[type];
}
