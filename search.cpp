#include "search.h"

#include <list>
#include <map>
#include <fstream>
#include <algorithm>

unsigned int last_id = 1;

void build_trie(Node & root, std::vector <std::pair <std::string, Node::Type> > const & patterns, int errors)
{
    for (auto it = patterns.begin(); it != patterns.end(); ++it) {
        Node * curr_node = &root;
        const std::string & pattern = it->first;
        size_t pattern_size = pattern.size();
        for (size_t j = 0; j < pattern_size; ++j) {
            Node * next = curr_node->next(pattern[j]);
            if (next == NULL) {
                Node * new_node = new Node(pattern[j]);
                curr_node->links.push_back(new_node);
                curr_node = new_node;
            } else {
                curr_node = next;
            }
            if (j == pattern_size - 1) {
                curr_node->type = it->second;
            }
            if (it->second == Node::Type::adapter && errors != 0) {
                if (errors == 1) {
                    if (j == pattern_size / 2) {
                        curr_node->update_node_stats(Node::Type::adapter, (size_t)(it - patterns.begin()), pattern_size/2);
                        curr_node = &root;
                    } if (j == pattern_size - 1) {
                        curr_node->update_node_stats(Node::Type::adapter, (size_t)(it - patterns.begin()), pattern_size - 1);
                    }
                } if (errors == 2) {
                    if (j == pattern_size / 3 || j == pattern_size * 2 / 3) {
                        if (j == pattern_size / 3) {
                            curr_node->update_node_stats(Node::Type::adapter, (size_t)(it - patterns.begin()), pattern_size / 3);
                        } else {
                            curr_node->update_node_stats(Node::Type::adapter, (size_t)(it - patterns.begin()), pattern_size * 2 / 3);
                        }
                        curr_node = &root;
                    } if (j == pattern_size - 1) {
                        curr_node->update_node_stats(Node::Type::adapter, (size_t)(it - patterns.begin()), pattern_size - 1);
                    }
                }
            }
        }
//		std::cout << patterns[i] << std::endl;
    }
}

void add_failures(Node & root)
{
    std::list <Node *> queue;
    queue.push_back(&root);
    root.fail = &root;
    do {
        Node * curr = queue.front();
        queue.pop_front();
        for (std::vector <Node *>::iterator it = curr->links.begin(); it != curr->links.end(); ++it) {
            if (curr != &root) {
                Node * parent = curr;
                do {
                    parent = parent->fail;
                    (*it)->fail = parent->next((*it)->label);
                } while(!(*it)->fail && parent != &root);
            }
            if (!(*it)->fail) {
                (*it)->fail = &root;
            }
        }
        queue.insert(queue.end(), curr->links.begin(), curr->links.end());
    } while(queue.size());
}

void go(Node * & curr, char c)
{
    while (!curr->next(c) && curr != curr->fail) {
        curr = curr->fail;
    }
    Node * next = curr->next(c);
    if (next) {
        curr = next;
    }
}

Node::Type find_all_matches(Node * node, size_t pos,
                      std::map <size_t, std::vector <std::pair<size_t, size_t> > > & matches,
                      size_t errors) {
    Node * curr = node;
    while (curr->fail != curr) {
        if (curr->type) {
            if (curr->type != Node::Type::adapter) {
                return curr->type;
            }
            for (auto it = curr->adapter_id_pos.begin(); it != curr->adapter_id_pos.end(); ++it) {
                std::vector <std::pair<size_t, size_t> > & curr_matches = matches[it->first];
                curr_matches.push_back(std::make_pair(pos, it->second));
                if (curr_matches.size() >= errors + 1){
                    size_t begin_pos = pos - it->second;
                    size_t pattern_size = it->second + 1;
                    if (errors == 1 &&
                        std::binary_search(curr_matches.begin(), curr_matches.end(), std::pair <size_t, size_t> (begin_pos + pattern_size / 2, pattern_size / 2))) {
                        return curr->type;
                    }
                    if (errors == 2 &&
                        std::binary_search(curr_matches.begin(), curr_matches.end(),
                                           std::pair <size_t, size_t> (begin_pos + pattern_size / 3, pattern_size / 3)) &&
                        std::binary_search(curr_matches.begin(), curr_matches.end(),
                                           std::pair <size_t, size_t> (begin_pos + pattern_size * 2 / 3, pattern_size * 2 / 3))) {
                        return curr->type;
                    }
                }
            }
        }
        curr = curr->fail;
    }
    return Node::Type::no_match;
}

Node::Type find_match(Node * node) {
    Node * curr = node;
    while (curr->fail != curr) {
        if (curr->type) {
            return curr->type;
        }
        curr = curr->fail;
    }
    return Node::Type::no_match;
}

bool cmp_upper(char i, char j)
{
    return std::toupper(i) == std::toupper(j);
}

int count_errors(std::string const & text, size_t text_pos,
                 std::string const & pattern, size_t pattern_pos,
                 size_t length, int err_max)
{
    if (text_pos < 0) {
        return err_max + 1;
    }
    auto text_it = text.begin() + text_pos;
    auto text_last = text_it + length;
    auto pattern_it = pattern.begin() + pattern_pos;
    int errors = 0;
    while (errors <= err_max) {
        auto mismatch = std::mismatch(text_it, text_last, pattern_it, cmp_upper);
        if (mismatch.first == text_last) {
            return errors;
        }
        ++errors;
        text_it = ++mismatch.first;
        pattern_it = ++mismatch.second;
    }
    return errors;
}

bool check_partial_matches(std::string const & text,
                           std::vector <std::pair<std::string, Node::Type> > const & patterns,
                           std::map <size_t, std::vector <std::pair <size_t, size_t> > > & matches, int errors)
{
    for (auto it = matches.begin(); it != matches.end(); ++it) {
        size_t pattern_size = patterns[it->first].first.size();
        while (!it->second.empty()) {
            auto start_match = it->second.begin();
            int res_errors = 0;
            size_t begin_pos = start_match->first;
            if (errors == 1) {
                if (start_match->second == pattern_size - 1) {
                    begin_pos -= pattern_size - 1;
                    res_errors = count_errors(text, begin_pos,
                                              patterns[it->first].first, 0, pattern_size/2, 1);
                } else {
                    begin_pos -= pattern_size/2;
                    res_errors = count_errors(text, start_match->first,
                                              patterns[it->first].first, start_match->second, pattern_size - start_match->second, 1);
                }
            } else {
                if (start_match->second == pattern_size - 1) {
                    begin_pos -= pattern_size - 1;
                    res_errors = count_errors(text, begin_pos,
                                              patterns[it->first].first, 0, pattern_size / 3, 1);
                    if (res_errors < 2) {
                        res_errors += count_errors(text, begin_pos + pattern_size / 3,
                                                  patterns[it->first].first, pattern_size / 3, pattern_size * 2/3 - pattern_size / 3, 1);
                    } else {
                        ++res_errors;
                    }
                } else if (start_match->second == pattern_size * 2/3) {
                    begin_pos -= pattern_size * 2/3;
                    auto found = std::lower_bound(it->second.begin(), it->second.end(),
                                                  std::pair <size_t, size_t> (begin_pos + pattern_size - 1, pattern_size - 1));
                    if (*found == std::pair <size_t, size_t> (begin_pos + pattern_size - 1, pattern_size - 1)){
                        res_errors = count_errors(text, begin_pos,
                                                  patterns[it->first].first, 0, pattern_size / 3, 2);
                        it->second.erase(found);
                    } else {
                        res_errors = count_errors(text, begin_pos,
                                                  patterns[it->first].first, 0, pattern_size / 3, 1);
                        if (res_errors < 2) {
                            res_errors += count_errors(text, begin_pos + pattern_size * 2/3,
                                                      patterns[it->first].first, pattern_size * 2/3, pattern_size - pattern_size * 2/3, 1);
                        } else {
                            ++res_errors;
                        }
                    }
                } else {
                    begin_pos -= pattern_size/3;
                    auto found = std::lower_bound(it->second.begin(), it->second.end(),
                                                  std::pair <size_t, size_t> (begin_pos + pattern_size * 2/3, pattern_size * 2/3));
                    if (*found == std::pair <size_t, size_t> (begin_pos + pattern_size * 2/3, pattern_size * 2/3)) {
                        res_errors = count_errors(text, begin_pos + pattern_size / 3,
                                                  patterns[it->first].first, pattern_size * 2/3, pattern_size - pattern_size * 2/3, 2);
                        it->second.erase(found);
                    } else {
                        found = std::lower_bound(it->second.begin(), it->second.end(),
                                                 std::pair <size_t, size_t> (begin_pos + pattern_size - 1, pattern_size - 1));
                        if (*found == std::pair <size_t, size_t> (begin_pos + pattern_size - 1, pattern_size - 1)) {
                            res_errors = count_errors(text, begin_pos + pattern_size / 3,
                                                      patterns[it->first].first, pattern_size / 3, pattern_size * 2/3 - pattern_size * 1/3, 2);
                            it->second.erase(found);
                        } else {
                            res_errors = count_errors(text, begin_pos + pattern_size / 3,
                                                      patterns[it->first].first, pattern_size / 3, pattern_size * 2/3 - pattern_size * 1/3, 1);
                            if (res_errors < 2) {
                                res_errors += count_errors(text, begin_pos + pattern_size * 2/3,
                                                          patterns[it->first].first, pattern_size * 2/3, pattern_size - pattern_size * 2/3, 1);
                            } else {
                                ++res_errors;
                            }
                        }
                    }
                }
            }
            it->second.erase(start_match);
            if (res_errors <= errors) {
                return true;
            }
        }
    }
    return false;
}

Node::Type search_inexact(const std::string & text, Node * root,
                          std::vector <std::pair<std::string, Node::Type> > const & patterns,
                          int errors)
{
    size_t text_len = text.size();
    std::map <size_t, std::vector <std::pair <size_t, size_t> > > matches; // value - <text_pos, adapter_pos>
    Node * curr = root;
    for (size_t i = 0; i < text_len; ++i) {
        char c = (text[i] > 96) ? text[i] - 32 : text[i];
        go(curr, c);
        Node::Type match_type = find_all_matches(curr, i, matches, errors);
        if(match_type) {
            return match_type;
        }
    }
    if (check_partial_matches(text, patterns, matches, errors)) {
        return Node::Type::adapter;
    }
    return Node::Type::no_match;
}

Node::Type search_any(const std::string & text, Node * root)
{
    size_t text_len = text.size();
    Node * curr = root;
    for (size_t i = 0; i < text_len; ++i) {
        char c = (text[i] > 96) ? text[i] - 32 : text[i];
        go(curr, c);
        Node::Type match_type = find_match(curr);
        if(match_type) {
            return match_type;
        }
    }
    return Node::Type::no_match;
}
