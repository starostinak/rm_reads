#include "search.h"

#include <list>
#include <fstream>
#include <algorithm>

unsigned int last_id = 1;

void build_trie(Node & root, std::vector <std::pair <std::string, Node::Type> > const & patterns)
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

Node::Type has_match(Node * node, int pos) {
    Node * curr = node;
    while (curr->fail != curr) {
        if (curr->type) {
            return curr->type;
        }
        curr = curr->fail;
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
        Node::Type match_type = has_match(curr, i);
        if(match_type) {
            return match_type;
        }
    }
    return Node::Type::no_match;
}
