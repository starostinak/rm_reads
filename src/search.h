#ifndef SEARCH_H
#define SEARCH_H

#include <vector>
#include <list>
#include <map>
#include <string>
#include <cstddef>

class Node
{
public:

    enum Type {
        no_match = 0,
        adapter = 1,
        n,
        polyG,
        polyC
    };

    Node() : fail(NULL), type(Type::no_match) {}
    Node(char label) :
        label(label), fail(NULL), type(Type::no_match)
    {}

    ~Node()
    {
        for (auto it = links.begin(); it != links.end(); ++it) {
            delete *it;
        }
    }

    Node * next(char c)
    {
        for (size_t i = 0; i < links.size(); ++i) {
            if (links[i]->label == c)
                return links[i];
        }
        return NULL;
    }

    void update_node_stats(Type t, size_t adapt_id, size_t adapt_pos)
    {
        type = t;
        adapter_id_pos.push_back(std::make_pair(adapt_id, adapt_pos));
    }

    char label;
    Node * fail;
    Type type;
    std::list <std::pair <size_t, size_t> > adapter_id_pos;
    std::vector <Node *> links;
};

void build_trie(Node & root,
                std::vector <std::pair <std::string, Node::Type> > const & patterns,
                int errors = 0);
void add_failures(Node & root);
void go(Node * & curr, char c);
Node::Type find_match(Node * node);
Node::Type find_all_matches(Node * node, size_t pos,
                      std::map <size_t, std::vector <std::pair<size_t, size_t> > > & matches,
                      size_t errors);
Node::Type search_inexact(const std::string & text, Node * root,
                          std::vector <std::pair<std::string, Node::Type> > const & patterns, int errors);
Node::Type search_any(const std::string & text, Node * root);

#endif // SEARCH_H
