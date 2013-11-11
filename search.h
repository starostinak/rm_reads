#ifndef SEARCH_H
#define SEARCH_H

#include <vector>
#include <string>
#include <cstddef>

extern unsigned int last_id;

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

    Node() : fail(NULL), id(last_id++), type(Type::no_match) {}
    Node(char label) :
        label(label), fail(NULL), id(last_id++), type(Type::no_match)
    {}

    Node * next(char c)
    {
        for (size_t i = 0; i < links.size(); ++i) {
            if (links[i]->label == c)
                return links[i];
        }
        return NULL;
    }

    char label;
    Node * fail;
    int id;
    Type type;
    std::vector <Node *> links;
};

void build_trie(Node & root,
                std::vector <std::pair <std::string, Node::Type> > const & patterns);
void add_failures(Node & root);
void go(Node * & curr, char c);
Node::Type has_match(Node * node, int pos);
Node::Type search_any(const std::string & text, Node * root);

#endif // SEARCH_H
