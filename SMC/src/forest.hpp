#pragma once    

#include <stack>    
#include <memory>
#include <iostream>
#include <boost/format.hpp>

#include "lot.hpp"
extern Lot rng;

#include "node.hpp"
#include "forest.hpp"
using namespace std;


class TreeManip;
class Likelihood;
class Updater;
class TreeUpdater;
class PolytomyUpdater;
class Particle;


class Forest {

        friend class TreeManip;
        friend class Likelihood;
        friend class Updater;
        friend class TreeUpdater;
        friend class PolytomyUpdater;
        friend class Particle;
        

    public:

                                    Forest();
                                    ~Forest();

        unsigned                    numLeaves() const;
        unsigned                    numInternals() const;
        unsigned                    numNodes() const;
        void                        showForest();
        static void                 setNumSpecies(unsigned n);
        
    
    private:

        void                        clear();

        Node *                      _root;
        unsigned                    _nleaves;
        unsigned                    _ninternals;
        std::vector<Node *>         _preorder;
        std::vector<Node *>         _levelorder;
        std::vector<Node>           _nodes;
        std::vector<Node *>         _unused_nodes;
        static unsigned             _nspecies;
        void                        refreshPreorder();
        Node *                      findNextPreorder(Node * nd);
        std::string                 makeNewick(unsigned precision, bool use_names) const;
        void                        nextStep();
        void                        detachSubtree(Node * s);
        void                        insertSubtreeOnLeft(Node * s, Node * u);
        Node *                      findLeftSib(Node * nd);
        Node *                      getSubtreeAt(unsigned i);
        unsigned                    getNumSubtrees();
        void                        setEdgeLength(Node * nd);
        void                        createNewSubtree(unsigned t1, unsigned t2, unsigned nsubtrees);
    
    public:

        typedef std::shared_ptr<Forest> SharedPtr;
};
    

inline Forest::Forest() {
    //std::cout << "Constructing a forest" << std::endl;
    clear();
}

inline Forest::~Forest() {
    //std::cout << "Destroying a Forest" << std::endl;
}

inline void Forest::clear() {
    _nodes.clear();
    _preorder.clear();
    _levelorder.clear();
    _nodes.resize(2*_nspecies);
    //creating root node
    _root = &_nodes[_nspecies];
    _root->_name="root";
    _root->_left_child=0;
    _root->_right_sib=0;
    _root->_parent=0;
    _root->_number=_nspecies;
    _root->_edge_length=0.0;
    
    //creating subroot node
    Node* subroot=&_nodes[_nspecies+1];
    subroot->_name="subroot";
    subroot->_left_child=0;
    subroot->_right_sib=0;
    subroot->_parent=_root;
    subroot->_number=_nspecies;
    subroot->_edge_length=0.0;
    _root->_left_child=subroot;
    
    //create species
    for (unsigned i = 0; i < _nspecies; i++) {
        Node* nd=&_nodes[i];
        if (i==0) {
            subroot->_left_child=nd;
        }
        else {
            _nodes[i-1]._right_sib=nd;
        }
        nd->_name=(char)('A'+i);
        nd->_left_child=0;
        nd->_right_sib=0;
        nd->_parent=subroot;
        nd->_number=i;
        nd->_edge_length=0.0;
        }
    _nleaves=_nspecies;
    _ninternals=2;
    refreshPreorder();
    }

inline unsigned Forest::numLeaves() const {
    return _nleaves;
}

inline unsigned Forest::numInternals() const {
    return _ninternals;
}

inline unsigned Forest::numNodes() const {
    return (unsigned)_nodes.size();
}

inline void Forest::refreshPreorder() {
    // Create vector of node pointers in preorder sequence
    _preorder.clear();
    _preorder.reserve(_nodes.size() - 1); // _preorder does not include root node

    if (!_root)
        return;

    Node * first_preorder = _root->_left_child;

    // sanity check: first preorder node should be the only child of the root node
    assert(first_preorder->_right_sib == 0);

    Node * nd = first_preorder;
    _preorder.push_back(nd);

    while (true) {
        nd = findNextPreorder(nd);
        if (nd)
            _preorder.push_back(nd);
        else
            break;
    }   // end while loop
}

inline Node * Forest::findNextPreorder(Node * nd) {
    assert(nd);
    Node * next = 0;
    if (!nd->_left_child && !nd->_right_sib) {
        // nd has no children and no siblings, so next preorder is the right sibling of
        // the first ancestral node that has a right sibling.
        Node * anc = nd->_parent;
        while (anc && !anc->_right_sib)
            anc = anc->_parent;
        if (anc) {
            // We found an ancestor with a right sibling
            next = anc->_right_sib;
        }
        else {
            // nd is last preorder node in the tree
            next = 0;
        }
    }
    else if (nd->_right_sib && !nd->_left_child) {
        // nd has no children (it is a tip), but does have a sibling on its right
        next = nd->_right_sib;
    }
    else if (nd->_left_child && !nd->_right_sib) {
        // nd has children (it is an internal node) but no siblings on its right
        next = nd->_left_child;
    }
    else {
        // nd has both children and siblings on its right
        next = nd->_left_child;
    }
    return next;
}

inline void Forest::showForest() {
    cout << " " << makeNewick(3, true) << "\n";
}

inline std::string Forest::makeNewick(unsigned precision, bool use_names) const {
    std::string newick;
    const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
    const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
    const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
    std::stack<Node *> node_stack;


    for (auto nd : _preorder) {
        if (nd->_left_child) {
            newick += "(";
            node_stack.push(nd);
        }
        else {
            if (use_names) {
                newick += boost::str(boost::format(tip_node_name_format)
                    % nd->_name
                    % nd->_edge_length);
            } else {
                newick += boost::str(boost::format(tip_node_number_format)
                    % (nd->_number + 1)
                    % nd->_edge_length);
            }
            if (nd->_right_sib)
                newick += ",";
            else {
                Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                while (popped && !popped->_right_sib) {
                    node_stack.pop();
                    if (node_stack.empty()) {
                        newick += ")";
                        popped = 0;
                    }
                    else {
                        newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                        popped = node_stack.top();
                    }
                }
                if (popped && popped->_right_sib) {
                    node_stack.pop();
                    newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                    newick += ",";
                }
            }
        }
    }

    return newick;
}

inline void Forest::setNumSpecies(unsigned n){
    _nspecies=n;
}

inline void Forest::nextStep(){
    unsigned t1=0;
    unsigned t2=1;
    unsigned nsubtrees = getNumSubtrees();
    //don't use this when there's only one choice
        //i.e. # of trees = 2
    if (nsubtrees > 2) {
        t1 = ::rng.randint(0, nsubtrees-1);
        t2 = ::rng.randint(0, nsubtrees-1);

        //keep calling t2 until it doesn't equal t1
        while (t2 == t1) {
            t2 = ::rng.randint(0, nsubtrees-1);
        }
    }

    cout << "joining taxon " << t1 << " with taxon " << t2 << endl;
    createNewSubtree(t1, t2, nsubtrees);

//    Node * subtree1=getSubtreeAt(t1);
//    Node * subtree2 = getSubtreeAt(t2);
//
//    detachSubtree(subtree1);
//    detachSubtree(subtree2);
//
//    //creating new node
//    Node* new_nd=&_nodes[_nleaves+_ninternals];
//    new_nd->_name=" ";
//    new_nd->_left_child=subtree1;
//    new_nd->_right_sib=0;
//    new_nd->_parent=_root->_left_child;
//    new_nd->_number=_nleaves+_ninternals;
//    new_nd->_edge_length=rng.gamma(1.0, 1.0/nsubtrees);
//
//    cout << "New node branch length is: " << new_nd->_edge_length << endl;
//
//    _ninternals++;
//    subtree1 -> _right_sib=subtree2;
//    subtree1->_parent=new_nd;
//    subtree2->_parent=new_nd;
//
//    insertSubtreeOnLeft(new_nd, _root->_left_child);
//
//    refreshPreorder();
}

inline void Forest::createNewSubtree(unsigned t1, unsigned t2, unsigned nsubtrees) {
    Node * subtree1=getSubtreeAt(t1);
    Node * subtree2 = getSubtreeAt(t2);

    detachSubtree(subtree1);
    detachSubtree(subtree2);
    //creating new node
    Node* new_nd=&_nodes[_nleaves+_ninternals];
    new_nd->_name=" ";
    new_nd->_left_child=subtree1;
    new_nd->_right_sib=0;
    new_nd->_parent=_root->_left_child;
    new_nd->_number=_nleaves+_ninternals;
    new_nd->_edge_length=rng.gamma(1.0, 1.0/nsubtrees);
//    new_nd->_edge_length=0;


    cout << "New node branch length is: " << new_nd->_edge_length << endl;

    _ninternals++;
    subtree1 -> _right_sib=subtree2;
    subtree1->_parent=new_nd;
    subtree2->_parent=new_nd;

    insertSubtreeOnLeft(new_nd, _root->_left_child);

    refreshPreorder();
}

inline void Forest::detachSubtree(Node * s) {
    assert(s);
    assert(s->_parent);
    
    // Save pointers to relevant nodes
    Node * s_leftsib  = findLeftSib(s);
    Node * s_rightsib = s->_right_sib;
    Node * s_parent   = s->_parent;

    // Completely detach s and seal up the wound
    s->_parent = 0;
    s->_right_sib = 0;
    if (s_leftsib)
        s_leftsib->_right_sib = s_rightsib;
    else
        s_parent->_left_child = s_rightsib;
}

inline void Forest::insertSubtreeOnLeft(Node * s, Node * u) {
    assert(u);
    assert(s);
    s->_right_sib  = u->_left_child;
    s->_parent     = u;
    u->_left_child = s;
}

inline Node * Forest::findLeftSib(Node * nd) {
    assert(nd);
    assert(nd->_parent);
    Node * child = nd->_parent->_left_child;
    while (child && child->_right_sib != nd)
        child = child->_right_sib;
    return child;
}

inline Node * Forest::getSubtreeAt(unsigned i) {
    unsigned nsubtrees = 0;
    Node * the_nd = 0;
    for (auto nd : _preorder) {
        if (nd->_parent == _root->_left_child) {
            if (nsubtrees == i) {
                the_nd = nd;
                break;
            }
            nsubtrees++;
        }
    }
    assert(the_nd);
    return the_nd;
}

inline unsigned Forest::getNumSubtrees() {
    unsigned nsubtrees = 0;
    for (auto nd : _preorder) {
        if (nd->_parent == _root->_left_child) {
            nsubtrees++;
        }
    }
    cout << "Number of subtrees is: " << nsubtrees << endl;
    return nsubtrees;
}
