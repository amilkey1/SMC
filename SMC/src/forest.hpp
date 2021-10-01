//
//  forest.hpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
// creates the forest class
//forest class contains a set of clades

#pragma once

#include <stack>
#include <memory>
#include <iostream>
#include "conditionals.hpp"


#include <list>
#include "tree_manip.hpp"
#include "tree.hpp"

#include "lot.hpp"
extern strom::Lot rng;

using namespace std;
//using namespace strom;

namespace strom {
class Forest {
    //set of clades
    //construct a vector of clades
    friend class TreeManip;
    
    public:
        Forest();
        void showForest();
        void chooseTrees();
    private:
        static unsigned _nspecies;
        int t1;
        int t2;
        std::list<strom::Tree::SharedPtr> _trees;
};

inline Forest::Forest() {
    for (unsigned i = 0; i < _nspecies; i++) {
        strom::Tree::SharedPtr tree = strom::Tree::SharedPtr(new strom::Tree());
        tree->_nodes.resize(2*_nspecies-1);
        for (unsigned j = 0; j < _nspecies; j++) {
            char x = (char)('A'+j);
            tree->_nodes[j]._name=x;
            tree->_nodes[j]._left_child=0;
            tree->_nodes[j]._right_sib=0;
            tree->_nodes[j]._parent=0;
            tree->_nodes[j]._number=j;
            tree->_nodes[j]._edge_length=0.0;
        }
        tree->_nleaves=_nspecies;
        tree->_is_rooted=true;
        tree->_root=&tree->_nodes[i];
        tree->_ninternals=0;
        tree->_preorder.push_back(tree->_root);
        _trees.push_back(tree);
    }
}

inline void Forest::showForest() {
    unsigned i = 0;
    for (auto t:_trees) {
        strom::TreeManip tm(t);
        cout << i << " " << tm.makeNewick(3, true) << "\n";
        i++;
    }
}

inline void Forest::chooseTrees() {
    //choose 2 trees to join
    t1 = rng.randint(0, (int) _trees.size()-1);
    t2 = rng.randint(0, (int) _trees.size()-1);
    
    //keep calling t2 until it doesn't equal t1
    while (t2 == t1) {
        t2 = rng.randint(0, (int) _trees.size()-1);
    }
    
    //don't use this when there's only one choice
        //i.e. # of trees = 2
    if (_trees.size() == 2) {
        //then choose t1 = 0 and t2 = 1
        t1 = 0;
        t2 = 1;
    }
    cout << "join taxon " << t1 << " with taxon " << t2 << endl;
//    cout <<"tree size is " << _trees.size() << endl;
}
}
