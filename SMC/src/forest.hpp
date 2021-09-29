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
    private:
        static unsigned _nspecies;
        std::list<strom::Tree::SharedPtr> _trees;
//        std::list<Clade> _clades;
        
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
        
        
//        TreeManip tm(tree);
//
//        clade._split.resize(_nspecies);
//        clade._split.setBitAt(i);
//        _clades.push_back(clade);
    }
    
//    createTrivialForest();
}


inline void Forest::showForest() {
    unsigned i = 0;
    for (auto t:_trees) {
        strom::TreeManip tm(t);
        cout << i << " " << tm.makeNewick(3, true) << "\n";
        i++;
    }
}
}
