#pragma once

#include <string>
#include <vector>
#include  <iostream>
#include "split.hpp"
#include "partial_store.hpp"
#include "conditionals.hpp"

namespace proj {

    class Likelihood;
    class Forest;
    class Particle;

    class Node {
        friend class Likelihood;
        friend class Forest;
        friend class Particle;

        public:
#if defined (FASTER_SECOND_LEVEL)
        typedef vector<Node *>  ptr_vect_t;
#endif

                                        Node();
                                        ~Node();

                    Node *              getParent()                 {return _parent;}
                    Node *              getLeftChild()              {return _left_child;}
                    Node *              getRightSib()               {return _right_sib;}
                    int                 getNumber()                 {return _number;}
                    std::string         getName()                   {return _name;}
                    Split               getSplit()                  {return _split;}

                    double              getEdgeLength()             {return _edge_length;}
                    void                setEdgeLength(double v);
        
                    unsigned            countChildren() const;

                    void                clearPointers()             {_left_child = _right_sib = _parent = 0;}
                    void                resetNode();
                                        
                    static const double _smallest_edge_length;

                    typedef std::vector<Node>    Vector;
                    typedef std::vector<Node *>  PtrVector;
        
#if defined (FASTER_SECOND_LEVEL)
            static string           taxonNameToSpeciesName(string taxon_name);
            static void             setSpeciesBits(G::species_t & to_species, const G::species_t & from_species, bool init_to_zero_first);
            static void             setSpeciesBit(G::species_t & to_species, unsigned i, bool init_to_zero_first);
            const G::species_t &    getSpecies() const {return _species;}
#endif
        
        private:
        
            enum Flag {
                Selected   = (1 << 0),
                SelPartial = (1 << 1),
                SelTMatrix = (1 << 2),
                AltPartial = (1 << 3),
                AltTMatrix = (1 << 4)
            };

            void                clear();

            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            std::string         _name;
            double              _edge_length;
            Split               _split;
            int                 _flags;
            PartialStore::partial_t _partial;
            int                 _position_in_lineages;
#if defined (FASTER_SECOND_LEVEL)
        // Bitset of species (indices) compatible with this node
            G::species_t    _species;
            double          _height;
#endif
    };
    
    
    inline Node::Node() {
        //std::cout << "Creating Node object" << std::endl;
        clear();
    }

    inline Node::~Node() {
        //std::cout << "Destroying Node object" << std::endl;
    }

    inline void Node::clear() {
        _flags = 0;
        clearPointers();
        _number = -1;
        _name = "";
        _edge_length = _smallest_edge_length;
        _partial.reset();
#if defined (FASTER_SECOND_LEVEL)
        _height = 0.0;
#endif
    }

    inline void Node::setEdgeLength(double v) {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

    inline unsigned Node::countChildren() const {
        unsigned n_children = 0;
        for (Node * child = _left_child; child; child=child->_right_sib) {
            n_children ++;
        }
        return n_children;
    }

#if defined (FASTER_SECOND_LEVEL)
    inline string Node::taxonNameToSpeciesName(string tname) {
        vector<string> before_after;
        split(before_after, tname, boost::is_any_of("^"));
        if (before_after.size() != 2)
            throw XProj(format("Expecting taxon name to conform to taxon^species pattern: %s") % tname);
        return before_after[1];
    }
#endif

#if defined (FASTER_SECOND_LEVEL)
    inline void Node::setSpeciesBits(G::species_t & to_species, const G::species_t & from_species, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (G::species_t)0;

        // Copy bits in from_species to to_species
        to_species |= from_species;
    }
#endif

#if defined (FASTER_SECOND_LEVEL)
    inline void Node::setSpeciesBit(G::species_t & to_species, unsigned i, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (G::species_t)0;
            
        // Set ith bit in to_species
        to_species |= ((G::species_t)1 << i);
    }
#endif

    inline void Node::resetNode() {
        _parent=0;
        _right_sib=0;
    }
}

