#pragma once

#include <stack>
#include <memory>
#include <iostream>
#include <boost/format.hpp>
#include <vector>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/math/tools/minima.hpp>
#include <thread>
#include <mutex>
#include <algorithm>
#include "conditionals.hpp"
#include <fstream>
#include <map>
#include <unordered_map>

#include "lot.hpp"
#include "g.hpp"
#include "stopwatch.hpp"
extern proj::Lot rng;

extern proj::StopWatch stopwatch;

#include "node.hpp"

namespace proj {

using namespace std;

class Likelihood;
class Particle;

class SpeciesForest {

        friend class Likelihood;
        friend class Particle;

    public:
                                    SpeciesForest();
                                    ~SpeciesForest();
        SpeciesForest(const SpeciesForest & other);
    
        typedef tuple<double, unsigned, vector<G::species_t> >  coalinfo_t;

        unsigned                        numInternals() const;
        unsigned                        numNodes() const;
        void                            showForest();
        void                            createDefaultTree(Lot::SharedPtr lot);
        void operator=(const SpeciesForest & other);
        void                            debugForest();

    private:

        void                            clear();
        Node *                          findNextPreorder(Node * nd);
        string                          makeNewick(unsigned precision, bool use_names);
        string                          makeAltNewick(unsigned precision, bool use_names);
        string                          makePartialNewick(unsigned precision, bool use_names);
 
        void                            setUpSpeciesForest();
        void                            setSpeciesFromNodeName(Node * nd);
        tuple<string,string, string>    speciesTreeProposalSim(Lot::SharedPtr lot);
#if defined (LAZY_COPYING)
        tuple<G::species_t,G::species_t, G::species_t>    speciesTreeProposal(Lot::SharedPtr lot);
#else
        tuple<string,string, string>    speciesTreeProposal(Lot::SharedPtr lot);
#endif
        void                            updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                            updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                            revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        void                            chooseSpeciesIncrement(Lot::SharedPtr lot);
        pair<double,double>             chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth);
        double                          calcTopologyPrior(unsigned nlineages);
    
        vector< pair<double, Node *>>   sortPreorder();
        void                            refreshPreorder();
        void                            buildFromNewick(const string newick, bool rooted, bool allow_polytomies);
        vector<pair<tuple<string, string, string>, double>> buildFromNewickMPI(const string newick, bool rooted, bool allow_polytomies, Lot::SharedPtr lot);
        void                            extractNodeNumberFromName(Node * nd, std::set<unsigned> & used);
        void                            stripOutNexusComments(std::string & newick);
        unsigned                        countNewickLeaves(const std::string newick);
        unsigned                        countNewickInternals(const std::string newick);
        void                            extractEdgeLen(Node * nd, std::string edge_length_string);
        void                            renumberInternals();
        bool                            canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);
        vector<tuple<string, string, string>>              buildFromNewickTopology(const string newick);
    
        vector<pair<tuple<string, string, string>, double>>                            resetLineages(Lot::SharedPtr lot);
        vector<pair<tuple<string, string, string>, double>> resetT();
    
        void                            saveCoalInfoInitial();
        void                            saveCoalInfoSpeciesTree(vector<SpeciesForest::coalinfo_t> & coalinfo_vect, bool cap);
        void                            addCoalInfoElem(const Node *, vector<coalinfo_t> & recipient);
        void                            buildCoalInfoVect();
        void                            fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const;
        static bool                     subsumed(G::species_t test_species, G::species_t subtending_species);
        void                            refreshAllPreorders() const;
        void                            refreshPreorderNew(vector<Node*> & preorder) const;
        Node *                          findNextPreorderNew(Node * nd) const;
        pair<double,double>             chooseSpeciesIncrementOnlySecondLevel(Lot::SharedPtr lot, double max_depth);
        void                            setTreeHeight();
        void                            storeSplits(set<Split> & internal_splits, set<Split> & leaf_splits);
    

        vector<coalinfo_t>                  _coalinfo;
        mutable vector<Node::ptr_vect_t>    _preorders;
        mutable unsigned                    _next_node_number;

        void                            setNodeHeights();
    
        std::vector<Node *>             _lineages;
        vector<Node>                    _nodes;

        unsigned                        _ninternals;
        unsigned                        _npatterns;
        double                          _last_edge_length;
    
        double                          _log_joining_prob;
        vector<pair<double, double>>    _increments_and_priors;
    
#if defined (DEBUG_MODE)
        pair<Node*, Node*>              _species_joined;
#endif
    
        vector<Node*>                   _preorder;
    
        double                          _forest_height;
        double                          _forest_length;
    
#if defined (DEBUG_MODE)
        void                            showSpeciesJoined();
#endif
        double                          getTreeLength();
        void                            calcTreeLength();
        double                          getLineageHeight(Node* nd);
        void                            addIncrement(double increment);
    
#if defined (LAZY_COPYING)
        void                            saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap = false) const;
#endif
            
    public:

        typedef std::shared_ptr<Forest> SharedPtr;
 };


    inline SpeciesForest::SpeciesForest() {
        //std::cout << "Constructing a forest" << std::endl;
        clear();
    }

    inline SpeciesForest::~SpeciesForest() {
        //std::cout << "Destroying a Forest" << std::endl;
    }

    inline void SpeciesForest::clear() {
        _nodes.clear();
        _lineages.clear();
        _npatterns = 0;
        _last_edge_length = 0.0;
        _lineages.clear();
        _log_joining_prob = 0.0;
        _ninternals=0;
        _preorder.clear();
        _forest_length = 0.0;
        _forest_height = 0.0;
        _coalinfo.clear();
        _preorders.clear();
    }

    inline SpeciesForest::SpeciesForest(const SpeciesForest & other) {
        clear();
        *this = other;
    }

    inline unsigned SpeciesForest::numInternals() const {
        return _ninternals;
    }

    inline unsigned SpeciesForest::numNodes() const {
        return (unsigned)_nodes.size();
    }

    inline Node * SpeciesForest::findNextPreorder(Node * nd) {
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

    inline void SpeciesForest::showForest() {
        cout << " species tree: " << endl;
        cout << " " << makeNewick(15, true) << "\n";
        cout << "\n";
    }

    inline string SpeciesForest::makePartialNewick(unsigned precision, bool use_names) {
        if (_lineages.size() == 1) {
            return makeNewick(precision, use_names);
        }
        else {
            // this function makes a newick string for a partially constructed tree
            string newick = "(";
            const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
            const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
            const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
            stack<Node *> node_stack;

            unsigned i = 0;
            unsigned a = 0;
            for (auto &lineage : _lineages) {
                Node * nd = lineage;
                while (nd) {
                    if (nd->_left_child) {
                        a++;
                        // internal node
                        newick += "(";
                        node_stack.push(nd);
                    }
                    else {
                        a++;
                        // leaf node
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
                                    //newick += ")";
                                    newick += boost::str(boost::format(internal_node_format) % lineage->_edge_length);
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
                    }   // leaf node
                    nd = findNextPreorder(nd);
                }   // while (subnd)...

                if (i < _lineages.size() - 1)
                    newick += ",";
                ++i;
            }
            newick += ")";

            return newick;
        }
    }

    inline string SpeciesForest::makeNewick(unsigned precision, bool use_names) {
        if (_lineages.size() > 1) {
            return makePartialNewick(precision, use_names);
        }

        else {
            string newick = "";
            const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
            const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
            const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
            stack<Node *> node_stack;

                unsigned i = 0;
                unsigned a = 0;
                for (auto &lineage : _lineages) {
                    Node * nd = lineage;
                    while (nd) {
                        if (nd->_left_child) {
                            a++;
                            // internal node
                            newick += "(";
                            node_stack.push(nd);
                        }
                        else {
                            a++;
                            // leaf node
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
                                            //newick += ")";
                                            if (lineage->_edge_length != 0.0) {
                                                newick += boost::str(boost::format(internal_node_format) % lineage->_edge_length);
                                            }
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
                            }   // leaf node
                        }
                        nd = findNextPreorder(nd);
                    }   // while (subnd)...

                    if (i < _lineages.size() - 1)
                        newick += ",";
                    ++i;
                }
                newick += ")";

                return newick;
            }
        }

    inline string SpeciesForest::makeAltNewick(unsigned precision, bool use_names) {
        use_names = false;
        if (_lineages.size() > 1) {
            return makePartialNewick(precision, use_names);
        }

        else {
            string newick = "";
            const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
            const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
            const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
            stack<Node *> node_stack;

                unsigned i = 0;
                unsigned a = 0;
                for (auto &lineage : _lineages) {
                    Node * nd = lineage;
                    while (nd) {
                        if (nd->_left_child) {
                            a++;
                            // internal node
                            newick += "(";
                            node_stack.push(nd);
                        }
                        else {
                            a++;
                            // leaf node
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
                                            //newick += ")";
                                            if (lineage->_edge_length != 0.0) {
                                                newick += boost::str(boost::format(internal_node_format) % lineage->_edge_length);
                                            }
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
                            }   // leaf node
                        }
                        nd = findNextPreorder(nd);
                    }   // while (subnd)...

                    if (i < _lineages.size() - 1)
                        newick += ",";
                    ++i;
                }
                newick += ")";

                return newick;
            }
        }

    inline void SpeciesForest::createDefaultTree(Lot::SharedPtr lot) {
        clear();
        //create taxa
        assert (lot != nullptr);
        double edge_length = lot->gamma(1.0, 1.0/G::_ntaxa);
        
        _lineages.reserve(_nodes.size());
        
        for (unsigned i = 0; i < G::_nspecies; i++) {
            Node * nd = &(_nodes[i]);
            nd->_right_sib=0;
            nd->_name="";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length = edge_length;
            nd->_position_in_lineages=i;
            }
        _ninternals=0;
        _last_edge_length = 0.0;
    }

    inline void SpeciesForest::operator=(const SpeciesForest & other) {
        _nodes.resize(other._nodes.size()); // don't need to clear these members because they will get overwritten
        _lineages.resize(other._lineages.size());
            
        _ninternals         = other._ninternals;
        _last_edge_length   = other._last_edge_length;
        _increments_and_priors = other._increments_and_priors;
        _preorder.resize(other._preorder.size());
        _forest_length = other._forest_length;
        _forest_height = other._forest_height;
            
        // the following data members apply only to the first round
//        if (!G::_in_second_level) {
            _log_joining_prob = other._log_joining_prob;
//        }
        
//        if (G::_in_second_level) {
//            _coalinfo = other._coalinfo; // TODO: don't copy this because it gets reset - but maybe it's faster to copy it?
//        }
            
    #if defined (DEBUG_MODE)
        _species_joined = other._species_joined;
    #endif

            // copy tree itself
            
            for (auto & othernd : other._nodes) {
                // get number of next node in preorder sequence (serves as index of node in _nodes vector)
                int k = othernd._number;

                if (k>-1) {
                    Node* nd = &_nodes[k];

                // copy parent
                    if (othernd._parent) {
                        unsigned parent_number = othernd._parent->_number;
                        Node * parent = &_nodes[parent_number];
                        nd->_parent = parent;
                    }
                    else {
                        nd->_parent = 0;
                    }

                // copy left child
                    if (othernd._left_child) {
                    unsigned left_child_number = othernd._left_child->_number;
                        Node * left_child = &_nodes[left_child_number];
                        nd->_left_child = left_child;
                }
                    else {
                        nd->_left_child = 0;
                    }

                // copy right sibling
                if (othernd._right_sib) {
                    unsigned right_sib_number = othernd._right_sib->_number;
                    Node * right_sib = &_nodes[right_sib_number];
                    nd->_right_sib = right_sib;
                }
                else {
                    nd->_right_sib = 0;
                }

                nd->_number = othernd._number;
                nd->_name = othernd._name;
                nd->_edge_length = othernd._edge_length;
                nd->_position_in_lineages = othernd._position_in_lineages;
                nd->_height = othernd._height;
                nd->_species = othernd._species;
                    
                    // don't need to copy these members
                    // nd->_split = othernd._split;
                    // nd->_species = othernd._species;
                }
            }

            unsigned j = 0;
            for (auto & othernd : other._lineages) {
                unsigned k = othernd->_number;
                Node * nd = &_nodes[k];
                _lineages[j] = nd;
                j++;
            }
            
            if (other._preorder.size() > 0) {
                unsigned m = 0;
                for (auto & othernd : other._preorder) {
                    unsigned n = othernd->_number;
                    Node * nd = &_nodes[n];
                    _preorder[m] = nd;
                    m++;
                }
            }
    }

    inline void SpeciesForest::setUpSpeciesForest() {
        //create species
        _nodes.resize(2*G::_nspecies - 1);
        _lineages.reserve(_nodes.size());
        
        //create taxa
        for (unsigned i = 0; i < G::_nspecies; i++) {
            Node * nd = &(_nodes[i]);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_height = 0.0;
            nd->_position_in_lineages=i;
            nd->_name=G::_species_names[i];
            _lineages.push_back(nd);
#if defined (LAZY_COPYING)
            setSpeciesFromNodeName(nd);
#endif
            }
        
        _ninternals=0;
    }

#if defined (LAZY_COPYING)
    inline void SpeciesForest::setSpeciesFromNodeName(Node * nd) {
        auto it = find(G::_species_names.begin(), G::_species_names.end(), nd->_name);
        if (it == G::_species_names.end())
            throw XProj(str(format("Could not find an index for the species name \"%s\"") % nd->_name));
        else {
            unsigned i = (unsigned)std::distance(G::_species_names.begin(), it);
            Node::setSpeciesBit(nd->_species, i, /*init_to_zero_first*/true);
        }
    }
#endif

    inline pair<double,double> SpeciesForest::chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth) {
        unsigned nlineages = (unsigned) _lineages.size();
        
        if (max_depth > 0.0) {
            double rate = (G::_lambda)*_lineages.size();
            
            double u = lot->uniform();
            double inner_term = 1-exp(-rate*max_depth);
            _last_edge_length = -log(1-u*inner_term)/rate;
            assert (_last_edge_length < max_depth);
            for (auto&  nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
             
#if !defined (HIERARCHICAL_FILTERING)
            // lorad only works if all topologies the same - then don't include the prior on joins b/c it is fixed
            double increment_prior = (log(rate)-_last_edge_length*rate);
            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior)); // do not include constrained factor in increment prior
#endif
        }
        else {
            double rate = G::_lambda*_lineages.size();
            
            assert (lot != nullptr);
            _last_edge_length = lot->gamma(1.0, 1.0/rate);

            for (auto & nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }

#if !defined (HIERARCHICAL_FILTERING)
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;
            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));
#endif
        }
        
        double constrained_factor = log(1 - exp(-1*nlineages*G::_lambda*max_depth));
        
        _forest_height += _last_edge_length;
        
        return make_pair(_last_edge_length, constrained_factor);

    }

    inline pair<double,double> SpeciesForest::chooseSpeciesIncrementOnlySecondLevel(Lot::SharedPtr lot, double max_depth) {
        double nlineages = (double) _lineages.size();
        
            max_depth = max_depth - _forest_height;
            
            assert (max_depth >= 0.0);
        
        if (max_depth > 0.0) {
            double rate = (G::_lambda)*_lineages.size();
            
            double u = lot->uniform();
            double inner_term = 1-exp(-rate*max_depth);
            _last_edge_length = -log(1-u*inner_term)/rate;
            assert (_last_edge_length < max_depth);
            for (auto& nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            // lorad only works if all topologies the same - then don't include the prior on joins b/c it is fixed
            double increment_prior = (log(rate)-_last_edge_length*rate);
                        
            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior)); // do not include constrained factor in increment prior
        }
        else {
            double rate = G::_lambda*_lineages.size();
            
            assert (lot != nullptr);
            _last_edge_length = lot->gamma(1.0, 1.0/rate);

            for (auto& nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));

        }
        
        double constrained_factor = log(1 - exp(-1*nlineages*G::_lambda*max_depth));
        
        _forest_height += _last_edge_length;
        
        return make_pair(_last_edge_length, constrained_factor);

    }

    inline void SpeciesForest::chooseSpeciesIncrement(Lot::SharedPtr lot) {
        double rate = G::_lambda*_lineages.size();
        
        assert (lot != nullptr);
        _last_edge_length = lot->gamma(1.0, 1.0/rate);

        for (auto & nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
        
        _forest_height += _last_edge_length;
    }


#if defined (LAZY_COPYING)
    inline tuple<G::species_t,G::species_t, G::species_t> SpeciesForest::speciesTreeProposal(Lot::SharedPtr lot) {
        // this function creates a new node and joins two species
        
        bool done = false;
        Node* subtree1 = nullptr;
        Node* subtree2 = nullptr;
        
        while (!done) {
            assert (lot != nullptr);
            pair<unsigned, unsigned> t = lot->nchoose2((unsigned) _lineages.size());
            assert (t.first != t.second);
            
            subtree1=_lineages[t.first];
            subtree2=_lineages[t.second]; // TODO: update species bits
            assert (t.first < _lineages.size());
            assert (t.second < _lineages.size());
            assert(!subtree1->_parent && !subtree2->_parent);
            
            if (G::_outgroup != "none") {
                if (subtree1->_name != G::_outgroup && subtree2->_name != G::_outgroup && _lineages.size() > 2) { // outgroup can only be chosen on the last step
                    done = true;
                }
                else if (_lineages.size() == 2) {
                    done = true;
                }
            }
            else {
                done = true;
            }
            if (G::_outgroup == "none") {
                assert (done == true);
            }
        }
        
        Node * new_nd = &_nodes[G::_nspecies + _ninternals];
        assert (new_nd->_parent==0);
       assert (new_nd->_number == -1);
       assert (new_nd->_right_sib == 0);
        new_nd->_number=G::_nspecies+_ninternals;
        new_nd->_name+=boost::str(boost::format("node-%d")%new_nd->_number);
        new_nd->_edge_length=0.0;
        _ninternals++;

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd);

    #if defined (DEBUG_MODE)
        if (_lineages.size() > 1) {
            _species_joined = make_pair(subtree1, subtree2); // last step just joins remaining two
        }
    #endif
        
        new_nd->_height = _forest_height;
        
        calcTopologyPrior((int) _lineages.size()+1);
        
#if defined (LAZY_COPYING)
        assert (new_nd->_left_child->_species == subtree1->_species);
        assert (new_nd->_left_child->_right_sib->_species == subtree2->_species);
        new_nd->_species = (new_nd->_left_child->_species | new_nd->_left_child->_right_sib->_species);
#endif
        
        return make_tuple(subtree1->_species, subtree2->_species, new_nd->_species);
    }
#else
    inline tuple<string,string, string> SpeciesForest::speciesTreeProposal(Lot::SharedPtr lot) {
        // this function creates a new node and joins two species
        
        bool done = false;
        Node* subtree1 = nullptr;
        Node* subtree2 = nullptr;
        
        while (!done) {
            assert (lot != nullptr);
            pair<unsigned, unsigned> t = lot->nchoose2((unsigned) _lineages.size());
            assert (t.first != t.second);
            
            subtree1=_lineages[t.first];
            subtree2=_lineages[t.second];
            assert (t.first < _lineages.size());
            assert (t.second < _lineages.size());
            assert(!subtree1->_parent && !subtree2->_parent);
            
            if (G::_outgroup != "none") {
                if (subtree1->_name != G::_outgroup && subtree2->_name != G::_outgroup && _lineages.size() > 2) { // outgroup can only be chosen on the last step
                    done = true;
                }
                else if (_lineages.size() == 2) {
                    done = true;
                }
            }
            else {
                done = true;
            }
            if (G::_outgroup == "none") {
                assert (done == true);
            }
        }
        
        Node * new_nd = &_nodes[G::_nspecies + _ninternals];
        assert (new_nd->_parent==0);
       assert (new_nd->_number == -1);
       assert (new_nd->_right_sib == 0);
//        _nodes.push_back(nd);
//        Node* new_nd = &_nodes.back();
//        new_nd->_parent=0;
        new_nd->_number=G::_nspecies+_ninternals;
        new_nd->_name+=boost::str(boost::format("node-%d")%new_nd->_number);
        new_nd->_edge_length=0.0;
        _ninternals++;
//        new_nd->_right_sib=0;

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd);

#if defined (DEBUG_MODE)
        if (_lineages.size() > 1) {
            _species_joined = make_pair(subtree1, subtree2); // last step just joins remaining two
        }
#endif
        
        new_nd->_height = _forest_height;
        
        calcTopologyPrior((int) _lineages.size()+1);
        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
    }
#endif

    inline tuple<string,string, string> SpeciesForest::speciesTreeProposalSim(Lot::SharedPtr lot) {
        // this function creates a new node and joins two species
        
        bool done = false;
        Node* subtree1 = nullptr;
        Node* subtree2 = nullptr;
        
        while (!done) {
            assert (lot != nullptr);
            pair<unsigned, unsigned> t = lot->nchoose2((unsigned) _lineages.size());
            assert (t.first != t.second);
            
            subtree1=_lineages[t.first];
            subtree2=_lineages[t.second];
            assert (t.first < _lineages.size());
            assert (t.second < _lineages.size());
            assert(!subtree1->_parent && !subtree2->_parent);
            
            if (G::_outgroup != "none") {
                if (subtree1->_name != G::_outgroup && subtree2->_name != G::_outgroup && _lineages.size() > 2) { // outgroup can only be chosen on the last step
                    done = true;
                }
                else if (_lineages.size() == 2) {
                    done = true;
                }
            }
            else {
                done = true;
            }
            if (G::_outgroup == "none") {
                assert (done == true);
            }
        }
        
        Node * new_nd = &_nodes[G::_nspecies + _ninternals];
        assert (new_nd->_parent==0);
       assert (new_nd->_number == -1);
       assert (new_nd->_right_sib == 0);
        new_nd->_number=G::_nspecies+_ninternals;
        new_nd->_name+=boost::str(boost::format("node-%d")%new_nd->_number);
        new_nd->_edge_length=0.0;
        _ninternals++;

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
        new_nd->_species = new_nd->_left_child->_species + new_nd->_left_child->_right_sib->_species;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd);
        
        new_nd->_height = _forest_height;
        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
    }

#if defined (DEBUG_MODE)
    inline void SpeciesForest::showSpeciesJoined() {
        assert (_index==0);
        if (_species_joined.first != NULL) {
            cout << "joining species " << _species_joined.first->_name << " and " << _species_joined.second->_name << endl;
        }
        else {
            cout << "no species joined" << endl;
        }
    }
#endif


    inline double SpeciesForest::calcTopologyPrior(unsigned nlineages) {
        _log_joining_prob += -log(0.5*nlineages*(nlineages-1));
        assert (!isinf(_log_joining_prob));
        return _log_joining_prob;
    }

    inline void SpeciesForest::debugForest() {
        cout << "debugging forest" << endl;
        for (auto &node : _nodes) {
            cout << "   node number " << node._number << " ";
            if (node._left_child) {
                cout << node._left_child->_number << " ";

                if (node._left_child->_right_sib) {
                    cout << node._left_child->_right_sib->_number << " ";
                }
            }
            else (cout << " - - ");

            if (node._partial!=nullptr) {
                cout << " * ";
            }
            cout << endl;
        }
        cout << "   _ninternals " << _ninternals << " ";
        cout << endl;
    }

    inline void SpeciesForest::updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode) {
        // Delete delnode1 from node_vector
        auto it1 = find(node_vector.begin(), node_vector.end(), delnode1);
        assert(it1 != node_vector.end());
        node_vector.erase(it1);

        // Delete delnode2 from node_vector
        auto it2 = find(node_vector.begin(), node_vector.end(), delnode2);
        assert(it2 != node_vector.end());
        node_vector.erase(it2);

        // Add addnode to node_vector
        node_vector.push_back(addnode);

        // reset _position_in_lineages
        for (int i=0; i < (int) _lineages.size(); i++) {
            _lineages[i] -> _position_in_lineages=i;
        }
    }

    inline void SpeciesForest::revertNodeVector(vector<Node *> &node_vector, Node *addnode1, Node *addnode2, Node *delnode1) {
        // Delete delnode1 from node_vector
        auto it = find(node_vector.begin(), node_vector.end(), delnode1);
        assert (it != node_vector.end());
        node_vector.erase(it);

        // find positions of nodes to insert
        auto position1 = addnode1->_position_in_lineages;
        auto iter1 = addnode1;

        auto position2 = addnode2->_position_in_lineages;
        auto iter2 = addnode2;

        // lower position must be inserted first
        if (position1 < position2) {
            node_vector.insert(node_vector.begin()+position1, iter1);
            node_vector.insert(node_vector.begin()+position2, iter2);
        }
        else {
            node_vector.insert(node_vector.begin()+position2, iter2);
            node_vector.insert(node_vector.begin()+position1, iter1);
        }

        assert(_lineages[addnode1->_position_in_lineages] == addnode1);
        assert(_lineages[addnode2->_position_in_lineages] == addnode2);

        // reset _position_in_lineages
        for (int i=0; i < (int) _lineages.size(); i++) {
            _lineages[i] -> _position_in_lineages=i;
        }
    }

    inline void SpeciesForest::updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode) {
        // Delete delnode1 from node_list
        auto it1 = find(node_list.begin(), node_list.end(), delnode1);
        assert(it1 != node_list.end());
        node_list.erase(it1);

        // Delete delnode2 from node_list
        auto it2 = find(node_list.begin(), node_list.end(), delnode2);
        assert(it2 != node_list.end());
        node_list.erase(it2);

        // Add addnode to node_list
        node_list.push_back(addnode);
    }

    inline void SpeciesForest::calcTreeLength() {
        // sum of all edge lengths in tree
        double sum_height = 0.0;
        
        for (auto &nd:_nodes) {
            // sum edge lengths from all nodes
            sum_height += nd._edge_length;
        }
        _forest_length = sum_height;
    }


    inline double SpeciesForest::getTreeLength() {
        // sum of all edge lengths in tree
        if (_forest_length == 0) { // check length is actually 0
            calcTreeLength();
        }
        return _forest_length;
    }

    inline void SpeciesForest::addIncrement(double increment) {
        for (auto &nd:_lineages) {
            nd->_edge_length += increment;
        }
        _last_edge_length = increment;
        
        _forest_height += _last_edge_length;
    }

    inline void SpeciesForest::refreshPreorder() {
        // this only works for complete trees - otherwise, need a separate preorder vector for each subtree - use refreshAllPreorders
       // Create vector of node pointers in preorder sequence
        Node *nd = &_nodes.back();
       _preorder.clear();
       _preorder.reserve(_nodes.size()); // _preorder must include root node

        _preorder.push_back(nd);
        
       while (true) {
           nd = findNextPreorder(nd);
           if (nd)
               _preorder.push_back(nd);
           else
               break;
       }   // end while loop
    }

    inline vector< pair<double, Node *>> SpeciesForest::sortPreorder() {
        vector< pair<double, Node *>> heights_and_nodes;
        for (auto it = _preorder.rbegin(); it != _preorder.rend(); it++) {
            Node * nd = *it;
            if (nd->_left_child) {
                // if internal node, store cumulative height in _height
                double height = nd->_height; //
                heights_and_nodes.push_back(make_pair(height, nd));
    //                nd->_height = nd->_left_child->_height + nd->_left_child->_edge_length;
    //                heights_and_nodes.push_back(make_pair(nd->_height, nd));
            }
            else {
                // if leaf node, initialize _height to zero
                nd->_height = 0.0;
            }
        }
         
        // sort heights_and_nodes so that smallest heights will be first
        sort(heights_and_nodes.begin(), heights_and_nodes.end());
        return(heights_and_nodes);
    }

    inline double SpeciesForest::getLineageHeight(Node* nd) {
        if (nd != nullptr) {
            double sum_height = 0.0;
            
            sum_height += nd->getEdgeLength();
            if (nd->_left_child) {
                for (Node* child = nd->_left_child; child; child=child->_left_child) {
                    sum_height += child->getEdgeLength();
                }
            }
            return sum_height;
        }
        else {
            return 0.0;
        }
    }

    inline void SpeciesForest::stripOutNexusComments(std::string & newick) {
        regex commentexpr("\\[.*?\\]");
        newick = std::regex_replace(newick, commentexpr, std::string(""));
    }

    inline unsigned SpeciesForest::countNewickInternals(const string newick) {
        size_t count = count_if(newick.begin(), newick.end(), []( char c ){return c =='(';});
        return count - 1;
    }

    inline unsigned SpeciesForest::countNewickLeaves(const std::string newick) {
        regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
        sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
        sregex_iterator m2;
        return (unsigned)std::distance(m1, m2);
    }

    inline bool SpeciesForest::canHaveSibling(Node * nd, bool rooted, bool allow_polytomies) {
        assert(nd);
        if (!nd->_parent) {
            // trying to give root node a sibling
            return false;
        }

        if (allow_polytomies)
            return true;

        bool nd_can_have_sibling = true;
        if (nd != nd->_parent->_left_child) {
            if (nd->_parent->_parent) {
                // trying to give a sibling to a sibling of nd, and nd's parent is not the root
                nd_can_have_sibling = false;
            }
            else {
                if (rooted) {
                    // root node has exactly 2 children in rooted trees
                    nd_can_have_sibling = false;
                }
                else if (nd != nd->_parent->_left_child->_right_sib) {
                    // trying to give root node more than 3 children
                    nd_can_have_sibling = false;
                }
            }
        }

        return nd_can_have_sibling;
    }

    inline void SpeciesForest::renumberInternals() {
        assert(_preorder.size() > 0);

        // Renumber internal nodes in postorder sequence
        unsigned curr_internal = G::_nspecies;
        for (auto & nd : boost::adaptors::reverse(_preorder)) {
            if (nd->_left_child) {
                // nd is an internal node
                nd->_number = curr_internal++;
            }
        }

        _ninternals = curr_internal - G::_nspecies;

        // If the tree has polytomies, then there are Node objects stored in
        // the _tree->_nodes vector that have not yet been numbered. These can
        // be identified because their _number is currently equal to -1.
        for (auto & nd : _nodes) {
            if (nd._number == -1)
                nd._number = curr_internal++;
        }
    }

    inline void SpeciesForest::extractEdgeLen(Node * nd, std::string edge_length_string) {
        assert(nd);
        bool success = true;
        double d = 0.0;
        try {
            d = std::stod(edge_length_string);
        }
        catch(std::invalid_argument &) {
            // edge_length_string could not be converted to a double value
            success = false;
        }

        if (success) {
            // conversion succeeded
            nd->setEdgeLength(d);
        }
        else
            throw XProj(boost::str(boost::format("%s is not interpretable as an edge length") % edge_length_string));
    }

    inline vector<tuple<string, string, string>> SpeciesForest::buildFromNewickTopology(string newick) {
        // assume tree is rooted
        // do not allow polytomies
        bool rooted = true;
        bool allow_polytomies = false;
        
        vector<tuple<string, string, string>> species_joined;
        species_joined.push_back(make_tuple("null", "null", "null"));
        
        set<unsigned> used; // used to ensure that no two leaf nodes have the same number
        unsigned curr_leaf = 0;
        unsigned num_edge_lengths = 0;
        unsigned curr_node_index = 0;

        // Remove comments from the supplied newick string
        string commentless_newick = newick;
        stripOutNexusComments(commentless_newick);

        // Resize the _nodes vector
//        _nleaves = countNewickLeaves(commentless_newick);
        
    //        if (G::_nspecies < 4) {
    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //        }
        unsigned max_nodes = 2*G::_nspecies - (rooted ? 0 : 2);
        _nodes.resize(max_nodes);
        for (auto & nd : _nodes ) {
            nd._name = "";
            nd._number = -1;
        }

        try {
            // Root node is the last node in _nodes
            auto l_front = _nodes.begin();
            std::advance(l_front, curr_node_index);
            Node *nd = &*l_front;

            if (rooted) {
                auto l_front = _nodes.begin();
                std::advance(l_front, ++curr_node_index);
                nd = &*l_front;

                auto parent = _nodes.begin();
                std::advance(parent, curr_node_index - 1);
                nd->_parent = &*parent;
                nd->_parent->_left_child = nd;
            }

            // Some flags to keep track of what we did last
            enum {
                Prev_Tok_LParen        = 0x01,    // previous token was a left parenthesis ('(')
                Prev_Tok_RParen        = 0x02,    // previous token was a right parenthesis (')')
                Prev_Tok_Colon        = 0x04,    // previous token was a colon (':')
                Prev_Tok_Comma        = 0x08,    // previous token was a comma (',')
                Prev_Tok_Name        = 0x10,    // previous token was a node name (e.g. '2', 'P._articulata')
                Prev_Tok_EdgeLen    = 0x20    // previous token was an edge length (e.g. '0.1', '1.7e-3')
            };
            unsigned previous = Prev_Tok_LParen;

            // Some useful flag combinations
            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

            // Set to true while reading an edge length
            bool inside_edge_length = false;
            std::string edge_length_str;
            unsigned edge_length_position = 0;

            // Set to true while reading a node name surrounded by (single) quotes
            bool inside_quoted_name = false;

            // Set to true while reading a node name not surrounded by (single) quotes
            bool inside_unquoted_name = false;

            // Set to start of each node name and used in case of error
            unsigned node_name_position = 0;

            // loop through the characters in newick, building up tree as we go
            unsigned position_in_string = 0;
            for (auto & ch : commentless_newick) {
                position_in_string++;

                if (inside_quoted_name) {
                    if (ch == '\'') {
                        inside_quoted_name = false;
                        node_name_position = 0;
                        if (!nd->_left_child) {
                            curr_leaf++;
                        }
                        previous = Prev_Tok_Name;
                    }
                    else if (iswspace(ch))
                        nd->_name += ' ';
                    else {
                        nd->_name += ch;
                    }

                    continue;
                }
                else if (inside_unquoted_name) {
                    if (ch == '(')
                        throw XProj(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));

                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                        inside_unquoted_name = false;

                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));

                        if (!nd->_left_child) {
                            curr_leaf++;
                        }

                        previous = Prev_Tok_Name;
                    }
                    else {
                        nd->_name += ch;
                        continue;
                    }
                }
                else if (inside_edge_length) {
                    if (ch == ',' || ch == ')' || iswspace(ch)) {
                        inside_edge_length = false;
                        edge_length_position = 0;
                        extractEdgeLen(nd, edge_length_str);
                        ++num_edge_lengths;
                        previous = Prev_Tok_EdgeLen;
                    }
                    else {
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid)
                            throw XProj(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
                        edge_length_str += ch;
                        continue;
                    }
                }

                if (iswspace(ch))
                    continue;

                switch(ch) {
                    case ';':
                        break;

                    case ')':
                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
                        if (!nd->_parent)
                            throw XProj(boost::str(boost::format("Too many right parentheses at position %d in tree description") % position_in_string));

                        // Expect right paren only after an edge length, a node name, or another right paren
                        if (!(previous & RParen_Valid))
                            throw XProj(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));

                        // Go down a level
                        nd = nd->_parent;
                        if (!nd->_left_child->_right_sib) {
                            cout << "newick is " << newick << endl;
                            throw XProj(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
                        }
                        previous = Prev_Tok_RParen;
                        break;

                    case ':':
                        // Expect colon only after a node name or another right paren
                        if (!(previous & Colon_Valid))
                            throw XProj(boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_Colon;
                        break;

                    case ',':
                    {
                        // Expect comma only after an edge length, a node name, or a right paren
                        if (!nd->_parent || !(previous & Comma_Valid))
                            throw XProj(boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));

                        // Check for polytomies
                        if (!canHaveSibling(nd, rooted, allow_polytomies)) {
                            throw XProj(boost::str(boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % newick));
                        }

                        // Create the sibling
                        curr_node_index++;
                        if (curr_node_index == _nodes.size())
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % G::_nspecies));

                        auto l_front = _nodes.begin();
                        std::advance(l_front, curr_node_index);
                        nd->_right_sib = &*l_front;

                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;
                        break;
                    }

                    case '(':
                    {
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw XProj(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % position_in_string));

                        // Create new node above and to the left of the current node
                        assert(!nd->_left_child);
                        curr_node_index++;
                        if (curr_node_index == _nodes.size())
                            throw XProj(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _nodes.size()));

                        auto l_front = _nodes.begin();
                        std::advance(l_front, curr_node_index);
                        nd->_left_child = &*l_front;

                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;
                        break;
                    }

                    case '\'':
                        // Encountered an apostrophe, which always indicates the start of a
                        // node name (but note that node names do not have to be quoted)

                        // Expect node name only after a left paren (child's name), a comma (sib's name)
                        // or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(boost::str(boost::format("Not expecting node name at position %d in tree description") % position_in_string));

                        // Get the rest of the name
                        nd->_name.clear();

                        inside_quoted_name = true;
                        node_name_position = position_in_string;

                        break;

                    default:
                        // Get here if ch is not one of ();:,'

                        // Expecting either an edge length or an unquoted node name
                        if (previous == Prev_Tok_Colon) {
                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                            inside_edge_length = true;
                            edge_length_position = position_in_string;
                            edge_length_str = ch;
                        }
                        else {
                            // Get the node name
                            nd->_name = ch;

                            inside_unquoted_name = true;
                            node_name_position = position_in_string;
                        }
                }   // end of switch statement
            }   // loop over characters in newick string

            if (inside_unquoted_name)
                throw XProj(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
            if (inside_edge_length)
                throw XProj(boost::str(boost::format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
            if (inside_quoted_name)
                throw XProj(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));

            if (rooted) {
                refreshPreorder();
            }
            renumberInternals();
        }
        catch(XProj &x) {
    //            if (_index == 0) {
    //                clearSpeciesForest();
    //            }
    //            else {
    //                clearGeneForest();
    //            }
    ////            clear();
            throw x;
        }

        // TODO: fix for _nodes as vector
        assert (1==2);
//        _nodes.pop_front(); // remove node at beginning of list because it's an extra root
//        // remove parent from new last node
//        _nodes.front()._parent = NULL;
//
//        _nodes.sort(
//             [this](Node& lhs, Node& rhs) {
//    //                 return lhs._left_child->_height < rhs._left_child->_height; } ); // TODO: is this just lhs->_height and rhs->_height?
//                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );

        // reset node numbers and names that are not tips
        int j = 0;
        for (auto &nd:_nodes) {
            nd._number = j;
    //            nd._edge_length = 0.0;
            if (nd._name == "") {
                nd._name=boost::str(boost::format("node-%d")%nd._number);
            }
            j++;
        }
        
        refreshPreorder();
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        for (auto &entry:heights_and_nodes) {
            species_joined.push_back(make_tuple(entry.second->_left_child->_name, entry.second->_left_child->_right_sib->_name, entry.second->_name));

        }
        
        unsigned count = 0;
        for (auto &nd:_nodes) {
            nd._edge_length = 0.0;
            if (count < G::_nspecies) {
                nd._name = "";
            }
            count++;
        }
        _preorder.clear();

        return species_joined;
    }

    inline void SpeciesForest::buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies) {
        list<Node> test_nodes;
        _nodes.clear();
        
        set<unsigned> used; // used to ensure that no two leaf nodes have the same number
        unsigned curr_leaf = 0;
        unsigned num_edge_lengths = 0;
        unsigned curr_node_index = 0;

        // Remove comments from the supplied newick string
        string commentless_newick = newick;
        stripOutNexusComments(commentless_newick);

        // Resize the _nodes vector
//        _nleaves = countNewickLeaves(commentless_newick);
    //        if (G::_nspecies < 4) {
    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //        }
        unsigned max_nodes = 2*G::_nspecies - (rooted ? 0 : 2);
        test_nodes.resize(max_nodes); // TODO: TEST
        
    //        int b=0;
        for (auto & nd : test_nodes ) {// TODO: TEST
            nd._name = "";
            nd._number = -1;
        }

        try {
            // Root node is the last node in _nodes
            auto l_front = test_nodes.begin(); // TODO: TEST
            std::advance(l_front, curr_node_index); // TODO: curr_node_index should be 0
            Node *nd = &*l_front;
            
//            Node *nd = &_nodes[curr_node_index];

            if (rooted) {
                auto l_front = test_nodes.begin(); // TODO: TEST
                std::advance(l_front, ++curr_node_index);
                nd = &*l_front;

                auto parent = test_nodes.begin(); // TODO: TEST
                std::advance(parent, curr_node_index - 1);
                nd->_parent = &*parent;
                nd->_parent->_left_child = nd;
            }

            // Some flags to keep track of what we did last
            enum {
                Prev_Tok_LParen        = 0x01,    // previous token was a left parenthesis ('(')
                Prev_Tok_RParen        = 0x02,    // previous token was a right parenthesis (')')
                Prev_Tok_Colon        = 0x04,    // previous token was a colon (':')
                Prev_Tok_Comma        = 0x08,    // previous token was a comma (',')
                Prev_Tok_Name        = 0x10,    // previous token was a node name (e.g. '2', 'P._articulata')
                Prev_Tok_EdgeLen    = 0x20    // previous token was an edge length (e.g. '0.1', '1.7e-3')
            };
            unsigned previous = Prev_Tok_LParen;

            // Some useful flag combinations
            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

            // Set to true while reading an edge length
            bool inside_edge_length = false;
            std::string edge_length_str;
            unsigned edge_length_position = 0;

            // Set to true while reading a node name surrounded by (single) quotes
            bool inside_quoted_name = false;

            // Set to true while reading a node name not surrounded by (single) quotes
            bool inside_unquoted_name = false;

            // Set to start of each node name and used in case of error
            unsigned node_name_position = 0;

            // loop through the characters in newick, building up tree as we go
            unsigned position_in_string = 0;
            for (auto ch : commentless_newick) {
                position_in_string++;

                if (inside_quoted_name) {
                    if (ch == '\'') {
                        inside_quoted_name = false;
                        node_name_position = 0;
                        if (!nd->_left_child) {
                            curr_leaf++;
                        }
                        previous = Prev_Tok_Name;
                    }
                    else if (iswspace(ch))
                        nd->_name += ' ';
                    else {
                        nd->_name += ch;
                    }

                    continue;
                }
                else if (inside_unquoted_name) {
                    if (ch == '(')
                        throw XProj(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));

                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                        inside_unquoted_name = false;

                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));

                        if (!nd->_left_child) {
                            curr_leaf++;
                        }

                        previous = Prev_Tok_Name;
                    }
                    else {
                        nd->_name += ch;
                        continue;
                    }
                }
                else if (inside_edge_length) {
                    if (ch == ',' || ch == ')' || iswspace(ch)) {
                        inside_edge_length = false;
                        edge_length_position = 0;
                        extractEdgeLen(nd, edge_length_str);
                        ++num_edge_lengths;
                        previous = Prev_Tok_EdgeLen;
                    }
                    else {
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid)
                            throw XProj(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
                        edge_length_str += ch;
                        continue;
                    }
                }

                if (iswspace(ch))
                    continue;

                switch(ch) {
                    case ';':
                        break;

                    case ')':
                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
                        if (!nd->_parent)
                            throw XProj(boost::str(boost::format("Too many right parentheses at position %d in tree description") % position_in_string));

                        // Expect right paren only after an edge length, a node name, or another right paren
                        if (!(previous & RParen_Valid))
                            throw XProj(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));

                        // Go down a level
                        nd = nd->_parent;
                        if (!nd->_left_child->_right_sib) {
                            cout << "newick is " << newick << endl;
                            throw XProj(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
                        }
                        previous = Prev_Tok_RParen;
                        break;

                    case ':':
                        // Expect colon only after a node name or another right paren
                        if (!(previous & Colon_Valid))
                            throw XProj(boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_Colon;
                        break;

                    case ',':
                    {
                        // Expect comma only after an edge length, a node name, or a right paren
                        if (!nd->_parent || !(previous & Comma_Valid))
                            throw XProj(boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));

                        // Check for polytomies
                        if (!canHaveSibling(nd, rooted, allow_polytomies)) {
                            throw XProj(boost::str(boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % newick));
                        }

                        // Create the sibling
                        curr_node_index++;
                        if (curr_node_index == _nodes.size())
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % G::_nspecies));

                        auto l_front = test_nodes.begin(); // TODO: TEST
                        std::advance(l_front, curr_node_index);
                        nd->_right_sib = &*l_front;

                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;
                        break;
                    }

                    case '(':
                    {
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw XProj(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % position_in_string));

                        // Create new node above and to the left of the current node
                        assert(!nd->_left_child);
                        curr_node_index++;
                        if (curr_node_index == _nodes.size())
                            throw XProj(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _nodes.size()));

                        auto l_front = test_nodes.begin(); // TODO: TEST
                        std::advance(l_front, curr_node_index);
                        nd->_left_child = &*l_front;

                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;
                        break;
                    }

                    case '\'':
                        // Encountered an apostrophe, which always indicates the start of a
                        // node name (but note that node names do not have to be quoted)

                        // Expect node name only after a left paren (child's name), a comma (sib's name)
                        // or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(boost::str(boost::format("Not expecting node name at position %d in tree description") % position_in_string));

                        // Get the rest of the name
                        nd->_name.clear();

                        inside_quoted_name = true;
                        node_name_position = position_in_string;

                        break;

                    default:
                        // Get here if ch is not one of ();:,'

                        // Expecting either an edge length or an unquoted node name
                        if (previous == Prev_Tok_Colon) {
                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                            inside_edge_length = true;
                            edge_length_position = position_in_string;
                            edge_length_str = ch;
                        }
                        else {
                            // Get the node name
                            nd->_name = ch;

                            inside_unquoted_name = true;
                            node_name_position = position_in_string;
                        }
                }   // end of switch statement
            }   // loop over characters in newick string

            if (inside_unquoted_name) {
                cout << "entering exception " << endl;
                throw XProj(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
            }
            if (inside_edge_length)
                throw XProj(boost::str(boost::format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
            if (inside_quoted_name)
                throw XProj(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));

            if (rooted) {
                for (auto &nd:test_nodes) {
                    _nodes.push_back(nd);
                }
                refreshPreorder(); // TODO: clearing for now
            }
            renumberInternals();
        }
        catch(XProj &x) {
            clear();
            throw x;
        }
        
        test_nodes.pop_front(); // remove node at beginning of list because it's an extra root
//        // remove parent from new last node
        test_nodes.front()._parent = NULL;
//
        test_nodes.sort(
             [this](Node& lhs, Node& rhs) {
                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );
        
        // first G::_nspecies nodes will be tips because the height will be 0
        // sort these nodes alphabetically so bit settings will be consistent
        
        // Extract the first G::_nspecies elements into a vector
//        list<Node> tip_nodes;
//        for (unsigned i = 0; i < G::_nspecies; ++i) {
//            auto it = test_nodes.begin();
//            std::advance(it, i);
//            tip_nodes.push_back(*it);
//        }
        
//        for (auto &nd:tip_nodes) {
//            if (nd._parent) {
//                Node parent = *nd._parent;
//                for (auto &nd_t:test_nodes) {
//                    if (nd_t._edge_length == parent._edge_length) {
//                        nd._parent = &nd_t;
//                    }
//                }
//            }
//            if (nd._right_sib) {
//                Node sib = *nd._right_sib;
//                for (auto &nd_t:test_nodes) {
//                    if (nd_t._edge_length == sib._edge_length) {
//                        nd._right_sib = &nd_t;
//                    }
//                }
//            }
//        }
        
//        tip_nodes.sort(
//             [this](Node& lhs, Node& rhs) {
//                 return lhs._name < rhs._name; } );
        
//        for (unsigned i=0; i<G::_nspecies; i++) {
//            auto it1 = test_nodes.begin();
//            auto it2 = tip_nodes.begin();
//            std::advance(it1, i);
//            std::advance(it2, i);
//            *it1 = *it2;
//        }
        
// //        test_nodes.splice(test_nodes.begin(), tip_nodes);
        
        _nodes.clear();
        
        _nodes.assign(test_nodes.begin(), test_nodes.end());
        
        refreshPreorder();
        
        _lineages.clear();
        
        _lineages.push_back(&_nodes.back());
        
        // reset node names
        int j = 0;
        for (auto &nd:_nodes) {
            nd._number = j;
            j++;
        }
        
        
        for (auto &nd:_nodes) {
            if (nd._name == "") {
                nd._name=boost::str(boost::format("node-%d")%nd._number);
            }
        }
        
        if (_lineages.size() == 1) {
            _nodes.back()._edge_length = 0.0;
        }
                
        // TODO: make sure tips are at front of list and have correct number associated with name
        unsigned number = 0;
        for (auto &nd:test_nodes) {
            nd._number = number;
            number++;
        }
        
        // reset pointers so they point to new vector, not old list
        for (auto &nd:_nodes) {
            if (nd._left_child) {
                nd._left_child = &_nodes[nd._left_child->_number];
            }
            if (nd._right_sib) {
                nd._right_sib = &_nodes[nd._right_sib->_number];
            }
            if (nd._parent) {
                nd._parent = &_nodes[nd._parent->_number];
            }
        }
        
//        test_nodes.clear();
    }

    inline void SpeciesForest::extractNodeNumberFromName(Node * nd, std::set<unsigned> & used) {
        assert(nd);
        bool success = true;
        unsigned x = 0;
        try {
            x = std::stoi(nd->_name);
        }
        catch(std::invalid_argument &) {
            // node name could not be converted to an integer value
            success = false;
        }

        if (success) {
            // conversion succeeded
            // attempt to insert x into the set of node numbers already used
            std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
            if (insert_result.second) {
                // insertion was made, so x has NOT already been used
                nd->_number = x - 1;
            }
            else {
                // insertion was not made, so set already contained x
                throw XProj(boost::str(boost::format("leaf number %d used more than once") % x));
            }
        }
        else
            throw XProj(boost::str(boost::format("node name (%s) not interpretable as a positive integer") % nd->_name));
    }

    inline vector<pair<tuple<string, string, string>, double>> SpeciesForest::buildFromNewickMPI(const std::string newick, bool rooted, bool allow_polytomies, Lot::SharedPtr lot) {
        vector<pair<tuple<string, string, string>, double>> new_t;
        // TODO: this doesn't seem to work if the tree is complete? leaves out the last node?
        _nodes.clear();
    //        cout << "in function buildfromnewickmpi" << endl;
    //        cout << "building newick " << newick << endl;
        std::set<unsigned> used; // used to ensure that no two leaf nodes have the same number
        unsigned curr_leaf = 0;
        unsigned num_edge_lengths = 0;
        unsigned curr_node_index = 0;

        // Remove comments from the supplied newick string
        std::string commentless_newick = newick;
        stripOutNexusComments(commentless_newick);

        // Resize the _nodes vector
        G::_nspecies = countNewickLeaves(commentless_newick);
        unsigned max_nodes = 2*G::_nspecies - (rooted ? 0 : 2);
        _nodes.resize(max_nodes);
        for (auto & nd : _nodes ) {
            nd._number = -1;
            nd._name = "";
        }

        try {
            // Root node
            auto it = _nodes.begin();
            std::advance(it, curr_node_index);
            Node * nd = &*it;

            if (rooted) {
                auto l_front = _nodes.begin();
                std::advance(l_front, ++curr_node_index);
                nd = &*l_front;

                auto parent = _nodes.begin();
                std::advance(parent, curr_node_index - 1);
                nd->_parent = &*parent;
                nd->_parent->_left_child = nd;
                
    //                std::list<Node>::iterator it1 = _nodes.begin();
    //                std::advance(it1, ++curr_node_index);
    //                nd = &*it1;
    //
    //                std::list<Node>::iterator it2 = _nodes.begin();
    //                std::advance(it2, curr_node_index - 1);
    //                Node parent = *it2;
    //
    //                nd->_parent = &parent;
    //                nd->_parent->_left_child = nd;
            }

            // Some flags to keep track of what we did last
            enum {
                Prev_Tok_LParen        = 0x01,    // previous token was a left parenthesis ('(')
                Prev_Tok_RParen        = 0x02,    // previous token was a right parenthesis (')')
                Prev_Tok_Colon        = 0x04,    // previous token was a colon (':')
                Prev_Tok_Comma        = 0x08,    // previous token was a comma (',')
                Prev_Tok_Name        = 0x10,    // previous token was a node name (e.g. '2', 'P._articulata')
                Prev_Tok_EdgeLen    = 0x20    // previous token was an edge length (e.g. '0.1', '1.7e-3')
            };
            unsigned previous = Prev_Tok_LParen;

            // Some useful flag combinations
            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

            // Set to true while reading an edge length
            bool inside_edge_length = false;
            std::string edge_length_str;
            unsigned edge_length_position = 0;

            // Set to true while reading a node name surrounded by (single) quotes
            bool inside_quoted_name = false;

            // Set to true while reading a node name not surrounded by (single) quotes
            bool inside_unquoted_name = false;

            // Set to start of each node name and used in case of error
            unsigned node_name_position = 0;

            // loop through the characters in newick, building up tree as we go
            unsigned position_in_string = 0;
            for (auto & ch : commentless_newick) {
                position_in_string++;

                if (inside_quoted_name) {
                    if (ch == '\'') {
                        inside_quoted_name = false;
                        node_name_position = 0;
                        if (!nd->_left_child) {
    //                            extractNodeNumberFromName(nd, used);
                            curr_leaf++;
                        }
                        previous = Prev_Tok_Name;
                    }
                    else if (iswspace(ch))
                        nd->_name += ' ';
                    else
                        nd->_name += ch;

                    continue;
                }
                else if (inside_unquoted_name) {
                    if (ch == '(')
                        throw XProj(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));

                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                        inside_unquoted_name = false;

                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));

                        if (!nd->_left_child) {
    //                            extractNodeNumberFromName(&nd, used);
                            curr_leaf++;
                        }

                        previous = Prev_Tok_Name;
                    }
                    else {
                        nd->_name += ch;
                        continue;
                    }
                }
                else if (inside_edge_length) {
                    if (ch == ',' || ch == ')' || iswspace(ch)) {
                        inside_edge_length = false;
                        edge_length_position = 0;
                        extractEdgeLen(nd, edge_length_str);
                        ++num_edge_lengths;
                        previous = Prev_Tok_EdgeLen;
                    }
                    else {
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid)
                            throw XProj(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
                        edge_length_str += ch;
                        continue;
                    }
                }

                if (iswspace(ch))
                    continue;

                switch(ch) {
                    case ';':
                        break;

                    case ')': {
                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
                        if (!nd->_parent)
                            throw XProj(boost::str(boost::format("Too many right parentheses at position %d in tree description") % position_in_string));

                        // Expect right paren only after an edge length, a node name, or another right paren
                        if (!(previous & RParen_Valid))
                            throw XProj(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));

                        // Go down a level
                        nd = nd->_parent;
                        if (!nd->_left_child->_right_sib) {
                            cout << "newick is " << newick << endl;
                            throw XProj(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
                        }
                        previous = Prev_Tok_RParen;
                        break;
                    }

                    case ':': {
                        // Expect colon only after a node name or another right paren
                        if (!(previous & Colon_Valid))
                            throw XProj(boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_Colon;
                        break;
                    }

                    case ',': { // TODO: comma doesn't necessarily mean siblings if there's a polytomy
                        // Expect comma only after an edge length, a node name, or a right paren
                        if (!nd->_parent || !(previous & Comma_Valid))
                            throw XProj(boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));

                        // Check for polytomies
                        if (!canHaveSibling(nd, rooted, allow_polytomies)) {
                            throw XProj(boost::str(boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % newick));
                        }

                        // Create the sibling
                        curr_node_index++;
                        if (curr_node_index == _nodes.size())
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % G::_nspecies));
                        auto l_front = _nodes.begin();
                        std::advance(l_front, curr_node_index);
                        nd->_right_sib = &*l_front;

                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;
                        break;
                        
    //                        std::list<Node>::iterator it_sib = _nodes.begin();
    //                        std::advance(it_sib, curr_node_index);
    //                        Node * nd_right_sib = &*it_sib;
    //                        nd->_right_sib = nd_right_sib;
    //
    ////                        nd._right_sib = &_nodes[curr_node_index];
    //                        nd->_right_sib->_parent = nd->_parent;
    //                        nd = nd->_right_sib;
    //                        previous = Prev_Tok_Comma;
    //                        break;
                    }

                    case '(': {
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw XProj(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % position_in_string));

                        // Create new node above and to the left of the current node
                        assert(!nd->_left_child);
                        curr_node_index++;
                        if (curr_node_index == _nodes.size())
                            throw XProj(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _nodes.size()));
                        auto l_front = _nodes.begin();
                        std::advance(l_front, curr_node_index);
                        nd->_left_child = &*l_front;

                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;
                        break;
                        
    //                        std::list<Node>::iterator it_child = _nodes.begin();
    //                        std::advance(it_child, curr_node_index);
    //                        Node * nd_left_child = &*it_child;
    //                        nd->_left_child = nd_left_child;
    //
    ////                        nd._left_child = &_nodes[curr_node_index];
    //                        nd->_left_child->_parent = nd;
    //                        nd = nd->_left_child;
    //                        previous = Prev_Tok_LParen;
    //                        break;
                    }

                    case '\'': {
                        // Encountered an apostrophe, which always indicates the start of a
                        // node name (but note that node names do not have to be quoted)

                        // Expect node name only after a left paren (child's name), a comma (sib's name)
                        // or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XProj(boost::str(boost::format("Not expecting node name at position %d in tree description") % position_in_string));

                        // Get the rest of the name
                        nd->_name.clear();

                        inside_quoted_name = true;
                        node_name_position = position_in_string;

                        break;
                    }

                    default: {
                        // Get here if ch is not one of ();:,'

                        // Expecting either an edge length or an unquoted node name
                        if (previous == Prev_Tok_Colon) {
                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                            inside_edge_length = true;
                            edge_length_position = position_in_string;
                            edge_length_str = ch;
                        }
                        else {
                            // Get the node name
                            nd->_name = ch;

                            inside_unquoted_name = true;
                            node_name_position = position_in_string;
                        }
                    }
                }   // end of switch statement
            }   // loop over characters in newick string

            if (inside_unquoted_name)
                throw XProj(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
            if (inside_edge_length)
                throw XProj(boost::str(boost::format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
            if (inside_quoted_name)
                throw XProj(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));

    //            unsigned max_nodes = countNewickInternals(newick) + G::_nspecies;
            unsigned ninternals = countNewickInternals(newick);
            unsigned max_internals = G::_nspecies-1;
            unsigned max_nodes = ninternals + G::_nspecies + 1;
            // TODO: fix for _nodes as vector
            assert (1 == 2);
//            if (ninternals != max_internals) {
//                _nodes.pop_front(); // if the tree is incomplete, delete both the root and subroot
//                _nodes.front()._parent = nullptr;
//                max_nodes--;
//            }
            
            _nodes.front()._name = "unused"; // break off root node since we are not using it
            
            for (auto &nd:_nodes) {
                if (nd._parent) {
                    if (nd._parent->_name == "unused") {
                        nd._parent = nullptr;
                        nd._right_sib = nullptr;
                    }
                }
            }
            // TODO: fix for _nodes as vector
            assert ( 1 == 2);
//            _nodes.pop_front();
//
//            unsigned nodes_to_delete = (unsigned)_nodes.size() - max_nodes;
//            for (unsigned i=0; i<nodes_to_delete; i++) {
//                _nodes.pop_back();
//            }
//            assert (_nodes.size() == max_nodes);
//
//            // sort nodes by height
//            _nodes.sort(
//                 [this](Node& lhs, Node& rhs) {
//                     return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );
                        
            // reset node numbers and names
            int j = 0;
            for (auto &nd:_nodes) {
                nd._number = j;
                if (nd._number > G::_ntaxa-1) {
                    nd._name = to_string(j);
                }
                j++;
            }
            
    //            for (auto &nd:_nodes) {
    //                cout << "name = " << nd._name << endl;
    //                cout << "number = " << nd._number << endl;
    //            }
            
            unsigned num = G::_nspecies;
            
            for (auto &nd:_nodes) {
                if (nd._name == "") {
                    nd._name=boost::str(boost::format("node-%d")%num);
                    num++;
//                        cout << "SETTING NODE NAME TO " << nd._name << endl;
                }
            }
            
           new_t = resetLineages(lot);
            
            if (_lineages.size() == 1) {
                _nodes.back()._edge_length = 0.0;
            }
            
    //            for (auto &nd:_nodes) {
    //                if (!nd._parent) {
    //                    _lineages.push_back(&nd);
    //                }
    //            }
    //
    //            unsigned count = 0;
    //            for (auto &nd:_lineages) {
    //                nd->_position_in_lineages = count;
    //                count++;
    //            }
            
            if (rooted) {
                refreshPreorder();
            }
            
            _forest_height = getLineageHeight(_lineages.back());
            for (auto &nd:_nodes) {
                if (nd._left_child) {
                    nd._height = getLineageHeight(nd._left_child);
                }
                else {
                    nd._height = 0.0;
                }
            }
    //            else {
    //                // Root at leaf whose _number = 0
    //                // refreshPreorder() and refreshLevelorder() called after rerooting
    //                rerootAtNodeNumber(0);
    //            }
    //            renumberInternals();
        }
        catch(XProj x) {
            clear();
            throw x;
        }
        return new_t;
    }

    inline vector<pair<tuple<string, string, string>, double>> SpeciesForest::resetLineages(Lot::SharedPtr lot) {
        // this function rebuilds the _lineages vector, setting _position_in_lineages
        // this function also rebuilds _t for use in reading in a species newick
        vector<pair<tuple<string, string, string>, double>> new_t;
        
        _ninternals = 0;
        _lineages.clear();
        unsigned i=0;
        for (auto &nd:_nodes) {
            if (!nd._left_child) {
                _lineages.push_back(&nd);
                nd._position_in_lineages = i;
                i++;
            }
            else {
                _ninternals++;
            }
        }
        
        bool draw_edge_lens_from_prior = false;
        if (G::_species_newick_name != "null") {
            draw_edge_lens_from_prior = true;
        }
        
        if (!draw_edge_lens_from_prior) {
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            double edge_len = _lineages.front()->_edge_length;
            new_t.push_back(make_pair(species_joined, edge_len));

            
            for (auto &nd : _nodes) {
                if (nd._left_child) {
                    assert (nd._left_child->_right_sib);
                    updateNodeVector(_lineages, nd._left_child, nd._left_child->_right_sib, &nd);
                    tuple<string, string, string> species_joined = make_tuple(nd._left_child->_name, nd._left_child->_right_sib->_name, nd._name);
                    double edge_len = nd._edge_length;
                    new_t.push_back(make_pair(species_joined, edge_len));
                }
            }
            
            // if there are 2 lineages left and neither has edge length 0, the tree is complete - add another node
            if (_lineages.size() == 2) {
                if (_lineages[0]->_edge_length > G::_small_enough && _lineages[1]->_edge_length > G::_small_enough) {
                    Node* subtree1 = _lineages[0];
                    Node* subtree2 = _lineages[1];
                    
//                    Node nd;
//                    _nodes.push_back(nd);
//                    Node* new_nd = &_nodes.back();
                    assert (1 == 2); // TODO: fix for _nodes as vector
                    Node * new_nd = &_nodes[G::_nspecies + _ninternals];
                    
                    new_nd->_parent=0;
                    new_nd->_number=G::_nspecies+_ninternals;
                    new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);
                    new_nd->_edge_length=0.0;
                    _ninternals++;
                    new_nd->_right_sib=0;

                    new_nd->_left_child=subtree1;
                    subtree1->_right_sib=subtree2;

                    subtree1->_parent=new_nd;
                    subtree2->_parent=new_nd;
                    
                    updateNodeVector (_lineages, subtree1, subtree2, new_nd);
                    assert (_lineages.size() == 1);
                    
                    tuple<string, string, string> species_joined = make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
                    double edge_len = new_nd->_edge_length;
                    new_t.push_back(make_pair(species_joined, edge_len));
                }
            }
        }
        else {
            _forest_height = 0.0;
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            for (auto &nd:_lineages) {
                nd->_edge_length = 0.0;
            }
            chooseSpeciesIncrementOnly(lot, 0.0);
    //            double edge_len = _lineages.front()->_edge_length
            new_t.push_back(make_pair(species_joined, _last_edge_length));

            
            for (auto &nd : _nodes) {
                if (nd._left_child) {
                    assert (nd._left_child->_right_sib);
                    
                    updateNodeVector(_lineages, nd._left_child, nd._left_child->_right_sib, &nd);
                    nd._edge_length = 0.0;
                    chooseSpeciesIncrementOnly(lot, 0.0);
                    
                    tuple<string, string, string> species_joined = make_tuple(nd._left_child->_name, nd._left_child->_right_sib->_name, nd._name);
                    double edge_len = nd._edge_length;
                    assert (edge_len == _last_edge_length);
                    new_t.push_back(make_pair(species_joined, _last_edge_length));
                }
            }
            
            // if there are 2 lineages left and neither has edge length 0, the tree is complete - add another node
            if (_lineages.size() == 2) {
                if ((_lineages[0]->_edge_length > G::_small_enough && _lineages[1]->_edge_length > G::_small_enough) || (G::_species_newick_name != "null")) {
                    Node* subtree1 = _lineages[0];
                    Node* subtree2 = _lineages[1];
                    
//                    Node nd;
//                    _nodes.push_back(nd);
//                    Node* new_nd = &_nodes.back();
                    assert (1 == 2); // TODO: fix for _nodes as vector
                    Node * new_nd = &_nodes[G::_nspecies + _ninternals];
                    new_nd->_parent=0;
                    new_nd->_number=G::_nspecies+_ninternals;
                    new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);
                    new_nd->_edge_length=0.0;
                    _ninternals++;
                    new_nd->_right_sib=0;

                    new_nd->_left_child=subtree1;
                    subtree1->_right_sib=subtree2;

                    subtree1->_parent=new_nd;
                    subtree2->_parent=new_nd;
                    
                    updateNodeVector (_lineages, subtree1, subtree2, new_nd);
                    
                    assert (_lineages.size() == 1);
                    
                    tuple<string, string, string> species_joined = make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
                    double edge_len = new_nd->_edge_length;
                    assert (edge_len == 0.0);
                    new_t.push_back(make_pair(species_joined, 0.0));
                }
            }
        }
        return new_t;
    }

    inline void SpeciesForest::addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient) {
        unsigned index = 0;
        // Assumes nd is an internal node
        assert(nd->_left_child);
        
        recipient.push_back(
            make_tuple(
                nd->_height,
                index,
                vector<G::species_t>({
                    nd->_left_child->_species,
                    nd->_left_child->_right_sib->_species,
                })
            )
        );
    }

    inline void SpeciesForest::buildCoalInfoVect() {
        
        refreshAllPreorders(); // TODO: don't always need to refresh preorders?
        
        _coalinfo.clear(); // TODO: why clear this and reconstruct each time?
        for (auto & preorder : _preorders) {
            for (auto & nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != G::_infinity);
                    assert(nd->_left_child->_right_sib);
                    assert(nd->_left_child->_right_sib->_right_sib == nullptr);
                    nd->_species = (nd->_left_child->_species | nd->_left_child->_right_sib->_species);
                    addCoalInfoElem(nd, _coalinfo);
                }
                else {
                    // nd is a leaf node
                    unsigned spp_index = G::_taxon_to_species.at(nd->_name);
                    nd->_species = (G::species_t)1 << spp_index;
                }
            }
        }
    }

    inline void SpeciesForest::saveCoalInfoInitial() {
        // coalinfo_t is a tuple with these elements:
        // - height of node
        // - 1-offset gene index (0 means speciation)
        // - vector of child species
        
        // Should only be called for complete gene trees
        assert (_lineages.size() == 1);
        
        _coalinfo.clear();
        
        refreshPreorder();
        
        // TODO: check node height - don't need to calculate, can just keep track as you go
            for (auto & nd : boost::adaptors::reverse(_preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
    //                    assert(nd->_height != _infinity);
                    assert(nd->_left_child->_right_sib);
                    assert(nd->_left_child->_right_sib->_right_sib == nullptr);
                    nd->_species = (nd->_left_child->_species | nd->_left_child->_right_sib->_species);
                    addCoalInfoElem(nd, _coalinfo);
                }
                else {
                    // nd is a leaf node
                    unsigned spp_index = G::_taxon_to_species.at(nd->_name);
                    nd->_species = (G::species_t)1 << spp_index;
                }
            }
        }

    inline bool SpeciesForest::subsumed(G::species_t test_species, G::species_t subtending_species) {
        bool not_equal = (test_species != subtending_species);
        bool is_subset = ((test_species & subtending_species) == test_species);
        if (not_equal && is_subset) {
            return true;
        }
        else {
            return false;
        }
    }

    inline void SpeciesForest::fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const {
        // No fixing up to do if there are no species tree joins
        if (sppinfo_vect.empty())
            return;
                    
    //        debugCheckCoalInfoSorted(coalinfo_vect);
    //        debugCheckCoalInfoSorted(sppinfo_vect);
        
        // Example:
        //          -- before --           -- after ---
        //  height   left  right  species    left   right
        // -------  -----  -----  -------  -----  -----
        // 0.01167      E      E              E       E
        // 0.01389      B      B              B       B
        // 0.02047      B      B              B       B
        // 0.02079      D      D              D       D
        // 0.02150      D      D              D       D
        // 0.03438      C      C              C       C
        // 0.08638  ------------       CD  ------------ <-- currently considering
        // 0.11725  ------------      ACD  ------------
        // 0.13886      C      A             CD       A
        // 0.13942      D      A             CD       A
        // 0.17349     AC     AD             AC      AD
        // 0.39425    ACD      E            ACD       E
        // 0.41254      B   ACDE              B    ACDE
        //
        //          -- before --           -- after ---
        //  height   left  right  species    left   right
        // -------  -----  -----  -------  -----  -----
        // 0.01167      E      E              E       E
        // 0.01389      B      B              B       B
        // 0.02047      B      B              B       B
        // 0.02079      D      D              D       D
        // 0.02150      D      D              D       D
        // 0.03438      C      C              C       C
        // 0.08638  ------------       CD  ------------
        // 0.11725  ------------      ACD  ------------ <-- currently considering
        // 0.13886     CD      A            ACD     ACD
        // 0.13942     CD      A            ACD     ACD
        // 0.17349     AC     AD            ACD     ACD
        // 0.39425    ACD      E            ACD       E
        // 0.41254      B   ACDE              B    ACDE

        unsigned jstart = 0;
        for (auto & si : sppinfo_vect) {
            // Get handles for current sppinfo_vect element
            double                 h0 = get<0>(si);
            vector<G::species_t> & v0 = get<2>(si);
            
            // Can't assume v0 has size 2 because this species tree join may be the fake
            // join that combines all remaining species into an ancestral species
            // for incomplete species forests
            G::species_t           s0 = accumulate(v0.begin(), v0.end(), (G::species_t)0, [](const G::species_t & next, const G::species_t & cum){return (cum | next);});
            
            // Advance through coalinfo_vect to the first event that might be
            // affected by the current species tree join
            unsigned j = jstart;
            assert(j < coalinfo_vect.size());
            double h = get<0>(coalinfo_vect[j]);
            while (h < h0) {
                j++;
                assert(j < coalinfo_vect.size());
                h = get<0>(coalinfo_vect[j]);
            }
            
            // Next time can start at this point in the coalinfo_vect
            jstart = j;
            
            // Perform replacements in all subsequent gene tree coalinfo elements
            map<G::species_t, G::species_t> replacements;
            while (j < coalinfo_vect.size()) {
                h = get<0>(coalinfo_vect[j]);
                vector<G::species_t> & v = get<2>(coalinfo_vect[j]);
                assert(v.size() == 2);
                for (auto & kv : replacements) {
                    if (subsumed(kv.first, v[0])) {
                        v[0] |= kv.second;
                    }
                    if (subsumed(kv.first, v[1])) {
                        v[1] |= kv.second;
                    }
                }
                if (subsumed(v[0], s0)) {
                    replacements[v[0]] = s0;
                    v[0] = s0;
                }
                if (subsumed(v[1], s0)) {
                    replacements[v[1]] = s0;
                    v[1] = s0;
                }
                j++;
            }
        }
    }

    inline void SpeciesForest::saveCoalInfoSpeciesTree(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap) {
        // Appends to coalinfo_vect; clear coalinfo_vect before calling if desired
        // Assumes heights and preorders are up-to-date
        
        // Dump _coalinfo into coalinfo_vect
        if (!_coalinfo.empty()) {
            coalinfo_vect.insert(coalinfo_vect.begin(), _coalinfo.begin(), _coalinfo.end());
        }
        
        if (cap) {
            // Create an entry pooling the remaining species
            if (_lineages.size() > 1) {
                vector<G::species_t> sppvect;
                for (auto & nd : _lineages) {
                    sppvect.push_back(nd->_species);
                }
                coalinfo_vect.push_back(make_tuple(
                    _forest_height,
                    0,
                    sppvect)
                );
            }
        }
    }

    inline void SpeciesForest::refreshAllPreorders() const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        _next_node_number = G::_nspecies;
        _preorders.clear();
        if (_lineages.size() == 0)
            return;
        
        for (auto & nd : _lineages) {
//            if (nd->_left_child) {
//                nd->_number = _next_node_number++;
//            }
            
            // lineage is a Node::ptr_vect_t (i.e. vector<Node *>)
            // lineage[0] is the first node pointer in the preorder sequence for this lineage
            // Add a new vector to _preorders containing, for now, just the root of the subtree
            _preorders.push_back({nd});
            
            // Now add the nodes above the root in preorder sequence
            Node::ptr_vect_t & preorder_vector = *_preorders.rbegin();
            refreshPreorderNew(preorder_vector);
        }
        
        unsigned number = 0;
        for (auto nd:_nodes) {
            nd._number = number;
            number++;
        }
    }

    inline void SpeciesForest::refreshPreorderNew(vector<Node*> & preorder) const {
        // Assumes preorder just contains the root node when this function is called
        // Also assumes that _next_node_number was initialized prior to calling this function
        assert(preorder.size() == 1);
        
        Node * nd = preorder[0];
        while (true) {
            nd = findNextPreorderNew(nd);
            if (nd) {
                preorder.push_back(nd);
//                if (nd->_left_child)
//                    nd->_number = _next_node_number++;
            }
            else
                break;
        }
    }

    inline void SpeciesForest::setTreeHeight() {
        _forest_height = _lineages.back()->_height;
    }

    inline void SpeciesForest::setNodeHeights() {
        assert (_preorder.size() > 0);
        for (auto &nd:_nodes) {
            if (!nd._left_child) {
                nd._height = 0.0;
            }
            else {
                nd._height = getLineageHeight(nd._left_child);
            }
        }
    }

    inline void SpeciesForest::storeSplits(set<Split> & internal_splits, set<Split> & leaf_splits) {
        // Start by clearing and resizing all splits
        
        // TODO: only works for 26 characters?
        unsigned count = 0;
        map<char, unsigned> names_and_bits;
        for (char letter = 'A'; letter < 'A' + G::_nspecies; letter++) {
            names_and_bits[letter] = count;
            count++;
         }
        
        refreshPreorder();
        for (auto & nd : _nodes) {
            nd._split.resize(G::_ntaxa);
        }

        // Now do a postorder traversal and add the bit corresponding
        // to the current node in its parent node's split
        for (auto nd : adaptors::reverse(_preorder)) {
            // Set split's edge length
            nd->_split.setEdgeLen(nd->_edge_length);

            if (nd->_left_child) {
                if (nd->_edge_length > 0.0) {// don't include root
                    // add this internal node's split to splitset
                    internal_splits.insert(nd->_split);
                }
            }
            else {
                // set bit corresponding to this leaf node's number
                nd->_split.setBitAt(names_and_bits[nd->_name[0]]);
//                nd->_split.setBitAt(nd->_number);
                leaf_splits.insert(nd->_split);
            }

            if (nd->_parent) {
                // parent's bits are the union of the bits set in all its children
                nd->_parent->_split.addSplit(nd->_split);
            }
        }
    }

    inline Node * SpeciesForest::findNextPreorderNew(Node * nd) const {
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

#if defined(LAZY_COPYING)
    inline void SpeciesForest::saveCoalInfo(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap) const {
        // Appends to coalinfo_vect; clear coalinfo_vect before calling if desired
        // Assumes heights and preorders are up-to-date

        // Dump _coalinfo into coalinfo_vect
        if (!_coalinfo.empty()) {
            coalinfo_vect.insert(coalinfo_vect.begin(), _coalinfo.begin(), _coalinfo.end());
        }
        
        if (cap) {
            // Create an entry pooling the remaining species
            if (_lineages.size() > 1) {
                vector<G::species_t> sppvect;
                for (auto nd : _lineages) {
                    sppvect.push_back(nd->_species);
                }
                coalinfo_vect.push_back(make_tuple(
                    _forest_height,
                    0,
                    sppvect)
                );
            }
        }
    }
#endif
}


