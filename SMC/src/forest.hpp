#pragma once

#include <stack>
#include <memory>
#include <iostream>
#include <boost/format.hpp>
#include <vector>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/algorithm/string/replace.hpp>


#include "lot.hpp"
extern proj::Lot rng;

#include "partial_store.hpp"
extern proj::PartialStore ps;

#include "node.hpp"

namespace proj {

using namespace std;

class Likelihood;
class Updater;
class TreeUpdater;
class PolytomyUpdater;
class Particle;

class Forest {

        friend class Likelihood;
        friend class Updater;
        friend class TreeUpdater;
        friend class PolytomyUpdater;
        friend class Particle;


    public:
                                    Forest();
                                    ~Forest();
        Forest(const Forest & other);

        unsigned                    numLeaves() const;
        unsigned                    numInternals() const;
        unsigned                    numNodes() const;
        void                        showForest();
        static void                 setNumSpecies(unsigned n);
        static void                 setNumTaxa(unsigned n);
        double                      calcLogLikelihood();
        void                        createDefaultTree();
        void operator=(const Forest & other);
        void                        debugForest();

    private:
    
        typedef std::vector <double> partial_array_t;
        void                        clear();
        void                        setData(Data::SharedPtr d, int index, map<string, string> &taxon_map);
        Node *                      findNextPreorder(Node * nd);
//        std::string                 makeNewick(unsigned precision, bool use_names) const;
        std::string                 makeNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>    chooseTaxaToJoin(double s);
        void                        detachSubtree(Node * s);
        void                        insertSubtreeOnLeft(Node * s, Node * u);
        Node *                      findLeftSib(Node * nd);
        Node *                      getSubtreeAt(unsigned i);
        unsigned                    getNumSubtrees();
        void                        setEdgeLength(Node * nd);
        void                        createNewSubtree(unsigned t1, unsigned t2);
        void                        drawCoalescenceTime();
        void                        calcPartialArray(Node* new_nd);
        Node*                       getNode(int i);
        void                        setUpGeneForest(map<string, string> &taxon_map);
        void                        setUpSpeciesForest(vector<string> &species_names);
        tuple<string,string, string> speciesTreeProposal();
        void                        geneTreeProposal(tuple<string, string, string> &species_merge_info, double time_increment);
//        void                        evolveSpeciesFor(list<Node*> &nodes, double time_increment) ;
        void                        evolveSpeciesFor(list <Node*> &nodes, double time_increment) ;
        void                        fullyCoalesceGeneTree(list<Node*> &nodes);
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        setUpLineages(vector<string> &species_names);

        std::vector<Node *>         _lineages;
        std::vector<Node>           _nodes;
    
        unsigned                    _nleaves;
        unsigned                    _ninternals;
        unsigned                    _npatterns;
        unsigned                    _nstates;
        double                      _last_edge_length;

        Data::SharedPtr             _data;
        static unsigned             _nspecies;
        static unsigned             _ntaxa;
        unsigned                    _first_pattern = 0;
        unsigned                    _index;
        map<string, list<Node*> > _species_partition;
        double                    _gene_tree_log_likelihood;

    public:
    
        typedef std::shared_ptr<Forest> SharedPtr;
        static double               _theta;
        static double               _speciation_rate;
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
        _lineages.clear();
        _nodes.resize(2*_ntaxa+2*_nspecies);
        _npatterns = 0;
        _nstates = 4;

        //create taxa
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd=&_nodes[i];
            if (i>0) {
                _nodes[i-1]._right_sib=nd;
            }
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            }
        _nleaves=_ntaxa;
        _ninternals=_nspecies;
        
        _lineages.reserve(_nodes.size()); //no root or subroot anymore

        for (int i=0; i < _ntaxa; i++) {
            _lineages.push_back(&_nodes[i]);
        }
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::setData(Data::SharedPtr d, int index, map<string, string> &taxon_map) {
        _data = d;
        _index = index;
        //index so it doesn't do this for species tree
//        _npatterns = d->getNumPatterns();
        
        if (index>0) {
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(index-1);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatternsInSubset(index-1);
            }
        
        const Data::taxon_names_t & taxon_names = _data->getTaxonNames();
        unsigned i = 0;
        auto data_matrix=_data->getDataMatrix();
        
        for (auto nd:_lineages) {
            if (!nd->_left_child) {
                // replace all spaces with underscores so that other programs do not have
                  // trouble parsing your tree descriptions
                  std::string name = taxon_names[i++];
                  boost::replace_all(name, " ", "_");
                nd->_name = name;
            
                if (index>0) {
                    nd->_partial=ps.getPartial(_npatterns*4);
                    for (unsigned p=0; p<_npatterns; p++) {
                        unsigned pp = _first_pattern+p;
                        for (unsigned s=0; s<_nstates; s++) {
                            Data::state_t state = (Data::state_t)1 << s;
                            Data::state_t d = data_matrix[nd->_number][pp];
                            double result = state & d;
                            (*nd->_partial)[p*_nstates+s]= (result == 0.0 ? 0.0:1.0);
                        }
                    }
                }
            }
        }

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
        if (_index > 0) {
            cout << " gene tree " << _index << ": " << endl;
            cout << " _gene_tree_log_likelihood: " << _gene_tree_log_likelihood << "\n";
        }
        else if (_index == 0) {
            cout << " species tree: " << endl;
        }
        cout << " " << makeNewick(9, true) << "\n";
        cout << "\n";
    }

//    inline std::string Forest::makeNewick(unsigned precision, bool use_names) const {
        inline string Forest::makeNewick(unsigned precision, bool use_names) {
                string newick = "(";
                const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
                const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
                const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
                stack<Node *> node_stack;
                
                unsigned i = 0;
                for (auto lineage : _lineages) {
                    Node * nd = lineage;
                    while (nd) {
                        if (nd->_left_child) {
                            // internal node
                            newick += "(";
                            node_stack.push(nd);
                        }
                        else {
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

    inline void Forest::setNumTaxa(unsigned n){
        _ntaxa=n;
    }

    inline void Forest::setNumSpecies(unsigned n) {
        _nspecies=n;
    }

    inline pair<unsigned, unsigned> Forest::chooseTaxaToJoin(double s){
        double nsubtrees = s;
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        if (nsubtrees > 2) {
            t1 = ::rng.randint(0, nsubtrees-1);
            t2 = ::rng.randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = ::rng.randint(0, nsubtrees-1);
            }
        }
        return make_pair(t1, t2);
    }

    inline void Forest::calcPartialArray(Node* new_nd) {
        auto & parent_partial_array = *(new_nd->_partial);
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
            double expterm = exp(-4.0*(child->_edge_length)/3.0);
            double prsame = 0.25+0.75*expterm;
            double prdif = 0.25 - 0.25*expterm;
            
            auto & child_partial_array = *(child->_partial);
            
            for (unsigned p = 0; p < _npatterns; p++) {
                for (unsigned s = 0; s <_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < _nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? prsame : prdif);
                        double child_partial = child_partial_array[p*_nstates + s_child];
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    if (child == new_nd->_left_child)
                        parent_partial_array[p*_nstates+s] = sum_over_child_states;
                    else
                        parent_partial_array[p*_nstates+s] *= sum_over_child_states;
                }   // parent state loop
            }   // pattern loop
        }   // child loop
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
        Node * the_nd = _lineages[i];
        assert (the_nd);
        return the_nd;
    }

    inline double Forest::calcLogLikelihood() {
        auto data_matrix=_data->getDataMatrix();
        
        //calc likelihood for each lineage separately
        auto &counts = _data->getPatternCounts();
        _gene_tree_log_likelihood = 0.0;
        
//        //if tree is fully resolved, calc likelihood from root
//        bool fully_resolved = false;
//        if (_lineages.size()==1) {
//            fully_resolved = true;
//        }
        
        for (auto nd:_lineages) {
            Node* base_node = nd;
//            if (fully_resolved){
//                Node* base_node = nd->_left_child;
//            }
            for (auto nd=base_node->_left_child; nd; nd=nd->_right_sib) {
                double log_like = 0.0;
                for (unsigned p=0; p<_npatterns; p++) {
                    double site_like = 0.0;
                    for (unsigned s=0; s<_nstates; s++) {
                        double partial = (*nd->_partial)[p*_nstates+s];
                        site_like += 0.25*partial;
                    }
                    assert(site_like>0);
                    log_like += log(site_like)*counts[_first_pattern+p];
                }
                _gene_tree_log_likelihood += log_like;
            }
        }
        return _gene_tree_log_likelihood;
    }

    inline void Forest::createDefaultTree() {
        clear();
        //create taxa
        double edge_length = rng.gamma(1.0, 1.0/_ntaxa);
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd=&_nodes[i];
            if (i>0) {
                _nodes[i-1]._right_sib=nd;
            }
            nd->_name="";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length = edge_length;
            }
        _nleaves=_ntaxa;
        _ninternals=_nspecies;
    }

    inline void Forest::operator=(const Forest & other) {
        _nstates = other._nstates;
        _npatterns = other._npatterns;
        _nodes.clear();
        _nodes.resize(other._nodes.size());
        _lineages.resize(other._lineages.size());
        _nleaves          = other._nleaves;
        _ninternals       = other._ninternals;
        _last_edge_length = other._last_edge_length;
        _index              = other._index;
        _first_pattern      = other._first_pattern;
        _gene_tree_log_likelihood = other._gene_tree_log_likelihood;
        _data                   = other._data;

        // copy tree itself
        
        _species_partition.clear();
        for (auto spiter : other._species_partition) {
            for (auto nd : spiter.second) {
                unsigned number = nd->_number;
                _species_partition[spiter.first].push_back(&_nodes[number]);
            }
        }
        
        for (auto othernd : other._nodes) {
            // get number of next node in preorder sequence (serves as index of node in _nodes vector)
            int k = othernd._number;
            
            if (k>-1) {

            // copy parent
//            assert(othernd._parent);
                if (othernd._parent) {
                    unsigned parent_number = othernd._parent->_number;
                    _nodes[k]._partial = othernd._partial;
                
                _nodes[k]._parent = &_nodes[parent_number];
            }

            // copy left child
                if (othernd._left_child) {
                unsigned left_child_number = othernd._left_child->_number;
                _nodes[k]._left_child = &_nodes[left_child_number];
            }
            else
                _nodes[k]._left_child = 0;

            // copy right sibling
            if (othernd._right_sib) {
                unsigned right_sib_number = othernd._right_sib->_number;
                _nodes[k]._right_sib = &_nodes[right_sib_number];
            }
            else
                _nodes[k]._right_sib = 0;

            _nodes[k]._number      = othernd._number;
            _nodes[k]._name        = othernd._name;
            _nodes[k]._edge_length = othernd._edge_length;
            }
        }
        
        unsigned i = 0;
        
        
        for (auto othernd : other._lineages) {
            unsigned k = othernd->_number;
            _lineages[i] = &_nodes[k];
            i++;
        }
    }

    inline void Forest::setUpSpeciesForest(vector<string> &species_names) {
        assert (_index==0);
        assert (_nspecies = (unsigned) species_names.size());
        clear();

        //create species
        double edge_length = 0.0;
        for (unsigned i = 0; i < _nspecies; i++) {
            Node* nd=&_nodes[i];
            nd->_name=species_names[i];
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length = edge_length;
            _lineages.push_back(nd);
            }
        _nleaves=_nspecies;
        _ninternals=0;
        _nodes.resize(_nspecies*2);
        _lineages.resize(_nspecies);
    }

    inline tuple<string,string, string> Forest::speciesTreeProposal() {
        if (_lineages.size() == 1) {
            _last_edge_length = -1.0;
            return make_tuple("null", "null", "null");
        }
        
        double rate = _speciation_rate*_lineages.size();
        _last_edge_length = rng.gamma(1.0, 1.0/rate);
        
        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
        
        pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size());
        Node* subtree1=_lineages[t.first];
        Node* subtree2=_lineages[t.second];
        
        Node* new_nd;

//        always need a new node now
        new_nd=&_nodes[_nleaves+_ninternals];
            new_nd->_parent=0;
            new_nd->_number=_nleaves+_ninternals;
            new_nd->_right_sib=0;
            
            _ninternals++;
            new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);
            new_nd->_left_child = _lineages[t.first];
            new_nd->_left_child->_right_sib = _lineages[t.second];
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
            child->_parent=new_nd;
        }
            updateNodeVector (_lineages, subtree1, subtree2, new_nd);
        
        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
    }


    inline void Forest::setUpGeneForest(map<string, string> &taxon_map) {
        assert (_index >0);
        _species_partition.clear();
        for (auto nd:_lineages) {
            if (!nd->_left_child) {
                string species_name = taxon_map[nd->_name];
                _species_partition[species_name].push_back(nd);
            }
        }
    }

    inline void Forest::setUpLineages(vector<string> &species_names) {
        assert (_index>0); //only set lineages for gene trees
        assert (_nspecies = (unsigned) species_names.size());
        _lineages.clear();

        //create species
        double edge_length = 0.0;
        for (unsigned i = 0; i < _nspecies; i++) {
            Node* nd=&_nodes[i+_ntaxa];
            nd->_name=species_names[i];
            nd->_number=i+_ntaxa;
            int j = 0;
            for (auto &s:_species_partition){
                if (j==i) {
                    //new node (species node) left child is first taxon in species list
                    nd->_left_child=s.second.front();
                    
                    //last taxon in lineage doesn't have a right sib
                    s.second.back()->_right_sib=0;
                    break;
                }
                j++;
            }
            
            //set parent to species node for each taxon in lineage
            //renumber each taxon in lineage
            int z = 0;
            for (Node* child=nd->_left_child; child; child=child->_right_sib) {
                child->_parent=nd;
//                child->_number=z;
                z++;
            }
            nd->_right_sib=0;
            nd->_parent=NULL;
            nd->_edge_length = edge_length;
            _lineages.push_back(nd);
            }
    }

    inline void Forest::evolveSpeciesFor(list<Node*> &nodes, double species_tree_increment) {
        bool done = false;
        double cum_time = 0.0;
        
        while (!done) {
            double s = nodes.size();
            double coalescence_rate = s*(s-1)/_theta;
            double increment = rng.gamma(1.0, 1.0/coalescence_rate);
            
            bool time_limited = species_tree_increment > 0.0;
            bool lineages_left_to_join = s > 1;
            bool time_limit_reached = cum_time+increment > species_tree_increment;
            if ((time_limited && time_limit_reached) || !lineages_left_to_join)  {
                done = true;
                increment = species_tree_increment - cum_time;
            }
            //add increment to each lineage

            for (auto nd:nodes) {
                //add new branch length to each lineage
                    nd->_edge_length += increment; //add most recently chosen branch length to each node whose parent is subroot
            }
            
            if (!done) {
                pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
                Node *subtree1 = nullptr;
                Node *subtree2 = nullptr;

                unsigned i = 0;
                for (auto iter=nodes.begin(); iter != nodes.end(); iter++){
                    if (i==t.first) {
                        subtree1 = *iter;
                    }
                    else if (i==t.second) {
                        subtree2 = *iter;
                    }
                    if (subtree1 && subtree2) {
                        break;
                    }
                    i++;
                }
                
                //new node is always needed now
                //set parent before subtree 1 is detached
                
                Node* new_nd=&_nodes[_nleaves+_ninternals];
                new_nd->_parent=subtree1->_parent;

                detachSubtree(subtree1);
                detachSubtree(subtree2);
                
//                new_nd->_parent->_left_child=new_nd;
//                for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
//                    child->_parent=new_nd;
//                }

                new_nd->_right_sib=0;
                new_nd->_number=_nleaves+_ninternals;

                _ninternals++;
                insertSubtreeOnLeft(new_nd, new_nd->_parent);

                new_nd->_left_child=subtree1;
                new_nd->_edge_length=0.0;

                subtree1->_right_sib=subtree2;
                subtree1->_parent=new_nd;
                subtree2->_parent=new_nd;

                //always calculating partials now
                assert (new_nd->_partial == nullptr);
                new_nd->_partial=ps.getPartial(_npatterns*4);
                assert(new_nd->_left_child->_right_sib);
                calcPartialArray(new_nd);

                //update species list
                updateNodeList(nodes, subtree1, subtree2, new_nd);
            }
            cum_time += increment;
        }
    }

inline void Forest::fullyCoalesceGeneTree(list<Node*> &nodes) {
    bool done = false;
    
    while (!done) {
        double s = nodes.size();
        double coalescence_rate = s*(s-1)/_theta;
        double increment = rng.gamma(1.0, 1.0/coalescence_rate);
        
        bool lineages_left_to_join = s > 1;
        if (!lineages_left_to_join)  {
            done = true;
        }
        //add increment to each lineage
        
        
        if (!done) {
            for (auto nd:nodes) {
                    nd->_edge_length += increment; //add most recently chosen branch length to each lineage
//                }
            }
            pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
            Node *subtree1 = nullptr;
            Node *subtree2 = nullptr;
            
            unsigned i = 0;
            for (auto iter=nodes.begin(); iter != nodes.end(); iter++){
                if (i==t.first) {
                    subtree1 = *iter;
                }
                else if (i==t.second) {
                    subtree2 = *iter;
                }
                if (subtree1 && subtree2) {
                    break;
                }
                i++;
            }
            
            //new node is always needed now
            //set parent before subtree 1 is detached
            
            Node* new_nd=&_nodes[_nleaves+_ninternals];
            new_nd->_parent=subtree1->_parent;

            detachSubtree(subtree1);
            detachSubtree(subtree2);

            new_nd->_number=_nleaves+_ninternals;

            _ninternals++;
            insertSubtreeOnLeft(new_nd, new_nd->_parent);

            new_nd->_left_child=subtree1;
            new_nd->_edge_length=0.0;

            subtree1->_right_sib=subtree2;
            subtree1->_parent=new_nd;
            subtree2->_parent=new_nd;
            

            //always calculating partials now
            assert (new_nd->_partial == nullptr);
            new_nd->_partial=ps.getPartial(_npatterns*4);
            assert(new_nd->_left_child->_right_sib);
            calcPartialArray(new_nd);

            //update species list
            updateNodeList(nodes, subtree1, subtree2, new_nd);
        }
    }
}

    inline void Forest::geneTreeProposal(tuple<string, string, string> &species_merge_info, double time_increment) {
        list<Node*> nodes_to_evolve ;
        for (auto &nd:_lineages){
            for (Node* child = nd->_left_child; child; child=child->_right_sib) {
                nodes_to_evolve.push_back(child);
            }
            if (_lineages.size()==1) {
                fullyCoalesceGeneTree(nodes_to_evolve);
                nodes_to_evolve.clear();
            }
            else {
                evolveSpeciesFor(nodes_to_evolve, time_increment);
                nodes_to_evolve.clear();
            }
        }

            //if not fully coalesced, join two lineages and update species partition
        if (_lineages.size() > 1) {
            string new_name = get<2>(species_merge_info);
            string species1 = get<0>(species_merge_info);
            string species2 = get<1>(species_merge_info);
            
            //look for nodes in updated _lineages list
            Node *lineage1;
            Node *lineage2;
            
            for (auto &nd:_lineages){
                if (nd->_name == species1) {
//                    lineage1=nd;
                    lineage1=nd;
                    break;
                }
            }
            
            for (auto &nd:_lineages){
                if (nd->_name == species2){
                    lineage2=nd;
                    break;
                }
            }
            
            //create new lineage
            Node* new_lineage=&_nodes[_nleaves+_ninternals];
            new_lineage->_name=new_name;
            new_lineage->_parent=lineage1->_parent;
            
            new_lineage->_right_sib=0;
            new_lineage->_number=_nleaves+_ninternals;
            _ninternals++;
            
            //set all parents of children of lineage1 and lineage2 to new_lineage
            for (Node*child=lineage1->_left_child; child; child=child->_right_sib) {
                child->_parent=new_lineage;
            }
            for (Node*child=lineage2->_left_child; child; child=child->_right_sib) {
                child->_parent=new_lineage;
            }
            
            new_lineage->_left_child=lineage1->_left_child;
            for (Node*child=lineage1->_left_child; child; child=child->_right_sib) {
                if (!child->_right_sib){
                    child->_right_sib=lineage2->_left_child;
                    break;
                }
            }
            
            //update _lineages list with species merge info
            updateNodeVector(_lineages, lineage1, lineage2, new_lineage);

            list<Node*> &nodes = _species_partition[new_name];
            nodes.push_back(lineage1->_left_child);
            nodes.push_back(lineage2->_left_child);
            
//            copy(_species_partition[species1].begin(), _species_partition[species1].end(), back_inserter(nodes));
//            copy(_species_partition[species2].begin(), _species_partition[species2].end(), back_inserter(nodes));
            _species_partition.erase(species1);
            _species_partition.erase(species2);
        }
    }

    
    inline void Forest::debugForest() {
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
        cout << "   _nleaves " << _nleaves << " ";
        cout << "   _ninternals " << _ninternals << " ";
        cout << endl;
    }

    inline void Forest::updateNodeVector(vector<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode) {
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

    inline void Forest::updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode) {
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


}
