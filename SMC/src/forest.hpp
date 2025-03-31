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

#include "lot.hpp"
#include "g.hpp"
#include "stopwatch.hpp"
extern proj::Lot rng;
std::mutex mtx;

extern proj::StopWatch stopwatch;

#include "partial_store.hpp"
extern proj::PartialStore ps;

#include "node.hpp"

namespace proj {

using namespace std;

class Likelihood;
class Particle;

class Forest {

        friend class Likelihood;
        friend class Particle;

    public:
                                    Forest();
                                    ~Forest();
        Forest(const Forest & other);
    
#if defined (FASTER_SECOND_LEVEL)
        typedef tuple<double, unsigned, vector<G::species_t> >  coalinfo_t;
#endif

        unsigned                        numLeaves() const;
        unsigned                        numInternals() const;
        unsigned                        numNodes() const;
        void                            showForest();
        double                          calcLogLikelihood();
        void                            createDefaultTree(Lot::SharedPtr lot);
        void operator=(const Forest & other);
        void                            debugForest();
        void                            debugLogLikelihood(Node* nd, double log_like);

    private:

        void                            clear();
        void                            setData(Data::SharedPtr d, int index, map<string, string> &taxon_map, bool partials);
        void                            setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa);
        Node *                          findNextPreorder(Node * nd);
        string                          makeNewick(unsigned precision, bool use_names);
        string                          makeAltNewick(unsigned precision, bool use_names);
        string                          makePartialNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>        chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        void                            calcPartialArrayJC(Node* new_nd);
        void                            calcPartialArrayHKY(Node* new_nd);

        void                            setUpGeneForest(map<string, string> &taxon_map);
        void                            setUpSpeciesForest(vector<string> &species_names);
        tuple<string,string, string>    speciesTreeProposal(Lot::SharedPtr lot);
        void                            updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                            updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                            revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                          getRunningSumChoices(vector<double> &log_weight_choices);
        vector<double>                  reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        int                             selectPair(vector<double> weight_vec, Lot::SharedPtr lot);
        void                            chooseSpeciesIncrement(Lot::SharedPtr lot);
        pair<double,double>             chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth);
        void                            allowCoalescence(string species_name, double increment, Lot::SharedPtr lot);
        vector<pair<double, string>>    calcForestRate(Lot::SharedPtr lot);
        void                            updateSpeciesPartition(tuple<string, string, string> species_info);
        double                          calcTopologyPrior(unsigned nlineages);
    
#if defined (UNUSED_FUNCTIONS)
        void                            calcIncrementPrior(double increment, string species_name, bool new_increment, bool coalesced_gene, bool gene_tree);
#endif
    
        void                            clearPartials();
        void                            setRelativeRate(double rel_rate) {_relative_rate = rel_rate;}
        unsigned                        getDeepCoal(tuple <string, string, string> species_joined);
        unsigned                        getMaxDeepCoal(tuple <string, string, string> species_joined);
        void                            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
#if !defined (FASTER_SECOND_LEVEL)
        void                            resetDepthVector(tuple<string, string, string> species_joined);
        vector<pair<double, pair<string, string>>>             getMinDepths();
        void                            calcMinDepth();
#endif
        vector< pair<double, Node *>>   sortPreorder();
        void                            refreshPreorder();
        void                            createThetaMap(Lot::SharedPtr lot);
        void                            createThetaMapFixedTheta(Lot::SharedPtr lot);
        void                            updateThetaMap(Lot::SharedPtr lot, string new_species_name);
        void                            updateThetaMapFixedTheta(Lot::SharedPtr lot, string new_species_name);
        void                            resetThetaMap(Lot::SharedPtr lot);
        void                            drawNewTheta(string new_species, Lot::SharedPtr lot);
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
        pair<Node*, Node*>              chooseAllPairs(list<Node*> &nodes, double increment, string species, Lot::SharedPtr lot);
        tuple<Node*, Node*, Node*>      createNewSubtree(pair<unsigned, unsigned> t, list<Node*> node_list, double increment, string species);
        pair<Node*, Node*>              getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list);
        void                            debugShowDistanceMatrix(const vector<double> & d) const;
    
        void                            constructUPGMA();
        void                            destroyUPGMA();
        void                            mergeDMatrixPair(vector<Split> & dmatrows, vector<double> & dmatrix, Split & s1, Split & s2);
        void                            buildStartingUPGMAMatrix();
#if defined (OLD_UPGMA)
        void                            buildStartingRow();
        void                            buildRestOfTreeFaster();
        void                            createSpeciesIndices();
#endif
        vector<pair<tuple<string, string, string>, double>>                            resetLineages(Lot::SharedPtr lot);
        vector<pair<tuple<string, string, string>, double>> resetT();
    
#if defined (FASTER_SECOND_LEVEL)
        void                            saveCoalInfoInitial();
        void                            saveCoalInfoGeneForest(vector<Forest::coalinfo_t> & coalinfo_vect) const;
        void                            saveCoalInfoSpeciesTree(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap);
        void                            addCoalInfoElem(const Node *, vector<coalinfo_t> & recipient);
        void                            buildCoalInfoVect();
        void                            fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const;
        static bool                     subsumed(G::species_t test_species, G::species_t subtending_species);
        void                            refreshAllPreorders() const;
        void                            refreshPreorderNew(vector<Node*> & preorder) const;
        Node *                          findNextPreorderNew(Node * nd) const;
        pair<double,double>             chooseSpeciesIncrementOnlySecondLevel(Lot::SharedPtr lot, double max_depth);
        void                            setTreeHeight();

        
        vector<coalinfo_t>              _coalinfo;
        mutable vector<Node::ptr_vect_t> _preorders;
        mutable unsigned                _next_node_number;
#endif

        void                            setNodeHeights();
        void                            resetSpeciesPartition(string species_partition);
        map<string, vector<string>>     saveSpeciesPartition();
        void                            setGeneUPGMAMatrices();
    
        map<string, double>             _theta_map;
        std::vector<Node *>             _lineages;
        std::list<Node>                 _nodes;
        std::vector<Node*>              _new_nodes;

        unsigned                        _nleaves;
        unsigned                        _ninternals;
        unsigned                        _npatterns;
        double                          _last_edge_length;

        Data::SharedPtr                 _data;
    
        unsigned                        _first_pattern = 0;
        unsigned                        _index;
        map<string, list<Node*> >       _species_partition;
        double                          _gene_tree_log_likelihood;
        double                          _log_joining_prob;
        vector<pair<double, double>>    _increments_and_priors;
        bool                            _done;
    
#if defined (DEBUG_MODE)
        pair<Node*, Node*>                  _species_joined;
#endif
    
#if !defined (FASTER_SECOND_LEVEL)
        vector<pair<double, pair<string, string>>>          _depths;
        vector<pair<tuple<string,string,string>, double>>   _species_build;
        unsigned                                            _nincrements = 0;
        double                          _log_coalescent_likelihood;
        double                          _panmictic_coalescent_likelihood;
        double                          _log_coalescent_likelihood_increment;
#endif
        vector<Node*>                   _preorder;
        map<string, string>             _taxon_map;
        double                          _log_weight;
        string                          _ancestral_species_name;
        vector<double>                  _vector_prior;
        double                          _theta_mean;
        vector<pair<Node*, Node*>>      _node_choices;
        vector<double>                  _log_likelihood_choices;
#if defined (OLD_UPGMA)
        map<Node*,  unsigned>           _starting_row;
        map<string, unsigned>           _species_indices;
#endif
        stack<Node *>                   _upgma_additions;
        map<Node *, double>             _upgma_starting_edgelen;
        
        // Local copy of the precalculated global JC pairwise distance matrix
        vector<double>                  _dmatrix;
        vector<Split>                   _dmatrix_rows;
    
        vector<string>                  _species_names;
        vector<double>                  _starting_dij;
        double                          _relative_rate;
        vector<pair<string, unsigned>>  _lineages_per_species;
        unsigned                        _partials_calculated_count;
        double                          _forest_height;
    
#if defined (DEBUG_MODE)
        void                            showSpeciesJoined();
#endif
        double                          calcTransitionProbabilityJC(double s, double s_child, double edge_length);
        double                          calcTransitionProbabilityHKY(double s, double s_child, double edge_length);
        double                          calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length);
        double                          getTreeLength();
        double                          getLineageHeight(Node* nd);
        void                            addIncrement(double increment);
        void                            simulateData(Lot::SharedPtr lot, unsigned starting_site, unsigned nsites);
    
#if !defined (FASTER_SECOND_LEVEL)
        double                          calcCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height);
        pair<vector<double>, vector<unsigned>>  calcCoalescentLikelihoodIntegratingOutTheta(vector<pair<tuple<string,string,string>, double>> species_build);
        pair<vector<double>, vector<unsigned>>  calcInitialCoalescentLikelihoodIntegratingOutTheta();
        pair<vector<double>, vector<unsigned>>  calcCoalescentLikelihoodIntegratingOutThetaLastStep(vector<pair<tuple<string,string,string>, double>> species_build);
#endif
        
        unsigned                        multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);

    public:

        typedef std::shared_ptr<Forest> SharedPtr;
        static double                   _kappa;
        static double                   _edge_rate_variance;
        static double                   _asrv_shape;
        static double                   _comphet;
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
        _npatterns = 0;
        _last_edge_length = 0.0;
        _lineages.clear();
        _log_weight = 0.0;
        _gene_tree_log_likelihood = 0.0;
        _log_joining_prob = 0.0;
        _nleaves=G::_ntaxa;
        _ninternals=0;
        _preorder.clear();
        _theta_mean = 0.0;
        _ancestral_species_name = "";
#if !defined (FASTER_SECOND_LEVEL)
        _species_build.clear();
        _depths.clear();
        _panmictic_coalescent_likelihood = 0.0;
        _log_coalescent_likelihood = 0.0;
        _log_coalescent_likelihood_increment = 0.0;
#endif
        _taxon_map.clear();
        _vector_prior.clear();
#if defined (OLD_UPGMA)
        _starting_row.clear();
        _nincrements = 0;
        _species_indices.clear();
#endif
        _upgma_additions = stack<Node*>();
        _upgma_starting_edgelen.clear();

        _species_names.clear();
        _starting_dij.clear();
        _lineages_per_species.clear();
        _partials_calculated_count = 0;
        _forest_height = 0.0;
#if defined (FASTER_SECOND_LEVEL)
        _coalinfo.clear();
        _preorders.clear();
        _next_node_number = 0;
#endif
        
#if defined (UNUSED_FUNCTIONS)
        _done = false;
#endif
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::setGeneUPGMAMatrices() {
        assert (_index > 0);
        _dmatrix      = G::_dmatrix[_index-1];
        _dmatrix_rows = G::_dmatrix_rows;
    }

    inline void Forest::setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa) {
        _index = index;
        assert (index > 0);         //don't set data for species tree, though it doesn't really matter for simulations
        G::_ntaxa = ntaxa;
        _nleaves = ntaxa;
        
        _data = d;
        
        _nodes.resize(G::_ntaxa);
        _lineages.reserve(_nodes.size());
        unsigned i= 0;
        
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_height = 0.0;
            nd->_position_in_lineages=i;
            _lineages.push_back(nd);
        }
        
        vector<string> taxon_names;
        for (auto &t:taxon_map) {
            taxon_names.push_back(t.first);
        }
        
        for (auto &nd:_lineages) {
            if (!nd->_left_child) {
                // replace all spaces with underscores so that other programs do not have
                  // trouble parsing your tree descriptions
                  std::string name = taxon_names[i++];
                  boost::replace_all(name, " ", "_");
                nd->_name = name;
            }
        }
    }

    inline void Forest::setData(Data::SharedPtr d, int index, map<string, string> &taxon_map, bool partials) {
        _data = d;
        _index = index;
        assert (index > 0);         //don't set data for species tree

        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(index-1);
        _first_pattern = gene_begin_end.first;
        _npatterns = _data->getNumPatternsInSubset(index-1);
        
        const Data::taxon_names_t & taxon_names = _data->getTaxonNames();
        unsigned i = 0;
        auto &data_matrix=_data->getDataMatrix();
        
        _nodes.resize(G::_ntaxa);
        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_height = 0.0;
            nd->_position_in_lineages=i;
            _lineages.push_back(nd);
            if (G::_upgma) {
                // set splits
                nd->_split.resize(G::_ntaxa);
                nd->_split.setBitAt(i);
            }
        }

        for (auto &nd:_lineages) {
            if (!nd->_left_child) {
                // replace all spaces with underscores so that other programs do not have
                  // trouble parsing your tree descriptions
                  std::string name = taxon_names[i++];
                  boost::replace_all(name, " ", "_");
                nd->_name = name;

                if (!G::_save_memory || (G::_save_memory && partials)) { // if save memory setting, don't set tip partials yet
                    nd->_partial=ps.getPartial(_npatterns*4);
                    for (unsigned p=0; p<_npatterns; p++) {
                        unsigned pp = _first_pattern+p;
                        for (unsigned s=0; s<G::_nstates; s++) {
                            Data::state_t state = (Data::state_t)1 << s;
                            Data::state_t d = data_matrix[nd->_number][pp];
                            double result = state & d;
                            (*nd->_partial)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
                        }
                    }
                }
            }
        }

# if defined  (LIKELIHOOD_TEST)
        for (auto &nd:_lineages) {
            nd->_edge_length += 1.4617;
        }
        calcLogLikelihood();
        
        Node* subtree1 = _lineages[5];
        Node* subtree2 = _lineages[4];
        
        Node nd;
        _nodes.push_back(nd);
        Node* new_nd = &_nodes.back();
        new_nd->_parent=0;
        new_nd->_number=_nleaves+_ninternals;
        new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);
        new_nd->_edge_length=0.0;
        _ninternals++;
        new_nd->_right_sib=0;

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd);
        
        assert (new_nd->_partial == nullptr);
        new_nd->_partial=ps.getPartial(_npatterns*4);
        assert(new_nd->_left_child->_right_sib);
        
        _relative_rate = 1.0;
        calcPartialArray(new_nd);
        
//        showForest();
        calcLogLikelihood();
        
        // next step
        for (auto &nd:_lineages) {
            nd->_edge_length += 4.1553;
//            nd->_edge_length += 0.01;
        }
        subtree1 = _lineages[1];
        subtree2 = _lineages[3];
        
        Node nd2;
        _nodes.push_back(nd2);
        Node* new_nd2 = &_nodes.back();
        new_nd2->_parent=0;
        new_nd2->_number=_nleaves+_ninternals;
        new_nd2->_name=boost::str(boost::format("node-%d")%new_nd2->_number);
        new_nd2->_edge_length=0.0;
        _ninternals++;
        new_nd2->_right_sib=0;

        new_nd2->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd2;
        subtree2->_parent=new_nd2;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd2);
        
        assert (new_nd2->_partial == nullptr);
        new_nd2->_partial=ps.getPartial(_npatterns*4);
        assert(new_nd2->_left_child->_right_sib);
        
        calcPartialArray(new_nd2);
        
//        showForest();
        calcLogLikelihood();
        
        // third step
        for (auto &nd:_lineages) {
            nd->_edge_length += 0.4611;
        }
        subtree1 = _lineages[2];
        subtree2 = _lineages[4];
        
        Node nd3;
        _nodes.push_back(nd3);
        Node* new_nd3 = &_nodes.back();
        new_nd3->_parent=0;
        new_nd3->_number=_nleaves+_ninternals;
        new_nd3->_name=boost::str(boost::format("node-%d")%new_nd3->_number);
        new_nd3->_edge_length=0.0;
        _ninternals++;
        new_nd3->_right_sib=0;

        new_nd3->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd3;
        subtree2->_parent=new_nd3;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd3);
        
        assert (new_nd3->_partial == nullptr);
        new_nd3->_partial=ps.getPartial(_npatterns*4);
        assert(new_nd3->_left_child->_right_sib);
        
        
//        showForest();
        
        calcPartialArray(new_nd3);
        
        calcLogLikelihood();
        
        // fourth step
        for (auto &nd:_lineages) {
            nd->_edge_length += 9.4541;
//            nd->_edge_length += 0.1;
        }
        subtree1 = _lineages[1];
        subtree2 = _lineages[2];
        
        Node nd4;
        _nodes.push_back(nd4);
        Node* new_nd4 = &_nodes.back();
        new_nd4->_parent=0;
        new_nd4->_number=_nleaves+_ninternals;
        new_nd4->_name=boost::str(boost::format("node-%d")%new_nd4->_number);
        new_nd4->_edge_length=0.0;
        _ninternals++;
        new_nd4->_right_sib=0;

        new_nd4->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd4;
        subtree2->_parent=new_nd4;
        
        updateNodeVector (_lineages, subtree1, subtree2, new_nd4);
        
        assert (new_nd4->_partial == nullptr);
        new_nd4->_partial=ps.getPartial(_npatterns*4);
        assert(new_nd4->_left_child->_right_sib);
        
        
//        showForest();
        
        calcPartialArray(new_nd4);
        
        calcLogLikelihood();
        
//        showForest();
#endif
        
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
//            cout << "log coalescent likelihood: " << _log_coalescent_likelihood << "\n";
        }
        else if (_index == 0) {
            cout << " species tree: " << endl;
        }
        cout << " " << makeNewick(15, true) << "\n";
        cout << "\n";
    }

    inline string Forest::makePartialNewick(unsigned precision, bool use_names) {
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

    inline string Forest::makeNewick(unsigned precision, bool use_names) {
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

    inline string Forest::makeAltNewick(unsigned precision, bool use_names) {
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

    inline pair<unsigned, unsigned> Forest::chooseTaxaToJoin(double s, Lot::SharedPtr lot){
        assert (s>1);
        double nsubtrees = s;
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        
        if (nsubtrees > 2) {
            t1 = lot->randint(0, nsubtrees-1);
            t2 = lot->randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = lot->randint(0, nsubtrees-1);
            }
        }
        assert(t1 < nsubtrees);
        assert (t2 < nsubtrees);

        return make_pair(t1, t2);
    }

    inline void Forest::calcPartialArrayJC(Node * new_nd) {
        _partials_calculated_count++;

        assert (_index > 0);
    
        if (!new_nd->_left_child) {
            auto &data_matrix=_data->getDataMatrix();
            assert (G::_save_memory || G::_start_mode == "sim");
//            assert (G::_save_memory || G::_start_mode_type == G::StartModeType::START_MODE_SIM);
            if (!new_nd->_left_child) {
                new_nd->_partial=ps.getPartial(_npatterns*4);
                for (unsigned p=0; p<_npatterns; p++) {
                    unsigned pp = _first_pattern+p;
                    for (unsigned s=0; s<G::_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[new_nd->_number][pp];
                        double result = state & d;
                        (*new_nd->_partial)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
                    }
                }
            }
        }
    
        auto & parent_partial_array = *(new_nd->_partial);
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {

            if (child->_partial == nullptr) {
                child->_partial = ps.getPartial(_npatterns*4);
                calcPartialArrayJC(child);
            }
            assert (child->_partial != nullptr);
            auto & child_partial_array = *(child->_partial);

            for (unsigned p = 0; p < _npatterns; p++) {
                for (unsigned s = 0; s <G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbabilityJC(s, s_child, child->_edge_length);
                        double child_partial = child_partial_array[p*G::_nstates + s_child];
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    if (child == new_nd->_left_child)
                        parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
                    else
                        parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
                }   // parent state loop
            }   // pattern loop
        }   // child loop
    }


    inline void Forest::calcPartialArrayHKY(Node * new_nd) {
        _partials_calculated_count++;

        assert (_index > 0);

        if (!new_nd->_left_child) {
            auto &data_matrix=_data->getDataMatrix();
            assert (G::_save_memory || G::_start_mode == "sim");
//            assert (G::_save_memory || G::_start_mode_type == G::StartModeType::START_MODE_SIM);
            if (!new_nd->_left_child) {
                new_nd->_partial=ps.getPartial(_npatterns*4);
                for (unsigned p=0; p<_npatterns; p++) {
                    unsigned pp = _first_pattern+p;
                    for (unsigned s=0; s<G::_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[new_nd->_number][pp];
                        double result = state & d;
                        (*new_nd->_partial)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
                    }
                }
            }
        }

        auto & parent_partial_array = *(new_nd->_partial);
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {

            if (child->_partial == nullptr) {
                child->_partial = ps.getPartial(_npatterns*4);
                calcPartialArrayHKY(child);
            }
            assert (child->_partial != nullptr);
            auto & child_partial_array = *(child->_partial);

            for (unsigned p = 0; p < _npatterns; p++) {
                for (unsigned s = 0; s <G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbabilityHKY(s, s_child, child->_edge_length);
                        double child_partial = child_partial_array[p*G::_nstates + s_child];
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    if (child == new_nd->_left_child)
                        parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
                    else
                        parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
                }   // parent state loop
            }   // pattern loop
        }   // child loop
    }

    inline double Forest::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) {
        assert(pi.size() == 4);
        assert(fabs(accumulate(pi.begin(), pi.end(), 0.0) - 1.0) < G::_small_enough);
        assert(_relative_rate > 0.0);
        double transition_prob = 0.0;
        
        // F81 transition probabilities
        double Pi[] = {pi[0] + pi[2], pi[1] + pi[3], pi[0] + pi[2], pi[1] + pi[3]};
        bool is_transition = (from == 0 && to == 2) || (from == 1 && to == 3) || (from == 2 && to == 0) || (from == 3 && to == 1);
        bool is_same = (from == 0 && to == 0) || (from == 1 && to == 1) | (from == 2 && to == 2) | (from == 3 && to == 3);
        bool is_transversion = !(is_same || is_transition);

        // HKY expected number of substitutions per site
        //  v = betat*(AC + AT + CA + CG + GC + GT + TA + TG) + kappa*betat*(AG + CT + GA + TC)
        //    = 2*betat*(AC + AT + CG + GT + kappa(AG + CT))
        //    = 2*betat*((A + G)*(C + T) + kappa(AG + CT))
        //  betat = v/[2*( (A + G)(C + T) + kappa*(AG + CT) )]
        double kappa = 1.0;
        double betat = 0.5*_relative_rate*edge_length/((pi[0] + pi[2])*(pi[1] + pi[3]) + kappa*(pi[0]*pi[2] + pi[1]*pi[3]));
        
        if (is_transition) {
            double pi_j = pi[to];
            double Pi_j = Pi[to];
            transition_prob = pi_j*(1.0 + (1.0 - Pi_j)*exp(-betat)/Pi_j - exp(-betat*(kappa*Pi_j + 1.0 - Pi_j))/Pi_j);
        }
        else if (is_transversion) {
            double pi_j = pi[to];
            transition_prob = pi_j*(1.0 - exp(-betat));
        }
        else {
            double pi_j = pi[to];
            double Pi_j = Pi[to];
            transition_prob = pi_j*(1.0 + (1.0 - Pi_j)*exp(-betat)/Pi_j) + (Pi_j - pi_j)*exp(-betat*(kappa*Pi_j + 1.0 - Pi_j))/Pi_j;
        }
        return transition_prob;
    }

    inline double Forest::calcTransitionProbabilityJC(double s, double s_child, double edge_length) {
        double child_transition_prob = 0.0;

            if (s == s_child) {
                child_transition_prob = 0.25 + 0.75*exp(-4.0*_relative_rate*edge_length/3.0);
            }
            
            else {
                child_transition_prob = 0.25 - 0.25*exp(-4.0*_relative_rate*edge_length/3.0);
            }
            return child_transition_prob;
    }

    inline double Forest::calcTransitionProbabilityHKY(double s, double s_child, double edge_length) {
    //        StopWatch sw;
    //        sw.start();
        double child_transition_prob = 0.0;
        
            double pi_A = G::_base_frequencies[0];
            double pi_C = G::_base_frequencies[1];
            double pi_G = G::_base_frequencies[2];
            double pi_T = G::_base_frequencies[3];

            double pi_j = 0.0;
            double PI_J = 0.0;

            double phi = (pi_A+pi_G)*(pi_C+pi_T)+_kappa*(pi_A*pi_G+pi_C*pi_T);
            double beta_t = 0.5*(edge_length * _relative_rate )/phi;

            // transition prob depends only on ending state
            if (s_child == 0) {
                // purine
                pi_j = pi_A;
                PI_J = pi_A + pi_G;
            }
            else if (s_child == 1) {
                // pyrimidine
                pi_j = pi_C;
                PI_J = pi_C + pi_T;
            }
            else if (s_child == 2) {
                // purine
                pi_j = pi_G;
                PI_J = pi_A + pi_G;
            }
            else if (s_child == 3) {
                // pyrimidine
                pi_j = pi_T;
                PI_J = pi_C + pi_T;
            }

            while (true) {
                if (s == s_child) {
                    // no transition or transversion
                    double first_term = 1+(1-PI_J)/PI_J*exp(-beta_t);
                    double second_term = (PI_J-pi_j)/PI_J*exp(-beta_t*(PI_J*_kappa+(1-PI_J)));
                    child_transition_prob = pi_j*first_term+second_term;
                    break;
                }

                else if ((s == 0 && s_child == 2) || (s == 2 && s_child == 0) || (s == 1 && s_child == 3) || (s == 3 && s_child==1)) {
                    // transition
                    double first_term = 1+(1-PI_J)/PI_J*exp(-beta_t);
                    double second_term = (1/PI_J)*exp(-beta_t*(PI_J*_kappa+(1-PI_J)));
                    child_transition_prob = pi_j*(first_term-second_term);
                    break;
                }

                else {
                    // transversion
                    child_transition_prob = pi_j*(1-exp(-beta_t));
                    break;
                }
            }
            assert (child_transition_prob > 0.0);
            return child_transition_prob;
    }

    inline double Forest::calcLogLikelihood() {
        //calc likelihood for each lineage separately
        auto &counts = _data->getPatternCounts();
        _gene_tree_log_likelihood = 0.0;
        
        if (_forest_height > 0 && G::_upgma) {
            constructUPGMA();
        } // TODO: be careful


        for (auto &nd:_lineages) {
            double log_like = 0.0;
            for (unsigned p=0; p<_npatterns; p++) {
                double site_like = 0.0;
                for (unsigned s=0; s<G::_nstates; s++) {
                    double partial = (*nd->_partial)[p*G::_nstates+s];
                    site_like += 0.25*partial;
                }
                assert(site_like>0);
                log_like += log(site_like)*counts[_first_pattern+p];
            }

            if (_forest_height > 0 && G::_upgma) {
                destroyUPGMA(); // TODO: be careful
            }

            _gene_tree_log_likelihood += log_like;
//            debugLogLikelihood(nd, log_like);
        }
        
        return _gene_tree_log_likelihood;
    }

    inline int Forest::selectPair(vector<double> weight_vec, Lot::SharedPtr lot) {
        // choose a random number [0,1]
        assert (lot != nullptr);
        double u = lot->uniform();
        
        double cum_prob = 0.0;
        int index = 0.0;
        for (int i=0; i < (int) weight_vec.size(); i++) {
            cum_prob += exp(weight_vec[i]);
            if (u <= cum_prob) {
                index = i;
                break;
            }
        }
        // return index of choice
        return index;
    }

    inline vector<double> Forest::reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood) {
        vector<double> weight_vec;
        for (int a = 0; a < (int) likelihood_vec.size(); a++) {
            weight_vec.push_back(likelihood_vec[a]-prev_log_likelihood);
        }
        return weight_vec;
    }

    inline double Forest::getRunningSumChoices(vector<double> &log_weight_choices) {
        double running_sum = 0.0;
        double log_weight_choices_sum = 0.0;
        double log_max_weight = *max_element(log_weight_choices.begin(), log_weight_choices.end());
        for (auto & i:log_weight_choices) {
            running_sum += exp(i - log_max_weight);
        }
        log_weight_choices_sum = log(running_sum) + log_max_weight;
        return log_weight_choices_sum;
    }

    inline unsigned Forest::getMaxDeepCoal(tuple <string, string, string> species_joined) {
        string species1 = get<0>(species_joined);
        string species2 = get<1>(species_joined);
        string species3 = get<2>(species_joined);
        
        unsigned lineages_to_coalesce = 0;
        for (auto &m:_lineages_per_species) {
            if (m.first == species1) {
                lineages_to_coalesce += m.second;
            }
            else if (m.first == species2) {
                lineages_to_coalesce += m.second;
            }
        }
        
        assert (lineages_to_coalesce > 0);
        
        unsigned max_deep_coal = lineages_to_coalesce - 1;
        
        // update _lineages_per_species
        for (auto it = _lineages_per_species.begin(); it != _lineages_per_species.end();) {
            if (it->first == species1) {
                it = _lineages_per_species.erase(it);
            }
            else if (it->first == species2) {
                it = _lineages_per_species.erase(it);
            }
            else {
                it++;
            }
        }
        
        _lineages_per_species.push_back(make_pair(species3, lineages_to_coalesce));
        
        
        return max_deep_coal;
    }

    inline void Forest::setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species) {
        unsigned i = 0;
        for (auto &s:_species_partition) {
            _lineages_per_species.push_back(make_pair(s.first, ntaxa_per_species[i]));
            i++;
        }
    }

    inline unsigned Forest::getDeepCoal(tuple <string, string, string> species_joined) {
        unsigned num_deep_coal = 0;
        
        string spp1 = get<0> (species_joined);
        string spp2 = get<1> (species_joined);
        
        unsigned nlineages1 = (unsigned) _species_partition[spp1].size();
        unsigned nlineages2 = (unsigned) _species_partition[spp2].size();
        
        num_deep_coal += nlineages1 - 1;

        num_deep_coal += nlineages2 - 1;
        
        return num_deep_coal;
    }

    inline void Forest::debugLogLikelihood(Node* nd, double log_like) {
        cout << "debugging likelihood" << endl;
        cout << "gene: " << _index << endl;
        cout << "calculating likelihood for node : " << endl;
        cout << "\t" << "name: " << nd->_name << endl;
        cout << "\t" << "number: " << nd->_number << endl;
        cout << "\t" << "edge length: " << nd->_edge_length << endl;
        if (nd->_left_child) {
            if (nd->_left_child->_name != "") {
                cout << "\t" << "left child name: " << nd->_left_child->_name << endl;
            }
            else
                cout << "\t" << "left child number: " << nd->_left_child->_number << endl;
        }
        else {
            cout << "\t" << "left child is null" << endl;
        }
        cout << "\t" << "node log likelihood: " << log_like << endl;
        cout << endl;
    }

    inline void Forest::createDefaultTree(Lot::SharedPtr lot) {
        clear();
        //create taxa
        assert (lot != nullptr);
        double edge_length = lot->gamma(1.0, 1.0/G::_ntaxa);
        
        _lineages.reserve(_nodes.size());
        
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name="";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length = edge_length;
            nd->_position_in_lineages=i;
            }
        _nleaves=G::_ntaxa;
        _ninternals=0;
        _last_edge_length = 0.0;
    }

    inline void Forest::operator=(const Forest & other) {
        _npatterns = other._npatterns;
        _edge_rate_variance = other._edge_rate_variance;
        _asrv_shape = other._asrv_shape;
        _comphet = other._comphet;
//        _nodes.clear(); // don't need to clear _nodes because they will get overwritten
        _nodes.resize(other._nodes.size());
        _lineages.resize(other._lineages.size());
        _new_nodes.resize(other._new_nodes.size());
        _nleaves            = other._nleaves;
        _ninternals         = other._ninternals;
        _last_edge_length   = other._last_edge_length;
        _index              = other._index;
        _first_pattern      = other._first_pattern;
        _gene_tree_log_likelihood = other._gene_tree_log_likelihood;
        _data               = other._data;
        _log_joining_prob = other._log_joining_prob;
        _increments_and_priors = other._increments_and_priors;
        _relative_rate = other._relative_rate;
        _preorder.resize(other._preorder.size());
        _theta_map = other._theta_map;
        _theta_mean = other._theta_mean;
        _ancestral_species_name = other._ancestral_species_name;
#if !defined (FASTER_SECOND_LEVEL)
        _species_build = other._species_build;
        _depths = other._depths;
        _panmictic_coalescent_likelihood = other._panmictic_coalescent_likelihood;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _log_coalescent_likelihood_increment = other._log_coalescent_likelihood_increment;
#endif
        _taxon_map = other._taxon_map;
        _vector_prior = other._vector_prior;
        _lineages_per_species = other._lineages_per_species;
#if defined (OLD_UPGMA)
        _starting_row = other._starting_row;
        _nincrements = other._nincrements;
        _species_indices = other._species_indices;
#endif
        _upgma_additions = other._upgma_additions;
        _upgma_starting_edgelen = other._upgma_starting_edgelen;

        _dmatrix = other._dmatrix;
        _dmatrix_rows = other._dmatrix_rows;
        
        _starting_dij = other._starting_dij;
        _species_names = other._species_names;
        _partials_calculated_count = other._partials_calculated_count;
        _forest_height = other._forest_height;
#if defined (FASTER_SECOND_LEVEL)
        _coalinfo = other._coalinfo;
        _preorders = other._preorders;
        _next_node_number = other._next_node_number;
#endif
        
#if defined (DEBUG_MODE)
        _species_joined = other._species_joined;
#endif
        
#if defined (UNUSED_FUNCTIONS)
        _done = other._done;
#endif

        // copy tree itself

        if (other._nodes.size() > 0) { // otherwise, there is no forest and nothing needs to be copied
            _species_partition.clear();
            for (auto spiter : other._species_partition) {
                for (auto s : spiter.second) {
                    unsigned number = s->_number;
                    Node* nd = &*next(_nodes.begin(), number);
                    _species_partition[spiter.first].push_back(nd);
                }
            }
            
#if defined (OLD_UPGMA)
            _starting_row.clear();
            for (auto strow : other._starting_row) {
                unsigned number = strow.first->_number;
                Node* nd = &*next(_nodes.begin(), number);
                _starting_row[nd] = strow.second;
            }
#endif

            for (auto othernd : other._nodes) {
                // get number of next node in preorder sequence (serves as index of node in _nodes vector)
                int k = othernd._number;

                if (k>-1) {
                    Node* nd = &*next(_nodes.begin(), k);

                // copy parent
                    if (othernd._parent) {
                        unsigned parent_number = othernd._parent->_number;
                        Node* parent = &*next(_nodes.begin(), parent_number);
                        nd->_parent = parent;
                    }

                // copy left child
                    if (othernd._left_child) {
                    unsigned left_child_number = othernd._left_child->_number;
                        Node* left_child = &*next(_nodes.begin(), left_child_number);
                        nd->_left_child = left_child;
                }
                    else {
                        nd->_left_child = 0;
                    }

                // copy right sibling
                if (othernd._right_sib) {
                    unsigned right_sib_number = othernd._right_sib->_number;
                    Node* right_sib = &*next(_nodes.begin(), right_sib_number);
                    nd->_right_sib = right_sib;
                }
                else
                    nd->_right_sib = 0;

                    nd->_number = othernd._number;
                    nd->_name = othernd._name;
                    nd->_edge_length = othernd._edge_length;
                    nd->_position_in_lineages = othernd._position_in_lineages;
                    nd->_partial = othernd._partial;
                    nd->_height = othernd._height;
                    nd->_split = othernd._split;
    #if defined (FASTER_SECOND_LEVEL)
                    nd->_species = othernd._species;
    #endif
                }
            }

            unsigned j = 0;
            for (auto othernd : other._lineages) {
                unsigned k = othernd->_number;
                Node* nd = &*next(_nodes.begin(), k);
                _lineages[j] = nd;
                j++;
            }
            
            if (other._preorder.size() > 0) {
                unsigned m = 0;
                for (auto othernd : other._preorder) {
                    unsigned n = othernd->_number;
                    Node* nd = &*next(_nodes.begin(), n);
                    _preorder[m] = nd;
                    m++;
                }
            }
        }
        
    }

    inline void Forest::setUpSpeciesForest(vector<string> &species_names) {
        _index = 0;
        assert (G::_nspecies = (unsigned) species_names.size());
        
        //create species
        _nodes.resize(G::_nspecies);
        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < G::_nspecies; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_height = 0.0;
            nd->_position_in_lineages=i;
            nd->_name=species_names[i];
            _lineages.push_back(nd);
#if defined (FASTER_SECOND_LEVEL)
            if (G::_start_mode != "sim") {
//            if (G::_start_mode_type != G::StartModeType::START_MODE_SIM) {
                if (G::_taxon_to_species.count(nd->_name) == 0) {
                    throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % nd->_name));
                }
                else {
                    Node::setSpeciesBit(nd->_species, G::_taxon_to_species.at(nd->_name), /*init_to_zero_first*/true);
                }
            }
#endif
            }
        
        _nleaves=G::_nspecies;
        _ninternals=0;
    }

    inline pair<double,double> Forest::chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth) {
        unsigned nlineages = (unsigned) _lineages.size();
        
        if (max_depth > 0.0) {
            double rate = (G::_lambda)*_lineages.size();
            
            double u = lot->uniform();
            double inner_term = 1-exp(-rate*max_depth);
            _last_edge_length = -log(1-u*inner_term)/rate;
            assert (_last_edge_length < max_depth);
            for (auto nd:_lineages) {
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

            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));

        }
        
#if !defined (FASTER_SECOND_LEVEL)
        if (_species_build.size() == 0) {
            _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
        }
        else {
            _species_build.back().second = _last_edge_length;
        }
#endif
        
        double constrained_factor = log(1 - exp(-1*nlineages*G::_lambda*max_depth));
        
        _forest_height += _last_edge_length;
        
        return make_pair(_last_edge_length, constrained_factor);

    }

#if defined (FASTER_SECOND_LEVEL)
    inline pair<double,double> Forest::chooseSpeciesIncrementOnlySecondLevel(Lot::SharedPtr lot, double max_depth) {
        double nlineages = (double) _lineages.size();
        
            max_depth = max_depth - _forest_height;
            
            assert (max_depth >= 0.0);
        
        if (max_depth > 0.0) {
            double rate = (G::_lambda)*_lineages.size();
            
            double u = lot->uniform();
            double inner_term = 1-exp(-rate*max_depth);
            _last_edge_length = -log(1-u*inner_term)/rate;
            assert (_last_edge_length < max_depth);
            for (auto nd:_lineages) {
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

            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));

        }
        
#if !defined (FASTER_SECOND_LEVEL)
        if (_species_build.size() == 0) {
            _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
        }
        else {
            _species_build.back().second = _last_edge_length;
        }
#endif
        
        double constrained_factor = log(1 - exp(-1*nlineages*G::_lambda*max_depth));
        
        _forest_height += _last_edge_length;
        
        return make_pair(_last_edge_length, constrained_factor);

    }
#endif


    inline void Forest::chooseSpeciesIncrement(Lot::SharedPtr lot) {
        double rate = G::_lambda*_lineages.size();
        
        assert (lot != nullptr);
        _last_edge_length = lot->gamma(1.0, 1.0/rate);

        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
        
#if defined (FASTER_SECOND_LEVEL)
        _forest_height += _last_edge_length;
#endif
    }


    inline tuple<string,string, string> Forest::speciesTreeProposal(Lot::SharedPtr lot) {
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
        
        Node nd;
        _nodes.push_back(nd);
        Node* new_nd = &_nodes.back();
        new_nd->_parent=0;
        new_nd->_number=_nleaves+_ninternals;
        new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);
        new_nd->_edge_length=0.0;
        _ninternals++;
        new_nd->_right_sib=0;

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

#if !defined (FASTER_SECOND_LEVEL)
        _species_build.push_back(make_pair(make_tuple(subtree1->_name, subtree2->_name, new_nd->_name), 0.0));
#endif
        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
    }

    inline void Forest::updateSpeciesPartition(tuple<string, string, string> species_info) {
        string spp1 = get<0>(species_info);
        string spp2 = get<1>(species_info);
        string new_spp = get<2>(species_info);
        
        unsigned before = (int) _species_partition.size();

        list<Node*> &nodes = _species_partition[new_spp];
        copy(_species_partition[spp1].begin(), _species_partition[spp1].end(), back_inserter(nodes));
        copy(_species_partition[spp2].begin(), _species_partition[spp2].end(), back_inserter(nodes));
        _species_partition.erase(spp1);
        _species_partition.erase(spp2);
        
        if (spp1 != "null") {
            assert (_species_partition.size() == before - 1);
        }
    }

#if defined (DEBUG_MODE)
    inline void Forest::showSpeciesJoined() {
        assert (_index==0);
        if (_species_joined.first != NULL) {
            cout << "joining species " << _species_joined.first->_name << " and " << _species_joined.second->_name << endl;
        }
        else {
            cout << "no species joined" << endl;
        }
    }
#endif

    inline void Forest::setUpGeneForest(map<string, string> &taxon_map) {
        _taxon_map = taxon_map;
        assert (_index >0);
        _species_partition.clear();
        
        unsigned count = 0;
        
        for (auto &nd:_nodes) {
            count++;
            assert (!nd._left_child);
            string species_name = taxon_map[nd._name];
            _species_partition[species_name].push_back(&nd);
            if (count == G::_ntaxa) {
                break;
            }
#if defined (FASTER_SECOND_LEVEL)
            if (G::_start_mode != "sim") {
//            if (G::_start_mode_type != G::StartModeType::START_MODE_SIM) {
                if (G::_taxon_to_species.count(nd._name) == 0) {
                    throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % nd._name));
                }
                else {
                    Node::setSpeciesBit(nd._species, G::_taxon_to_species.at(nd._name), /*init_to_zero_first*/true);
                }
            }
#endif
        }
        
        assert (_species_partition.size() == G::_nspecies);
    }
    inline double Forest::calcTopologyPrior(unsigned nlineages) {
        _log_joining_prob += -log(0.5*nlineages*(nlineages-1));
        assert (!isinf(_log_joining_prob));
        return _log_joining_prob;
    }

    inline void Forest::clearPartials() {
        for (auto &nd:_nodes) {
            nd._partial = nullptr;
        }
    }

#if defined (UNUSED_FUNCTIONS)
    inline void Forest::calcIncrementPrior(double increment, string species_name, bool new_increment, bool coalesced_gene, bool gene_tree) {
#if defined (DRAW_NEW_THETA)
        double log_increment_prior = 0.0;
        
        if (!_done) {
            if (gene_tree) {
                if (coalesced_gene) {
                    // calculate increment prior
                    for (auto &s:_species_partition) {
                        bool coalescence = false;
                        if (s.first == species_name) {
                            coalescence = true;
                        }
                        else {
                            coalescence = false;
                        }

                        if (coalescence) {
                            assert (s.second.size() > 1);
                            // if there is coalescence, need to use number of lineages before the join
                            double population_theta = _theta_map[s.first];
                            double coalescence_rate = (s.second.size())*(s.second.size()-1) / population_theta;
                            assert (coalescence_rate > 0.0); // rate should be >0 if there is coalescence
                            double nChooseTwo = (s.second.size())*(s.second.size()-1);
                            double log_prob_join = log(2/nChooseTwo);
                            log_increment_prior += log(coalescence_rate) - (increment*coalescence_rate) + log_prob_join;
                        }
                        else {
                            // no coalescence
                            double population_theta = _theta_map[s.first];
                            double coalescence_rate = s.second.size()*(s.second.size() - 1) / population_theta;
                            log_increment_prior -= increment*coalescence_rate;
                        }
                    }
                }
            
                else if (!coalesced_gene) {
                    // no coalescence
                    for (auto &s:_species_partition) {
                        double population_theta = _theta_map[s.first];
                        double coalescence_rate = s.second.size() * (s.second.size() - 1) / population_theta;
                        log_increment_prior -= increment*coalescence_rate;
                    }
                }
            
                if (new_increment) { // add a new increment to the list
                    _increments_and_priors.push_back(make_pair(increment, log_increment_prior));
                }
            
                else if (!new_increment) { // add to existing increment and prior
                    assert (_increments_and_priors.size() > 0);
                    _increments_and_priors.back().first += increment;
                    _increments_and_priors.back().second += log_increment_prior;
                }
                _log_coalescent_likelihood += log_increment_prior;
            }
        
            else {
                // species tree
                if (coalesced_gene) {
                    double rate = G::_lambda*(_lineages.size());
                    // calculate increment prior
                    double nChooseTwo = (_lineages.size())*(_lineages.size()-1);
                    double log_prob_join = log(2/nChooseTwo);
//                    log_increment_prior = log(_lambda) - (increment * rate) + log_prob_join;
//                    log_increment_prior = log(_lambda) - (increment * rate);
                    log_increment_prior = log(rate) - (increment*rate) + log_prob_join;
                }
                else {
                    double rate = G::_lambda*(_lineages.size());
                    // calculate increment prior
                    log_increment_prior = - (increment*rate);
                }
                
                if (new_increment) {
                    _increments_and_priors.push_back(make_pair(increment, log_increment_prior));
                }
                else {
                    _increments_and_priors.back().first += increment;
                    _increments_and_priors.back().second += log_increment_prior;
                }
                _log_coalescent_likelihood += log_increment_prior;
            }
            _log_coalescent_likelihood_increment = log_increment_prior;
        }
#else
        double log_increment_prior = 0.0;
        
        if (!_done) {
            if (gene_tree) {
                if (coalesced_gene) {
                    // calculate increment prior
                    for (auto &s:_species_partition) {
                        bool coalescence = false;
                        if (s.first == species_name) {
                            coalescence = true;
                        }
                        else {
                            coalescence = false;
                        }

                        if (coalescence) {
                            assert (s.second.size() > 1);
                            // if there is coalescence, need to use number of lineages before the join
                            double coalescence_rate = (s.second.size())*(s.second.size()-1) / G::_theta;
                            assert (coalescence_rate > 0.0); // rate should be >0 if there is coalescence
                            double nChooseTwo = (s.second.size())*(s.second.size()-1);
                            double log_prob_join = log(2/nChooseTwo);
                            log_increment_prior += log(coalescence_rate) - (increment*coalescence_rate) + log_prob_join;
                        }
                        else {
                            // no coalescence
                            double coalescence_rate = s.second.size()*(s.second.size() - 1) / G::_theta;
                            log_increment_prior -= increment*coalescence_rate;
                        }
                    }
                }
            
                else if (!coalesced_gene) {
                    // no coalescence
                    for (auto &s:_species_partition) {
                        double coalescence_rate = s.second.size() * (s.second.size() - 1) / G::_theta;
                        log_increment_prior -= increment*coalescence_rate;
                    }
                }
            
                if (new_increment) { // add a new increment to the list
                    _increments_and_priors.push_back(make_pair(increment, log_increment_prior));
                }
            
                else if (!new_increment) { // add to existing increment and prior
                    assert (_increments_and_priors.size() > 0);
                    _increments_and_priors.back().first += increment;
                    _increments_and_priors.back().second += log_increment_prior;
                }
                _log_coalescent_likelihood += log_increment_prior;
            }
        
            else {
                // species tree
                if (coalesced_gene) {
                    double rate = G::_lambda*(_lineages.size());
                    // calculate increment prior
                    double nChooseTwo = (_lineages.size())*(_lineages.size()-1);
                    double log_prob_join = log(2/nChooseTwo);
                    log_increment_prior = log(rate) - (increment*rate) + log_prob_join;
                }
                else {
                    double rate = G::_lambda*(_lineages.size());
                    // calculate increment prior
                    log_increment_prior = - (increment*rate);
                }
                
                if (new_increment) {
                    _increments_and_priors.push_back(make_pair(increment, log_increment_prior));
                }
                else {
                    _increments_and_priors.back().first += increment;
                    _increments_and_priors.back().second += log_increment_prior;
                }
                _log_coalescent_likelihood += log_increment_prior;
            }
            _log_coalescent_likelihood_increment = log_increment_prior;
        }
#endif
    }
#endif

#if defined (OLD_UPGMA)
    struct negLogLikeDist {
        negLogLikeDist(unsigned npatterns, unsigned first, const Data::pattern_counts_t & counts, const vector<double> & same, const vector<double> & diff, double v0)
            : _npatterns(npatterns), _first(first), _counts(counts), _same(same), _diff(diff), _v0(v0) {}
        
        double operator()(double const & v) {
            double edgelen = v + _v0;
            double tprob_same = 0.25 + 0.75*exp(-4.0*edgelen/3.0);
            double tprob_diff = 0.25 - 0.25*exp(-4.0*edgelen/3.0);

            double log_like = 0.0;
            for (unsigned p = 0; p < _npatterns; p++) {
                double site_like = 0.25 * (tprob_same * _same[p] + tprob_diff * _diff[p]);
                log_like += log(site_like) * _counts[_first + p];
            }
            
            return -log_like;
        }
        
        private:
            unsigned _npatterns;
            unsigned _first;
            const Data::pattern_counts_t & _counts;
            const vector<double> & _same;
            const vector<double> & _diff;
            double _v0;
    };
#endif
        
#if defined (OLD_UPGMA)
    inline void Forest::buildStartingUPGMAMatrix() {
        bool use_minimizer = true;
        
        if (!use_minimizer) {
            // Get the number of patterns
            unsigned npatterns = _data->getNumPatternsInSubset(_index - 1); // forest index starts at 1 but subsets start at 0

            // Get the first and last pattern index for this gene's data

            // Get pattern counts
            auto counts = _data->getPatternCounts();
            
            // Create vectors to store products of same-state and different-state partials
            vector<double> same_state(npatterns, 0.0);
            vector<double> diff_state(npatterns, 0.0);
            
            // Create a map relating position in dij vector to row,col in distance matrix
            map<unsigned, pair<unsigned, unsigned>> dij_row_col;
            
            // Create distance matrix dij and workspace dij2 used to build next dij
            // Both dij and dij2 are 1-dimensional vectors that store only the
            // lower diagonal of the distance matrix (excluding diagonal elements)
            assert (_lineages.size() == G::_ntaxa);
            unsigned n = G::_ntaxa;
            vector<double> dij(n*(n-1)/2, G::_infinity);
            vector<double> dij2;
            
            vector<tuple<unsigned, unsigned, unsigned, unsigned>> sites_tuples = _data->_partition->getSubsetRangeVect();
            
            for (unsigned i = 1; i < n; i++) {
                for (unsigned j = 0; j < i; j++) {
                    double ndiff = 0;
                    double ntotal = 0;
                    unsigned start_index = _index - 1;
                    unsigned start = get<0>(sites_tuples[start_index]) - 1;
                    assert (start >= 0);
                    assert (_data->_original_data_matrix.size() > 0);
                    unsigned end = get<1>(sites_tuples[start_index]); // include last site
                    for (unsigned m = start; m<end; m++) {
                        if (_data->_original_data_matrix[i][m] < 15 && _data->_original_data_matrix[j][m] < 15) {// 15 is ambiguity?
                            if (_data->_original_data_matrix[i][m] != _data->_original_data_matrix[j][m]) {
                                ndiff++;
                            }
                            ntotal++;
                        }
                    }
                        
                    assert (ntotal > 0);
                    
                    double p = ndiff / ntotal;
                    
                    if (p >= 0.75) {
                        p = 0.7499;
                    }
                    
                    // TODO: if p > 0.75, this will cause a crash - v will be NaN
                    // TODO: for now, just reset p to 0.7499
                    
                    double v = -0.75 * log(1 - 4.0/3.0 * p);
                    
                    assert (v != G::_infinity);
                    
                    unsigned k = i*(i-1)/2 + j;
                    dij[k] = v;
                    dij_row_col[k] = make_pair(i,j);
                    
                    assert (v == v);
                }
            }
            
            _starting_dij = dij;
            
            for (auto &d:_starting_dij) {
                assert (d == d);
            }
            
//            debugShowDistanceMatrix(_starting_dij);
        }
        
        else {
            _data->_original_data_matrix.clear();
            // Get the number of patterns
            unsigned npatterns = _data->getNumPatternsInSubset(_index - 1); // forest index starts at 1 but subsets start at 0

            // Get the first and last pattern index for this gene's data
            Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_index - 1);
            unsigned first_pattern = be.first;
            
            // Get the name of the gene (data subset)
            //string gene_name = _data->getSubsetName(_gene_index);

            // Get pattern counts
            auto counts = _data->getPatternCounts();
            
            // Create vectors to store products of same-state and different-state partials
            vector<double> same_state(npatterns, 0.0);
            vector<double> diff_state(npatterns, 0.0);
            
            // Create a map relating position in dij vector to row,col in distance matrix
            map<unsigned, pair<unsigned, unsigned>> dij_row_col;
            
            // Create distance matrix dij and workspace dij2 used to build next dij
            // Both dij and dij2 are 1-dimensional vectors that store only the
            // lower diagonal of the distance matrix (excluding diagonal elements)
            unsigned n = (unsigned)_lineages.size();
            vector<double> dij(n*(n-1)/2, G::_infinity);
            vector<double> dij2;

            // Calculate distances between all pairs of lineages
            
            for (unsigned i = 1; i < n; i++) {
                for (unsigned j = 0; j < i; j++) {
                    Node * lnode = _lineages[i];
                    Node * rnode = _lineages[j];
                    
                    // Fill same_state and diff_state vectors
                    same_state.assign(npatterns, 0.0);
                    diff_state.assign(npatterns, 0.0);
                    for (unsigned p = 0; p < npatterns; p++) {
                        for (unsigned lstate = 0; lstate < G::_nstates; lstate++) {
                            auto & l_partial_array = *(lnode->_partial);
                            double lpartial = l_partial_array[p*G::_nstates + lstate];
                            for (unsigned rstate = 0; rstate < G::_nstates; rstate++) {
                                auto & r_partial_array = *(rnode->_partial);
                                double rpartial = r_partial_array[p*G::_nstates + rstate];
                                if (lstate == rstate)
                                    same_state[p] += lpartial*rpartial;
                                else
                                    diff_state[p] += lpartial*rpartial;
                            }
                        }
                    }
                    
                    double min_dist = 0.0;
                    double max_dist = min_dist + 5.0; //TODO: replace arbitrary value 5.0
                    
                    double v0 = 0.0; // don't need to get edge lengths since we are starting from the trivial forest
                    
                    // TODO: what is this doing with missing data?
                    
                    negLogLikeDist f(npatterns, first_pattern, counts, same_state, diff_state, v0);
                    auto r = boost::math::tools::brent_find_minima(f, min_dist, max_dist, std::numeric_limits<double>::digits);
    //                double maximized_log_likelihood = -r.second;
                    unsigned k = i*(i-1)/2 + j;
                    dij[k] = r.first;
                    dij_row_col[k] = make_pair(i,j);
                                    
        //                output(format("d[%d, %d] = %.7f (logL = %.5f") % i % j % d[ij] % maximized_log_likelihood, 1);
                
                }
            }
            _starting_dij = dij;
            
//                debugShowDistanceMatrix(_starting_dij);
            
            for (auto &d:_starting_dij) {
                assert (d == d);
            }
        }
        
    }
#endif
        
#if defined (OLD_UPGMA)
    inline void Forest::buildStartingRow() {
        unsigned n = G::_ntaxa;
        map<Node *, unsigned> row;
        _upgma_starting_edgelen.clear();
        for (unsigned i = 0; i < n; i++) {
            Node * nd = _lineages[i];
            _upgma_starting_edgelen[nd] = nd->_edge_length;
            row[nd] = i;
        }
        _starting_row = row;
    }
#endif

    inline void Forest::destroyUPGMA() {
        if (G::_save_memory) {
             for (auto &nd:_nodes) {
                 nd._partial=nullptr;
             }
         }
        
        while (!_upgma_additions.empty()) {
            Node * parent = _upgma_additions.top();
            Node * child1 = parent->_left_child;
            Node * child2 = parent->_left_child->_right_sib;
            
            assert(child1);
            assert(child2);
            
            
            revertNodeVector(_lineages, child1, child2, parent);

            //reset siblings and parents of original nodes back to 0
            child1->resetNode(); //subtree1
            child2->resetNode(); //subtree2

            if (G::_save_memory) {
                child1->_partial = nullptr;
                child2->_partial = nullptr;
            }

            // clear new node from _nodes
            //clear new node that was just created
            parent->clear(); //new_nd
            
            _upgma_additions.pop();
            _nodes.pop_back(); // remove unused node from node list
        
            _ninternals--;
        }
        
        // Restore starting edge lengths
        for (auto nd : _lineages) {
            nd->_edge_length = _upgma_starting_edgelen.at(nd);
        }
        
        _upgma_starting_edgelen.clear();
        
    #if defined(DEBUG_UPGMA)
        output("\nIn GeneForest::destroyUPGMA:\n");
        output(format("  Height before refreshAllHeightsAndPreorders = %g\n") % _forest_height);
        refreshAllHeightsAndPreorders();
        output(format("  newick = %s\n") % makeNewick(9, /*use_names*/true, /*coalunits*/false));
        output(format("  Height after refreshAllHeightsAndPreorders = %g\n") % _forest_height);
        output("\n");
    #endif
//        showForest();
        
        if (_lineages.size() == 1) {
            _lineages.back()->_partial = nullptr; // last step, only node with partials should be the new node, and partials are no longer needed if likelihood has been calculated
        }
    }

    inline void Forest::constructUPGMA() {
//        unsigned partial_arrays = 0;

//        double prev_log_likelihood = _gene_tree_log_likelihood;
//        showForest();
        assert(_index >= 0);

        // Create two maps:
        //   "row_of_node" (key = node in _lineages, value = index of row of _dmatrix)
        //   "node_for_row" (key = index of row of _dmatrix, value = node in _lineages)
        // Also save starting edge lengths so they can be restored in destroyUPGMA()
        map<Node *, unsigned> row_of_node;
        map<unsigned, Node *> node_for_row;
        _upgma_starting_edgelen.clear();
        unsigned n = (unsigned)_lineages.size();
        for (unsigned i = 0; i < n; i++) {
            Node * nd = _lineages[i];
            
            _upgma_starting_edgelen[nd] = nd->_edge_length;

            // Find index of row corresponding to each lineage
            auto it1 = find(_dmatrix_rows.begin(), _dmatrix_rows.end(), nd->_split);
            assert(it1 != _dmatrix_rows.end());
            unsigned row_index = (unsigned)std::distance(_dmatrix_rows.begin(), it1);
            row_of_node[nd] = row_index;
            node_for_row[row_index] = nd;
        }
        
        // Build UPGMA tree on top of existing forest
        
        // Begin by making copies of the original distance matrix and row splits vector
        vector<double> dij     = _dmatrix;
        vector<Split>  dijrows = _dmatrix_rows;
        
        assert(_upgma_additions.empty());
        unsigned nsteps = n - 1;
        double h0 = _forest_height;
        
        while (nsteps > 0) {
            // Find smallest entry in dij
            
            auto it = min_element(dij.begin(), dij.end());
            unsigned k = (unsigned)std::distance(dij.begin(), it);
            
            // Use quadratic formula to figure out i. Note that k = i*(i-1)/2 + j, so solve
            // for a*i^2 + b*i + c = 0, where a=1, b=-1, and c=-2k, ignoring the "minus" root.
            //   i = floor{[-b + sqrt(b^2 - 4ac)]/(2a)}
            //     = floor{[1 + sqrt(1 + 8*k)]/2}
            //   j = k - i*(i-1)/2
            unsigned i = (unsigned)(floor((1.0 + sqrt(1.0 + 8.0*k))/2.0));
            unsigned j = k - i*(i-1)/2;
            
            // Update all leading edge lengths
            // If current gene forest looks like this:
            //
            //  A   B   C   D   E   F   G   H   I   J
            //  |   |   |   |   |   |   |   |   |   |
            //  |   +-+-+   |   |   |   |   |   |   |
            //  |     |     |   |   |   +-+-+   |   |
            //  |     +--+--+   |   |     |     |   | <-- height = h0
            //
            // and we've decided to join E and F next, with v equal to the distance
            // between E and F, then the amount to add to each leading edge is
            //   dh = max(0.0, v/2 - h0)
            //
            double v = *it;
            double dh = 0.5*v - h0;
            if (dh < 0.0)
                dh = 0.0;
            h0 += dh;
            for (auto nd : _lineages) {
                nd->_edge_length += dh;
            }
            
            //debugShowLineages();

            // Join lineages i and j
            Node nd;
            _nodes.push_back(nd);
            Node* new_nd = &_nodes.back();

            Node * subtree1 = node_for_row[i];
            Node * subtree2 = node_for_row[j];
            
            new_nd->_parent=0;
            new_nd->_number=_nleaves+_ninternals;
            new_nd->_right_sib=0;

            new_nd->_left_child=subtree1;
            subtree1->_right_sib=subtree2;

            subtree1->_parent=new_nd;
            subtree2->_parent=new_nd;
            
            new_nd->_split.resize(G::_ntaxa);
            new_nd->_split = subtree1->_split + subtree2->_split;
            
            _ninternals++;
            
            // Nodes added to _upgma_additions will be removed in destroyUPGMA()
            _upgma_additions.push(new_nd);
            
            updateNodeVector(_lineages, subtree1, subtree2, new_nd);
            
            row_of_node[new_nd] = i;

            assert (new_nd->_partial == nullptr);
            new_nd->_partial=ps.getPartial(_npatterns*4);
            assert(new_nd->_left_child->_right_sib);
            
//            StopWatch sw;
//            sw.start();
            G::_partial_arrays++;
            calcPartialArrayJC(new_nd); // TODO: use JC model for all UPGMA
//            double test = sw.stop();
//            G::_test += test;
            
            // Update distance matrix
            mergeDMatrixPair(dijrows, dij, subtree1->_split, subtree2->_split);
                     
            // Reset maps for next round
            unsigned n = (unsigned)_lineages.size();
            row_of_node.clear();
            node_for_row.clear();
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                                
                // Find index of row corresponding to each lineage
                auto it1 = find(dijrows.begin(), dijrows.end(), nd->_split);
                assert(it1 != dijrows.end());
                unsigned row_index = (unsigned)std::distance(dijrows.begin(), it1);
                row_of_node[nd] = row_index;
                node_for_row[row_index] = nd;
            }
            
            --nsteps;
        }
//        cout << "partial arrays = " << partial_arrays << endl;

        
//        _gene_tree_log_likelihood = calcLogLikelihood();
//        _log_weight = _gene_tree_log_likelihood - prev_log_likelihood; // previous likelihood is the entire tree
        
        // destroy upgma
//        destroyUPGMA();
//        showForest();
    }

    inline void Forest::mergeDMatrixPair(vector<Split> & dmatrows, vector<double> & dmatrix, Split & s1, Split & s2) {
        // dmatrix = [d10, d20, d21, d31, d31, d32]
        //   index =   0    1    2    3    4    5
        //
        // where the actual 2-dimensional matrix looks like this:
        //
        //   d00  d01  d02  d03
        //   d10  d11  d12  d13
        //   d20  d21  d22  d23
        //   d30  d31  d32  d33
        //
        // Only the lower triangle (not including diagonals) is used,
        // and elements are storwd by row.
        //
        // The index k of the (i,j)th element, where i > j, can be
        // obtained as follows:
        //
        // k = i*(i-1)/2 + j
        //
        // For example, i = 3, j = 1, k = 3*2/2 + 1 = 4
        //
        
        unsigned n = (unsigned)dmatrows.size();
        
        // Find index of s1 in dmatrows
        auto it1 = find(dmatrows.begin(), dmatrows.end(), s1);
        assert(it1 != dmatrows.end());
        unsigned i = (unsigned)std::distance(dmatrows.begin(), it1);
        
        // Find index of s2 in dmatrows
        auto it2 = find(dmatrows.begin(), dmatrows.end(), s2);
        assert(it2 != dmatrows.end());
        unsigned j = (unsigned)std::distance(dmatrows.begin(), it2);
        
        // Update distance matrix
        unsigned ij = (i > j) ? (i*(i-1)/2 + j) : (j*(j-1)/2 + i);
        dmatrix[ij] = G::_infinity;
        for (unsigned k = 0; k < n; k++) {
            if (k != i && k != j) {
                unsigned ik = (i > k) ? (i*(i-1)/2 + k) : (k*(k-1)/2 + i);
                unsigned jk = (j > k) ? (j*(j-1)/2 + k) : (k*(k-1)/2 + j);
                double a = dmatrix[ik];
                double b = dmatrix[jk];
                
                // Put average in cell with the smaller index
                if (i < j) {
                    dmatrix[ik] = 0.5*(a + b);
                    dmatrix[jk] = G::_infinity;
                }
                else {
                    dmatrix[ik] = G::_infinity;
                    dmatrix[jk] = 0.5*(a + b);
                }
            }
        }

        // Example (5x5 distance matrix):
        //   dmatrix [d10, d20, d21, d30, d31, d32, d40, d41, d42, d43]
        //   number of elements: 5*4/2 = 10
        //
        // Original (5x5):
        //
        //              0       1       2       3       4
        //            ----*   ---*-   --*--   -*---   *----
        //          +-------+-------+-------+-------+-------+
        // 0 ----*  |       |       |       |       |       |
        //          +-------+-------+-------+-------+-------+
        // 1 ---*-  |  d10  |       |       |       |       |
        //          +-------+-------+-------+-------+-------+
        // 2 --*--  |  d20  |  d21  |       |       |       |
        //          +-------+-------+-------+-------+-------+
        // 3 -*---  |  d30  |  d31  |  d32  |       |       |
        //          +-------+-------+-------+-------+-------+
        // 4 *----  |  d40  |  d41  |  d42  |  d43  |       |
        //          +-------+-------+-------+-------+-------+
        //
        // Intermediate during merging of 0 and 2 (4x4)
        //
        // Set d20 = inf and make these changes as well
        //
        //      k     ik     jk  dmatrix[ik]   dmatrix[jk]
        //  -----  -----  -----  -----------  ------------
        //  i = 0
        //      1    d10    d21  (d10+d21)/2           inf
        //  j = 2
        //      3    d30    d32  (d30+d32)/2           inf
        //      4    d40    d42  (d40+d42)/2           inf
        //
        //                  0           1       2       3       4
        //                ----*       ---*-   --*--   -*---   *----
        //          +---------------+-------+-------+-------+-------+
        // 0 ----*  |               |       |       |       |       | <-- kept
        //          +---------------+-------+-------+-------+-------+
        // 1 ---*-  |  (d10+d21)/2  |       |       |       |       |
        //          +---------------+-------+-------+-------+-------+
        // 2 --*--  |  inf          |  inf  |       |       |       | <-- elim
        //          +---------------+-------+-------+-------+-------+
        // 3 -*---  |  (d30+d32)/2  |  d31  |  inf  |       |       |
        //          +---------------+-------+-------+-------+-------+
        // 4 *----  |  (d40+d42)/2  |  d41  |  inf  |  d43  |       |
        //          +---------------+-------+-------+-------+-------+
        //                 kept                elim
        //
        // New matrix after eliminating row 2 and column 2 (3x3):
        //
        //                  0           1       2       3
        //                --*-*       ---*-   -*---   *----
        //          +---------------+-------+-------+-------+
        // 0 --*-*  |               |       |       |       |
        //          +---------------+-------+-------+-------+
        // 1 ---*-  |  (d10+d21)/2  |       |       |       |
        //          +---------------+-------+-------+-------+
        // 2 -*---  |  (d30+d32)/2  |  d31  |       |       |
        //          +---------------+-------+-------+-------+
        // 3 *----  |  (d40+d42)/2  |  d41  |  d43  |       |
        //          +---------------+-------+-------+-------+

        // Build new distance matrix
        unsigned n2 = n - 1;
        unsigned dim2 = n2*(n2-1)/2;
        vector<double> dmat2(dim2, G::_infinity);
        vector<Split> drows2(n2);
        
        // Eliminated row and column is the larger of i and j
        unsigned elim = i > j ? i : j;
        unsigned kept = i > j ? j : i;
        unsigned newrow = 0;
        unsigned newcol = 0;
        
        for (unsigned oldrow = 0; oldrow < n; oldrow++) {
            if (oldrow != elim) {
                if (oldrow == kept) {
                    // Set new row split to union of splits for rows i and j
                    drows2[newrow] = dmatrows[i] + dmatrows[j];
                }
                else {
                    drows2[newrow] = dmatrows[oldrow];
                }
                
                // Copy matrix row
                for (unsigned oldcol = 0; oldcol < oldrow; oldcol++) {
                    if (oldcol != elim) {
                        unsigned oldk = oldrow*(oldrow-1)/2 + oldcol;
                        unsigned newk = newrow*(newrow-1)/2 + newcol;
                        dmat2[newk] = dmatrix[oldk];
                        newcol++;
                    }
                }
                newrow++;
                newcol = 0;
            }
        }

        // Replace distance matrix with new one
        dmatrix = dmat2;
        dmatrows = drows2;
    }
    
#if defined (OLD_UPGMA)
    inline void Forest::buildRestOfTreeFaster() {
        
        if (_data->_original_data_matrix.size() > 0) {
            _data->_original_data_matrix.clear(); // if not using minimizer, can't clear this until all genes have gone through once
        }

    //            debugShowDistanceMatrix(_starting_dij);
        
        double prev_log_likelihood = _gene_tree_log_likelihood;
        
        vector<double> dij = _starting_dij;
        vector<double> dij2;
        
    //            debugShowDistanceMatrix(dij);
        
        // Create a map relating position in dij vector to row,col in distance matrix
        map<unsigned, pair<unsigned, unsigned>> dij_row_col;
        
                // Create distance matrix dij and workspace dij2 used to build next dij
                // Both dij and dij2 are 1-dimensional vectors that store only the
                // lower diagonal of the distance matrix (excluding diagonal elements)
        
        unsigned temp1 = _lineages.back()->_left_child->_right_sib->_position_in_lineages;
        unsigned temp2 = _lineages.back()->_left_child->_position_in_lineages;
        unsigned i_to_delete = temp1;
        unsigned j_to_delete = temp2; // i is larger number

        if (temp2 > temp1) {
            i_to_delete = temp2;
            j_to_delete = temp1;
        }

        
        Node* parent = _lineages.back();
        
        unsigned n = (unsigned) _lineages.size() + 1;
        
        _starting_row[parent] = i_to_delete;
        
            for (unsigned k = 0; k < n; k++) {
                if (k != i_to_delete && k != j_to_delete) {
                    unsigned ik = (i_to_delete > k) ? (i_to_delete*(i_to_delete-1)/2 + k) : (k*(k-1)/2 + i_to_delete);
                    unsigned jk = (j_to_delete > k) ? (j_to_delete*(j_to_delete-1)/2 + k) : (k*(k-1)/2 + j_to_delete);
                    double a = dij[ik];
                    double b = dij[jk];
                    dij[ik] = 0.5*(a + b);
                    dij[jk] = G::_infinity;
                }
            }
        
            // Build new distance matrix
             unsigned n2 = (unsigned)_lineages.size();
            assert(n2 == n - 1);
            unsigned dim2 = n2*(n2-1)/2;
            dij2.resize(dim2);
            dij2.assign(dim2, G::_infinity);
        
            // Calculate distances between all pairs of lineages
            dij_row_col.clear();
            for (unsigned i2 = 1; i2 < n2; i2++) {
                for (unsigned j2 = 0; j2 < i2; j2++) {
                    Node * lnode2 = _lineages[i2];
                    Node * rnode2 = _lineages[j2];
                    assert(_starting_row.find(lnode2) != _starting_row.end());
                    assert(_starting_row.find(rnode2) != _starting_row.end());
                    unsigned i = _starting_row[lnode2];
                    unsigned j = _starting_row[rnode2];
                    unsigned k2 = i2*(i2-1)/2 + j2;
                    unsigned k = i*(i-1)/2 + j;
                    if (j > i) {
                        k = j*(j-1)/2 + i;
                    }
                    dij2[k2] = dij[k];
                    dij_row_col[k2] = make_pair(i2,j2);
                    
                    assert (dij[k] == dij[k]);
                }
            }
        
        map<Node*, unsigned> row;
            dij = dij2;
            n = n2;
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                row[nd] = i;
            }

                // save starting distance matrix to reuse in next step
        _starting_dij = dij;
        
        // Create a map relating nodes in _lineages to rows of dij
        // Also save starting edge lengths so they can be restored in destroyUPGMA()
            row.clear();
        _upgma_starting_edgelen.clear();
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                _upgma_starting_edgelen[nd] = nd->_edge_length;
                row[nd] = i;
            }
            
            // Build UPGMA tree on top of existing forest
            assert(_upgma_additions.empty());
            
            double upgma_height = _lineages.back()->_height;
                        
        unsigned nsteps = n - 1;
    //    for (auto &d:dij) {
    //        cout << d << endl;
    //    }
        while (nsteps > 0) {
            // Find smallest entry in d
            auto it = min_element(dij.begin(), dij.end());
            unsigned offset = (unsigned)std::distance(dij.begin(), it);
            auto p = dij_row_col.at(offset);
            unsigned i = p.first;
            unsigned j = p.second;
            
            // Update all leading edge lengths
            double v = *it; // TODO: need to subtract existing node height?
            
            assert (v != G::_infinity);
            
            double edge_len_to_add = 0.5*v - upgma_height;
            if (edge_len_to_add <= 0.0) {
                edge_len_to_add = _small_enough; // TODO: to avoid likelihood issues, set v to very small if 0
                v = _small_enough;
            }
            
            assert (edge_len_to_add > 0.0);
            assert (v == v); // check v is not NaN
            for (auto nd : _lineages) {
        //                nd->_edge_length += (0.5*v - upgma_height);
                nd->_edge_length += edge_len_to_add;
            }
            
            upgma_height += edge_len_to_add;
            
            //debugShowLineages();
            
            // Join lineages i and j
            Node nd;
            _nodes.push_back(nd);
            Node* new_nd = &_nodes.back();

            Node * subtree1 = _lineages[i];
            Node * subtree2 = _lineages[j];
            
            new_nd->_parent=0;
            new_nd->_number=_nleaves+_ninternals;
            new_nd->_right_sib=0;

            new_nd->_left_child=subtree1;
            subtree1->_right_sib=subtree2;

            subtree1->_parent=new_nd;
            subtree2->_parent=new_nd;
            
            _ninternals++;
            
            // Nodes added to _upgma_additions will be removed in destroyUPGMA()
            _upgma_additions.push(new_nd);
            
            // Remove lnode and rnode from _lineages and add anc at the end
            updateNodeVector(_lineages, subtree1, subtree2, new_nd);
            row[new_nd] = i;
                        
            //debugShowLineages();
            // output(format("\nJoining lineages %d and %d\n") % i % j, 0);

            assert (new_nd->_partial == nullptr);
            new_nd->_partial=ps.getPartial(_npatterns*4);
            assert(new_nd->_left_child->_right_sib);
    //        sw.start();
            calcPartialArrayJC(new_nd); // TODO: for now, always complete UPGMA with JC
    //        double test = sw.stop();
    //        G::_test += test;
                        
            // Update distance matrix
            for (unsigned k = 0; k < n; k++) {
                if (k != i && k != j) {
                    unsigned ik = (i > k) ? (i*(i-1)/2 + k) : (k*(k-1)/2 + i);
                    unsigned jk = (j > k) ? (j*(j-1)/2 + k) : (k*(k-1)/2 + j);
                    double a = dij[ik];
                    double b = dij[jk];
                    dij[ik] = 0.5*(a + b);
                    dij[jk] = G::_infinity;
                }
            }
            
            // Sanity check
            for (auto nd : _lineages) {
                assert(!nd->_right_sib);
                assert(!nd->_parent);
            }
            
            // Build new distance matrix
            unsigned n2 = (unsigned)_lineages.size();
            assert(n2 == n - 1);
            unsigned dim2 = n2*(n2-1)/2;
            dij2.resize(dim2);
            dij2.assign(dim2, G::_infinity);
            
            // Calculate distances between all pairs of lineages
            dij_row_col.clear();
            for (unsigned i2 = 1; i2 < n2; i2++) {
                for (unsigned j2 = 0; j2 < i2; j2++) {
                    Node * lnode2 = _lineages[i2];
                    Node * rnode2 = _lineages[j2];
                    unsigned i = row[lnode2];
                    unsigned j = row[rnode2];
                    unsigned k2 = i2*(i2-1)/2 + j2;
                    unsigned k = i*(i-1)/2 + j;
                    if (j > i) {
                        k = j*(j-1)/2 + i;
                    }
                    dij2[k2] = dij[k];
                    dij_row_col[k2] = make_pair(i2,j2);
                }
            }
                
            // Set up for next iteration
            dij = dij2;
            n = n2;
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                row[nd] = i;
            }
            
            --nsteps;
        }

        // debugging output
        // output(format("\nGene forest for locus \"%s\" after UPGMA:\n%s\n") % gene_name % makeNewick(9, /*use_names*/true, /*coalunits*/false), 0);
        // output(format("  Height after UPGMA = %g\n") % _forest_height, 0);
            
        _gene_tree_log_likelihood = calcLogLikelihood();
        _log_weight = _gene_tree_log_likelihood - prev_log_likelihood; // previous likelihood is the entire tree
            
    //                showForest();
            
        if (G::_save_memory) {
            for (auto &nd:_nodes) {
                nd._partial=nullptr;
            }
        }

        // destroy upgma
        while (!_upgma_additions.empty()) {
            Node * parent = _upgma_additions.top();
            Node* child1 = parent->_left_child;
            Node* child2 = parent->_left_child->_right_sib;
            
            assert(child1);
            assert(child2);
            
            revertNodeVector(_lineages, child1, child2, parent);

            //reset siblings and parents of original nodes back to 0
            child1->resetNode(); //subtree1
            child2->resetNode(); //subtree2
            
            if (G::_save_memory) {
                child1->_partial = nullptr;
                child2->_partial = nullptr;
            }

            // clear new node from _nodes
            //clear new node that was just created
            parent->clear(); //new_nd

            _upgma_additions.pop();
            _nodes.pop_back(); // remove unused node from node list
            
            _ninternals--;
        }

        // Restore starting edge lengths
        for (auto nd : _lineages) {
            nd->_edge_length = _upgma_starting_edgelen.at(nd);
        }

        _upgma_starting_edgelen.clear();
            
            n = (unsigned) _lineages.size();
            
            _starting_row.clear();
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                _starting_row[nd] = i; // TODO: can save this earlier to not remake it
            }
            
            if (_lineages.size() == 1) {
                _lineages.back()->_partial = nullptr; // last step, only node with partials should be the new node, and partials are no longer needed if likelihood has been calculated
            }
        }
#endif
        
    inline void Forest::debugShowDistanceMatrix(const vector<double> & d) const {
        // d is a 1-dimensional vector that stores the lower triangle of a square matrix
        // (not including diagonals) in row order
        //
        // For example, for a 4x4 matrix (- means non-applicable):
        //
        //       0  1  2  3
        //     +-----------
        //  0  | -  -  -  -
        //  1  | 0  -  -  -
        //  2  | 1  2  -  -
        //  3  | 3  4  5  -
        //
        // For this example, d = {0, 1, 2, 3, 4 ,5}
        //
        // See this explanation for how to index d:
        //   https://math.stackexchange.com/questions/646117/how-to-find-a-function-mapping-matrix-indices
        //
        // In short, d[k] is the (i,j)th element, where k = i(i-1)/2 + j
        //       i   j   k = i*(i-1)/2 + j
        //       1   0   0 = 1*0/2 + 0
        //       2   0   1 = 2*1/2 + 0
        //       2   1   2 = 2*1/2 + 1
        //       3   0   3 = 3*2/2 + 0
        //       3   1   4 = 3*2/2 + 1
        //       3   2   5 = 3*2/2 + 2
        //
        // Number of elements in d is n(n-1)/2
        // Solving for n, and letting x = d.size(),
        //  x = n(n-1)/2
        //  2x = n^2 - n
        //  0 = a n^2 + b n + c, where a = 1, b = -1, c = -2x
        //  n = (-b += sqrt(b^2 - 4ac))/(2a)
        //    = (1 + sqrt(1 + 8x))/2
        double x = (double)d.size();
        double dbln = (1.0 + sqrt(1.0 + 8.0*x))/2.0;
        unsigned n = (unsigned)dbln;
        
        cout << format("\nDistance matrix (%d x %d):\n") % n % n;

        // Column headers
        cout << format("%12d") % " ";
        for (unsigned j = 0; j < n; j++) {
            cout << format("%12d") % j;
        }
        cout << "\n";
        
        unsigned k = 0;
        for (unsigned i = 0; i < n; i++) {
            cout << format("%12d") % i;
            for (unsigned j = 0; j < n; j++) {
                if (j < i) {
                    double v = d[k++];
                    if (v == G::_infinity)
                        cout << "         inf";
                    else
                        cout << format("%12.5f") % v;
                }
                else {
                    cout << "         inf";
                }
            }
            cout << "\n";
        }
        cout << "\n";
    }

    inline pair<Node*, Node*> Forest::chooseAllPairs(list<Node*> &node_list, double increment, string species, Lot::SharedPtr lot) {
          double prev_log_likelihood = _gene_tree_log_likelihood;
    //         _node_choices.clear();
          assert (_node_choices.size() == 0);
    //         _log_likelihood_choices.clear();
          assert (_log_likelihood_choices.size() == 0);
           _log_weight = 0.0;
          
           // choose pair of nodes to try
           for (unsigned i = 0; i < node_list.size()-1; i++) {
               for (unsigned j = i+1; j < node_list.size(); j++) {
                   // createNewSubtree returns subtree1, subtree2, new_nd
                   
                   tuple<Node*, Node*, Node*> t = createNewSubtree(make_pair(i,j), node_list, increment, species);
                   
                   _log_likelihood_choices.push_back(calcLogLikelihood());
                   // gene tree log coalescent likelihood is the same for every possible join

                   // revert _lineages if > 1 choice
                   revertNodeVector(_lineages, get<0>(t), get<1>(t), get<2>(t));

                   //reset siblings and parents of original nodes back to 0
                   get<0>(t)->resetNode(); //subtree1
                   get<1>(t)->resetNode(); //subtree2

                   // clear new node from _nodes
                   //clear new node that was just created
                   get<2>(t)->clear(); //new_nd
    //                 _nodes.pop_back();
               }
           }
           
           // reweight each choice of pairs
          vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, prev_log_likelihood);


           // sum unnormalized weights before choosing the pair
           // must include the likelihoods of all pairs in the final particle weight
           double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
           _log_weight = log_weight_choices_sum;
           for (unsigned b=0; b < log_weight_choices.size(); b++) {
               log_weight_choices[b] -= log_weight_choices_sum;
           }
          
           // randomly select a pair
           unsigned index_of_choice = selectPair(log_weight_choices, lot);

           // find nodes to join in node_list
           Node* subtree1 = _node_choices[index_of_choice].first;
           Node* subtree2 = _node_choices[index_of_choice].second;
       
           _gene_tree_log_likelihood = _log_likelihood_choices[index_of_choice]; // reset the log likelihood
          
           // erase extra nodes created from node list
           for (unsigned i = 0; i < _node_choices.size(); i++) {
               _nodes.pop_back();
           }
          
          _node_choices.clear();
          _log_likelihood_choices.clear();
           return make_pair(subtree1, subtree2);
       }

    inline pair<Node*, Node*> Forest::getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list) {
          Node *subtree1 = nullptr;
          Node *subtree2 = nullptr;

          unsigned a = 0;
          for (auto iter=node_list.begin(); iter != node_list.end(); iter++){
              if (a==t.first) {
                  subtree1 = *iter;
              }
              else if (a==t.second) {
                  subtree2 = *iter;
              }
              if (subtree1 && subtree2) {
                  break;
              }
              a++;
          }

          pair<Node*, Node*> s = make_pair(subtree1, subtree2);
          _node_choices.push_back(make_pair(subtree1, subtree2));

          return s;
      }


    inline tuple<Node*, Node*, Node*> Forest::createNewSubtree(pair<unsigned, unsigned> t, list<Node*> node_list, double increment, string species) {
         pair<Node*, Node*> p = getSubtreeAt(t, node_list);

         Node* subtree1 = p.first;
         Node* subtree2 = p.second;

    //        new node is always needed
         Node nd;
         _nodes.push_back(nd);
         Node* new_nd = &_nodes.back();
         new_nd->_parent=0;
         new_nd->_number=_nleaves+_ninternals;
         new_nd->_edge_length=0.0;
         new_nd->_right_sib=0;

         new_nd->_left_child=subtree1;
         subtree1->_right_sib=subtree2;

         subtree1->_parent=new_nd;
         subtree2->_parent=new_nd;

         //always calculating partials now
         assert (new_nd->_partial == nullptr);
         new_nd->_partial=ps.getPartial(_npatterns*4);
         assert(new_nd->_left_child->_right_sib);
        if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
            calcPartialArrayJC(new_nd);
        }
        else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
            calcPartialArrayHKY(new_nd);
        }
        else {
            throw XProj("model must be either JC or HKY");
        }
        new_nd->_height = _forest_height;

         // don't update the species list
         updateNodeVector(_lineages, subtree1, subtree2, new_nd);
                  
         return make_tuple(subtree1, subtree2, new_nd);
     }

    inline void Forest::allowCoalescence(string species_name, double increment, Lot::SharedPtr lot) {
         double prev_log_likelihood = _gene_tree_log_likelihood;

        // TODO: BE CAREFUL
//        for (auto &nd:_lineages) {
//            nd->_edge_length = 0.003343;
//        }
         Node *subtree1 = nullptr;
         Node *subtree2 = nullptr;
         list<Node*> nodes = _species_partition[species_name];
        
        assert (nodes.size() > 0);

         unsigned s = (unsigned) nodes.size();
         calcTopologyPrior(s);

         assert (s > 1);
         bool one_choice = false;
         if (nodes.size() == 2) {
             one_choice = true;
         }

         if (G::_proposal == "prior-post" && (!one_choice)) {
             if (G::_save_memory) {
                 for (auto &nd:_lineages) {
                     if (nd->_partial == nullptr) {
                         nd->_partial = ps.getPartial(_npatterns*4);
                         if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                             calcPartialArrayJC(nd);
                         }
                         else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
                             calcPartialArrayHKY(nd);
                         }
                         else {
                             throw XProj("model must be either HKY or JC");
                         }
                     }
                 }
             }
             
             pair<Node*, Node*> t = chooseAllPairs(nodes, increment, species_name, lot);
             
             subtree1 = t.first;
             subtree2 = t.second;
         }
         
         else {
             assert (G::_proposal == "prior-prior" || one_choice);
             // prior-prior proposal
             pair<unsigned, unsigned> t = chooseTaxaToJoin(s, lot);
             auto it1 = std::next(nodes.begin(), t.first);
             subtree1 = *it1;

             auto it2 = std::next(nodes.begin(), t.second);
             subtree2 = *it2;
             assert (t.first < nodes.size());
             assert (t.second < nodes.size());

//             subtree1 = _lineages[7];
//             subtree2 = _lineages[2]; // TODO: BE CAREFUL
             // TODO: BE CAREFUL
//             _forest_height = _lineages[0]->_edge_length;
             assert (subtree1 != subtree2);
         }

         //new node is always needed
         Node nd;
         _nodes.push_back(nd);
         Node* new_nd = &_nodes.back();

         new_nd->_parent=0;
         new_nd->_number=_nleaves+_ninternals;
         new_nd->_edge_length=0.0;
         _ninternals++;
         new_nd->_right_sib=0;
        new_nd->_name = to_string(new_nd->_number);

         new_nd->_left_child=subtree1;
         subtree1->_right_sib=subtree2;

         subtree1->_parent=new_nd;
         subtree2->_parent=new_nd;

         if (!G::_run_on_empty) {
             //always calculating partials now
             assert (new_nd->_partial == nullptr);
             new_nd->_partial=ps.getPartial(_npatterns*4);
             assert(new_nd->_left_child->_right_sib);

             if (G::_save_memory) {
                 for (auto &nd:_lineages) {
                     if (nd->_partial == nullptr) {
                         nd->_partial = ps.getPartial(_npatterns*4);
                         if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                             calcPartialArrayJC(nd);
                         }
                         else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
                             calcPartialArrayHKY(nd);
                         }
                         else {
                             throw XProj("model specified must be JC or HKY");
                         }
                     }
                 }
             }
             
             if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                 calcPartialArrayJC(new_nd);
             }
             else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
                 calcPartialArrayHKY(new_nd);
             }
             else {
                 throw XProj("model specified must be JC or HKY");
             }
             new_nd->_height = _forest_height;
             
    #if defined (FASTER_SECOND_LEVEL)
             new_nd->_species = 0;
    #endif

             subtree1->_partial=nullptr; // throw away subtree partials now, no longer needed
             subtree2->_partial=nullptr;
             
             if (G::_upgma) {
                // Update _dmatrix and _dmatrix_rows
                Split s1 = subtree1->getSplit();
                Split s2 = subtree2->getSplit();
                
                //output(format("s1 = %s\n") % s1.createPatternRepresentation());
                //output(format("s2 = %s\n") % s2.createPatternRepresentation());
                //output("_dmatrix_rows:\n");
                //unsigned z = 0;
                //for (auto s : _dmatrix_rows) {
                //    output(format("%6d = %s\n") % z++ % s.createPatternRepresentation());
                //}
                
                mergeDMatrixPair(_dmatrix_rows, _dmatrix, s1, s2);
             }

         }
        
        //update species list
         updateNodeList(nodes, subtree1, subtree2, new_nd); // TODO: BE CAREFUL
         updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        if (G::_upgma) {
            new_nd->_split.resize(G::_ntaxa);
            new_nd->_split += subtree1->_split;
            new_nd->_split += subtree2->_split;
        }

        _species_partition[species_name] = nodes;

//        if (!G::_upgma) { // TODO: be careful
             if ((G::_proposal == "prior-prior" || one_choice) && (!G::_run_on_empty) ) {
                 _gene_tree_log_likelihood = calcLogLikelihood();
                 _log_weight = _gene_tree_log_likelihood - prev_log_likelihood;
             }
        
            if (G::_save_memory) {
                for (auto &nd:_nodes) {
                    nd._partial=nullptr;
                }
            }
//        }
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

    inline void Forest::updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode) {
    //        showForest();
        // Delete delnode1 from node_vector
        auto it1 = find(node_vector.begin(), node_vector.end(), delnode1);
    //        if (it1 == node_vector.end()) {
    //            cout << "node to delete 1 has position " << delnode1->_position_in_lineages << " and name " << delnode1->_name << endl;
    //            cout << delnode1->_name << " position " << delnode1->_position_in_lineages << endl;
    //            cout << "\tedge len " << delnode1->_edge_length << endl;
    //            if (delnode1->_parent) {
    //                cout << "\tparent " << delnode1->_parent << endl;
    //            }
    //            cout << "node to delete 2 has position " << delnode2->_position_in_lineages << " and name " << delnode2->_name << endl;
    //            cout << delnode2->_name << " position " << delnode2->_position_in_lineages << endl;
    //            cout << "\tedge len " << delnode2->_edge_length << endl;
    //            if (delnode2->_parent) {
    //                cout << "\tparent " << delnode2->_parent << endl;
    //            }
    //            cout << "node to add has position " << addnode->_position_in_lineages << " and name " << addnode->_name << endl;
    //            cout << "lineages is: ";
    //            for (auto &nd:_lineages) {
    //                cout << nd->_name << " position " << nd->_position_in_lineages << endl;
    //                cout << "\tedge len " << nd->_edge_length << endl;
    //                if (nd->_parent) {
    //                    cout << "\tparent " << nd->_parent << endl;
    //                }
    //            }
    //
    //            cout << "nodes is: ";
    //            for (auto &nd:_nodes) {
    //                cout << nd._name << " position " << nd._position_in_lineages << endl;
    //                cout << "\tedge len " << nd._edge_length << endl;
    //                if (nd._parent) {
    //                    cout << "\tparent " << nd._parent << endl;
    //                }
    //            }
    //        }
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

    inline void Forest::revertNodeVector(vector<Node *> &node_vector, Node *addnode1, Node *addnode2, Node *delnode1) {
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

    inline double Forest::getTreeLength() {
        // sum of all edge lengths in tree
        double sum_height = 0.0;
        
        for (auto &nd:_nodes) {
            // sum edge lengths from all nodes
            sum_height += nd._edge_length;
        }
        return sum_height;
    }

    inline vector<pair<double, string>> Forest::calcForestRate(Lot::SharedPtr lot) {
        vector<pair<double, string>> rates;
        pair<double, string> rate_and_name;

        for (auto &s:_species_partition) {
            if (s.second.size() > 1) { // if size == 0, no possibility of coalescence and rate is 0
                double population_coalescence_rate = 0.0;
    #if defined (DRAW_NEW_THETA)
                    double population_theta = _theta_map[s.first];
                population_coalescence_rate = s.second.size()*(s.second.size()-1)/population_theta;
    #else
                population_coalescence_rate = s.second.size()*(s.second.size()-1)/G::_theta;
    #endif
                string name = s.first;
                rate_and_name = make_pair(population_coalescence_rate, name);
                rates.push_back(rate_and_name);
            }
        }
        return rates;
    }

    inline void Forest::addIncrement(double increment) {
        for (auto &nd:_lineages) {
            nd->_edge_length += increment;
        }
        _last_edge_length = increment;
        
        _forest_height += _last_edge_length;
    }

    inline void Forest::resetThetaMap(Lot::SharedPtr lot) { // TODO: not sure if this works if not doing jones coalescent likelihood - double check
        assert (_theta_map.size() == 0);
        // map should be 2*nspecies - 1 size
        unsigned number = 0;
        vector<string> species_names;
        
        for (auto &s:_species_partition) {
            species_names.push_back(s.first);
            number++;
        }
        
        for (int i=0; i<G::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            species_names.push_back(name);
        }
        
        _ancestral_species_name = species_names.back();
        
        assert (species_names.size() == 2*G::_nspecies - 1);
        
        // draw thetas for tips of species trees and ancestral population
        // for all other populations, theta = -1
        
        if (G::_theta_proposal_mean == 0.0) {
            assert (G::_theta > 0.0);
            G::_theta_proposal_mean = G::_theta;
        }
        double scale = 1 / G::_theta_proposal_mean;
        
        unsigned count = 0;
        for (auto &name:species_names) {
            if (count < G::_nspecies || count == 2*G::_nspecies-2) {
                double new_theta = 0.0;
                if (new_theta < G::_small_enough) {
                    new_theta = 1 / (lot->gamma(2.0, scale));
                    assert (new_theta > 0.0);
                    _theta_map[name] = new_theta;
                }
                // pop mean = theta / 4
                double a = 2.0;
                double b = scale;
                double x = new_theta;
                double log_inv_gamma_prior = (a*log(b) - lgamma(a) - (a+1)*log(x) - b/x);
                _vector_prior.push_back(log_inv_gamma_prior);
            }
            else {
                _theta_map[name] = -1;
            }
            count++;
        }
    }

    inline void Forest::drawNewTheta(string new_species, Lot::SharedPtr lot) {
        // draw a new theta for the newest species population
        double scale = 1 / _theta_mean;
        double new_theta = 0.0;
        if (new_theta < G::_small_enough) {
            new_theta = 1 / lot->gamma(2.0, scale);
            _theta_map[new_species] = new_theta;
        }
        // pop mean = theta / 4
        double a = 2.0;
        double b = scale;
        double x = new_theta;
        double log_inv_gamma_prior = (a*log(b) - lgamma(a) - (a+1)*log(x) - b/x);
        _vector_prior.push_back(log_inv_gamma_prior);
    }

    inline void Forest::updateThetaMap(Lot::SharedPtr lot, string new_species_name) {
        // add a new theta for the most recently drawn species
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        double new_theta = 0.0;
        if (new_theta < G::_small_enough) {
    //            new_theta = 1 / (lot->gamma(2.01, scale));
            new_theta = 1 / (lot->gamma(2.0, scale));
            assert (new_theta > 0.0);
            _theta_map[new_species_name] = new_theta;
        }
        // pop mean = theta / 4
        double a = 2.0;
        double b = scale;
        double x = new_theta; //  x is theta, not theta / 4 like it is for starbeast3
        double log_inv_gamma_prior = - 1 / (b*x) - (a + 1) * log(x) - a*log(b) - lgamma(a);
        _vector_prior.push_back(log_inv_gamma_prior);
    }
        
    inline void Forest::updateThetaMapFixedTheta(Lot::SharedPtr lot, string new_species_name) {
        _theta_map[new_species_name] = G::_theta;
    }
        
#if defined (OLD_UPGMA)
    inline void Forest::createSpeciesIndices() {
        unsigned number = 0;
        for (auto &s:_species_partition) {
            number++;
            _species_indices[s.first] = number - 1;
        }
        
        for (int i=0; i<G::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            _species_indices[name] = number - 1;
        }
    }
#endif
        
    inline void Forest::createThetaMapFixedTheta(Lot::SharedPtr lot) {
        // map should be 2*nspecies - 1 size
        unsigned number = 0;
        _species_names.clear();
        
        for (auto &s:_species_partition) {
            _species_names.push_back(s.first);
            number++;
#if defined (OLD_UPGMA)
            _species_indices[s.first] = number - 1;
#endif
        }
        assert (_species_names.size() == G::_nspecies);
        
        for (int i=0; i<G::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            _species_names.push_back(name);
#if defined (OLD_UPGMA)
            _species_indices[name] = number - 1;
#endif
        }
        
        _theta_mean = G::_theta;
        
        for (auto &name:_species_names) {
            _theta_map[name] = G::_theta;
        }
    }

    inline void Forest::createThetaMap(Lot::SharedPtr lot) {
        // map should be 2*nspecies - 1 size
        unsigned number = 0;
    //        vector<string> species_names;
        _species_names.clear();
        
        for (auto &s:_species_partition) {
            _species_names.push_back(s.first);
            number++;
#if defined (OLD_UPGMA)
            _species_indices[s.first] = number - 1;
#endif
        }
        assert (_species_names.size() == G::_nspecies);
        
        for (int i=0; i<G::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            _species_names.push_back(name);
#if defined (OLD_UPGMA)
            _species_indices[name] = number - 1;
#endif
        }
        
        // gamma mean = shape * scale
        // draw mean from lognormal distribution
        // shape = 2.0 to be consistent with starbeast3
        // scale = 1 / mean;
        
        if (G::_theta_proposal_mean > 0.0) {
            assert (_theta_mean == 0.0);
            _theta_mean = lot->gamma(1, G::_theta_proposal_mean); // equivalent to exponential(exponential_rate)
        }
        else {
            _theta_mean = G::_theta; // if no proposal distribution specified, use one theta mean for all particles
        }
        
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        for (auto &name:_species_names) {
            double new_theta = 0.0;
            if (new_theta < G::_small_enough) {
                new_theta = 1 / (lot->gamma(2.0, scale));
                assert (new_theta > 0.0);
                _theta_map[name] = new_theta;
            }
            // pop mean = theta / 4
            double a = 2.0;
            double b = scale;
            double x = new_theta; //  x is theta, not theta / 4 like it is for starbeast3
            double log_inv_gamma_prior = - 1 / (b*x) - (a + 1) * log(x) - a*log(b) - lgamma(a);
            _vector_prior.push_back(log_inv_gamma_prior);

        }
    }

#if !defined (FASTER_SECOND_LEVEL)
    inline double Forest::calcCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height) {
        _panmictic_coalescent_likelihood = 0.0;

        double neg_inf = -1*numeric_limits<double>::infinity();
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        double log_coalescent_likelihood = 0.0;
        double cum_time = 0.0;
        int a = 0;

        assert (heights_and_nodes.size() > 0);
        if (species_increment > 0) {
            for (unsigned i=_nincrements; i < heights_and_nodes.size(); i++) {
                Node* node = nullptr;
                if (heights_and_nodes[i].first < species_tree_height) {
                    // calc coalescent prob and update species partition
                    double increment = heights_and_nodes[i].first - (species_tree_height - species_increment) - cum_time;
                    cum_time += increment;
                    assert (increment > 0.0);
                    // find node in species partition and calculate nlineages
                    string species;
                    for (auto &s:_species_partition) {
                        for (auto &nd:s.second) {
                            if (nd == heights_and_nodes[i].second->_left_child) {
                                species = s.first; // TODO: can break out of this to be faster
                                node = nd;
                            }
                        }
                    }
                    for (auto &s:_species_partition) {
                        if (s.first == species) {
                            // coalescence
                            unsigned nlineages = (int) s.second.size();

                            if (nlineages == 1) {
                                log_coalescent_likelihood = neg_inf; // if there is only 1 lineage left, there cannot be coalescence
                                return log_coalescent_likelihood;
                            }
    #if defined (DRAW_NEW_THETA)
                            double population_theta = _theta_map[s.first];
                            assert (population_theta > 0.0);
                            double coalescence_rate = nlineages*(nlineages-1) / population_theta;
    #else
                            double coalescence_rate = nlineages*(nlineages-1) / G::_theta;
    #endif
                            double nChooseTwo = nlineages*(nlineages-1);
                            double log_prob_join = log(2/nChooseTwo);
                            log_coalescent_likelihood += log_prob_join + log(coalescence_rate) - (increment * coalescence_rate);

                            bool found = (find(s.second.begin(), s.second.end(), node->_right_sib) != s.second.end());

                            if (found) {
                                updateNodeList(s.second, node, node->_right_sib, node->_parent); // update the species lineage
                            }
                            else {
                                log_coalescent_likelihood = neg_inf;
                            }
                            a++;
                        }
                        else {
                            // no coalescence
                            unsigned nlineages = (int) s.second.size();

    #if defined (DRAW_NEW_THETA)
                            double population_theta = _theta_map[s.first];
                            assert (population_theta > 0.0);
                            double coalescence_rate = nlineages*(nlineages-1) / population_theta;
    #else
                            double coalescence_rate = nlineages*(nlineages-1) / G::_theta;
    #endif
                            log_coalescent_likelihood -= increment * coalescence_rate;
                        }
                    }

                }
            }
            double remaining_chunk_of_branch = species_increment - cum_time;
            // no coalescence
            for (auto &s:_species_partition) {
                unsigned nlineages = (int) s.second.size();
    #if defined (DRAW_NEW_THETA)
                double population_theta = _theta_map[s.first];
                assert (population_theta > 0.0);
                double coalescence_rate = nlineages*(nlineages-1) / population_theta;
    #else
                double coalescence_rate = nlineages*(nlineages-1) / G::_theta;
    #endif
                log_coalescent_likelihood -= remaining_chunk_of_branch * coalescence_rate;

    //                cout << "pr no coalescence = " << -1 * remaining_chunk_of_branch * coalescence_rate << endl;
            }
            _nincrements += a;
        }

        // calculate coalescent likelihood for the rest of the panmictic tree
        if (log_coalescent_likelihood != neg_inf) {
        if (species_increment > 0.0) {
            double cum_time = 0.0;
            // calculate height of each join in the gene tree, minus the species tree height
            // get number of lineages at beginning of that step and then decrement with each join
            // don't increment _nincrements?

            // start at _nincrements and go through heights and nodes list
            int panmictic_nlineages = 0;
            for (auto &s:_species_partition) {
                panmictic_nlineages += s.second.size();
            }

            for (unsigned i=_nincrements; i < heights_and_nodes.size(); i++) {
                // calculate increment (should be heights_and_nodes[i].first - species_tree_height - cum_time (no need to worry about species increment now)
                // calculate prob of coalescence, ln(rate) - rate*increment
                double increment = heights_and_nodes[i].first - species_tree_height - cum_time;
    #if defined (DRAW_NEW_THETA)
                double population_theta = _theta_map[_ancestral_species_name];
                assert (population_theta > 0.0);
                double coalescence_rate = panmictic_nlineages*(panmictic_nlineages-1) / population_theta;
    #else
                double coalescence_rate = panmictic_nlineages*(panmictic_nlineages-1) / G::_theta;
    #endif
                double nChooseTwo = panmictic_nlineages*(panmictic_nlineages-1);
                double log_prob_join = log(2/nChooseTwo);
                _panmictic_coalescent_likelihood += log_prob_join + log(coalescence_rate) - (increment * coalescence_rate);

                cum_time += increment;
                panmictic_nlineages--;
            }
            assert (panmictic_nlineages == 1);
        }

        else {
            // final step; no deep coalescence; one species
            assert (species_increment == 0.0);
            assert (_species_partition.size() == 1);
            for (unsigned i=_nincrements; i < heights_and_nodes.size(); i++) {
                // there must be coalescence at this point
                Node* node = nullptr;
                string species;
                for (auto &s:_species_partition) {
                    for (auto &nd:s.second) {
                        if (nd == heights_and_nodes[i].second->_left_child) {
                            node = nd;
                        }
                    }
                    for (auto &s:_species_partition) {
                        // coalescence
                        unsigned nlineages = (int) s.second.size();

    #if defined (DRAW_NEW_THETA)
                        double population_theta = _theta_map[s.first];
                        double coalescence_rate = nlineages*(nlineages-1) / population_theta;
    #else
                        double coalescence_rate = nlineages*(nlineages-1) / G::_theta;
    #endif
                        double nChooseTwo = nlineages*(nlineages-1);
                        double log_prob_join = log(2/nChooseTwo);
                        double increment = heights_and_nodes[i].first - species_tree_height - cum_time;
                        log_coalescent_likelihood += log_prob_join + log(coalescence_rate) - (increment * coalescence_rate);

                        updateNodeList(s.second, node, node->_right_sib, node->_parent); // update the species lineage
                        a++;
                        cum_time += increment;
                    }
                        // no deep coalescence to deal with
                    }
                }
        }
        }
        _log_coalescent_likelihood += log_coalescent_likelihood;
        return log_coalescent_likelihood;
    }
#endif

#if !defined (FASTER_SECOND_LEVEL)
    inline pair<vector<double>, vector<unsigned>> Forest::calcCoalescentLikelihoodIntegratingOutThetaLastStep(vector<pair<tuple<string,string,string>, double>> species_build) {
        vector<double> gamma_b; // contains gamma_b by species
        vector<unsigned> q_b; // contains q_b by species
        
        gamma_b.resize(G::_nspecies + species_build.size());
        q_b.resize(G::_nspecies + species_build.size());
        
        _panmictic_coalescent_likelihood = 0.0;
        _log_coalescent_likelihood = 0.0;
        
        double neg_inf = -1*numeric_limits<double>::infinity();
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        assert (heights_and_nodes.size() > 0);
        
        double log_coalescent_likelihood = 0.0;
        double cum_time = 0.0;
        string species;
        bool found = false;
        bool coalescence = false;
        unsigned nlineages = 0;
        unsigned a = 0;
        double species_tree_height = 0.0;
        
        unsigned count = -1;
        
        assert (_taxon_map.size() > 0);
        setUpGeneForest(_taxon_map); // TODO: don't have to redo this because it's the same at every step?
        
        for (unsigned gen=0; gen < species_build.size() - 1; gen++) { // don't go into last step because panmictic part will take care of that
            count = -1;
            tuple<string, string, string> species_joined = species_build[gen].first;
            double species_increment = species_build[gen].second;
            species_tree_height += species_increment;
            
            updateSpeciesPartition(species_joined);
        
            for (auto &s:_species_partition) { // TODO: how to check if gene tree violates species tree?
                if (gamma_b.size() > 0 && gamma_b[0] != neg_inf) {
                    count = _species_indices[s.first];
                    
                    bool done = false;
                    cum_time = 0.0;
                    species = s.first;
                    for (unsigned i=0; i < heights_and_nodes.size(); i++) {
                        if (!done) {
                            found = false;
                            coalescence = false;
                            Node*search_nd = heights_and_nodes[i].second->_left_child;
                            for (auto &nd:s.second) {
                                if (nd == search_nd) {
                                    search_nd = nd;
                                    nlineages = (unsigned) s.second.size();
                                    if (heights_and_nodes[i].first < species_tree_height) {
                                        coalescence = true;
                                        a++;
                                    }
                                    else {
                                        done = true; // stop once gene increment goes below species tree height to avoid double counting some deep coalescences
                                    }
                                    found = true;
                                    break;
                                }
                            }
                            if (found) {
                                if (species_increment == 0.0) {
                                    assert (coalescence); // no deep coalescence if last step
                                }
                                
                                if (coalescence) {
                                    if (s.second.size() == 1) {
                                        q_b.clear(); // don't push back to q_b because it is full of ints
                                        gamma_b.clear();
                                        gamma_b.push_back(neg_inf);
                                        done = true;
                                        break;
                                    }
                                    
                                    double increment = heights_and_nodes[i].first - (species_tree_height - species_increment) - cum_time; // find increment height
                                    cum_time += increment;
                                    assert (increment > 0.0);
                                    assert (q_b.size() > 0);
                                    assert (gamma_b.size() > 0);
                                    
                                    q_b[count] += 1;
                                    gamma_b[count] += 4* increment * nlineages *(nlineages-1) / (2*G::_ploidy);
                                    for (auto &s:_species_partition) {
                                        if (s.first == species) {
                                            bool found = (find(s.second.begin(), s.second.end(), search_nd->_right_sib) != s.second.end());
                                            if (found ) {
                                                updateNodeList(s.second, search_nd, search_nd->_right_sib, search_nd->_parent); // update the species lineage
                                            }
                                            else {
                                                q_b.clear(); // don't push back to q_b because it is full of ints
                                                gamma_b.clear();
                                                gamma_b.push_back(neg_inf);
                                                done = true;
                                                break;
                                            }
                                            break;
                                        }
                                    }
                                }
                                else {
                                    double increment = species_increment - cum_time; // find increment height - can't be any more deep coalescence, so gene increment is whatever remains of the species increment
                                    cum_time += increment;
                                    assert (increment > 0.0);
                                    
                                    assert (gamma_b.size() > 0);
                                    gamma_b[count] += 4*increment * nlineages *(nlineages-1) / (2*G::_ploidy);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if (gamma_b.back() != neg_inf) {
            assert (q_b.back() == 0.0);
            assert (gamma_b.back() == 0);
            
            // calculate coalescent likelihood for the rest of the panmictic tree
            double species_increment = species_build.back().second;
            if (log_coalescent_likelihood != neg_inf) {
                if (species_increment > 0.0) {
                    double cum_time = 0.0;
                    // calculate height of each join in the gene tree, minus the species tree height
                    // get number of lineages at beginning of that step and then decrement with each join

                    // start at _nincrements and go through heights and nodes list
                    int panmictic_nlineages = 0;
                    for (auto &s:_species_partition) {
                        panmictic_nlineages += s.second.size();
                    }

                    for (unsigned i=a; i < heights_and_nodes.size(); i++) {
                        // calculate increment (should be heights_and_nodes[i].first - species_tree_height - cum_time (no need to worry about species increment now)
                        // calculate prob of coalescence
                        double increment = heights_and_nodes[i].first - species_tree_height - cum_time;
                        q_b.back() += 1;
                        gamma_b.back() += 4*increment * panmictic_nlineages *(panmictic_nlineages-1) / (2*G::_ploidy);

                        cum_time += increment;
                        panmictic_nlineages--;
                    }
                    assert (panmictic_nlineages == 1);
                }
            }
            assert (gamma_b.size() == q_b.size());
        }
        
        return make_pair(gamma_b, q_b);
    }
#endif

#if !defined (FASTER_SECOND_LEVEL)
    inline pair<vector<double>, vector<unsigned>> Forest::calcInitialCoalescentLikelihoodIntegratingOutTheta() {
        vector<double> gamma_b; // contains gamma_b by species
        vector<unsigned> q_b; // contains q_b by species
        
        gamma_b.resize(1);
        q_b.resize(1);
        
        _panmictic_coalescent_likelihood = 0.0;
        _log_coalescent_likelihood = 0.0;
        
        double neg_inf = -1*numeric_limits<double>::infinity();
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        assert (heights_and_nodes.size() > 0);
        
        unsigned a = 0;
        double species_tree_height = 0.0;
        
        assert (_taxon_map.size() > 0);
        setUpGeneForest(_taxon_map);
        
        if (gamma_b.back() != neg_inf) {
            assert (q_b.back() == 0.0);
            assert (gamma_b.back() == 0);
            
        // calculate coalescent likelihood for the rest of the panmictic tree
        double cum_time = 0.0;
        // calculate height of each join in the gene tree, minus the species tree height
        // get number of lineages at beginning of that step and then decrement with each join

        // start at _nincrements and go through heights and nodes list
        int panmictic_nlineages = 0;
        for (auto &s:_species_partition) {
            panmictic_nlineages += s.second.size();
        }

        for (unsigned i=a; i < heights_and_nodes.size(); i++) {
            // calculate increment (should be heights_and_nodes[i].first - species_tree_height - cum_time (no need to worry about species increment now)
            // calculate prob of coalescence
            double increment = heights_and_nodes[i].first - species_tree_height - cum_time;
            q_b.back() += 1;
            gamma_b.back() += 4*increment * panmictic_nlineages *(panmictic_nlineages-1) / (2*G::_ploidy);

            cum_time += increment;
            panmictic_nlineages--;
        }
        assert (panmictic_nlineages == 1);

        assert (gamma_b.size() == q_b.size());
        }
        
        return make_pair(gamma_b, q_b);
    }
#endif

#if !defined (FASTER_SECOND_LEVEL)
    inline pair<vector<double>, vector<unsigned>> Forest::calcCoalescentLikelihoodIntegratingOutTheta(vector<pair<tuple<string,string,string>, double>> species_build) { // TODO: issues when newicks read in have very small branch lengths
         vector<double> gamma_b; // contains gamma_b by species
        vector<unsigned> q_b; // contains q_b by species
        
        gamma_b.resize(G::_nspecies + species_build.size());
        q_b.resize(G::_nspecies + species_build.size());
        
        _panmictic_coalescent_likelihood = 0.0;
        _log_coalescent_likelihood = 0.0;
        
        double neg_inf = -1*numeric_limits<double>::infinity();
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        assert (heights_and_nodes.size() > 0);
        
        double log_coalescent_likelihood = 0.0;
        double cum_time = 0.0;
        string species;
        bool found = false;
        bool coalescence = false;
        unsigned nlineages = 0;
        unsigned a = 0;
        double species_tree_height = 0.0;
        
        unsigned count = -1;
        
        assert (_taxon_map.size() > 0);
        
        setUpGeneForest(_taxon_map);
        
        for (unsigned gen=0; gen < species_build.size(); gen++) {
            count = -1;
            tuple<string, string, string> species_joined = species_build[gen].first;
            double species_increment = species_build[gen].second;
            species_tree_height += species_increment;
            
            updateSpeciesPartition(species_joined);
        
            for (auto &s:_species_partition) { // TODO: how to check if gene tree violates species tree?
                if (gamma_b.size() > 0 && gamma_b[0] != neg_inf) {
    //                    count = _species_indices[s.first];
                    assert (_species_indices.size() > 0);
                    count = _species_indices.at(s.first); // TODO: need to update species indices with species merges?
                    
                    bool done = false;
                    cum_time = 0.0;
                    species = s.first;
                    for (unsigned i=0; i < heights_and_nodes.size(); i++) {
                        if (!done) {
                            found = false;
                            coalescence = false;
                            Node*search_nd = heights_and_nodes[i].second->_left_child;
                            for (auto &nd:s.second) {
                                if (nd == search_nd) {
                                    search_nd = nd;
                                    nlineages = (unsigned) s.second.size();
                                    if (heights_and_nodes[i].first < species_tree_height) {
                                        coalescence = true;
                                        a++;
                                    }
                                    else {
                                        done = true; // stop once gene increment goes below species tree height to avoid double counting some deep coalescences
                                    }
                                    found = true;
                                    break;
                                }
                            }
                            if (found) {
                                if (species_increment == 0.0) {
                                    assert (coalescence); // no deep coalescence if last step
                                }
                                
                                if (coalescence) {
                                    if (s.second.size() == 1) {
                                        q_b.clear(); // don't push back to q_b because it is full of doubles
                                        gamma_b.clear();
                                        gamma_b.push_back(neg_inf);
                                        done = true;
                                        break;
                                    }
                                    
                                    double increment = heights_and_nodes[i].first - (species_tree_height - species_increment) - cum_time; // find increment height
                                    cum_time += increment;
                                    assert (increment > 0.0);
                                    assert (q_b.size() > 0);
                                    assert (gamma_b.size() > 0);
                                    
                                    q_b[count] += 1;
                                    gamma_b[count] += 4* increment * nlineages *(nlineages-1) / (2*G::_ploidy);
                                    for (auto &s:_species_partition) {
                                        if (s.first == species) {
                                            bool found = (find(s.second.begin(), s.second.end(), search_nd->_right_sib) != s.second.end());
                                            if (found ) {
                                                updateNodeList(s.second, search_nd, search_nd->_right_sib, search_nd->_parent); // update the species lineage
                                            }
                                            else {
                                                q_b.clear(); // don't push back to q_b because it is full of doubles
                                                gamma_b.clear();
                                                gamma_b.push_back(neg_inf);
                                                done = true;
                                                break;
                                            }
                                            break;
                                        }
                                    }
                                }
                                else {
                                    double increment = species_increment - cum_time; // find increment height - can't be any more deep coalescence, so gene increment is whatever remains of the species increment
                                    cum_time += increment;
                                    assert (increment > 0.0);
                                    
                                    assert (gamma_b.size() > 0);
                                    gamma_b[count] += 4*increment * nlineages *(nlineages-1) / (2*G::_ploidy);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if (gamma_b.back() != neg_inf) {
            assert (q_b.back() == 0.0);
            assert (gamma_b.back() == 0);
            
            // calculate coalescent likelihood for the rest of the panmictic tree
            double species_increment = species_build.back().second;
            if (log_coalescent_likelihood != neg_inf) {
                if (species_increment > 0.0) {
                    double cum_time = 0.0;
                    // calculate height of each join in the gene tree, minus the species tree height
                    // get number of lineages at beginning of that step and then decrement with each join

                    // start at _nincrements and go through heights and nodes list
                    int panmictic_nlineages = 0;
                    for (auto &s:_species_partition) {
                        panmictic_nlineages += s.second.size();
                    }

                    for (unsigned i=a; i < heights_and_nodes.size(); i++) {
                        // calculate increment (should be heights_and_nodes[i].first - species_tree_height - cum_time (no need to worry about species increment now)
                        // calculate prob of coalescence
                        double increment = heights_and_nodes[i].first - species_tree_height - cum_time;
                        q_b.back() += 1;
                        gamma_b.back() += 4*increment * panmictic_nlineages *(panmictic_nlineages-1) / (2*G::_ploidy);

                        cum_time += increment;
                        panmictic_nlineages--;
                    }
                    assert (panmictic_nlineages == 1);
                }
            }
            
            assert (gamma_b.size() == q_b.size());
        }
        
        _species_partition.clear(); // avoid having to copy extra nodes by clearing this since it will be reset
        return make_pair(gamma_b, q_b);
    }
#endif

    #if !defined (FASTER_SECOND_LEVEL)
    inline void Forest::resetDepthVector(tuple<string, string, string> species_joined) {
        // this function replaces species names with new species name after a join
        for (auto &d:_depths) {
            bool match1 = false;
            bool match2 = false;
            if (d.second.first == get<0>(species_joined) || d.second.first == get<1>(species_joined)) {
                match1 = true;
            }
            if (d.second.second == get<0>(species_joined) || d.second.second == get<1>(species_joined)) {
                match2 = true;
            }
            if (match1) {
                d.second.first = get<2>(species_joined);
            }
            if (match2) {
                d.second.second = get<2>(species_joined);
            }
        }
        
        vector<unsigned> indices_to_erase;
        for (int i=0; i<_depths.size(); i++) {
            if (_depths[i].second.first == _depths[i].second.second) {
    //                // species have already been joined in the species tree, so they are no longer a constraint
                indices_to_erase.push_back(i);
            }
        }
        
        for (unsigned i= (unsigned) indices_to_erase.size() - 1; i != -1; i--) {
            unsigned num = indices_to_erase[i];
            _depths.erase(_depths.begin() + num);
        }
        assert (_depths.size() > 0);
    }
    #endif

    #if !defined (FASTER_SECOND_LEVEL)
    inline vector<pair<double, pair<string, string>>> Forest::getMinDepths() {
        return _depths;
    }
    #endif

    #if !defined (FASTER_SECOND_LEVEL)
    inline void Forest::calcMinDepth() {
        assert (_index != 0);
        _depths.clear();
        // walk through nodes
        // find children of node and determine if children are in the same species
        // if children are in different species, calculate the node height
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
    // post order traversal is the reverse of pre order, now sorted by node heights
        string spp_left_child;
        string spp_right_child;
        bool done = false;
    //        while (spp_left_child == "" || spp_right_child == "") {
        for (int i=0; i<heights_and_nodes.size(); i++) {
            Node* nd = heights_and_nodes[i].second;
            done = false;
                // figure out if descendants of internal node are in the same species
            if (nd->_left_child) {
                while (!done) {
                    Node* left_child = nd->_left_child;
                    Node* right_child = nd->_left_child->_right_sib;
                    while (spp_left_child == "" || spp_right_child == "") {
                        for (auto &s:_species_partition) {
                            for (auto &nd:s.second) {
                                if (nd == left_child && spp_left_child == "") {
                                    spp_left_child = s.first;
                                }
                                else if (nd == right_child && spp_right_child == "") {
                                    spp_right_child = s.first;
                                }
                            }
                            if (spp_left_child != "" && spp_right_child != "") {
                                break;
                            }
                        }
                        
                        if (spp_left_child == "") {
                            if (left_child->_left_child) {
                                left_child = left_child->_left_child;
                                // this is saying if the species of the node is not found, the parent will be dealt with later & ignore the node for now
                                done = true;
                            }
                        }
                        if (spp_right_child == "") {
                            if (right_child->_left_child) {
                                right_child = right_child->_left_child;
                                done = true;
                            }
                        }
                    }
                    if (spp_left_child != "" && spp_right_child != "") {
                        done = true;
                    }
                    if (spp_left_child != spp_right_child) {
                        double height = nd->_height;
                        _depths.push_back(make_pair(height, make_pair(spp_left_child, spp_right_child)));
                    }
                    spp_left_child = "";
                    spp_right_child = "";
                }
            }
        }
        assert (_depths.size() > 0);
    }
    #endif

    inline void Forest::refreshPreorder() {
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

    inline vector< pair<double, Node *>> Forest::sortPreorder() {
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

    inline double Forest::getLineageHeight(Node* nd) {
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
        
    inline unsigned Forest::multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs) {
        // Compute cumulative probababilities
        vector<double> cumprobs(probs.size());
        partial_sum(probs.begin(), probs.end(), cumprobs.begin());
        assert(fabs(*(cumprobs.rbegin()) - 1.0) < 0.0001);

        // Draw a Uniform(0,1) random deviate
        double u = lot->uniform();

        // Find first element in cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6 <-- u
        //
        // cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(cumprobs.begin(), cumprobs.end(), [u](double cumpr){return cumpr > u;});
        if (it == cumprobs.end()) {
            double last_cumprob = *(cumprobs.rbegin());
            throw XProj(format("G::multinomialDraw failed: u = %.9f, last cumprob = %.9f") % u % last_cumprob);
        }

        auto d = std::distance(cumprobs.begin(), it);
        assert(d >= 0);
        assert(d < probs.size());
        return (unsigned)d;
    }

    inline void Forest::simulateData(Lot::SharedPtr lot, unsigned starting_site, unsigned nsites) {
        if (G::_model != "JC") {
            throw XProj("must use JC model for data simulations");
        }
        
        // Create vector of states for each node in the tree
        unsigned nnodes = (unsigned)_nodes.size();
        vector< vector<unsigned> > sequences(nnodes);
        for (unsigned i = 0; i < nnodes; i++) {
            sequences[i].resize(nsites, 4);
        }
        
        // Walk through tree in preorder sequence, simulating all sites as we go
        //    DNA   state      state
        //         (binary)  (decimal)
        //    A      0001        1
        //    C      0010        2
        //    G      0100        4
        //    T      1000        8
        //    ?      1111       15
        //    R      0101        5
        //    Y      1010       10
        
        // Draw equilibrium base frequencies from Dirichlet
        // having parameter G::_comphet
        vector<double> basefreq = {0.25, 0.25, 0.25, 0.25};
        if (Forest::_comphet != G::_infinity) {
            // Draw 4 Gamma(G::_comphet, 1) variates
            double A = lot->gamma(Forest::_comphet, 1.0);
            double C = lot->gamma(Forest::_comphet, 1.0);
            double G = lot->gamma(Forest::_comphet, 1.0);
            double T = lot->gamma(Forest::_comphet, 1.0);
            double total = A + C + G + T;
            basefreq[0] = A/total;
            basefreq[1] = C/total;
            basefreq[2] = G/total;
            basefreq[3] = T/total;
        }
        
        // Simulate starting sequence at the root node

        Node * nd = *(_lineages.begin());
        unsigned ndnum = nd->_number;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = multinomialDraw(lot, basefreq);
        }
        
        nd = findNextPreorder(nd);
        while (nd) {
            ndnum = nd->_number;
            assert(ndnum < nnodes);

            // Get reference to parent sequence
            assert(nd->_parent);
            unsigned parnum = nd->_parent->_number;
            assert(parnum < nnodes);
            
            // Choose relative rate for this node's edge
            // Lognormal m=mean v=variance
            //   sigma^2 = log(v/m^2 + 1)
            //   mu = log m - sigma^2/2
            // Mean m=1 because these are relative rates, so
            //   sigma^2 = log(v + 1)
            //   mu = -sigma^2/2
            double edge_relrate = 1.0;
            if (Forest::_edge_rate_variance > 0.0) {
                double sigma2 = log(1.0 + Forest::_edge_rate_variance);
                double sigma = sqrt(sigma2);
                double mu = -0.5*sigma2;
                double normal_variate = sigma*lot->normal() + mu;
                edge_relrate = exp(normal_variate);
            }

            // Evolve nd's sequence given parent's sequence and edge length
            for (unsigned i = 0; i < nsites; i++) {
                // Choose relative rate for this site
                double site_relrate = 1.0;
                if (Forest::_asrv_shape != G::_infinity)
                    site_relrate = lot->gamma(Forest::_asrv_shape, 1.0/Forest::_asrv_shape);
                unsigned from_state = sequences[parnum][i];
                double cum_prob = 0.0;
                double u = lot->uniform();
                for (unsigned to_state = 0; to_state < 4; to_state++) {
                    cum_prob += calcSimTransitionProbability(from_state, to_state, basefreq, site_relrate*edge_relrate*nd->_edge_length);
                    if (u < cum_prob) {
                        sequences[ndnum][i] = to_state;
                        break;
                    }
                }
                assert(sequences[ndnum][i] < 4);
            }

            // Move to next node in preorder sequence
            nd = findNextPreorder(nd);
        }

        assert(_data);
        Data::data_matrix_t & dm = _data->getDataMatrixNonConst();

        // Copy sequences to data object
    //        for (unsigned t = 0; t < _ntaxa; t++) {
        unsigned t = 0;
        for (auto &nd:_nodes) {
            if (t < G::_ntaxa) {
                // Allocate row t of _data's _data_matrix data member
                dm[t].resize(starting_site + nsites);

                // Get reference to nd's sequence
                unsigned ndnum = nd._number;

                // Translate to state codes and copy
                for (unsigned i = 0; i < nsites; i++) {
                    dm[t][starting_site + i] = (Data::state_t)1 << sequences[ndnum][i];
                }
            }
            t++;
        }
    }

    inline void Forest::stripOutNexusComments(std::string & newick) {
        regex commentexpr("\\[.*?\\]");
        newick = std::regex_replace(newick, commentexpr, std::string(""));
    }

    inline unsigned Forest::countNewickInternals(const string newick) {
        size_t count = count_if(newick.begin(), newick.end(), []( char c ){return c =='(';});
        return count - 1;
    }

    inline unsigned Forest::countNewickLeaves(const std::string newick) {
        regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
        sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
        sregex_iterator m2;
        return (unsigned)std::distance(m1, m2);
    }

    inline bool Forest::canHaveSibling(Node * nd, bool rooted, bool allow_polytomies) {
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

    inline void Forest::renumberInternals() {
        assert(_preorder.size() > 0);

        // Renumber internal nodes in postorder sequence
        unsigned curr_internal = _nleaves;
        for (auto nd : boost::adaptors::reverse(_preorder)) {
            if (nd->_left_child) {
                // nd is an internal node
                nd->_number = curr_internal++;
            }
        }

        _ninternals = curr_internal - _nleaves;

        // If the tree has polytomies, then there are Node objects stored in
        // the _tree->_nodes vector that have not yet been numbered. These can
        // be identified because their _number is currently equal to -1.
        for (auto & nd : _nodes) {
            if (nd._number == -1)
                nd._number = curr_internal++;
        }
    }

    inline void Forest::extractEdgeLen(Node * nd, std::string edge_length_string) {
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

    inline vector<tuple<string, string, string>> Forest::buildFromNewickTopology(string newick) {
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
        _nleaves = countNewickLeaves(commentless_newick);
        
    //        if (_nleaves < 4) {
    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //        }
        unsigned max_nodes = 2*_nleaves - (rooted ? 0 : 2);
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
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % _nleaves));

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

        _nodes.pop_front(); // remove node at beginning of list because it's an extra root
        // remove parent from new last node
        _nodes.front()._parent = NULL;

        _nodes.sort(
             [this](Node& lhs, Node& rhs) {
    //                 return lhs._left_child->_height < rhs._left_child->_height; } ); // TODO: is this just lhs->_height and rhs->_height?
                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );

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

    inline void Forest::buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies) {
        
        set<unsigned> used; // used to ensure that no two leaf nodes have the same number
        unsigned curr_leaf = 0;
        unsigned num_edge_lengths = 0;
        unsigned curr_node_index = 0;

        // Remove comments from the supplied newick string
        string commentless_newick = newick;
        stripOutNexusComments(commentless_newick);

        // Resize the _nodes vector
        _nleaves = countNewickLeaves(commentless_newick);
    //        if (_nleaves < 4) {
    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //        }
        unsigned max_nodes = 2*_nleaves - (rooted ? 0 : 2);
        _nodes.resize(max_nodes);
    //        int b=0;
        for (auto & nd : _nodes ) {
            nd._name = "";
            nd._number = -1;
        }

        try {
            // Root node is the last node in _nodes
            auto l_front = _nodes.begin();
            std::advance(l_front, curr_node_index); // TODO: curr_node_index should be 0
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
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % _nleaves));

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

            if (inside_unquoted_name) {
                cout << "entering exception " << endl;
                throw XProj(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
            }
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
            clear();
            throw x;
        }

        _nodes.pop_front(); // remove node at beginning of list because it's an extra root
        // remove parent from new last node
        _nodes.front()._parent = NULL;
        
        _nodes.sort(
             [this](Node& lhs, Node& rhs) {
    //                 return lhs._left_child->_height < rhs._left_child->_height; } );  // TODO: is this just lhs->_height and rhs->_height?
                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );
        _lineages.clear();
        
        _lineages.push_back(&_nodes.back());
        
        // reset node names
        int j = 0;
        for (auto &nd:_nodes) {
            nd._number = j;
            j++;
        }
        
        if (_index == 0) { // rename internal nodes for species tree only
            for (auto &nd:_nodes) {
                if (nd._name == "") {
                    nd._name=boost::str(boost::format("node-%d")%nd._number);
                }
            }
        }
        
        if (_index == 0) {
            if (_lineages.size() == 1) {
                _nodes.back()._edge_length = 0.0;
            }
        }
    }

    //    inline void Forest::buildFromNewickMPI(const std::string newick, bool rooted, bool allow_polytomies) {
    //        unsigned node_count = 0;
    //        set<unsigned> used; // used to ensure that no two leaf nodes have the same number
    //        unsigned curr_leaf = 0;
    //        unsigned num_edge_lengths = 0;
    //        int curr_node_index = 0;
    //
    //        // Remove comments from the supplied newick string
    //        string commentless_newick = newick;
    //        stripOutNexusComments(commentless_newick);
    //
    //        // Resize the _nodes vector
    //        _nleaves = countNewickLeaves(commentless_newick);
    //        unsigned ninternals = countNewickInternals(commentless_newick);
    //    //        if (_nleaves < 4) {
    //    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //    //        }
    //        unsigned max_nodes = _nleaves + ninternals;
    //        // TODO: if polytomies, max nodes is not this
    ////        unsigned max_nodes = 2*_nleaves - (rooted ? 0 : 2);
    //        _nodes.resize(max_nodes);
    //    //        int b=0;
    //        for (auto & nd : _nodes ) {
    //            nd._name = "";
    //            nd._number = -1;
    //        }
    //
    //        try {
    //            // no root node
    //
    //            // Root node is the last node in _nodes
    //            auto l_front = _nodes.begin();
    //            // TODO: adding
    ////            l_front->_name = "unused";
    //            std::advance(l_front, curr_node_index); // TODO: curr_node_index should be 0
    //            Node *nd = &*l_front;
    //            curr_node_index = -1;
    //
    ////            if (rooted) {
    ////                auto l_front = _nodes.begin();
    ////                std::advance(l_front, ++curr_node_index);
    ////                nd = &*l_front;
    ////
    ////                auto parent = _nodes.begin();
    ////                std::advance(parent, curr_node_index - 1);
    ////                nd->_parent = &*parent;
    ////                nd->_parent->_left_child = nd;
    ////                nd->_name = "unused";
    ////            }
    //
    //            // Some flags to keep track of what we did last
    //            enum {
    //                Prev_Tok_LParen        = 0x01,    // previous token was a left parenthesis ('(')
    //                Prev_Tok_RParen        = 0x02,    // previous token was a right parenthesis (')')
    //                Prev_Tok_Colon        = 0x04,    // previous token was a colon (':')
    //                Prev_Tok_Comma        = 0x08,    // previous token was a comma (',')
    //                Prev_Tok_Name        = 0x10,    // previous token was a node name (e.g. '2', 'P._articulata')
    //                Prev_Tok_EdgeLen    = 0x20    // previous token was an edge length (e.g. '0.1', '1.7e-3')
    //            };
    //            unsigned previous = Prev_Tok_LParen;
    //
    //            // Some useful flag combinations
    //            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
    //            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    //            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    //            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
    //            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);
    //
    //            // Set to true while reading an edge length
    //            bool inside_edge_length = false;
    //            std::string edge_length_str;
    //            unsigned edge_length_position = 0;
    //
    //            // Set to true while reading a node name surrounded by (single) quotes
    //            bool inside_quoted_name = false;
    //
    //            // Set to true while reading a node name not surrounded by (single) quotes
    //            bool inside_unquoted_name = false;
    //
    //            // Set to start of each node name and used in case of error
    //            unsigned node_name_position = 0;
    //
    //            // loop through the characters in newick, building up tree as we go
    //            unsigned position_in_string = 0;
    //            for (auto ch : commentless_newick) {
    //                position_in_string++;
    //
    //                if (inside_quoted_name) {
    //                    if (ch == '\'') {
    //                        inside_quoted_name = false;
    //                        node_name_position = 0;
    //                        if (!nd->_left_child) {
    //                            curr_leaf++;
    //                        }
    //                        previous = Prev_Tok_Name;
    //                    }
    //                    else if (iswspace(ch))
    //                        nd->_name += ' ';
    //                    else {
    //                        nd->_name += ch;
    //                    }
    //
    //                    continue;
    //                }
    //                else if (inside_unquoted_name) {
    //                    if (ch == '(')
    //                        throw XProj(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));
    //
    //                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
    //                        inside_unquoted_name = false;
    //
    //                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
    //                        if (!(previous & Name_Valid))
    //                            throw XProj(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));
    //
    //                        if (!nd->_left_child) {
    //                            curr_leaf++;
    //                        }
    //
    //                        previous = Prev_Tok_Name;
    //                    }
    //                    else {
    //                        nd->_name += ch;
    //                        continue;
    //                    }
    //                }
    //                else if (inside_edge_length) {
    //                    if (ch == ',' || ch == ')' || iswspace(ch)) {
    //                        inside_edge_length = false;
    //                        edge_length_position = 0;
    //                        extractEdgeLen(nd, edge_length_str);
    //                        ++num_edge_lengths;
    //                        previous = Prev_Tok_EdgeLen;
    //                        node_count++;
    //                    }
    //                    else {
    //                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
    //                        if (!valid)
    //                            throw XProj(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
    //                        edge_length_str += ch;
    //                        continue;
    //                    }
    //                }
    //
    //                if (iswspace(ch))
    //                    continue;
    //
    //                switch(ch) {
    //                    case ';':
    //                        break;
    //
    //                    case ')':
    //                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
    //                        if (!nd->_parent)
    //                            throw XProj(boost::str(boost::format("Too many right parentheses at position %d in tree description") % position_in_string));
    //
    //                        // Expect right paren only after an edge length, a node name, or another right paren
    //                        if (!(previous & RParen_Valid))
    //                            throw XProj(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));
    //
    //                        // Go down a level
    //                        nd = nd->_parent;
    //                        if (!nd->_left_child->_right_sib)
    //                            throw XProj(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
    //                        previous = Prev_Tok_RParen;
    //                        break;
    //
    //                    case ':':
    //                        // Expect colon only after a node name or another right paren
    //                        if (!(previous & Colon_Valid))
    //                            throw XProj(boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));
    //                        previous = Prev_Tok_Colon;
    //                        break;
    //
    //                    case ',':
    //                    {
    //                        // Expect comma only after an edge length, a node name, or a right paren
    //                        if (!nd->_parent || !(previous & Comma_Valid))
    //                            throw XProj(boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));
    //
    //                        // Check for polytomies
    //                        if (!canHaveSibling(nd, rooted, allow_polytomies)) {
    //                            throw XProj(boost::str(boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % newick));
    //                        }
    //
    //                        // Create the sibling
    //                        curr_node_index++;
    //                        if (curr_node_index == _nodes.size())
    //                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % _nleaves));
    //
    //                        auto l_front = _nodes.begin();
    //                        std::advance(l_front, curr_node_index);
    //                        nd->_right_sib = &*l_front;
    //
    //                        nd->_right_sib->_parent = nd->_parent;
    //                        nd = nd->_right_sib;
    //                        previous = Prev_Tok_Comma;
    //                        break;
    //                    }
    //
    //                    case '(':
    //                    {
    //                        // Expect left paren only after a comma or another left paren
    //                        if (!(previous & LParen_Valid))
    //                            throw XProj(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % position_in_string));
    //
    //                        // Create new node above and to the left of the current node
    //                        assert(!nd->_left_child);
    //                        curr_node_index++;
    //                        if (curr_node_index == _nodes.size())
    //                            throw XProj(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _nodes.size()));
    //
    //                        auto l_front = _nodes.begin();
    //                        std::advance(l_front, curr_node_index);
    //                        // TODO: don't need next part because no root?
    ////                        nd->_left_child = &*l_front;
    //
    ////                        nd->_left_child->_parent = nd;
    ////                        nd = nd->_left_child;
    //                        previous = Prev_Tok_LParen;
    //                        break;
    //                    }
    //
    //                    case '\'':
    //                        // Encountered an apostrophe, which always indicates the start of a
    //                        // node name (but note that node names do not have to be quoted)
    //
    //                        // Expect node name only after a left paren (child's name), a comma (sib's name)
    //                        // or a right paren (parent's name)
    //                        if (!(previous & Name_Valid))
    //                            throw XProj(boost::str(boost::format("Not expecting node name at position %d in tree description") % position_in_string));
    //
    //                        // Get the rest of the name
    //                        nd->_name.clear();
    //
    //                        inside_quoted_name = true;
    //                        node_name_position = position_in_string;
    //
    //                        break;
    //
    //                    default:
    //                        // Get here if ch is not one of ();:,'
    //
    //                        // Expecting either an edge length or an unquoted node name
    //                        if (previous == Prev_Tok_Colon) {
    //                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
    //                            inside_edge_length = true;
    //                            edge_length_position = position_in_string;
    //                            edge_length_str = ch;
    //                        }
    //                        else {
    //                            // Get the node name
    //                            nd->_name = ch;
    //
    //                            inside_unquoted_name = true;
    //                            node_name_position = position_in_string;
    //                        }
    //                }   // end of switch statement
    //            }   // loop over characters in newick string
    //
    //            if (inside_unquoted_name) {
    //                cout << "entering exception " << endl;
    //                throw XProj(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
    //            }
    //            if (inside_edge_length)
    //                throw XProj(boost::str(boost::format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
    //            if (inside_quoted_name)
    //                throw XProj(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));
    //
    ////            // remove extra fake nodes
    ////            unsigned count = (unsigned) _nodes.size()-1;
    ////            for (auto iter=_nodes.rbegin(); iter != _nodes.rend(); iter++) {
    ////                if (iter->_edge_length == Node::_smallest_edge_length) {
    ////                    _nodes.erase(iter);
    ////                    count--;
    ////                }
    ////            }
    //
    //            for (auto &nd:_nodes) {
    //                if (nd._parent) {
    //                    if (nd._parent->_name == "unused") {
    //                        nd._parent = nullptr;
    //                    }
    //                }
    //            }
    //            _nodes.pop_front(); // remove node at beginning of list because it's an extra root
    //            // remove parent from new last node
    //            _nodes.front()._parent = NULL;
    //
    //            if (allow_polytomies) {
    //                _nodes.pop_front(); // TODO: need to remove the extra 0.0 at the end of the partial newick
    //            }
    //
    //            unsigned nnodes_to_remove = (unsigned) _nodes.size() - node_count;
    //            if (!allow_polytomies) {
    //                assert (nnodes_to_remove == 0);
    //            }
    //            for (unsigned n=0; n<nnodes_to_remove; n++) {
    //                _nodes.pop_back();
    //            }
    //
    //            if (rooted) {
    //                refreshPreorder();
    //            }
    //            renumberInternals();
    //        }
    //
    //        catch(XProj &x) {
    //            clear();
    //            throw x;
    //        }
    //
    //        // TODO: moved this up
    ////        _nodes.pop_front(); // remove node at beginning of list because it's an extra root
    ////        // remove parent from new last node
    ////        _nodes.front()._parent = NULL;
    ////
    //        _nodes.sort(
    //             [this](Node& lhs, Node& rhs) {
    ////                 return lhs._left_child->_height < rhs._left_child->_height; } );  // TODO: is this just lhs->_height and rhs->_height?
    //                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );
    //        _lineages.clear();
    //
    //        for (auto &nd:_nodes) {
    //            if (!nd._parent) {
    //                _lineages.push_back(&nd);
    //            }
    //        }
    ////        _lineages.push_back(&_nodes.back());
    //
    //        // reset node names
    //        int j = 0;
    //        for (auto &nd:_nodes) {
    //            nd._number = j;
    //            j++;
    //        }
    //    }

    inline void Forest::extractNodeNumberFromName(Node * nd, std::set<unsigned> & used) {
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

    inline vector<pair<tuple<string, string, string>, double>> Forest::buildFromNewickMPI(const std::string newick, bool rooted, bool allow_polytomies, Lot::SharedPtr lot) {
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
        _nleaves = countNewickLeaves(commentless_newick);
        unsigned max_nodes = 2*_nleaves - (rooted ? 0 : 2);
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
            for (auto ch : commentless_newick) {
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
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % _nleaves));
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

    //            unsigned max_nodes = countNewickInternals(newick) + _nleaves;
            unsigned ninternals = countNewickInternals(newick);
            unsigned max_internals = _nleaves-1;
            unsigned max_nodes = ninternals + _nleaves + 1;
            if (ninternals != max_internals) {
                _nodes.pop_front(); // if the tree is incomplete, delete both the root and subroot
                _nodes.front()._parent = nullptr;
                max_nodes--;
            }
            
            _nodes.front()._name = "unused"; // break off root node since we are not using it
            
            for (auto &nd:_nodes) {
                if (nd._parent) {
                    if (nd._parent->_name == "unused") {
                        nd._parent = nullptr;
                        nd._right_sib = nullptr;
                    }
                }
            }
            _nodes.pop_front();
            
            unsigned nodes_to_delete = (unsigned)_nodes.size() - max_nodes;
            for (unsigned i=0; i<nodes_to_delete; i++) {
                _nodes.pop_back();
            }
            assert (_nodes.size() == max_nodes);
            
            // sort nodes by height
            _nodes.sort(
                 [this](Node& lhs, Node& rhs) {
                     return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );
                        
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
            
            unsigned num = _nleaves;
            
            if (_index == 0) { // rename internal nodes for species tree only
                for (auto &nd:_nodes) {
                    if (nd._name == "") {
                        nd._name=boost::str(boost::format("node-%d")%num);
                        num++;
    //                        cout << "SETTING NODE NAME TO " << nd._name << endl;
                    }
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

    inline vector<pair<tuple<string, string, string>, double>> Forest::resetLineages(Lot::SharedPtr lot) {
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
                    
                    Node nd;
                    _nodes.push_back(nd);
                    Node* new_nd = &_nodes.back();
                    new_nd->_parent=0;
                    new_nd->_number=_nleaves+_ninternals;
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
                    
                    Node nd;
                    _nodes.push_back(nd);
                    Node* new_nd = &_nodes.back();
                    new_nd->_parent=0;
                    new_nd->_number=_nleaves+_ninternals;
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

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::addCoalInfoElem(const Node * nd, vector<coalinfo_t> & recipient) {
        assert(_index >= 0);

        // Assumes nd is an internal node
        assert(nd->_left_child);
        
        recipient.push_back(
            make_tuple(
                nd->_height,
                _index,
                vector<G::species_t>({
                    nd->_left_child->_species,
                    nd->_left_child->_right_sib->_species,
                })
            )
        );
    }
    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::saveCoalInfoGeneForest(vector<Forest::coalinfo_t> & coalinfo_vect) const {
        // Appends to coalinfo_vect; clear before calling if desired
        // GeneForest version ignores cap argument.
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this
        
        // coalinfo_t is a tuple with these elements:
        // - height of node
        // - 1-offset gene index (0 means speciation)
        // - vector of child species
        
        // Should only be called for complete gene trees
    //        assert(_lineages.size() == 1);

        // Copy tuples stored in _coalinfo to end of coalinfo_vect
        coalinfo_vect.insert(coalinfo_vect.end(), _coalinfo.begin(), _coalinfo.end());
    }
    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::buildCoalInfoVect() {
        assert (_index == 0);
        
        refreshAllPreorders(); // TODO: don't always need to refresh preorders?
        
        _coalinfo.clear();
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
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

    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::saveCoalInfoInitial() {
        // coalinfo_t is a tuple with these elements:
        // - height of node
        // - 1-offset gene index (0 means speciation)
        // - vector of child species
        
        // Should only be called for complete gene trees
        assert (_lineages.size() == 1);
        
        _coalinfo.clear();
        
        refreshPreorder();
        
        // TODO: check node height - don't need to calculate, can just keep track as you go
            for (auto nd : boost::adaptors::reverse(_preorder)) {
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
    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline bool Forest::subsumed(G::species_t test_species, G::species_t subtending_species) {
        bool not_equal = (test_species != subtending_species);
        bool is_subset = ((test_species & subtending_species) == test_species);
        if (not_equal && is_subset)
            return true;
        else
            return false;
    }

    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const {
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
    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::saveCoalInfoSpeciesTree(vector<Forest::coalinfo_t> & coalinfo_vect, bool cap) {
        // Appends to coalinfo_vect; clear coalinfo_vect before calling if desired
        // Assumes heights and preorders are up-to-date; call
        //   heightsInternalsPreorders() beforehand to ensure this
        
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

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::refreshAllPreorders() const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        assert (_index == 0); // only need to call this for a species tree
        _next_node_number = G::_nspecies;
        _preorders.clear();
        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
            if (nd->_left_child) {
                nd->_number = _next_node_number++;
            }
            
            // lineage is a Node::ptr_vect_t (i.e. vector<Node *>)
            // lineage[0] is the first node pointer in the preorder sequence for this lineage
            // Add a new vector to _preorders containing, for now, just the root of the subtree
            _preorders.push_back({nd});
            
            // Now add the nodes above the root in preorder sequence
            Node::ptr_vect_t & preorder_vector = *_preorders.rbegin();
            refreshPreorderNew(preorder_vector);
        }
    }

    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::refreshPreorderNew(vector<Node*> & preorder) const {
        // Assumes preorder just contains the root node when this function is called
        // Also assumes that _next_node_number was initialized prior to calling this function
        assert(preorder.size() == 1);
        
        Node * nd = preorder[0];
        while (true) {
            nd = findNextPreorderNew(nd);
            if (nd) {
                preorder.push_back(nd);
                if (nd->_left_child)
                    nd->_number = _next_node_number++;
            }
            else
                break;
        }
    }
    #endif

    #if defined (FASTER_SECOND_LEVEL)
    inline void Forest::setTreeHeight() {
        _forest_height = _lineages.back()->_height;
    }
    #endif

    inline void Forest::setNodeHeights() {
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

    inline void Forest::resetSpeciesPartition(string species_partition_string) {
        _species_partition.clear();
        cout << species_partition_string << endl;
        for (auto s:species_partition_string) {
            cout << s << endl;
        }
        
    //        for (auto &m:species_partition) {
    //            _species_partition[m.first];
    //            for (auto &s:m.second) {
    //                for (auto &nd:_nodes) {
    //                    if (nd._name == s) {
    //                        _species_partition[m.first].push_back(&nd);
    //                        break;
    //                    }
    //                }
    //            }
    //        }
    }

    inline map<string, vector<string>> Forest::saveSpeciesPartition() {
        map<string, vector<string>> species_partition_strings;
        for (auto &s:_species_partition) {
            species_partition_strings[s.first];
            for (auto &nd:s.second) {
    //                species_partition_strings[s.first].push_back(to_string(nd->_number));
                species_partition_strings[s.first].push_back(nd->_name);
            }
        }
        return species_partition_strings;
    }

    #if defined (FASTER_SECOND_LEVEL)
    inline Node * Forest::findNextPreorderNew(Node * nd) const {
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

    #endif
    }

