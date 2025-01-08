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
extern proj::Lot rng;
std::mutex mtx;

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

        unsigned                    numLeaves() const;
        unsigned                    numInternals() const;
        unsigned                    numNodes() const;
        void                        showForest();
        static void                 setNumSpecies(unsigned n);
        static void                 setNumTaxa(unsigned n);
        double                      calcLogLikelihood();
        void                        createDefaultTree(Lot::SharedPtr lot);
        void operator=(const Forest & other);
        void                        debugForest();
        void                        debugLogLikelihood(Node* nd, double log_like);

    private:

        void                        clear();
        void                        setData(Data::SharedPtr d, int index, map<string, string> &taxon_map, bool partials);
        void                        setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa);
        Node *                      findNextPreorder(Node * nd);
        string                      makeNewick(unsigned precision, bool use_names);
        string                      makeAltNewick(unsigned precision, bool use_names);
        string                      makePartialNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>    chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        void                        calcPartialArray(Node* new_nd);
        void                        setUpGeneForest(map<string, string> &taxon_map);
        void                        setUpSpeciesForest(vector<string> &species_names);
        tuple<string,string, string> speciesTreeProposal(Lot::SharedPtr lot);
#if defined (FOSSILS)
        tuple<string,string, string> speciesTreeProposalFossils(Lot::SharedPtr lot);
#endif
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                        revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                      getRunningSumChoices(vector<double> &log_weight_choices);
        vector<double>              reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        int                         selectPair(vector<double> weight_vec, Lot::SharedPtr lot);
        void                        chooseSpeciesIncrement(Lot::SharedPtr lot);
        void                        chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth);
#if defined (FOSSILS)
        bool                        chooseSpeciesIncrementFossil(Lot::SharedPtr lot, double next_fossil_time, string fossil_name);
        bool                        chooseSpeciesIncrementFossilFirstStep(Lot::SharedPtr lot, double next_fossil_time, string fossil_name);
#endif
        void                        addSpeciesIncrement();
        void                        allowCoalescence(string species_name, double increment, Lot::SharedPtr lot);
        vector<pair<double, string>>      calcForestRate(Lot::SharedPtr lot);
        void                        updateSpeciesPartition(tuple<string, string, string> species_info);
        double                      calcTopologyPrior(unsigned nlineages);
        void                        calcIncrementPrior(double increment, string species_name, bool new_increment, bool coalesced_gene, bool gene_tree);
        void                        clearPartials();
        void                        setStartMode(string mode) {_start_mode = mode;}
        void                        setRelativeRate(double rel_rate) {_relative_rate = rel_rate;}
        unsigned                    getDeepCoal(tuple <string, string, string> species_joined);
        unsigned                    getMaxDeepCoal(tuple <string, string, string> species_joined);
        void                        setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
        void                        resetDepthVector(tuple<string, string, string> species_joined);
        vector<pair<double, pair<string, string>>>             getMinDepths();
        void                        calcMinDepth();
        vector< pair<double, Node *>>      sortPreorder();
        void                        refreshPreorder();
        void                        createThetaMap(Lot::SharedPtr lot);
        void                        createThetaMapFixedTheta(Lot::SharedPtr lot);
        void                        createSpeciesIndices();
        void                        updateThetaMap(Lot::SharedPtr lot, string new_species_name);
        void                        updateThetaMapFixedTheta(Lot::SharedPtr lot, string new_species_name);
        void                        resetThetaMap(Lot::SharedPtr lot);
        void                        drawNewTheta(string new_species, Lot::SharedPtr lot);
        void                        buildFromNewick(const string newick, bool rooted, bool allow_polytomies);
        void                        stripOutNexusComments(std::string & newick);
        unsigned                    countNewickLeaves(const std::string newick);
        void                        extractEdgeLen(Node * nd, std::string edge_length_string);
        void                        renumberInternals();
        bool                        canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);
        vector<tuple<string, string, string>>              buildFromNewickTopology(const string newick);
        pair<Node*, Node*>          chooseAllPairs(list<Node*> &nodes, double increment, string species, Lot::SharedPtr lot);
        tuple<Node*, Node*, Node*>  createNewSubtree(pair<unsigned, unsigned> t, list<Node*> node_list, double increment, string species);
        pair<Node*, Node*>          getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list);
#if defined (BUILD_UPGMA_TREE)
# if defined (BUILD_UPGMA_TREE_CONSTRAINED)
    void                        buildRestOfTree(Lot::SharedPtr lot, vector<pair<tuple<string, string, string>, double>> species_info);
# else
        void                        buildRestOfTree(Lot::SharedPtr lot);
#endif
#endif
        void                        debugShowDistanceMatrix(const vector<double> & d) const;
    
#if defined (FASTER_UPGMA_TREE)
        void                        buildRestOfTreeFaster();
        void                        buildStartingUPGMAMatrix();
        void                        buildStartingRow();
#endif

    
        map<string, double>         _theta_map;

        std::vector<Node *>         _lineages;
    
        std::list<Node>             _nodes;
        std::vector<Node*>          _new_nodes;

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
        map<string, list<Node*> >   _species_partition;
        double                      _gene_tree_log_likelihood;
        pair<Node*, Node*>          _species_joined;
        double                      _log_joining_prob;
        vector<pair<double, double>> _increments_and_priors;
        bool                        _done;
        double                      _log_coalescent_likelihood;
        double                      _panmictic_coalescent_likelihood;
        double                      _log_coalescent_likelihood_increment;
        double                      _cum_height;
        vector<pair<double, pair<string, string>>>              _depths;
        unsigned                    _nincrements = 0;
        vector<Node*>               _preorder;
        double                      _small_enough;
        vector<pair<tuple<string,string,string>, double>>    _species_build;
        map<string, string>         _taxon_map;
        map<string, unsigned>       _species_indices;
        double                      _log_weight;
        string                      _ancestral_species_name;
        vector<double>              _vector_prior;
        double                      _theta_mean;
        vector<pair<Node*, Node*>>  _node_choices;
        vector<double>              _log_likelihood_choices;
     // TODO: only define these as needed with conditionals
        stack<Node *>               _upgma_additions;
        map<Node *, double>         _upgma_starting_edgelen;
        vector<string>              _species_names;
        vector<double>              _starting_dij;
        map<Node*,  unsigned>       _starting_row;
        double                      _relative_rate;
        vector<pair<string, unsigned>>  _lineages_per_species;
    
        void                        showSpeciesJoined();
        double                      calcTransitionProbability(Node* child, double s, double s_child);
        double                      calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length);
        double                      getTreeHeight();
        double                      getTreeLength();
        double                      getSpeciesTreeIncrement();
        double                      getLineageHeight(Node* nd);
        void                        addIncrement(double increment);
        void                        simulateData(Lot::SharedPtr lot, unsigned starting_site, unsigned nsites);
        double                      calcCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height);
        pair<vector<double>, vector<unsigned>>        calcCoalescentLikelihoodIntegratingOutTheta(vector<pair<tuple<string,string,string>, double>> species_build);
        pair<vector<double>, vector<unsigned>> calcInitialCoalescentLikelihoodIntegratingOutTheta();
        pair<vector<double>, vector<unsigned>>        calcCoalescentLikelihoodIntegratingOutThetaLastStep(vector<pair<tuple<string,string,string>, double>> species_build);
        
        unsigned                    multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);

    public:

        typedef std::shared_ptr<Forest> SharedPtr;
        static double               _theta;
        static double               _theta_proposal_mean;
        static double               _theta_prior_mean;
        double                      _lambda;
        double                      _extinction_rate;
        static string               _proposal;
        static string               _model;
        static double               _kappa;
        static vector<double>       _base_frequencies;
        static string               _string_base_frequencies;
        static bool                 _save_memory;
        static string               _outgroup;
        static bool                 _run_on_empty;
        static double               _ploidy;
        static string               _start_mode;
        static double               _edge_rate_variance;
        static double               _asrv_shape;
        static double               _comphet;
        static double               _infinity;
        static double               _clock_rate;
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
        _nstates = 4;
        _last_edge_length = 0.0;
        _lineages.clear();
        _log_weight = 0.0;
        _gene_tree_log_likelihood = 0.0;
        _log_joining_prob = 0.0;
        _done = false;
        _log_coalescent_likelihood = 0.0;
        _log_coalescent_likelihood_increment = 0.0;
        _cum_height = 0.0;
        _nleaves=_ntaxa;
        _ninternals=0;
        _depths.clear();
        _nincrements = 0;
        _preorder.clear();
        _panmictic_coalescent_likelihood = 0.0;
        _small_enough = 0.0000001;
        _theta_mean = 0.0;
        _ancestral_species_name = "";
        _ploidy = 2;
        _species_build.clear();
        _taxon_map.clear();
        _species_indices.clear();
        _vector_prior.clear();
        _upgma_additions = stack<Node*>();
        _upgma_starting_edgelen.clear();
        _species_names.clear();
        _starting_dij.clear();
        _starting_row.clear();
        _lineages_per_species.clear();
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa) {
        _index = index;
        assert (index > 0);         //don't set data for species tree, though it doesn't really matter for simulations
//        const Data::taxon_names_t & taxon_names = _data->getTaxonNames();
        _ntaxa = ntaxa;
        _nleaves = ntaxa;
        
        _data = d;

//        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(index-1);
//        _first_pattern = gene_begin_end.first;
//        _npatterns = _data->getNumPatternsInSubset(index-1);
        
        _nodes.resize(_ntaxa);
        _lineages.reserve(_nodes.size());
        unsigned i= 0;
        
        //create taxa
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_position_in_lineages=i;
            _lineages.push_back(nd);
        }
        
        vector<string> taxon_names;
        for (auto &t:taxon_map) {
            taxon_names.push_back(t.first);
        }

//        if (!_save_memory || (_save_memory && partials)) { // if save memory setting, don't set tip partials yet
//            nd->_partial=ps.getPartial(_npatterns*4);
//            for (unsigned p=0; p<_npatterns; p++) {
//                unsigned pp = _first_pattern+p;
//                for (unsigned s=0; s<_nstates; s++) {
//                    Data::state_t state = (Data::state_t)1 << s;
//                    Data::state_t d = data_matrix[nd->_number][pp];
//                    double result = state & d;
//                    (*nd->_partial)[p*_nstates+s]= (result == 0.0 ? 0.0:1.0);
//                }
//            }
//        }
        
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
        
        _nodes.resize(_ntaxa);
        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_position_in_lineages=i;
            _lineages.push_back(nd);
        }

        for (auto &nd:_lineages) {
            if (!nd->_left_child) {
                // replace all spaces with underscores so that other programs do not have
                  // trouble parsing your tree descriptions
                  std::string name = taxon_names[i++];
                  boost::replace_all(name, " ", "_");
                nd->_name = name;

                if (!_save_memory || (_save_memory && partials)) { // if save memory setting, don't set tip partials yet
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
//            cout << "log coalescent likelihood: " << _log_coalescent_likelihood << "\n";
        }
        else if (_index == 0) {
            cout << " species tree: " << endl;
        }
        cout << " " << makeNewick(15, true) << "\n";
        cout << "\n";
    }

    inline string Forest::makePartialNewick(unsigned precision, bool use_names) {
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

    inline void Forest::setNumTaxa(unsigned n){
        _ntaxa=n;
    }

    inline void Forest::setNumSpecies(unsigned n) {
        _nspecies=n;
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
                                                        
    inline void Forest::calcPartialArray(Node* new_nd) {
        assert (_index > 0);
        
        if (!new_nd->_left_child) {
            auto &data_matrix=_data->getDataMatrix();
            assert (_save_memory || _start_mode == "sim");
            if (!new_nd->_left_child) {
                new_nd->_partial=ps.getPartial(_npatterns*4);
                for (unsigned p=0; p<_npatterns; p++) {
                    unsigned pp = _first_pattern+p;
                    for (unsigned s=0; s<_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[new_nd->_number][pp];
                        double result = state & d;
                        (*new_nd->_partial)[p*_nstates+s]= (result == 0.0 ? 0.0:1.0);
                    }
                }
            }
        }
        
//        for (auto &nd:_lineages) {
//            assert (nd->_right_sib != nd);
//        }
        
        auto & parent_partial_array = *(new_nd->_partial);
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
            
            if (child->_partial == nullptr) {
                child->_partial = ps.getPartial(_npatterns*4);
                calcPartialArray(child);
            }
            assert (child->_partial != nullptr);
            auto & child_partial_array = *(child->_partial);

            for (unsigned p = 0; p < _npatterns; p++) {
                for (unsigned s = 0; s <_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < _nstates; s_child++) {
                        double child_transition_prob = calcTransitionProbability(child, s, s_child);
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

    inline double Forest::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) {
        assert(pi.size() == 4);
        assert(fabs(accumulate(pi.begin(), pi.end(), 0.0) - 1.0) < Forest::_small_enough);
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

    inline double Forest::calcTransitionProbability(Node* child, double s, double s_child) {
        double child_transition_prob = 0.0;

        if (_model == "JC" ) {
            double expterm = exp(-4.0*(child->_edge_length)*_relative_rate/3.0); // TODO: double check relative rate
            double prsame = 0.25+0.75*expterm;
            double prdif = 0.25 - 0.25*expterm;

            child_transition_prob = (s == s_child ? prsame : prdif);
            assert (child_transition_prob > 0.0);
            return child_transition_prob;
        }

        if (_model == "HKY") {
            double pi_A = _base_frequencies[0];
            double pi_C = _base_frequencies[1];
            double pi_G = _base_frequencies[2];
            double pi_T = _base_frequencies[3];

            double pi_j = 0.0;
            double PI_J = 0.0;

            double phi = (pi_A+pi_G)*(pi_C+pi_T)+_kappa*(pi_A*pi_G+pi_C*pi_T);
            double beta_t = 0.5*(child->_edge_length * _relative_rate )/phi; // TODO: double check relative rate

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
        }
        assert (child_transition_prob > 0.0);
        return child_transition_prob;
    }

    inline double Forest::calcLogLikelihood() { // TODO: if relative rates specified, multiply branch lengths by rel rate for the locus when calculating the likelihood
        //calc likelihood for each lineage separately
        auto &counts = _data->getPatternCounts();
        _gene_tree_log_likelihood = 0.0;

        for (auto &nd:_lineages) {
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
        // TODO: need to know the number of lineages in each tip species
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
        // TODO: add new element to _lienages_per_species
        
        
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
        
//        if (nlineages1 > 1) {
            num_deep_coal += nlineages1 - 1;
//        }
//        if (nlineages2 > 1) {
            num_deep_coal += nlineages2 - 1;
//        }
        
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
        double edge_length = lot->gamma(1.0, 1.0/_ntaxa);
        
        _lineages.reserve(_nodes.size());
        
        for (unsigned i = 0; i < _ntaxa; i++) {
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
        _nleaves=_ntaxa;
        _ninternals=0;
        _last_edge_length = 0.0;
    }

    inline void Forest::operator=(const Forest & other) {
        _nstates = other._nstates;
        _npatterns = other._npatterns;
        _edge_rate_variance = other._edge_rate_variance;
        _asrv_shape = other._asrv_shape;
        _comphet = other._comphet;
        _nodes.clear();
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
        _nspecies           = other._nspecies;
        _ntaxa              = other._ntaxa;
        _species_joined = other._species_joined;
        _log_joining_prob = other._log_joining_prob;
        _increments_and_priors = other._increments_and_priors;
        _done = other._done;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _log_coalescent_likelihood_increment = other._log_coalescent_likelihood_increment;
        _cum_height = other._cum_height;
        _outgroup = other._outgroup;
        _run_on_empty = other._run_on_empty;
        _start_mode = other._start_mode;
        _relative_rate = other._relative_rate;
        _depths = other._depths;
        _nincrements = other._nincrements;
        _preorder.resize(other._preorder.size());
        _panmictic_coalescent_likelihood = other._panmictic_coalescent_likelihood;
        _theta_map = other._theta_map; // TODO: unsure if this will copy correctly
        _small_enough = other._small_enough;
        _theta_proposal_mean = other._theta_proposal_mean;
        _theta_mean = other._theta_mean;
        _theta_prior_mean = other._theta_prior_mean;
        _ancestral_species_name = other._ancestral_species_name;
        _ploidy = other._ploidy;
        _species_build = other._species_build;
        _taxon_map = other._taxon_map;
        _species_indices = other._species_indices;
        _vector_prior = other._vector_prior;
        _infinity = other._infinity;
        _lambda = other._lambda;
        _lineages_per_species = other._lineages_per_species;
#if defined(BUILD_UPGMA_TREE)
        _upgma_additions = other._upgma_additions;
        _upgma_starting_edgelen = other._upgma_starting_edgelen;
        _starting_dij = other._starting_dij;
        _starting_row = other._starting_row;
#endif
        _species_names = other._species_names;

        // copy tree itself

        _species_partition.clear();
        for (auto spiter : other._species_partition) {
            for (auto s : spiter.second) {
                unsigned number = s->_number;
                Node* nd = &*next(_nodes.begin(), number);
                _species_partition[spiter.first].push_back(nd);
            }
        }
        
#if defined(BUILD_UPGMA_TREE)
        _starting_row.clear();
        for (auto strow : other._starting_row) {
            unsigned number = strow.first->_number;
            Node* nd = &*next(_nodes.begin(), number);
            _starting_row[nd] = strow.second;
        } // TODO: check this is copying correctly
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

    inline void Forest::setUpSpeciesForest(vector<string> &species_names) {
        _index = 0;
        assert (_nspecies = (unsigned) species_names.size());
        
        //create species
        _nodes.resize(_nspecies);
        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < _nspecies; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_position_in_lineages=i;
            nd->_name=species_names[i];
            _lineages.push_back(nd);
            }
        
        _nleaves=_nspecies;
        _ninternals=0;
    }

    inline void Forest::chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth) {
        assert (max_depth >= 0.0);
        if (max_depth > 0.0) {
            double rate = (_lambda)*_lineages.size();
            
#if defined (FOSSILS)
        rate = (_lambda - _extinction_rate) * _lineages.size(); // TODO: unsure, also number of lineages changes based on whether fossil has been added
#endif
            
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
            double rate = _lambda*_lineages.size();
            
#if defined (FOSSILS)
        rate = (_lambda - _extinction_rate) * _lineages.size(); // TODO: unsure, also number of lineages changes based on whether fossil has been added
#endif
            
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
        
        if (_species_build.size() == 0) {
            _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
        }
        else {
            _species_build.back().second = _last_edge_length;
        }
    }

#if defined (FOSSILS)
    inline bool Forest::chooseSpeciesIncrementFossil(Lot::SharedPtr lot, double next_fossil_time, string fossil_name) {
        bool fossil_join = false;
        
        double rate = (_lambda - _extinction_rate) * _lineages.size(); // TODO: unsure, also number of lineages changes based on whether fossil has been added
            
        assert (lot != nullptr);
        double edge_len = lot->gamma(1.0, 1.0/rate) * _clock_rate; // TODO: unsure
        double species_tree_height = getTreeHeight() + edge_len;
        
//        _last_edge_length = lot->gamma(1.0, 1.0/rate);
        
        if (species_tree_height < next_fossil_time || next_fossil_time == 0.0) {
            _last_edge_length = edge_len;
            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));
        
            if (_species_build.size() == 0) {
                _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
            }
            else {
                _species_build.back().second = _last_edge_length;
            }
        }
        
        else {
            // time chosen exceed fossil; add fossil instead
            species_tree_height -= edge_len;
            double edge_len = next_fossil_time - species_tree_height;
            assert (edge_len > 0.0);
            _last_edge_length = edge_len;
            
            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            // add the fossil to _lineages and _nodes in the species tree
            Node nd;
            _nodes.push_back(nd);
            Node* new_nd = &_nodes.back();
            new_nd->_parent=0;
            new_nd->_number=_nleaves+_ninternals;
            new_nd->_name=fossil_name + "_FOSSIL";
            new_nd->_edge_length=edge_len;
            _ninternals++;
            new_nd->_right_sib=0;
            _lineages.push_back(new_nd);
                        
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));

        if (_species_build.size() == 0) {
                _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
            }
            else {
                _species_build.back().second = _last_edge_length;
            }
            fossil_join = true;
        }
        return fossil_join;
    }
#endif

#if defined (FOSSILS)
    inline bool Forest::chooseSpeciesIncrementFossilFirstStep(Lot::SharedPtr lot, double next_fossil_time, string fossil_name) {
        bool fossil_join = false;
        double rate = (_lambda - _extinction_rate) * _lineages.size(); // TODO: unsure, also number of lineages changes based on whether fossil has been added
            
        assert (lot != nullptr);
        double edge_len = lot->gamma(1.0, 1.0/rate);
        edge_len = 100.0; // TODO: BE CAREFUL
        double species_tree_height = getTreeHeight() + edge_len;
                
        if (species_tree_height < next_fossil_time || next_fossil_time == 0.0) {
            _last_edge_length = edge_len;
            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));
        
            if (_species_build.size() == 0) {
                _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
            }
            else {
                _species_build.back().second = _last_edge_length;
            }
            fossil_join = false;
        }
        
        else {
            // time chosen exceed fossil; add fossil increment instead; next step will need to join that fossil immediately
            species_tree_height -= edge_len;
            double edge_len = next_fossil_time - species_tree_height;
            assert (edge_len > 0.0);
            _last_edge_length = edge_len;
            
            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
                    
            // add the fossil to _lineages and _nodes in the species tree
            Node nd;
            _nodes.push_back(nd);
            Node* new_nd = &_nodes.back();
            new_nd->_parent=0;
            new_nd->_number=_nleaves+_ninternals;
            new_nd->_name=fossil_name + "_FOSSIL";
            new_nd->_edge_length=edge_len;
            _ninternals++;
            new_nd->_right_sib=0;
            _lineages.push_back(new_nd);
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));

        if (_species_build.size() == 0) {
                _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
            }
            else {
                _species_build.back().second = _last_edge_length;
            }
            fossil_join = true;;
        }
        return fossil_join;
    }
#endif

    inline void Forest::chooseSpeciesIncrement(Lot::SharedPtr lot) {
        double rate = _lambda*_lineages.size();
        
#if defined (FOSSILS)
        rate = (_lambda - _extinction_rate) * _lineages.size(); // TODO: unsure, also number of lineages changes based on whether fossil has been added
#endif
        
        assert (lot != nullptr);
        _last_edge_length = lot->gamma(1.0, 1.0/rate);

        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
    }


    inline tuple<string,string, string> Forest::speciesTreeProposal(Lot::SharedPtr lot) {
        _cum_height = 0.0; // reset cum height for a new lineage
        // this function creates a new node and joins two species
        
        bool done = false;
        Node* subtree1 = nullptr;
        Node* subtree2 = nullptr;
        
        while (!done) {
        
    //        pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size(), lot);
            assert (lot != nullptr);
            pair<unsigned, unsigned> t = lot->nchoose2((unsigned) _lineages.size());
            
            assert (t.first != t.second);
            subtree1=_lineages[t.first];
            subtree2=_lineages[t.second];
            assert (t.first < _lineages.size());
            assert (t.second < _lineages.size());
            assert(!subtree1->_parent && !subtree2->_parent);
            
            if (_outgroup != "none") {
                if (subtree1->_name != _outgroup && subtree2->_name != _outgroup && _lineages.size() > 2) { // outgroup can only be chosen on the last step
                    done = true;
                }
                else if (_lineages.size() == 2) {
                    done = true;
                }
            }
            else {
                done = true;
            }
            if (_outgroup == "none") {
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
        
        calcTopologyPrior((int) _lineages.size()+1);

        _species_build.push_back(make_pair(make_tuple(subtree1->_name, subtree2->_name, new_nd->_name), 0.0));
        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
    }

    inline tuple<string,string, string> Forest::speciesTreeProposalFossils(Lot::SharedPtr lot) {
        _cum_height = 0.0; // reset cum height for a new lineage
        // this function creates a new node and joins two species
        
        bool done = false;
        Node* subtree1 = nullptr;
        Node* subtree2 = nullptr;
        
        while (!done) {
        
    //        pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size(), lot);
            assert (lot != nullptr);
            unsigned t1 = (unsigned) (lot->randint(0, _lineages.size() - 2)); // last element of _lineages is the fossil
//            pair<unsigned, unsigned> t = lot->nchoose2((unsigned) _lineages.size());
            pair<unsigned, unsigned> t = make_pair (t1, _lineages.size() - 1);
            
            assert (t.first != t.second);
            subtree1=_lineages[t.first];
            subtree2=_lineages[t.second];
            assert (t.first < _lineages.size());
            assert (t.second < _lineages.size());
            assert(!subtree1->_parent && !subtree2->_parent);
            
            if (_outgroup != "none") {
                if (subtree1->_name != _outgroup && subtree2->_name != _outgroup && _lineages.size() > 2) { // outgroup can only be chosen on the last step
                    done = true;
                }
                else if (_lineages.size() == 2) {
                    done = true;
                }
            }
            else {
                done = true;
            }
            if (_outgroup == "none") {
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
        
        calcTopologyPrior((int) _lineages.size()+1);

        _species_build.push_back(make_pair(make_tuple(subtree1->_name, subtree2->_name, new_nd->_name), 0.0));
        
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

    inline void Forest::showSpeciesJoined() {
        assert (_index==0);
        if (_species_joined.first != NULL) {
            cout << "joining species " << _species_joined.first->_name << " and " << _species_joined.second->_name << endl;
        }
        else {
            cout << "no species joined" << endl;
        }
    }

    inline void Forest::setUpGeneForest(map<string, string> &taxon_map) {
        _taxon_map = taxon_map;
        assert (_index >0);
        _species_partition.clear();
        
        unsigned count = 0;
        
        for (auto &nd:_nodes) {
            count++;
            assert (!nd._left_child);
//            if (!nd._left_child) {
            string species_name = taxon_map[nd._name];
            _species_partition[species_name].push_back(&nd);
            if (count == Forest::_ntaxa) {
                break;
            }
        }
        
        assert (_species_partition.size() == Forest::_nspecies);
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
                    double rate = _lambda*(_lineages.size());
                    // calculate increment prior
                    double nChooseTwo = (_lineages.size())*(_lineages.size()-1);
                    double log_prob_join = log(2/nChooseTwo);
//                    log_increment_prior = log(_lambda) - (increment * rate) + log_prob_join;
//                    log_increment_prior = log(_lambda) - (increment * rate);
                    log_increment_prior = log(rate) - (increment*rate) + log_prob_join;
                }
                else {
                    double rate = _lambda*(_lineages.size());
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
                            double coalescence_rate = (s.second.size())*(s.second.size()-1) / _theta;
                            assert (coalescence_rate > 0.0); // rate should be >0 if there is coalescence
                            double nChooseTwo = (s.second.size())*(s.second.size()-1);
                            double log_prob_join = log(2/nChooseTwo);
                            log_increment_prior += log(coalescence_rate) - (increment*coalescence_rate) + log_prob_join;
                        }
                        else {
                            // no coalescence
                            double coalescence_rate = s.second.size()*(s.second.size() - 1) / _theta;
                            log_increment_prior -= increment*coalescence_rate;
                        }
                    }
                }
            
                else if (!coalesced_gene) {
                    // no coalescence
                    for (auto &s:_species_partition) {
                        double coalescence_rate = s.second.size() * (s.second.size() - 1) / _theta;
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
                    double rate = _lambda*(_lineages.size());
                    // calculate increment prior
                    double nChooseTwo = (_lineages.size())*(_lineages.size()-1);
                    double log_prob_join = log(2/nChooseTwo);
                    log_increment_prior = log(rate) - (increment*rate) + log_prob_join;
                }
                else {
                    double rate = _lambda*(_lineages.size());
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


# if defined (BUILD_UPGMA_TREE)
# if defined (BUILD_UPGMA_TREE_CONSTRAINED)
    inline void Forest::buildRestOfTree(Lot::SharedPtr lot, vector<pair<tuple<string, string, string>, double>> species_info) {
#else
    inline void Forest::buildRestOfTree(Lot::SharedPtr lot) {
#endif
        
        // TODO: try using starting distances rather than likelihood minimizer
# if defined (FASTER_UPGMA_TREE)
        buildRestOfTreeFaster();
#else
        double prev_log_likelihood = _gene_tree_log_likelihood;
        
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
        vector<double> dij(n*(n-1)/2, _infinity);
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
                    for (unsigned lstate = 0; lstate < _nstates; lstate++) {
                        auto & l_partial_array = *(lnode->_partial);
                        double lpartial = l_partial_array[p*_nstates + lstate];
                        for (unsigned rstate = 0; rstate < _nstates; rstate++) {
                            auto & r_partial_array = *(rnode->_partial);
                            double rpartial = r_partial_array[p*_nstates + rstate];
                            if (lstate == rstate)
                                same_state[p] += lpartial*rpartial;
                            else
                                diff_state[p] += lpartial*rpartial;
                        }
                    }
                }
                
                double min_dist = 0.0;
                double max_dist = min_dist + 5.0; //TODO: replace arbitrary value 5.0
                
                double v0 = 0.0;
                
#if defined(BUILD_UPGMA_TREE_CONSTRAINED)
                
//                // TODO: for now, walk through species partition to find them, but can make this faster
//                // find the distance from the tips to the deeper of the two species, then subtract the height of the existing gene tree
                string lspp = "";
                string rspp = "";

                for (auto &s:_species_partition) {
                    for (auto &nd:s.second) {
                        if (nd == lnode) {
                            lspp = s.first;
                        }
                        else if (nd == rnode) {
                            rspp = s.first;
                        }
                    }
                    if (lspp != "" && rspp != "") {
                        break;
                    }
                }
                
                if (lspp != rspp) {
                    
                    double mrca_height = 0.0;
                
                    for (auto &s:species_info) {
                        if (get<0>(s.first) == lspp) {
                            lspp = get<2>(s.first);
                        }
                        else if (get<1>(s.first) == lspp) {
                            lspp = get<2>(s.first);
                        }
                        
                        if (get<0>(s.first) == rspp) {
                            rspp = get<2>(s.first);
                        }
                        else if (get<1>(s.first) == rspp) {
                            rspp = get<2>(s.first);
                        }
                        
                        if (lspp != rspp) {
                            mrca_height += s.second;
                            // don't include any increment for the population where the two species have joined for the first time
                        }
                        
                        if (lspp == rspp) {
                            break;
                        }
                    }
                    
                        double gene_tree_height = getTreeHeight();
                    
                    mrca_height -= gene_tree_height;
                    
                    min_dist = 2.0 * mrca_height;

                    assert (min_dist > 0.0);
                }
#endif
                
                v0 = lnode->getEdgeLength() + rnode->getEdgeLength();
                
                negLogLikeDist f(npatterns, first_pattern, counts, same_state, diff_state, v0);
                auto r = boost::math::tools::brent_find_minima(f, min_dist, max_dist, std::numeric_limits<double>::digits);
//                double maximized_log_likelihood = -r.second;
                unsigned k = i*(i-1)/2 + j;
                dij[k] = r.first;
                dij_row_col[k] = make_pair(i,j);
                                
    //                output(format("d[%d, %d] = %.7f (logL = %.5f") % i % j % d[ij] % maximized_log_likelihood, 1);
            
            }
        }
        
        // Create a map relating nodes in _lineages to rows of dij
        // Also save starting edge lengths so they can be restored in destroyUPGMA()
        map<Node *, unsigned> row;
        _upgma_starting_edgelen.clear();
        for (unsigned i = 0; i < n; i++) {
            Node * nd = _lineages[i];
            _upgma_starting_edgelen[nd] = nd->_edge_length;
            row[nd] = i;
        }

        
//        debugShowDistanceMatrix(dij);
        
        double upgma_height = 0.0;
        
        // Build UPGMA tree on top of existing forest
        assert(_upgma_additions.empty());
        unsigned nsteps = n - 1;
        while (nsteps > 0) {
            // Find smallest entry in d
            auto it = min_element(dij.begin(), dij.end());
            unsigned offset = (unsigned)std::distance(dij.begin(), it);
            auto p = dij_row_col.at(offset);
            unsigned i = p.first;
            unsigned j = p.second;
            
            // Update all leading edge lengths
            double v = *it;
            for (auto nd : _lineages) {
                nd->_edge_length += (0.5*v - upgma_height);
            }
            
            upgma_height = v / 2.0;
            
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
            calcPartialArray(new_nd);
                        
            // Update distance matrix
            for (unsigned k = 0; k < n; k++) {
                if (k != i && k != j) {
                    unsigned ik = (i > k) ? (i*(i-1)/2 + k) : (k*(k-1)/2 + i);
                    unsigned jk = (j > k) ? (j*(j-1)/2 + k) : (k*(k-1)/2 + j);
                    double a = dij[ik];
                    double b = dij[jk];
                    dij[ik] = 0.5*(a + b);
                    dij[jk] = _infinity;
                }
//                debugShowDistanceMatrix(dij);
            }
            
//            debugShowDistanceMatrix(dij);
            
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
            dij2.assign(dim2, _infinity);
            
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
                    
//                    debugShowDistanceMatrix(dij);
//                    debugShowDistanceMatrix(dij2);
                }
            }
            
//            debugShowDistanceMatrix(dij);
//             debugShowDistanceMatrix(dij2);
            
            // Set up for next iteration
            dij = dij2;
            n = n2;
            for (unsigned i = 0; i < n; i++) {
                Node * nd = _lineages[i];
                row[nd] = i;
            }
            
//            debugShowDistanceMatrix(dij);
            
            --nsteps;
        }
        
        // debugging output
        // output(format("\nGene forest for locus \"%s\" after UPGMA:\n%s\n") % gene_name % makeNewick(9, /*use_names*/true, /*coalunits*/false), 0);
        // output(format("  Height after UPGMA = %g\n") % _forest_height, 0);
        
        _gene_tree_log_likelihood = calcLogLikelihood();
        _log_weight = _gene_tree_log_likelihood - prev_log_likelihood; // previous likelihood is the entire tree
                        
        if (_save_memory) {
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
            
            if (_save_memory) {
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
        
        // output("\nIn GeneForest::destroyUPGMA:\n", 0);
        // output(format("  Height before refreshAllHeightsAndPreorders = %g\n") % _forest_height, 0);
        // refreshAllHeightsAndPreorders();
        // output(format("  newick = %s\n") % makeNewick(9, /*use_names*/true, /*coalunits*/false), 0);
        // output(format("  Height after refreshAllHeightsAndPreorders = %g\n") % _forest_height, 0);
        // output("\n", 0);
#endif
    }
#endif
        
#if defined (FASTER_UPGMA_TREE)
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
                assert (_lineages.size() == _ntaxa);
                unsigned n = _ntaxa;
                vector<double> dij(n*(n-1)/2, _infinity);
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
                        
                        assert (v != _infinity);
                        
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
                vector<double> dij(n*(n-1)/2, _infinity);
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
                            for (unsigned lstate = 0; lstate < _nstates; lstate++) {
                                auto & l_partial_array = *(lnode->_partial);
                                double lpartial = l_partial_array[p*_nstates + lstate];
                                for (unsigned rstate = 0; rstate < _nstates; rstate++) {
                                    auto & r_partial_array = *(rnode->_partial);
                                    double rpartial = r_partial_array[p*_nstates + rstate];
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
        
#if defined (FASTER_UPGMA_TREE)
        inline void Forest::buildStartingRow() {
            unsigned n = _ntaxa;
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
        
#if defined (FASTER_UPGMA_TREE)
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
            
//            // Create distance matrix dij and workspace dij2 used to build next dij
//            // Both dij and dij2 are 1-dimensional vectors that store only the
//            // lower diagonal of the distance matrix (excluding diagonal elements)
            
            unsigned temp1 = _lineages.back()->_left_child->_right_sib->_position_in_lineages;
            unsigned temp2 = _lineages.back()->_left_child->_position_in_lineages;
            unsigned i_to_delete = temp1;
            unsigned j_to_delete = temp2; // TODO: does i need to be the larger number?

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
                        dij[jk] = _infinity;
                    }
                }
            
                // Build new distance matrix
                 unsigned n2 = (unsigned)_lineages.size();
                assert(n2 == n - 1);
                unsigned dim2 = n2*(n2-1)/2;
                dij2.resize(dim2);
                dij2.assign(dim2, _infinity);
            
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
            
//        double upgma_height = 0.0;
        
        // Build UPGMA tree on top of existing forest
        assert(_upgma_additions.empty());
            
            double upgma_height = getLineageHeight(_lineages.back());
                        
        unsigned nsteps = n - 1;
        while (nsteps > 0) {
            // Find smallest entry in d
            auto it = min_element(dij.begin(), dij.end());
            unsigned offset = (unsigned)std::distance(dij.begin(), it);
            auto p = dij_row_col.at(offset);
            unsigned i = p.first;
            unsigned j = p.second;
            
            // Update all leading edge lengths
            double v = *it; // TODO: need to subtract existing node height?
            
            assert (v != _infinity);
            
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
            calcPartialArray(new_nd);
                        
            // Update distance matrix
            for (unsigned k = 0; k < n; k++) {
                if (k != i && k != j) {
                    unsigned ik = (i > k) ? (i*(i-1)/2 + k) : (k*(k-1)/2 + i);
                    unsigned jk = (j > k) ? (j*(j-1)/2 + k) : (k*(k-1)/2 + j);
                    double a = dij[ik];
                    double b = dij[jk];
                    dij[ik] = 0.5*(a + b);
                    dij[jk] = _infinity;
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
            dij2.assign(dim2, _infinity);
            
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
                        
        if (_save_memory) {
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
            
            if (_save_memory) {
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
                    if (v == Forest::_infinity)
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
           double log_weight_choices_sum = getRunningSumChoices(log_weight_choices); // TODO: just push back node numbers to _node_choices?
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
     calcPartialArray(new_nd);

     // don't update the species list
     updateNodeVector(_lineages, subtree1, subtree2, new_nd);
              
     return make_tuple(subtree1, subtree2, new_nd);
 }

    inline void Forest::allowCoalescence(string species_name, double increment, Lot::SharedPtr lot) {
         double prev_log_likelihood = _gene_tree_log_likelihood;

         Node *subtree1 = nullptr;
         Node *subtree2 = nullptr;
         list<Node*> nodes;

        nodes = _species_partition[species_name];
        
        assert (nodes.size() > 0);

         unsigned s = (unsigned) nodes.size();
         calcTopologyPrior(s);

         assert (s > 1);
         bool one_choice = false;
         if (nodes.size() == 2) {
             one_choice = true;
         }

         if (Forest::_proposal == "prior-post" && (!one_choice)) {
             if (_save_memory) {
                 for (auto &nd:_lineages) {
    //                for (auto &nd:nodes) {
                     if (nd->_partial == nullptr) {
                         nd->_partial = ps.getPartial(_npatterns*4);
                         calcPartialArray(nd);
                     }
                 }
             }
             
             pair<Node*, Node*> t = chooseAllPairs(nodes, increment, species_name, lot);
             
             subtree1 = t.first;
             subtree2 = t.second;
         }
         
         else {
             assert (Forest::_proposal == "prior-prior" || one_choice);
             // prior-prior proposal
             pair<unsigned, unsigned> t = chooseTaxaToJoin(s, lot);
             auto it1 = std::next(nodes.begin(), t.first);
             subtree1 = *it1;

             auto it2 = std::next(nodes.begin(), t.second);
             subtree2 = *it2;
             assert (t.first < nodes.size());
             assert (t.second < nodes.size());

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

         new_nd->_left_child=subtree1;
         subtree1->_right_sib=subtree2;

         subtree1->_parent=new_nd;
         subtree2->_parent=new_nd;

         if (!_run_on_empty) {
             //always calculating partials now
             assert (new_nd->_partial == nullptr);
             new_nd->_partial=ps.getPartial(_npatterns*4);
             assert(new_nd->_left_child->_right_sib);

             if (_save_memory) {
                 for (auto &nd:_lineages) {
                     if (nd->_partial == nullptr) {
                         nd->_partial = ps.getPartial(_npatterns*4);
                         calcPartialArray(nd);
                     }
                 }
             }
             calcPartialArray(new_nd);

             subtree1->_partial=nullptr; // throw away subtree partials now, no longer needed
             subtree2->_partial=nullptr;
         }

             //update species list
             updateNodeList(nodes, subtree1, subtree2, new_nd);
             updateNodeVector(_lineages, subtree1, subtree2, new_nd);

        _species_partition[species_name] = nodes;

# if !defined (BUILD_UPGMA_TREE)
             if ((_proposal == "prior-prior" || one_choice) && (!_run_on_empty) ) {
                 _gene_tree_log_likelihood = calcLogLikelihood();
                 _log_weight = _gene_tree_log_likelihood - prev_log_likelihood;
             }
        
        if (_save_memory) {
            for (auto &nd:_nodes) {
                nd._partial=nullptr;
            }
        }
        
#endif
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

    inline void Forest::addSpeciesIncrement() {
        // add the previously chosen edge length
        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
    }

    inline double Forest::getTreeHeight() {
        double sum_height = 0.0;

        // calculate height of lineage
        Node* base_node = _lineages[0];
        sum_height += base_node->getEdgeLength();
        for (Node* child=base_node->_left_child; child; child=child->_left_child) {
            sum_height += child->getEdgeLength();
        }
        return sum_height;
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

    inline double Forest::getSpeciesTreeIncrement() {
        assert (_index == 0);
        return _cum_height;
    }

    inline vector<pair<double, string>> Forest::calcForestRate(Lot::SharedPtr lot) {
        vector<pair<double, string>> rates;
        pair<double, string> rate_and_name;

        for (auto &s:_species_partition) {
            if (s.second.size() > 1) { // if size == 0, no possibility of coalescence and rate is 0
                double population_coalescence_rate = 0.0;
#if defined (DRAW_NEW_THETA)
                    double population_theta = _theta_map[s.first];
//#if defined (RATE_HET_SIM)
//                if (_index == 1) {
//                    population_theta *= 10.0;
//                }
//#endif
                population_coalescence_rate = s.second.size()*(s.second.size()-1)/population_theta;
#else
                population_coalescence_rate = s.second.size()*(s.second.size()-1)/_theta;
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
        if (_index == 0) {
            _last_edge_length = increment;
            _cum_height += increment;
        }
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
        
        // set ancestral species name for use in calculating panmictic coalescent likelihood
//        number = _nspecies - 2;
//        _ancestral_species_name = boost::str(boost::format("node-%d")%number);
        
        for (int i=0; i<_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            species_names.push_back(name);
        }
        
        _ancestral_species_name = species_names.back();
        
        assert (species_names.size() == 2*_nspecies - 1);
        
        // draw thetas for tips of species trees and ancestral population
        // for all other populations, theta = -1
        
//        assert (_theta_proposal_mean > 0.0);
        if (_theta_proposal_mean == 0.0) {
            assert (_theta > 0.0);
            _theta_proposal_mean = _theta;
        }
        double scale = 1 / _theta_proposal_mean;
        
        unsigned count = 0;
        for (auto &name:species_names) {
            if (count < _nspecies || count == 2*_nspecies-2) {
                double new_theta = 0.0;
                if (new_theta < _small_enough) {
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
        if (new_theta < _small_enough) {
//            new_theta = 1 / rng.gamma(2.0, scale);
            new_theta = 1 / lot->gamma(2.0, scale);
            _theta_map[new_species] = new_theta;
        }
        // pop mean = theta / 4
        double a = 2.0;
        double b = scale;
//            double x = new_theta / 4.0;
        double x = new_theta;
        double log_inv_gamma_prior = (a*log(b) - lgamma(a) - (a+1)*log(x) - b/x);
        _vector_prior.push_back(log_inv_gamma_prior);
    }

    inline void Forest::updateThetaMap(Lot::SharedPtr lot, string new_species_name) {
        // add a new theta for the most recently drawn species
//        double scale = (2.01 - 1.0) / (_theta_mean);
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        double new_theta = 0.0;
        if (new_theta < _small_enough) {
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
        _theta_map[new_species_name] = Forest::_theta;
    }
        
    inline void Forest::createSpeciesIndices() {
        unsigned number = 0;
        for (auto &s:_species_partition) {
            number++;
            _species_indices[s.first] = number - 1;
        }
        
        for (int i=0; i<_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            _species_indices[name] = number - 1;
        }
    }
        
    inline void Forest::createThetaMapFixedTheta(Lot::SharedPtr lot) {
        // map should be 2*nspecies - 1 size
        unsigned number = 0;
        _species_names.clear();
        
        for (auto &s:_species_partition) {
            _species_names.push_back(s.first);
            number++;
            _species_indices[s.first] = number - 1;
        }
        assert (_species_names.size() == _nspecies);
        
        for (int i=0; i<_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            _species_names.push_back(name);
            _species_indices[name] = number - 1;
        }
        
        _theta_mean = Forest::_theta;
        
        for (auto &name:_species_names) {
            _theta_map[name] = Forest::_theta;
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
            _species_indices[s.first] = number - 1;
        }
        assert (_species_names.size() == _nspecies);
        
        for (int i=0; i<_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            _species_names.push_back(name);
            _species_indices[name] = number - 1;
        }
        
        // gamma mean = shape * scale
        // draw mean from lognormal distribution
        // shape = 2.0 to be consistent with starbeast3
        // scale = 1 / mean;
        
        if (_theta_proposal_mean > 0.0) {
            assert (_theta_mean == 0.0);
            _theta_mean = lot->gamma(1, _theta_proposal_mean); // equivalent to exponential(exponential_rate)
        }
        else {
            _theta_mean = Forest::_theta; // if no proposal distribution specified, use one theta mean for all particles
        }
        
//        double scale = (2.01 - 1.0) / (_theta_mean);
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        for (auto &name:_species_names) {
            double new_theta = 0.0;
            if (new_theta < _small_enough) {
                new_theta = 1 / (lot->gamma(2.0, scale));
//                new_theta = 1 / (lot->gamma(2.01, scale));
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
                            double coalescence_rate = nlineages*(nlineages-1) / Forest::_theta;
#endif
                            double nChooseTwo = nlineages*(nlineages-1);
                            double log_prob_join = log(2/nChooseTwo);
                            log_coalescent_likelihood += log_prob_join + log(coalescence_rate) - (increment * coalescence_rate);

    //                            cout << "pr coalescence = " << log_prob_join + log(coalescence_rate) - (increment * coalescence_rate) << endl;

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
                            double coalescence_rate = nlineages*(nlineages-1) / _theta;
#endif
                            log_coalescent_likelihood -= increment * coalescence_rate;

    //                            cout << "pr no coalescence = " << -1 * increment * coalescence_rate << endl;
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
                double coalescence_rate = nlineages*(nlineages-1) / _theta;
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
                double coalescence_rate = panmictic_nlineages*(panmictic_nlineages-1) / _theta;
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
                        double coalescence_rate = nlineages*(nlineages-1) / _theta;
#endif
                        double nChooseTwo = nlineages*(nlineages-1);
                        double log_prob_join = log(2/nChooseTwo);
                        double increment = heights_and_nodes[i].first - species_tree_height - cum_time;
                        log_coalescent_likelihood += log_prob_join + log(coalescence_rate) - (increment * coalescence_rate);

    //                        cout << "pr coalescence = " << log_prob_join + log(coalescence_rate) - (increment * coalescence_rate) << endl;

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
//#endif
    }

    inline pair<vector<double>, vector<unsigned>> Forest::calcCoalescentLikelihoodIntegratingOutThetaLastStep(vector<pair<tuple<string,string,string>, double>> species_build) {
        vector<double> gamma_b; // contains gamma_b by species
        vector<unsigned> q_b; // contains q_b by species
        
        gamma_b.resize(_nspecies + species_build.size());
        q_b.resize(_nspecies + species_build.size());
        
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
                                    gamma_b[count] += 4* increment * nlineages *(nlineages-1) / (2*_ploidy);
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
                                    gamma_b[count] += 4*increment * nlineages *(nlineages-1) / (2*_ploidy);
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
                        gamma_b.back() += 4*increment * panmictic_nlineages *(panmictic_nlineages-1) / (2*_ploidy);

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
            gamma_b.back() += 4*increment * panmictic_nlineages *(panmictic_nlineages-1) / (2*_ploidy);

            cum_time += increment;
            panmictic_nlineages--;
        }
        assert (panmictic_nlineages == 1);
    
        assert (gamma_b.size() == q_b.size());
        }
        
        return make_pair(gamma_b, q_b);
    }

    inline pair<vector<double>, vector<unsigned>> Forest::calcCoalescentLikelihoodIntegratingOutTheta(vector<pair<tuple<string,string,string>, double>> species_build) { // TODO: issues when newicks read in have very small branch lengths
         vector<double> gamma_b; // contains gamma_b by species
        vector<unsigned> q_b; // contains q_b by species
        
        gamma_b.resize(_nspecies + species_build.size());
        q_b.resize(_nspecies + species_build.size());
        
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
                                    gamma_b[count] += 4* increment * nlineages *(nlineages-1) / (2*_ploidy);
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
                                    gamma_b[count] += 4*increment * nlineages *(nlineages-1) / (2*_ploidy);
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
                        gamma_b.back() += 4*increment * panmictic_nlineages *(panmictic_nlineages-1) / (2*_ploidy);

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

    inline vector<pair<double, pair<string, string>>> Forest::getMinDepths() {
        return _depths;
    }

    inline void Forest::calcMinDepth() {
        assert (_index != 0);
        _depths.clear();
        // walk through nodes
        // find children of node and determine if children are in the same species
        // if children are in different species, calculate the node height
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
    //        vector<double> depths;
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
                                left_child = left_child->_left_child; // TODO: not sure about these; will the parent get accounted for later?
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
                        double height = getLineageHeight(nd->_left_child);
                        _depths.push_back(make_pair(height, make_pair(spp_left_child, spp_right_child)));
    //                        _depths.push_back(height);
                    }
                    spp_left_child = "";
                    spp_right_child = "";
                }
            }
        }
        assert (_depths.size() > 0);
    //        return depths[0];
    }

    inline void Forest::refreshPreorder() {
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
                double height = getLineageHeight(nd->_left_child); //
                heights_and_nodes.push_back(make_pair(height, nd));
//                nd->_height = nd->_left_child->_height + nd->_left_child->_edge_length;
//                heights_and_nodes.push_back(make_pair(nd->_height, nd));
            }
            else {
                // if leaf node, initialize _height to zero
//                nd->_height = 0.0;
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
        if (_model != "JC") {
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
        if (Forest::_comphet != Forest::_infinity) {
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
                if (Forest::_asrv_shape != _infinity)
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
            if (t < _ntaxa) {
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
                        if (!nd->_left_child->_right_sib)
                            throw XProj(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
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
                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );

//        _lineages.clear();
//        _lineages.push_back(&_nodes.back());

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
            if (count < _nspecies) {
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
                        if (!nd->_left_child->_right_sib)
                            throw XProj(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
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
                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } );
        _lineages.clear();
        
        _lineages.push_back(&_nodes.back());
        
        // reset node names
        int j = 0;
        for (auto &nd:_nodes) {
            nd._number = j;
            j++;
        }
    }

}

