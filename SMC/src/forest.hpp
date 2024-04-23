#pragma once

#include <stack>
#include <memory>
#include <iostream>
#include <boost/format.hpp>
#include <vector>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <thread>
#include <mutex>
#include <algorithm>
#include "conditionals.hpp"
#include <fstream>

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

        typedef std::vector <double> partial_array_t;
        void                        clear();
        void                        setData(Data::SharedPtr d, int index, map<string, string> &taxon_map, bool partials);
        void                        setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa);
        Node *                      findNextPreorder(Node * nd);
        string                      makeNewick(unsigned precision, bool use_names);
        string                      makePartialNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>    chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        tuple<Node*, Node*, Node*>  createNewSubtree(pair<unsigned, unsigned> p, list<Node*> node_list, string species_name, double increment);
        void                        calcPartialArray(Node* new_nd);
        void                        setUpGeneForest(map<string, string> &taxon_map);
        void                        setUpSpeciesForest(vector<string> &species_names);
        tuple<string,string, string> speciesTreeProposal(Lot::SharedPtr lot);
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                        hybridizeNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node* delnode3, Node * addnode1);
        void                        revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                      getRunningSumChoices(vector<double> &log_weight_choices);
        double                      getRunningSumHybridChoices(vector<double> &log_weight_choices);
        vector<double>              reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        pair<Node*, Node*>          getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list);
        int                         selectPair(vector<double> weight_vec, Lot::SharedPtr lot);
        void                        chooseSpeciesIncrement(Lot::SharedPtr lot);
        void                        chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth);
        void                        addSpeciesIncrement();
        string                      chooseEvent();
        void                        allowMigration(list<Node*> &nodes);
        double                      chooseTaxonToMigrate(double s);
        string                      findKeyToDel(Node* taxon_to_migrate);
        void                        migrateTaxon(unsigned taxon_choice, string key_to_del, Node* taxon_to_migrate);
        string                      chooseLineage(Node* taxon_to_migrate, string key_to_del);
        void                        addMigratingTaxon(string key_to_add, string key_to_del, Node* taxon_to_migrate);
        void                        deleteTaxon(string key_to_del, unsigned taxon_choice);
        void                        allowCoalescence(string species_name, double increment, Lot::SharedPtr lot);
        tuple<unsigned, unsigned, unsigned> chooseTaxaToHybridize();
        vector<string>              hybridizeSpecies();
        void                        moveGene(string new_nd, string parent, string hybrid);
        void                        rebuildSpeciesPartition(vector<string> names, vector<list<Node*>> nodes);
        void                        switchParents(string parent, string parent2);
        void                        resetLineages(vector<double> branch_lengths);
        vector<double>              saveBranchLengths();
        int                         chooseDirectionOfHybridization(vector<double> likelihood_vec, Lot::SharedPtr lot);
        void                        hybridGeneTreeProposal(double species_tree_increment, string species_name);
        vector<pair<double, string>>      calcForestRate(Lot::SharedPtr lot);
        void                        updateSpeciesPartition(tuple<string, string, string> species_info);
        double                      calcTopologyPrior(unsigned nlineages);
        void                        calcIncrementPrior(double increment, string species_name, bool new_increment, bool coalesced_gene, bool gene_tree);
        void                        clearPartials();
        void                        setStartMode(string mode) {_start_mode = mode;}
        unsigned                    getDeepCoal(tuple <string, string, string> species_joined);
        void                        resetDepthVector(tuple<string, string, string> species_joined);
        vector<pair<double, pair<string, string>>>             getMinDepths();
        void                        calcMinDepth();
        vector< pair<double, Node *>>      sortPreorder();
        void                        refreshPreorder();
        void                        createThetaMap(Lot::SharedPtr lot);
        void                        resetThetaMap(Lot::SharedPtr lot);
        void                        drawNewTheta(string new_species, Lot::SharedPtr lot);
        void                        buildFromNewick(const string newick, bool rooted, bool allow_polytomies);
        void                        stripOutNexusComments(std::string & newick);
        unsigned                    countNewickLeaves(const std::string newick);
        void                        extractEdgeLen(Node * nd, std::string edge_length_string);
        void                        renumberInternals();
        bool                        canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);
    
        map<string, double>              _theta_map;

        std::vector<Node *>         _lineages;
    
        std::list<Node>             _nodes;
        std::vector<Node*>          _new_nodes;

        unsigned                    _nleaves;
        unsigned                    _ninternals;
        unsigned                    _npatterns;
        unsigned                    _nstates;
        double                      _last_edge_length;
        vector<double>              _gamma;

        Data::SharedPtr             _data;
        static unsigned             _nspecies;
        static unsigned             _ntaxa;
        unsigned                    _first_pattern = 0;
        unsigned                    _index;
        map<string, list<Node*> >   _species_partition;
        double                      _gene_tree_log_likelihood;
        vector<pair<Node*, Node*>>  _node_choices;
        vector<double>              _log_likelihood_choices;
        int                         _index_of_choice;
        pair<Node*, Node*>          _species_joined;
        tuple<Node*, Node*, Node*>  _hybrid_species_joined;
        string                      _last_direction;
        double                      _log_joining_prob;
        vector<pair<double, double>> _increments_and_priors;
        bool                        _done;
        double                      _log_coalescent_likelihood;
        double                      _panmictic_coalescent_likelihood;
        double                      _log_coalescent_likelihood_increment;
        double                      _cum_height;
        vector<string>              _species_for_coalescent_events;
        string                      _start_mode;
        vector<pair<double, double>>         _increments;
        vector<pair<double, pair<string, string>>>              _depths;
        unsigned                    _nincrements = 0;
        vector<Node*>               _preorder;
        double                      _small_enough;
        vector<pair<tuple<string,string,string>, double>>    _species_build;
        map<string, string>         _taxon_map;
        map<string, unsigned>       _species_indices;
    
        void                        showSpeciesJoined();
        double                      calcTransitionProbability(Node* child, double s, double s_child);
        double                      calcSimTransitionProbability(unsigned from, unsigned to, double edge_length);
        double                      calculateNewEdgeLength(string key_to_add, Node* taxon_to_migrate);
        void                        setNewEdgeLength(double difference, Node* taxon_to_migrate, string key_to_add);
        void                        hybridizeGene(vector<string> hybridized_nodes, double species_tree_increment, string species_name, Lot::SharedPtr lot);
        void                        resetToMinor(vector<Node*> minor_nodes, vector<Node*> minor_left_children, vector<Node*> minor_right_children, vector<double> minor_left_edge_lengths, vector<double> minor_right_edge_lengths);
        double                      getTreeHeight();
        double                      getTreeLength();
        double                      getSpeciesTreeIncrement();
        double                      getLineageHeight(Node* nd);
        double                      getTotalLineageHeight(Node* nd);
        double                      _log_weight;
        double                      _other_log_weight;
        void                        addIncrement(double increment);
        void                        simulateData(Lot::SharedPtr lot, unsigned starting_site, unsigned nsites);
        double                      calcCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height);
        pair<vector<double>, vector<unsigned>>        calcCoalescentLikelihoodIntegratingOutTheta(vector<pair<tuple<string,string,string>, double>> species_build);
        inline pair<vector<double>, vector<unsigned>> calcInitialCoalescentLikelihoodIntegratingOutTheta();
        pair<vector<double>, vector<unsigned>>        calcCoalescentLikelihoodIntegratingOutThetaLastStep(vector<pair<tuple<string,string,string>, double>> species_build);
        string                      _ancestral_species_name;
        vector<double>              _vector_prior;
        double                      _scale;
        double                      _theta_mean;

    public:

        typedef std::shared_ptr<Forest> SharedPtr;
        static double               _theta;
        static double               _theta_proposal_mean;
        static double               _theta_constant_mean;
        static double               _theta_prior_mean;
        static double               _lambda;
        static string               _proposal;
        static string               _model;
        static double               _kappa;
        static vector<double>       _base_frequencies;
        static string               _string_base_frequencies;
        static double               _migration_rate;
        static double               _hybridization_rate;
        static bool                 _save_memory;
        static string               _outgroup;
        static bool                 _run_on_empty;
        static double               _ploidy;
        static double               _gamma_scale;
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
        _other_log_weight = 0.0;
        _cum_height = 0.0;
        _nleaves=_ntaxa;
        _ninternals=0;
        _increments.clear();
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
        _scale = 0.0;
        
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
        if (nd->_major_parent) { // TODO: not sure
            next = nd->_parent->_right_sib;
        }
        else if (!nd->_left_child && !nd->_right_sib) {
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
        cout << " " << makeNewick(9, true) << "\n";
        cout << "\n";
    }

    inline string Forest::makePartialNewick(unsigned precision, bool use_names) {
            // this function makes a newick string for a partially constructed tree
            string newick = "(";
            const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
            const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
            const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
            stack<Node *> node_stack; // TODO: this function doesn't work with >1 cycles

        // find hybrid nodes, don't visit minor or major parent until hybrid node is visited
        vector<Node*> hybrid_nodes;
        int n = _nspecies + 1;
        for (auto &node:_nodes) {
            if (node._parent2) {
                hybrid_nodes.push_back(&node);
                string name = "#H" + to_string(n);
                node._hybrid_newick_name = name;
                n++;
            }
        }

            unsigned i = 0;
            unsigned a = 0;
            for (auto &lineage : _lineages) {
                Node * nd = lineage;
                while (nd) {
                    bool skip = false;
                    for (auto &i:hybrid_nodes) {
                        if (nd == i->_major_parent || nd == i->_minor_parent) {
                            skip = true;
                        }
                    }
                    if (nd->_minor_parent && !nd->_visited && !skip) {
                        a++;
                        // hybrid node with minor parent
                        if (use_names) {
                            newick += "(";
                            newick += boost::str(boost::format(tip_node_name_format)
                                % nd->_minor_parent->_name
                                % nd->_minor_parent->_edge_length);
                            newick += ",";
                            newick += boost::str(boost::format(tip_node_name_format)
                                % nd->_hybrid_newick_name
                                % nd->_edge_length);
                                nd->_minor_parent->_visited = true;
                            newick += ")";
                        }

                        // hybrid node with major parent
                        if (use_names) {
                            newick += ",(";
                            newick += boost::str(boost::format(tip_node_name_format)
                                % nd->_major_parent->_name
                                % nd->_major_parent->_edge_length);
                            newick += ",(";
    //                            newick += "#H_";
                            newick += boost::str(boost::format(tip_node_name_format)
                                % nd->_name
                                % nd->_edge_length);
                            newick += "),";
                            newick += boost::str(boost::format(tip_node_name_format)
                                % nd->_hybrid_newick_name
                                % nd->_edge_length);
                                nd->_major_parent->_visited = true; // TODO: I think this only works if major and minor parents are tip nodes
                        }

                        else {

                        }

                    }
    //                    else if (nd->_left_child && !nd->_visited && !skip) {
                    if (nd->_left_child) {
                        a++;
                        nd->_visited = true;
                        // internal node
                        newick += "(";
                        node_stack.push(nd);
                    }
    //                    else if (!nd->_left_child && !nd->_visited && !skip) {
                    else {
                        a++;
                        nd->_visited = true;
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
                    if (a >= _ninternals + _nleaves - 1 && hybrid_nodes.size()>0) {
                        break;
                    }
                    nd->_visited = true;
                    nd = findNextPreorder(nd);
                }   // while (subnd)...

                if (i < _lineages.size() - 1)
                    newick += ",";
                ++i;
            }
        if (hybrid_nodes.size() >0) {
            newick.pop_back();
            newick.pop_back();
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
                stack<Node *> node_stack; // TODO: this function doesn't work with >1 cycles

                // find hybrid nodes, don't visit minor or major parent until hybrid node is visited
                vector<Node*> hybrid_nodes;
                int n = _nspecies + 1;
                for (auto &node:_nodes) {
                    if (node._parent2) {
                        hybrid_nodes.push_back(&node);
                        string name = "#H" + to_string(n);
                        node._hybrid_newick_name = name;
                        n++;
                    }
                }

                    unsigned i = 0;
                    unsigned a = 0;
                    for (auto &lineage : _lineages) {
                        Node * nd = lineage;
                        while (nd) {
                            bool skip = false;
                            for (auto &i:hybrid_nodes) {
                                if (nd == i->_major_parent || nd == i->_minor_parent) {
                                    skip = true;
                                }
                            }
                            if (nd->_minor_parent && !nd->_visited && !skip) {
                                a++;
                                // hybrid node with minor parent
                                if (use_names) {
                                    newick += "(";
                                    newick += boost::str(boost::format(tip_node_name_format)
                                        % nd->_minor_parent->_name
                                        % nd->_minor_parent->_edge_length);
                                    newick += ",";
                                    newick += boost::str(boost::format(tip_node_name_format)
                                        % nd->_hybrid_newick_name
                                        % nd->_edge_length);
            //                                nd->_minor_parent->_visited = true;
                                    newick += ")";
                                }

                                // hybrid node with major parent
                                if (use_names) {
                                    newick += ",(";
                                    newick += boost::str(boost::format(tip_node_name_format)
                                        % nd->_major_parent->_name
                                        % nd->_major_parent->_edge_length);
                                    newick += ",(";
            //                            newick += "#H_";
                                    newick += boost::str(boost::format(tip_node_name_format)
                                        % nd->_name
                                        % nd->_edge_length);
                                    newick += "),";
                                    newick += boost::str(boost::format(tip_node_name_format)
                                        % nd->_hybrid_newick_name
                                        % nd->_edge_length);
            //                                nd->_major_parent->_visited = true; // TODO: I think this only works if major and minor parents are tip nodes
                                }

                                else {

                                }

                            }
            //                    else if (nd->_left_child && !nd->_visited && !skip) {
                            if (nd->_left_child) {
                                a++;
            //                        nd->_visited = true;
                                // internal node
                                newick += "(";
                                node_stack.push(nd);
                            }
            //                    else if (!nd->_left_child && !nd->_visited && !skip) {
                            else {
                                a++;
            //                        nd->_visited = true;
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
                            if (a >= _ninternals + _nleaves - 1 && hybrid_nodes.size()>0) {
                                break;
                            }
            //                    nd->_visited = true;
                            nd = findNextPreorder(nd);
                        }   // while (subnd)...

                        if (i < _lineages.size() - 1)
                            newick += ",";
                        ++i;
                    }
                if (hybrid_nodes.size() >0) {
                    newick.pop_back();
                    newick.pop_back();
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
                                                        
    inline tuple<unsigned, unsigned, unsigned> Forest::chooseTaxaToHybridize(){
        double nsubtrees = _lineages.size();
        unsigned t1;
        unsigned t2;
        unsigned t3;
        //don't use this when there's only one choice (2 subtrees)
        // thread safe random number generator with mutex
        mtx.lock();
//        if (nsubtrees > 3) {
            t1 = ::rng.randint(0, nsubtrees-1);
            t2 = ::rng.randint(0, nsubtrees-1);
            t3 = ::rng.randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1 or t3
            while (t2 == t1 || t2 == t3) {
                t2 = ::rng.randint(0, nsubtrees-1);
            }
            // keep calling t3 until it doesn't equal t1 or t2
            while (t3 == t1 || t3 == t2) {
                t3 = ::rng.randint(0, nsubtrees-1);
            }
//        }
        mtx.unlock();
        return make_tuple(t1, t2, t3);
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
        
        for (auto &nd:_lineages) {
            assert (nd->_right_sib != nd);
        }
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

        inline double Forest::calcSimTransitionProbability(unsigned from, unsigned to, double edge_length) {
            double transition_prob = 0.0;
            assert (_model == "JC");
            if (from == to) {
                transition_prob = 0.25 + 0.75*exp(-4.0*edge_length/3.0);
            }

            else {
                transition_prob = 0.25 - 0.25*exp(-4.0*edge_length/3.0);
            }
            return transition_prob;
        }

    inline double Forest::calcTransitionProbability(Node* child, double s, double s_child) {
        double child_transition_prob = 0.0;

        if (_model == "JC" ) {
            double expterm = exp(-4.0*(child->_edge_length)/3.0);
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
            double beta_t = 0.5*(child->_edge_length)/phi;

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

    inline double Forest::calcLogLikelihood() {
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

    inline double Forest::getRunningSumHybridChoices(vector<double> &log_weight_choices) {
        double running_sum = 0.0;
        double log_weight_choices_sum = 0.0;
        vector<double> adjustedLogLikelihood;
        
        // multiply major likelihood * gamma, multiply minor likelihood * (1 - gamma)
        adjustedLogLikelihood.push_back(log_weight_choices[0]+0.15);
        adjustedLogLikelihood.push_back(log_weight_choices[1]+0.85);
        
        double log_max_weight = *max_element(adjustedLogLikelihood.begin(), adjustedLogLikelihood.end());
        for (auto & i:adjustedLogLikelihood) {
            running_sum += exp(i - log_max_weight);
        }
        
        log_weight_choices_sum = log(running_sum) + log_max_weight;
        return log_weight_choices_sum;
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
        _node_choices.push_back(s);

        return s;
    }

    inline unsigned Forest::getDeepCoal(tuple <string, string, string> species_joined) {
        unsigned num_deep_coal = 0;
//        if (_species_partition.size() > 2) { // don't count ancestral population as deep coalescence
            unsigned count1 = 0;
            unsigned count2 = 0;
            
            // first lineage
            for (auto &nd:_species_partition[get<0>(species_joined)]) {
    //            if (nd->_deep_coalescence_counted) {
    //                count1 = 1; // avoid double counting
    //                break;
    //            }
    //            else {
                    count1 += 1;
                    nd->_deep_coalescence_counted = true;
    //            }
            }
            
            // ensure all nodes are accounted for
            for (auto &nd:_species_partition[get<0>(species_joined)]) {
                nd->_deep_coalescence_counted = true;
            }
            
            // second lineage
            for (auto &nd:_species_partition[get<1>(species_joined)]) {
    //            if (nd->_deep_coalescence_counted) {
    //                count2 = 1; // avoid double counting
    //                break;
    //            }
    //            else {
                    count2 += 1;
                    nd->_deep_coalescence_counted = true;
    //            }
            }
            
            // ensure all nodes are accounted for
            for (auto &nd:_species_partition[get<1>(species_joined)]) {
                nd->_deep_coalescence_counted = true;
            }
            
            num_deep_coal = (count1 + count2) - 1;
            
    //        for (auto &nd:_species_partition[new_spp]) {
    //            if (!nd->_deep_coalescence_counted) {
    //                num_deep_coal += 1;
    //                nd->_deep_coalescence_counted = true;
    //            }
    //        }
            //                _num_deep_coalescences += _forests[i]._species_partition[new_spp].size() - 1;
//        }
        return num_deep_coal;
    }

    inline tuple<Node*, Node*, Node*> Forest::createNewSubtree(pair<unsigned, unsigned> t, list<Node*> node_list, string species_name, double increment) {
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

        // update node vector
        // don't update the species list
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        return make_tuple(subtree1, subtree2, new_nd);
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
        _index_of_choice    = other._index_of_choice;
        _node_choices = other._node_choices;
        _log_likelihood_choices = other._log_likelihood_choices;
        _species_joined = other._species_joined;
        _hybrid_species_joined = other._hybrid_species_joined;
        _migration_rate = other._migration_rate;
        _hybridization_rate = other._hybridization_rate;
        _last_direction = other._last_direction;
        _gamma = other._gamma;
        _log_weight = other._log_weight;
        _log_joining_prob = other._log_joining_prob;
        _increments_and_priors = other._increments_and_priors;
        _done = other._done;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _log_coalescent_likelihood_increment = other._log_coalescent_likelihood_increment;
        _other_log_weight = other._other_log_weight;
        _cum_height = other._cum_height;
        _species_for_coalescent_events = other. _species_for_coalescent_events;
        _outgroup = other._outgroup;
        _run_on_empty = other._run_on_empty;
        _start_mode = other._start_mode;
        _increments = other._increments;
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
        _scale = other._scale;
        _gamma_scale = other._gamma_scale;

        // copy tree itself

        _species_partition.clear();
        for (auto spiter : other._species_partition) {
            for (auto s : spiter.second) {
                unsigned number = s->_number;
                Node* nd = &*next(_nodes.begin(), number);
                _species_partition[spiter.first].push_back(nd);
            }
        }

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

            // copy parent2
                if (othernd._parent2) {
                    unsigned parent2_number = othernd._parent2->_number;
                    Node* parent2 = &*next(_nodes.begin(), parent2_number);
                    nd->_parent2 = parent2;
                }

            // copy major parent
                if (othernd._major_parent) {
                    unsigned major_parent_number = othernd._major_parent->_number;
                    Node* major_parent = &*next(_nodes.begin(), major_parent_number);
                    nd->_major_parent = major_parent;
                }

                if (othernd._minor_parent) {
                    unsigned minor_parent_number = othernd._minor_parent->_number;
                    Node* minor_parent = &*next(_nodes.begin(), minor_parent_number);
                    nd->_minor_parent = minor_parent;
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
                nd->_visited = othernd._visited;
                nd->_hybrid_newick_name = othernd._hybrid_newick_name;
                nd->_deep_coalescence_counted = othernd._deep_coalescence_counted;
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

    inline string Forest::chooseEvent() {
        string event;
        // hybridization prior
        double rate = (_lambda+_hybridization_rate)*_lineages.size();
        
        double hybridization_prob = _hybridization_rate/(_hybridization_rate+_lambda);
        
        double u = rng.uniform();
        if (u<hybridization_prob && _lineages.size()>2) {
            event = "hybridization";
        }
        else if (_lineages.size() == 1) {
            event = "null";
        }
        else {
            event = "speciation";
        }
        // choose edge length but don't add it yet
        _last_edge_length = rng.gamma(1.0, 1.0/rate);
        
        return event;
    }

    inline void Forest::chooseSpeciesIncrementOnly(Lot::SharedPtr lot, double max_depth) {
        assert (max_depth >= 0.0);
        if (max_depth > 0.0) {
            // hybridization prior
            double rate = (_lambda)*_lineages.size();
            
            double u = lot->uniform();
            double inner_term = 1-exp(-rate*max_depth);
            _last_edge_length = -log(1-u*inner_term)/rate;
            assert (_last_edge_length < max_depth);

            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            // lorad only works if all topologies the same - then don't include the prior on joins b/c it is fixed
            double increment_prior = (log(rate)-_last_edge_length*rate);
                        
            _increments.push_back(make_pair(_last_edge_length, increment_prior));
            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior)); // do not include constrained factor in increment prior
        }
        else {
            double rate = _lambda*_lineages.size();
            
            assert (lot != nullptr);
            _last_edge_length = lot->gamma(1.0, 1.0/rate);

            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size() - 1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;

            _increments.push_back(make_pair(_last_edge_length, increment_prior));
            _increments_and_priors.push_back(make_pair(_last_edge_length, increment_prior));

        }
        
        if (_species_build.size() == 0) {
            _species_build.push_back(make_pair(make_tuple("null", "null", "null"), _last_edge_length));
        }
        else {
            _species_build.back().second = _last_edge_length;
        }
    }


    inline void Forest::chooseSpeciesIncrement(Lot::SharedPtr lot) {
        // hybridization prior
        double rate = (_lambda+_hybridization_rate)*_lineages.size();
        
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

        if (_lineages.size() > 1) {
            _species_joined = make_pair(subtree1, subtree2); // last step just joins remaining two
        }
        
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
        
        else if (get<0>(_hybrid_species_joined) != NULL) {
            cout << "hybridizing species " << get<0>(_hybrid_species_joined)->_name << " (hybrid) and " << get<1>(_hybrid_species_joined)->_name << " (parent) and " << get<2>(_hybrid_species_joined)->_name <<  " (parent2) " << endl;
        }
        
        else {
            cout << "no species joined" << endl;
        }
    }

    inline void Forest::setUpGeneForest(map<string, string> &taxon_map) {
        _taxon_map = taxon_map;
        assert (_index >0);
        _species_partition.clear();
//        for (auto nd:_lineages) {
        for (auto &nd:_nodes) {
            if (!nd._left_child) {
                string species_name = taxon_map[nd._name];
                _species_partition[species_name].push_back(&nd);
            }
        }
        assert (_species_partition.size() > 0);
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

    inline void Forest::allowCoalescence(string species_name, double increment, Lot::SharedPtr lot) {
        _log_likelihood_choices.clear();
        _node_choices.clear();
        _other_log_weight = 0.0;
        double prev_log_likelihood = _gene_tree_log_likelihood;
            
        Node *subtree1 = nullptr;
        Node *subtree2 = nullptr;
        list<Node*> nodes;
        
        for (auto &s:_species_partition) {
            if (s.first == species_name) {
                nodes = s.second;
                break;
            }
        }
        
        unsigned s = (unsigned) nodes.size();
        calcTopologyPrior(s);
        
        assert (s > 1);
        bool one_choice = false;

            // prior-prior proposal
            pair<unsigned, unsigned> t = chooseTaxaToJoin(s, lot);
            auto it1 = std::next(nodes.begin(), t.first);
            subtree1 = *it1;

            auto it2 = std::next(nodes.begin(), t.second);
            subtree2 = *it2;
            assert (t.first < nodes.size());
            assert (t.second < nodes.size());
        
            assert (subtree1 != subtree2);
            
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
            
            for (auto &s:_species_partition) {
                if (s.first == species_name) {
                    s.second = nodes; // TODO: this should happen automatically
                    break;
                }
            }
            
            if ((_proposal == "prior-prior" || one_choice) && (!_run_on_empty) ) {
                _gene_tree_log_likelihood = calcLogLikelihood();
                _log_weight = _gene_tree_log_likelihood - prev_log_likelihood;
            }
        if (_save_memory) {
            for (auto &nd:_nodes) {
                nd._partial=nullptr;
            }
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

    inline void Forest::hybridizeNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * delnode3, Node * addnode1) {
        // Delete delnode1 from node_vector
        auto it1 = find(node_vector.begin(), node_vector.end(), delnode1);
        assert(it1 != node_vector.end());
        node_vector.erase(it1);

        // Delete delnode2 from node_vector
        auto it2 = find(node_vector.begin(), node_vector.end(), delnode2);
        assert(it2 != node_vector.end());
        node_vector.erase(it2);
        
        // Delete delnode3 from node_vector
        auto it3 = find(node_vector.begin(), node_vector.end(), delnode3);
        assert(it3 != node_vector.end());
        node_vector.erase(it3);

        // Add addnode1 to node_vector
        node_vector.push_back(addnode1);

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

    inline void Forest::allowMigration(list<Node*> &nodes) {
        unsigned taxon_choice;
        Node* taxon_to_migrate;
        taxon_choice = chooseTaxonToMigrate(nodes.size());
        auto iter = std::next(nodes.begin(), taxon_choice);
        taxon_to_migrate = *iter;

        // find taxa to delete and lineage to add to
        string key_to_del = findKeyToDel(taxon_to_migrate);

        // delete migrating taxon from its original lineage and add to new lineage
        migrateTaxon(taxon_choice, key_to_del, taxon_to_migrate);
    }

    inline double Forest::chooseTaxonToMigrate(double s) {
        mtx.lock();
        unsigned taxon_choice = ::rng.randint(0, s-1);
        mtx.unlock();
        return taxon_choice;
    }

    inline string Forest::findKeyToDel(Node* taxon_to_migrate) {
        // find lineage to move taxon from in species partition
        string key_to_del;
        // TODO: is there a better way to do this?
        for (auto &s:_species_partition) {
            for (auto &t:s.second) {
                if (t == taxon_to_migrate) {
                    key_to_del = s.first;
                    break;
                    }
                }
            }
        return key_to_del;
    }

    inline void Forest::migrateTaxon(unsigned taxon_choice, string key_to_del, Node* taxon_to_migrate) {
        string key_to_add;
        // choose lineage to migrate into and add migrating taxon to chosen lineage
        bool choose_lineage = true;
        while (choose_lineage) {
            key_to_add = chooseLineage(taxon_to_migrate, key_to_del);
            // make sure chosen lineage is different from original lineage (otherwise, migration isn't really happening)
            if (key_to_add != key_to_del) {
                addMigratingTaxon(key_to_add, key_to_del, taxon_to_migrate);
                choose_lineage = false;
            }
        }

        deleteTaxon(key_to_del, taxon_choice);
        
        // calculate difference between lineage length of migrating taxon and target edge length (lineage taxon is migrating into)
        double difference = calculateNewEdgeLength(key_to_add, taxon_to_migrate);
        
        // set edge length of migrating taxon or target edge length (lineage taxon is migrating into)
        setNewEdgeLength(difference, taxon_to_migrate, key_to_add);
    }

    inline string Forest::chooseLineage (Node* taxon_to_migrate, string key_to_del) {
        // make vector of species names
        vector<string> species_names;
        
        for (auto & s:_species_partition) {
            species_names.push_back(s.first);
        }
        
        // choose lineage to migrate into
        string key_to_add;
        mtx.lock();
        unsigned lineage_choice = ::rng.randint(0, (unsigned) _species_partition.size()-1);
        mtx.unlock();

        // find lineage to migrate to in species partition
        key_to_add = species_names[lineage_choice];

        return key_to_add;
    }

    inline void Forest::addMigratingTaxon(string key_to_add, string key_to_del, Node* taxon_to_migrate) {
         // add migrating taxon to the chosen lineage
        _species_partition[key_to_add].push_back(taxon_to_migrate);
    }

    inline void Forest::deleteTaxon(string key_to_del, unsigned taxon_choice) {
        // delete migrating taxon from its original lineage
        bool done = false;
        for (auto &s:_species_partition) {
            if (done) {
                break;
            }
            if (s.first == key_to_del) {
                unsigned i = -1;
                for (auto itr = s.second.begin(); itr != s.second.end(); itr++) {
                    i++;
                    if (i == taxon_choice) {
                        s.second.erase(itr);
                        done = true;
                        break;
                    }
                }
            }
        }
    }

    inline double Forest::calculateNewEdgeLength(string key_to_add, Node* taxon_to_migrate) {
        // find target lineage
        Node* taxon_in_target_lineage = nullptr;
        double target_edge_length = 0.0;
        for (auto &s:_species_partition) {
            if (s.first == key_to_add) {
                taxon_in_target_lineage = s.second.front();
                break;
            }
        }
        
        // get target lineage edge length
        for (Node* child = taxon_in_target_lineage; child; child=child->_left_child) {
            target_edge_length += child->_edge_length;
        }
        
        // get edge length of migrating lineage
        double migrating_edge_length = 0.0;
        for (Node* child = taxon_to_migrate; child; child = child->_left_child) {
            migrating_edge_length += child->_edge_length;
        }
            
        // compare target edge length to migrating edge length
        double difference = target_edge_length - migrating_edge_length;
        return difference;
    }

    inline void Forest::setNewEdgeLength(double difference, Node* taxon_to_migrate, string key_to_add) {
        // if migrating lineage is shorter than target lineage, extend migrating taxon edge length
        if (difference > 0.0) {
            taxon_to_migrate->_edge_length += difference;
        }
        
        // if migrating lineage is longer than target lineage, extend target edge length
        // TODO: does this work if migrating taxon has edge length 0 but its children do not?
        // TODO: can this be simplified?
        else if (difference < 0.0) {
            for (auto &s:_species_partition) {
                // find target species and extend each lineage
                if (s.first == key_to_add) {
                    for (auto iter = s.second.begin(); iter != s.second.end(); iter++) {
                        Node* taxon = *iter;
                        if (taxon == taxon_to_migrate) {
                            break;
                        }
                        taxon->_edge_length += -1.0*difference;
                    }
                    break;
                }
            }
        }
    }

    inline void Forest::hybridizeGene(vector<string> hybridized_nodes, double species_tree_increment, string species_name, Lot::SharedPtr lot) {
        // parent, parent2, hybrid_node, new_nd
        string parent = hybridized_nodes[0];
        string parent2 = hybridized_nodes[1];
        string hybrid = hybridized_nodes[2];
        string new_nd = hybridized_nodes[3];
        string new_nd2 = hybridized_nodes[4];

        // find hybridizing lineage
        // move the gene in the hybrid node left or right
        
        // prior-prior
# if false
        if (u < gamma) {
            double gamma = 0.85;
            double u = rng.uniform();
            // move gene in direction of major parent
            _last_direction = "major";
            moveGene(new_nd, parent, hybrid);
        }
        else {
            // move gene in direction of minor parent
            _last_direction = "minor";
            moveGene(new_nd2, parent2, hybrid);
        }

        for (auto & s:_species_partition) {
            assert (s.second.size()>0);
//            evolveSpeciesFor(s.second, species_tree_increment, s.first);
        }
# endif
        // prior-post
        // move gene both ways, save likelihood of each, reweight, then draw random number to decide which way to move
        
        vector<list<Node*>> original_nodes;
        vector<string> original_names;
        for (auto &s:_species_partition) {
            original_names.push_back(s.first);
            original_nodes.push_back(_species_partition[s.first]);
        }
        
        // save branch lengths of original _lineages vector
        vector<double> branch_lengths = saveBranchLengths();
        
        _new_nodes.clear();
        // move towards minor parent
        moveGene(new_nd, parent2, hybrid);
        
        // go through coalescence
        hybridGeneTreeProposal(species_tree_increment, species_name);
        
        // save likelihood
        vector<double> likelihood_vec;
        likelihood_vec.reserve(2);
        likelihood_vec.push_back(calcLogLikelihood());
        
        // save minor move information
        vector<Node*> minor_nodes = _new_nodes;
        vector<double> minor_branch_lengths = saveBranchLengths();
        
        // save children of minor nodes
        vector<Node*> minor_left_children;
        vector<Node*> minor_right_children;
        vector<double> minor_left_edge_lengths;
        vector<double> minor_right_edge_lengths;
        
        for (auto &nd:minor_nodes) {
            minor_left_children.push_back(nd->_left_child);
            minor_right_children.push_back(nd->_left_child->_right_sib);
            minor_left_edge_lengths.push_back(nd->_left_child->_edge_length);
            minor_right_edge_lengths.push_back(nd->_left_child->_right_sib->_edge_length);
        }
        
        // save minor spp partition
        vector<list<Node*>> minor_partition;
        vector<string> minor_names;
        for (auto &s:_species_partition) {
            minor_names.push_back(s.first);
            minor_partition.push_back(_species_partition[s.first]);
        }
        
        resetLineages(branch_lengths);
        
        // rebuild species partition
        _species_partition.clear();
        rebuildSpeciesPartition(original_names, original_nodes);
        
        _new_nodes.clear();
        // move towards major parent
        moveGene(new_nd, parent, hybrid);
        
        // go through coalescence
        hybridGeneTreeProposal(species_tree_increment, species_name);
        
        // save likelihood
        likelihood_vec.push_back(calcLogLikelihood());
        
        // choose major or minor path
        int index_of_choice = chooseDirectionOfHybridization(likelihood_vec, lot);
        if (index_of_choice == 0) {
            // minor choice
            _last_direction = "minor";
            
            // rebuild species partition to original
            _species_partition.clear();
            rebuildSpeciesPartition(original_names, original_nodes);
            
            // revert _lineages to original
            resetLineages(branch_lengths);
            
            // clear major nodes
            for (auto &nd:_new_nodes) {
                nd->_left_child = 0;
                nd->_name = "unused";
                nd->_edge_length = 0;
                nd->_partial->clear();
                nd->_position_in_lineages = -1;
            }
//            _new_nodes.clear();
            
            resetToMinor(minor_nodes, minor_left_children, minor_right_children, minor_left_edge_lengths, minor_right_edge_lengths);

            // revert all _lineages edge lengths to minor nodes
            for (int i=0; i < (int) _lineages.size(); i++) {
                _lineages[i]->_edge_length = minor_branch_lengths[i];
            }
            // rebuild species partition
            _species_partition.clear();
            rebuildSpeciesPartition(minor_names, minor_partition);

            // switch parent and parent2
            switchParents(parent, parent2);
        }
        else {
            // major choice
            _last_direction = "major";
            for (auto &nd:minor_nodes) {
                nd->_left_child = 0;
                nd->_name = "unused";
                nd->_edge_length = 0;
                nd->_partial->clear();
                nd->_position_in_lineages = -1;
            }
        }
            // remove unused nodes
        for (auto iter = _nodes.begin(); iter != _nodes.end(); iter++) {
            if (iter->_name == "unused") {
                iter = _nodes.erase(iter);
                --iter;
                _ninternals--;
            }
        }
        // reset node numbers
        int n = 0;
        for (auto &nd:_nodes) {
            nd._number = n;
            n++;
        }
        assert(_nodes.size() > 0);
        }

    inline void Forest::hybridGeneTreeProposal(double species_tree_increment, string species_name) {
        if (_species_partition.size() == 1) {
//            fullyCoalesceGeneTree(_species_partition.begin()->second);
        }

        else {
            for (auto &s:_species_partition) {
                assert (s.second.size()>0);
//                evolveSpeciesFor(s.second, species_tree_increment, s.first);
            }
        }
    }

    inline void Forest::resetToMinor(vector<Node*> minor_nodes, vector<Node*>minor_left_children, vector<Node*>minor_right_children, vector<double> minor_left_edge_lengths, vector<double> minor_right_edge_lengths) {
        // find new nodes
        int k = -1;
        // reset _lineages to minor coalescence
        for (auto &minor_node:minor_nodes) {
            k++;
            for (auto &nd:_lineages) {
                if (minor_node->_left_child == nd) {
                    updateNodeVector(_lineages, minor_left_children[k], minor_right_children[k], minor_nodes[k]);
                    minor_nodes[k]->_parent = 0;
                    minor_nodes[k]->_left_child = minor_left_children[k];
                    minor_nodes[k]->_left_child->_right_sib = minor_right_children[k];
                    minor_left_children[k]->_parent = minor_nodes[k];
                    minor_right_children[k]->_parent = minor_nodes[k];
                    minor_nodes[k]->_right_sib = 0;
                    minor_nodes[k]->_left_child->_edge_length = minor_left_edge_lengths[k];
                    minor_nodes[k]->_left_child->_right_sib->_edge_length = minor_right_edge_lengths[k];
                }
            }
        }
    }

    inline void Forest::moveGene(string new_nd, string parent, string hybrid) {
        // update species partition
        list<Node*> &nodes = _species_partition[new_nd];
        copy(_species_partition[parent].begin(), _species_partition[parent].end(), back_inserter(nodes));
        copy(_species_partition[hybrid].begin(), _species_partition[hybrid].end(), back_inserter(nodes));
        _species_partition.erase(parent);
        _species_partition.erase(hybrid);
        assert(nodes.size()>0);
    }

    inline void Forest::rebuildSpeciesPartition(vector<string> names, vector<list<Node*>> nodes) {
        assert(_species_partition.size() == 0);
        int i = 0;
        for (auto &name:names) {
            _species_partition[name] = nodes[i];
            i++;
        }
    }

    inline int Forest::chooseDirectionOfHybridization(vector<double> likelihood_vec, Lot::SharedPtr lot) {
        // choose a direction
        vector<double> log_weight_choices;
        log_weight_choices.reserve(2);
        
        log_weight_choices.push_back(likelihood_vec[0]+log(.15)); // multiply minor likelihood by (1-gamma)
        log_weight_choices.push_back(likelihood_vec[1]+log(.85)); // multiply major likelihood by (gamma)
        
        // normalize weights
//        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        double log_weight_choices_sum = getRunningSumHybridChoices(log_weight_choices);
        for (int b=0; b < (int) log_weight_choices.size(); b++) {
            log_weight_choices[b] -= log_weight_choices_sum;
        }
        
        // select a direction
        int index_of_choice = selectPair(log_weight_choices, lot);
        return index_of_choice;
    }

    inline vector<double> Forest::saveBranchLengths() {
        vector<double> branch_lengths;
        for (auto &nd:_lineages) {
            branch_lengths.push_back(nd->_edge_length);
        }
        return branch_lengths;
    }

    inline void Forest::switchParents(string parent, string parent2) {
        list<Node*> &nodes = _species_partition[parent2];
        copy(_species_partition[parent].begin(), _species_partition[parent].end(), back_inserter(nodes));
        _species_partition.erase(parent);
        assert(nodes.size()>0);
    }

    inline void Forest::resetLineages(vector<double> branch_lengths) {
        for (int a = (int) _new_nodes.size()-1; a>=0; a--) {
            revertNodeVector(_lineages, _new_nodes[a]->_left_child, _new_nodes[a]->_left_child->_right_sib, _new_nodes[a]);
        }
        
        // reset _lineages edge lengths
        assert (_lineages.size() == branch_lengths.size());
        for (int i=0; i < (int) _lineages.size(); i++) {
            _lineages[i]->_edge_length = branch_lengths[i];
            _lineages[i]->_parent = 0;
            _lineages[i]->_right_sib = 0;
        }
    }

    inline vector<string> Forest::hybridizeSpecies() {
        tuple<unsigned, unsigned, unsigned> t = chooseTaxaToHybridize();
        Node* parent = _lineages[get<0>(t)];
        Node* parent2 = _lineages[get<1>(t)];
        Node* hybrid_node = _lineages[get<2>(t)];
        
        _hybrid_species_joined = make_tuple(hybrid_node, parent, parent2);

        assert (!parent->_parent && !hybrid_node->_parent && !parent2->_parent);

//        create a new node
        Node nd;
        _nodes.push_back(nd);
        Node* new_nd = &_nodes.back();
//        Node* new_nd = &_nodes[_nleaves+_ninternals];
        new_nd->_parent=0;
        new_nd->_number=_nleaves+_ninternals;
        new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);
        new_nd->_edge_length=0.0;
        _ninternals++;
        new_nd->_right_sib=0;
        new_nd->_left_child=parent;
        parent->_right_sib=hybrid_node;
        parent->_parent=new_nd;
        hybrid_node->_parent=new_nd;
        
//        create another new node
        Node nd2;
        _nodes.push_back(nd2);
        Node* new_nd2 = &_nodes.back();
//        Node* new_nd2 = &_nodes[_nleaves+_ninternals];
        new_nd2->_parent=0;
        new_nd2->_number=_nleaves+_ninternals;
        new_nd2->_name=boost::str(boost::format("node-%d")%new_nd2->_number);
        new_nd2->_edge_length=0.0;
        _ninternals++;
        new_nd2->_left_child=parent2;
        hybrid_node->_right_sib = parent2;
        hybrid_node->_parent2=new_nd2;
        new_nd2->_left_child->_right_sib = hybrid_node;
        
        updateNodeVector(_lineages, parent, hybrid_node, new_nd);
        
        vector<string> hybridized_nodes;
        
        hybridized_nodes.push_back(parent->_name);
        hybridized_nodes.push_back(parent2->_name);
        hybridized_nodes.push_back(hybrid_node->_name);
        hybridized_nodes.push_back(new_nd->_name);
        hybridized_nodes.push_back(new_nd2->_name);
        
        hybrid_node->_major_parent = parent;
        hybrid_node->_minor_parent = parent2;
        
        // update _lineages vector with major new_nd
        return hybridized_nodes;
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
        
        for (int i=0; i<_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            species_names.push_back(name);
        }
        
        assert (species_names.size() == 2*_nspecies - 1);
        
        // draw thetas for tips of species trees and ancestral population
        // for all other populations, theta = -1
        
        assert (_theta_proposal_mean > 0.0);
        double scale = 1 / _theta_proposal_mean;
        _scale = scale;
        
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
        scale = _scale;
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

    inline void Forest::createThetaMap(Lot::SharedPtr lot) {
        // map should be 2*nspecies - 1 size
        unsigned number = 0;
        vector<string> species_names;
        
        for (auto &s:_species_partition) {
            species_names.push_back(s.first);
            number++;
            _species_indices[s.first] = number - 1;
        }
        for (int i=0; i<_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            species_names.push_back(name);
            _species_indices[name] = number - 1;
            if (i == _nspecies-2) {
                _ancestral_species_name = name;
            }
        }
        
        assert (species_names.size() == 2*_nspecies - 1);
        
        // gamma mean = shape * scale
        // draw mean from lognormal distribution
        // shape = 2.0 to be consistent with starbeast3
        // scale = 1 / mean;
        
        if (_theta_proposal_mean > 0.0) {
            assert (_theta_mean == 0.0);
//            double exponential_rate = 1 / _theta_proposal_mean;
//            _theta_mean = -log(1.0 - rng.uniform()) / exponential_rate;
//            double mean = 1 / exponential_rate;
//            _theta_mean = lot->gamma(1, mean); // mean = 10, equivalent to exponential(exponential_rate)
            _theta_mean = lot->gamma(1, _theta_proposal_mean);
        }
        
//        if (_theta_mean == 0.0) { // TODO: cannot specify theta_mean and theta_proposal_mean / theta_prior_mean
            // Draw _theta_mean from Exponential prior
//            _theta_mean = -log(1.0 - rng.uniform())/exponential_prior_rate;
//            _theta_proposal_mean = lot->logNormal(-4.6, 2.14); // TODO: mean = 0.1, sd = 1 --> make sd 0.1? - or try gamma
//        }
//        _theta_mean = lot->gamma(1, 10); // mean = 10, equivalent to exponential(0.1)
//        _theta_mean = lot->gamma(2, _gamma_scale); // mean = 2 * scale, var = 2 * scale^2 // TODO: variance is very small when scale small

        
        double scale = 1 / _theta_mean;
        _scale = scale;
        assert (scale > 0.0);
        for (auto &name:species_names) {
            double new_theta = 0.0;
            if (new_theta < _small_enough) {
//                new_theta = 1 / (lot->gamma(2.00001, scale));
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

    inline pair<vector<double>, vector<unsigned>> Forest::calcCoalescentLikelihoodIntegratingOutTheta(vector<pair<tuple<string,string,string>, double>> species_build) {
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
        
//        showForest();
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
        for (int i=0; i<_depths.size(); i++) { // TODO: can't iterate through a vecotr and change it at the same time - keep track of which ones to erase and then erase them
            if (_depths[i].second.first == _depths[i].second.second) {
//                // species have already been joined in the species tree, so they are no longer a constraint
//                _depths.erase(_depths.begin()+i);
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
    //                if (nd->_left_child && !nd->_visited) {
            if (nd->_left_child) {
                while (!done) {
                    Node* left_child = nd->_left_child;
                    Node* right_child = nd->_left_child->_right_sib;
                    while (spp_left_child == "" || spp_right_child == "") {
                        for (auto &s:_species_partition) {
                            for (auto &nd:s.second) {
                                if (nd == left_child) {
                                    spp_left_child = s.first;
                                }
                                else if (nd == right_child) {
                                    spp_right_child = s.first;
                                }
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
    //                    nd->_visited = true;
                    if (spp_left_child != "" && spp_right_child != "") {
                        left_child->_parent->_visited = true;
                        done = true;
                    }
                    if (spp_left_child != spp_right_child) {
                        double height = getTotalLineageHeight(nd->_left_child);
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
                double height = getTotalLineageHeight(nd->_left_child); //
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

    inline double Forest::getTotalLineageHeight(Node* nd) {
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
        
        // Simulate starting sequence at the root node

        Node * nd = *(_lineages.begin());
        unsigned ndnum = nd->_number;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = lot->randint(0,3);
        }
        
        nd = findNextPreorder(nd);
        while (nd) {
            ndnum = nd->_number;
            assert(ndnum < nnodes);

            // Get reference to parent sequence
            assert(nd->_parent);
            unsigned parnum = nd->_parent->_number;
            assert(parnum < nnodes);

            // Evolve nd's sequence given parent's sequence and edge length
            for (unsigned i = 0; i < nsites; i++) {
                unsigned from_state = sequences[parnum][i];
                double cum_prob = 0.0;
                double u = lot->uniform();
                for (unsigned to_state = 0; to_state < 4; to_state++) {
                    cum_prob += calcSimTransitionProbability(from_state, to_state, nd->_edge_length);
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
//                unsigned ndnum = _nodes[t]._number;
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
                 return getLineageHeight(lhs._left_child) < getLineageHeight(rhs._left_child); } ); // TODO: this isn't working?
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

