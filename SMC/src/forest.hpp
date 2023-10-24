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

#include "lot.hpp"
extern proj::Lot rng;
std::mutex mtx;

#include "partial_store.hpp"
extern proj::PartialStore ps;
extern int my_rank;

#include "node.hpp"

namespace proj {

using namespace std;

class Likelihood;
class Particle;

class Forest {

        friend class Likelihood;
        friend class Particle;

    public:
                                    Forest(string type);
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
        void                        debugLogLikelihood(Node* nd, double log_like);
        double                      calcTopologyPrior(unsigned nlineages);

    private:

        typedef std::vector <double> partial_array_t;
//        void                        clear();
        void                        clearGeneForest();
        void                        clearSpeciesForest();
        void                        clearForestForCopying();
        void                        setData(Data::SharedPtr d, int index, map<string, string> &taxon_map);
        Node *                      findNextPreorder(Node * nd);
        std::string                 makeNewick(unsigned precision, bool use_names);
        std::string                 makePartialNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>    chooseTaxaToJoin(double s);
        tuple<Node*, Node*, Node*>  createNewSubtree(pair<unsigned, unsigned> p, list<Node*> node_list, double increment, string species);
        tuple<Node*, Node*, Node*>  createNewSubtreeFromPrior(pair<unsigned, unsigned> p, double increment);
        void                        calcPartialArray(Node* new_nd);
        void                        setUpGeneForest(map<string, string> &taxon_map);
        void                        setUpSpeciesForest(vector<string> &species_names);
        tuple<string,string, string> speciesTreeProposal();
        void                        geneTreeProposal(pair<double, string> species_info, vector<pair<tuple<string, string, string>, double>> _t);
        void                        evolveSpeciesFor(list <Node*> &nodes, double increment, string species);
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                        hybridizeNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node* delnode3, Node * addnode1);
        void                        revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        void                        revertNodeList(list<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                      getRunningSumChoices(vector<double> &log_weight_choices);
        double                      getRunningSumHybridChoices(vector<double> &log_weight_choices);
        vector<double>              reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        pair<Node*, Node*>          getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list);
        int                         selectPair(vector<double> weight_vec);
        void                        chooseSpeciesIncrement(double max_depth);
        void                        addSpeciesIncrement();
        string                      chooseEvent();
        void                        allowMigration(list<Node*> &nodes);
        double                      chooseTaxonToMigrate(double s);
        string                      findKeyToDel(Node* taxon_to_migrate);
        void                        migrateTaxon(unsigned taxon_choice, string key_to_del, Node* taxon_to_migrate);
        string                      chooseLineage(Node* taxon_to_migrate, string key_to_del);
        void                        addMigratingTaxon(string key_to_add, string key_to_del, Node* taxon_to_migrate);
        void                        deleteTaxon(string key_to_del, unsigned taxon_choice);
        void                        allowCoalescence(list<Node*> &nodes, double increment, string species);
        tuple<unsigned, unsigned, unsigned> chooseTaxaToHybridize();
        vector<string>              hybridizeSpecies();
        void                        moveGene(string new_nd, string parent, string hybrid);
        void                        rebuildSpeciesPartition(vector<string> names, vector<list<Node*>> nodes);
        void                        switchParents(string parent, string parent2);
        void                        resetLineages(vector<double> branch_lengths);
        vector<double>              saveBranchLengths();
        int                         chooseDirectionOfHybridization(vector<double> likelihood_vec);
        void                        hybridGeneTreeProposal(double species_tree_increment);
        unsigned                    countDescendants(Node* nd, unsigned count);
        double                      getTreeHeight();
        double                      getLineageHeight(Node* nd); // getLineageHeight(const Node* nd) const;
        void                        extendGeneTreeLineages(double species_tree_increment, bool deep_coalescence);
        double                      extendSingleLineage(double gene_increment, Node* nd, string spp_right_child, string spp_left_child);
        void                        updateSpeciesPartition(tuple<string, string, string> species_info);
        pair<double, string>        chooseDelta(vector<pair<tuple<string, string, string>, double>> species_info);
        pair<Node*, Node*>          chooseAllPairs(list<Node*> &nodes, double increment, string species);
        pair<Node*, Node*>          chooseAllPairsFromPrior(double increment);
        double                      calcCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height);
        vector< pair<double, Node *>>      sortPreorder();
        void                        calcMinDepth();
        vector<pair<double, pair<string, string>>>             getMinDepths();
        void                        resetDepthVector(tuple<string, string, string> species_joined);
        void                        buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies);
    vector<pair<tuple<string, string, string>, double>>                         buildFromNewickTopology(const std::string newick, bool topology_only);
        void                        stripOutNexusComments(std::string & newick);
        unsigned                    countNewickLeaves(const std::string newick);
        void                        extractEdgeLen(Node * nd, std::string edge_length_string);
        bool                        canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);
        void                        renumberInternals();
        void                        remakeGeneTree(map<string, string> &taxon_map);
        void                        resetIncrements();
        vector<string>              updateExistingLineagesVector(vector<string> existing_lineages, tuple<string, string, string> species_joined);
        vector<string>              setUpExistingLineagesVector();
        void                        chooseSpeciesIncrementFromNewick(vector<string> existing_lineages);
        void                        drawFromGeneTreePrior();
        double                      calcLogSpeciesTreeDensity(double lambda);
        void                        setIndex(int n) {_index = n;}
        void                        resetLineages();
//        void                        stowForestPartials();

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
        double                      _gene_tree_log_weight;
        vector<double>              _log_weight_vec;
        vector<pair<Node*, Node*>>  _node_choices;
        vector<double>              _log_likelihood_choices;
        int                         _index_of_choice;
        tuple<Node*, Node*, Node*>  _hybrid_species_joined;
        string                      _last_direction;
        vector<pair<double, double>>    _increments;
        double                      _topology_prior;
        double                      _prev_log_likelihood;
        double                      _log_joining_prob;
        vector<double>              _increment_choices;
        double                      _extended_increment;
        unsigned                    _species_join_number;
        double                      _prev_gene_tree_log_likelihood;
        bool                        _ready_to_join_species;
        vector<Node*>               _preorder;
        vector<pair<double, pair<string, string>>>              _depths;
        double                      _gene_tree_log_coalescent_likelihood;
        double                      _panmictic_coalescent_likelihood;
        int                         _nincrements = 0;
        double                      _count = 0;

        double                      calcTransitionProbability(Node* child, double s, double s_child);
        double                      calculateNewEdgeLength(string key_to_add, Node* taxon_to_migrate);
        void                        setNewEdgeLength(double difference, Node* taxon_to_migrate, string key_to_add);
        void                        hybridizeGene(vector<string> hybridized_nodes, double species_tree_increment);
        void                        resetToMinor(vector<Node*> minor_nodes, vector<Node*> minor_left_children, vector<Node*> minor_right_children, vector<double> minor_left_edge_lengths, vector<double> minor_right_edge_lengths);
        vector<pair<double, double>>              _deep_coalescent_increments; // increment, prior
        void                        refreshPreorder();
        double                      calcLogCoalLikeGivenTheta(double proposed_theta, vector<pair<tuple<string, string, string>, double>> species_info, bool both);

    public:

        typedef std::shared_ptr<Forest> SharedPtr;
        static double               _theta;
        static double               _lambda;
        static string               _proposal;
        static string               _model;
        static double               _kappa;
        static vector<double>       _base_frequencies;
        static string               _string_base_frequencies;
        static double               _migration_rate;
        static double               _hybridization_rate;
        static string               _outgroup;
        static double               _theta_prior_mean;
        static double               _lambda_prior_mean;
};


    inline Forest::Forest(string type) {
        //std::cout << "Constructing a forest" << std::endl;
        if (type == "species") {
            clearSpeciesForest();
        }
        else if (type == "new") {
//            // don't need to create nodes because they will be overwritten when particle is copying
            clearForestForCopying();
        }
        else {
            clearGeneForest();
        }
    }

    inline Forest::~Forest() {
        //std::cout << "Destroying a Forest" << std::endl;
    }

    inline void Forest::clearForestForCopying() {
        _nodes.clear();
        _lineages.clear();
        _npatterns = 0;
        _nstates = 4;
        _last_edge_length = 0.0;
        _lineages.clear();
        _log_joining_prob = 0.0;
        _extended_increment = 0.0;
        _species_join_number = 0;
        _ready_to_join_species = false;
        _gene_tree_log_likelihood = 0.0;
        _gene_tree_log_weight = 0.0;
        _gene_tree_log_coalescent_likelihood = 0.0;
        _panmictic_coalescent_likelihood = 0.0;
        _log_weight_vec.clear();
        _node_choices.clear();
        _prev_log_likelihood = 0.0;
        _increment_choices.clear();
        _extended_increment = 0.0;
        _species_join_number = 0;
        _deep_coalescent_increments.clear();
        _prev_gene_tree_log_likelihood = 0.0;
        _preorder.clear();
        _depths.clear();
        _nleaves=0;
        _ninternals=0;
        _topology_prior = 0.0;
        _index_of_choice = 0;
    }

    inline void Forest::clearGeneForest() {
        _nodes.clear();
        _lineages.clear();
        _nodes.resize(_ntaxa);
        _npatterns = 0;
        _nstates = 4;
        _last_edge_length = 0.0;
        _lineages.reserve(_nodes.size());
        _log_joining_prob = 0.0;
        _extended_increment = 0.0;
        _species_join_number = 0;
        _ready_to_join_species = false;
        _preorder.clear();
        _gene_tree_log_likelihood = 0.0;
        _gene_tree_log_weight = 0.0;
        _gene_tree_log_coalescent_likelihood = 0.0;
        _panmictic_coalescent_likelihood = 0.0;
        _log_weight_vec.clear();
        _node_choices.clear();
        _prev_log_likelihood = 0.0;
        _increment_choices.clear();
        _extended_increment = 0.0;
        _species_join_number = 0;
        _deep_coalescent_increments.clear();
        _prev_gene_tree_log_likelihood = 0.0;
        _depths.clear();
        _topology_prior = 0.0;
        _index_of_choice = 0;

        //create taxa
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd = &*next(_nodes.begin(), i);
//            Node new_nd;
//            _nodes.push_back(new_nd);
//            Node* nd = &_nodes.back();
            
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_position_in_lineages=i;
            _lineages.push_back(nd);
            nd->_partial.reset();
            }
        _nleaves=_ntaxa;
        _ninternals=0;

    }

    inline void Forest::clearSpeciesForest() {
        _nodes.clear();
        _lineages.clear();
        _nodes.resize(_nspecies);
        _npatterns = 0;
        _nstates = 4;
        _last_edge_length = 0.0;
        _lineages.reserve(_nodes.size());
        _log_joining_prob = 0.0;
        _extended_increment = 0.0;
        _species_join_number = 0;
        _ready_to_join_species = false;
        _preorder.clear();
        _gene_tree_log_likelihood = 0.0;
        _gene_tree_log_weight = 0.0;
        _gene_tree_log_coalescent_likelihood = 0.0;
        _panmictic_coalescent_likelihood = 0.0;
        _log_weight_vec.clear();
        _node_choices.clear();
        _prev_log_likelihood = 0.0;
        _increment_choices.clear();
        _extended_increment = 0.0;
        _species_join_number = 0;
        _deep_coalescent_increments.clear();
        _prev_gene_tree_log_likelihood = 0.0;
        _depths.clear();
        _topology_prior = 0.0;
        _index_of_choice = 0;

        //create taxa
        for (unsigned i = 0; i < _nspecies; i++) {
//            Node new_nd;
//            _nodes.push_back(new_nd);
//            Node* nd = &_nodes.back();
            
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
        _nleaves=_nspecies;
        _ninternals=0;

    }

    inline Forest::Forest(const Forest & other) {
        if (_index == 0) {
            clearSpeciesForest();
        }
        else {
            clearGeneForest();
        }
        *this = other;
    }

    inline void Forest::setData(Data::SharedPtr d, int index, map<string, string> &taxon_map) {
        assert (index > 0); // don't set data for species tree
        _data = d;
        _index = index;

        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(index-1);
        _first_pattern = gene_begin_end.first;
        _npatterns = _data->getNumPatternsInSubset(index-1);
        ps.setNElements(_npatterns, index-1);

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
                    nd->_partial=ps.getPartial(_npatterns*4, index-1);
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
        setUpGeneForest(taxon_map);
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

    inline vector< pair<double, Node *>> Forest::sortPreorder() {
            vector< pair<double, Node *>> heights_and_nodes;
            for (auto it = _preorder.rbegin(); it != _preorder.rend(); it++) {
                Node * nd = *it;
                if (nd->_left_child) {
                    // if internal node, store cumulative height in _height
                    double height = getLineageHeight(nd->_left_child); //
                    heights_and_nodes.push_back(make_pair(height, nd));
                }
            }
             
            // sort heights_and_nodes so that smallest heights will be first
            sort(heights_and_nodes.begin(), heights_and_nodes.end());
            return(heights_and_nodes);
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
        for (auto lineage : _lineages) {
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
                for (auto lineage : _lineages) {
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

    inline pair<unsigned, unsigned> Forest::chooseTaxaToJoin(double s){
        // save random numbers
        assert (s>1);
        double nsubtrees = s;
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        // thread safe random number generator with mutex
        mtx.lock();
        if (nsubtrees > 2) {
            t1 = ::rng.randint(0, nsubtrees-1);
            t2 = ::rng.randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = ::rng.randint(0, nsubtrees-1);
            }
        }
        mtx.unlock();
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
        mtx.unlock();
        return make_tuple(t1, t2, t3);
    }

    inline void Forest::calcPartialArray(Node* new_nd) {
        for (auto &nd:_lineages) {
            assert (nd->_right_sib != nd);
        }
        auto & parent_partial_array = *(new_nd->_partial);
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
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
//        auto data_matrix=_data->getDataMatrix();

        //calc likelihood for each lineage separately
        auto &counts = _data->getPatternCounts();
        _gene_tree_log_likelihood = 0.0;

        for (auto nd:_lineages) {
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

    inline pair<Node*, Node*> Forest::chooseAllPairsFromPrior(double increment) {
        _node_choices.clear();
        _log_likelihood_choices.clear();
        _gene_tree_log_weight = 0.0;
       
        // choose pair of nodes to try
        for (unsigned i = 0; i < _lineages.size()-1; i++) {
            for (unsigned j = i+1; j < _lineages.size(); j++) {
                // createNewSubtree returns subtree1, subtree2, new_nd
                
                tuple<Node*, Node*, Node*> t = createNewSubtreeFromPrior(make_pair(i,j), increment);
                
                _log_likelihood_choices.push_back(calcLogLikelihood()); // no coalescent likelihood
                // gene tree log coalescent likelihood is the same for every possible join

                // revert _lineages
                revertNodeVector(_lineages, get<0>(t), get<1>(t), get<2>(t));

                //reset siblings and parents of original nodes back to 0
                get<0>(t)->resetNode(); //subtree1
                get<1>(t)->resetNode(); //subtree2

                // clear new node from _nodes
                //clear new node that was just created
                get<2>(t)->clear(); //new_nd
            }
        }
        
        // reweight each choice of pairs
       double prev_log_likelihood = _prev_gene_tree_log_likelihood; // no coalescent likelihood
       vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, prev_log_likelihood);


        // sum unnormalized weights before choosing the pair
        // must include the likelihoods of all pairs in the final particle weight
        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        _gene_tree_log_weight = log_weight_choices_sum;
        for (unsigned b=0; b < log_weight_choices.size(); b++) {
            log_weight_choices[b] -= log_weight_choices_sum;
        }
        
        // randomly select a pair
        _index_of_choice = selectPair(log_weight_choices);

        // find nodes to join in node_list
        Node *subtree1 = _node_choices[_index_of_choice].first;
        Node *subtree2 = _node_choices[_index_of_choice].second;
    
        _gene_tree_log_likelihood = _log_likelihood_choices[_index_of_choice]; // no coalescent likelihood
       
        // erase extra nodes created from node list
        for (unsigned i = 0; i < _node_choices.size(); i++) {
            _nodes.pop_back();
        }
       
       // no coalescent likelihood
        return make_pair(subtree1, subtree2);
    }

    inline pair<Node*, Node*> Forest::chooseAllPairs(list<Node*> &node_list, double increment, string species) {
         _node_choices.clear();
         _log_likelihood_choices.clear();
         _gene_tree_log_weight = 0.0;
         
# if defined(GENE_TREE_COALESCENT_LIKELIHOOD)
        double starting_gene_tree_log_coalescent_likelihood = _gene_tree_log_coalescent_likelihood;
# else
        double starting_gene_tree_log_coalescent_likelihood = 0.0;
# endif
        
         // choose pair of nodes to try
         for (unsigned i = 0; i < node_list.size()-1; i++) {
             for (unsigned j = i+1; j < node_list.size(); j++) {
                 // createNewSubtree returns subtree1, subtree2, new_nd
                 
                 tuple<Node*, Node*, Node*> t = createNewSubtree(make_pair(i,j), node_list, increment, species);
                 
#if !defined(GENE_TREE_COALESCENT_LIKELIHOOD)
                 assert (_gene_tree_log_coalescent_likelihood == 0.0);
                 assert (starting_gene_tree_log_coalescent_likelihood == 0.0);
#endif
                 _log_likelihood_choices.push_back(calcLogLikelihood()+_gene_tree_log_coalescent_likelihood);
                 // gene tree log coalescent likelihood is the same for every possible join

                 // revert _lineages if > 1 choice
                 revertNodeVector(_lineages, get<0>(t), get<1>(t), get<2>(t));

                 //reset siblings and parents of original nodes back to 0
                 get<0>(t)->resetNode(); //subtree1
                 get<1>(t)->resetNode(); //subtree2

                 // clear new node from _nodes
                 //clear new node that was just created
                 
                 // stow partial before removing it
//                 ps.stowPartial(get<2>(t)->_partial, _index-1);
                 get<2>(t)->clear(); //new_nd
             }
         }
         
         // reweight each choice of pairs
        double prev_log_likelihood = _prev_gene_tree_log_likelihood + starting_gene_tree_log_coalescent_likelihood;
        vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, prev_log_likelihood);

         // sum unnormalized weights before choosing the pair
         // must include the likelihoods of all pairs in the final particle weight
         double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
         _gene_tree_log_weight = log_weight_choices_sum;
         for (unsigned b=0; b < log_weight_choices.size(); b++) {
             log_weight_choices[b] -= log_weight_choices_sum;
         }
        
         // randomly select a pair
         _index_of_choice = selectPair(log_weight_choices);

         // find nodes to join in node_list
         Node* subtree1 = _node_choices[_index_of_choice].first;
         Node* subtree2 = _node_choices[_index_of_choice].second;
     
         _gene_tree_log_likelihood = _log_likelihood_choices[_index_of_choice] - _gene_tree_log_coalescent_likelihood; // remove the coalescent likelihood to get the ordinary gene tree log likelihood
        
         // erase extra nodes created from node list
         for (unsigned i = 0; i < _node_choices.size(); i++) {
             _nodes.pop_back();
         }
        
        // reset gene tree log coalescent likelihood
        _gene_tree_log_coalescent_likelihood = starting_gene_tree_log_coalescent_likelihood;
        
         return make_pair(subtree1, subtree2);
     }

    inline vector<pair<double, pair<string, string>>> Forest::getMinDepths() {
        return _depths;
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
            advance(l_front, curr_node_index);
            Node *nd = &*l_front;

            if (rooted) {
                auto l_front = _nodes.begin();
                advance(l_front, ++curr_node_index);
                nd = &*l_front;

                auto parent = _nodes.begin();
                advance(parent, curr_node_index - 1);
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
                        advance(l_front, curr_node_index);
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
                        advance(l_front, curr_node_index);
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
            if (_index == 0) {
                clearSpeciesForest();
            }
            else {
                clearGeneForest();
            }
//            clear();
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

    inline void Forest::resetIncrements() {
        // clear all branch lengths from forest
        for (auto &nd:_nodes) {
            nd._edge_length = 0.0;
        }
        _increments.clear();
    }

    inline vector<string> Forest::setUpExistingLineagesVector() {
        vector<string> existing_lineages;
        
        unsigned b=0;
        for (auto &nd:_nodes) {
            if (b < Forest::_nspecies) {
                existing_lineages.push_back(nd._name);
                b++;
            }
        }
        
        return existing_lineages;
    }

    inline void Forest::chooseSpeciesIncrementFromNewick(vector<string> existing_lineages) {
        // not conditioning on anything
        unsigned nlineages = (int) existing_lineages.size();
        
        // hybridization prior
        double rate = (_lambda + _hybridization_rate)*nlineages;

        _last_edge_length = rng.gamma(1.0, 1.0/rate);

        for (auto &nd:_nodes) {
            // add increment to tip nodes not already involved in a join
            if (count(existing_lineages.begin(), existing_lineages.end(), nd._name)) {
                // add increment to nodes in existing lineages
                nd._edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            // add increment to tip nodes that have already joined
        }
        double nChooseTwo = _lineages.size()*(_lineages.size()-1);
        double log_prob_join = log(2/nChooseTwo);
        
        double log_increment_prior = log(rate)-(_last_edge_length*rate) + log_prob_join;
        
        _increments.push_back(make_pair(_last_edge_length, log_increment_prior));
    }

    inline vector<string> Forest::updateExistingLineagesVector(vector<string> existing_lineages, tuple<string, string, string> species_joined) {
        string species1 = get<0>(species_joined);
        string species2 = get<1>(species_joined);
        string new_spp = get<2>(species_joined);
        
        // do not need to update for first generation since nothing has been joined
        if (species1 != "null" && species2 != "null") {
            // remove species1 and species2
            existing_lineages.erase(remove(existing_lineages.begin(), existing_lineages.end(), species1), existing_lineages.end());
            existing_lineages.erase(remove(existing_lineages.begin(), existing_lineages.end(), species2), existing_lineages.end());
            
            // ad new species
            existing_lineages.push_back(new_spp);
        }
        
        return existing_lineages;
    }

    inline vector<pair<tuple<string, string, string>, double>>  Forest::buildFromNewickTopology(const std::string newick, bool topology_only) {
        // assume tree is rooted
        // do not allow polytomies
        bool rooted = true;
        bool allow_polytomies = false;
        
        vector<pair<tuple<string, string, string>, double>> species_joined;
        species_joined.push_back(make_pair(make_tuple("null", "null", "null"), 0.0));
        
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
            advance(l_front, curr_node_index);
            Node *nd = &*l_front;

            if (rooted) {
                auto l_front = _nodes.begin();
                advance(l_front, ++curr_node_index);
                nd = &*l_front;

                auto parent = _nodes.begin();
                advance(parent, curr_node_index - 1);
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
                        advance(l_front, curr_node_index);
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
                        advance(l_front, curr_node_index);
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
            if (_index == 0) {
                clearSpeciesForest();
            }
            else {
                clearGeneForest();
            }
//            clear();
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
        
        // reset node numbers and names that are not tips
        int j = 0;
        for (auto &nd:_nodes) {
            nd._number = j;
            if (nd._name == "") {
                nd._name=boost::str(boost::format("node-%d")%nd._number);
            }
            j++;
        }
        
        refreshPreorder();
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        for (auto &entry:heights_and_nodes) {
            species_joined.push_back(make_pair(make_tuple(entry.second->_left_child->_name, entry.second->_left_child->_right_sib->_name, entry.second->_name), 0.0));

        }
        
        if (!topology_only) {
            double height = 0.0;
            // include branch lengths
            for (unsigned a=0; a < heights_and_nodes.size(); a++) {
                species_joined[a].second = heights_and_nodes[a].first - height;
                height += species_joined[a].second;
            }
        }
        double height = 0.0;
        for (auto &s:species_joined) {
            height += s.second;
        }
        return species_joined;
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
        for (unsigned i=0; i < _depths.size(); i++) {
            if (_depths[i].second.first == _depths[i].second.second) {
                // species have already been joined in the speceis tree, so they are no longer a constraint
                _depths.erase(_depths.begin()+i);
            }
        }
        assert (_depths.size() > 0);
    }

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

        for (unsigned i=0; i < heights_and_nodes.size(); i++) {
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
                                left_child = left_child->_left_child;
                                // if the species of the node is not found, the parent will be dealt with later; ignore the node for now
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
                        left_child->_parent->_visited = true;
                        done = true;
                    }
                    if (spp_left_child != spp_right_child) {
                        double height = getLineageHeight(nd->_left_child);
                        _depths.push_back(make_pair(height, make_pair(spp_left_child, spp_right_child)));
                    }
                    spp_left_child = "";
                    spp_right_child = "";
                }
            }
        }
        assert (_depths.size() > 0);

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

    inline double Forest::calcLogCoalLikeGivenTheta(double proposed_theta, vector<pair<tuple<string, string, string>, double>> species_info, bool both) {
        // walk through species increments and calculate coalescent likelihood
        int nincrements = 0.0;
        double log_coalescent_likelihood = 0.0;
        
        for (unsigned species_join_number = 0; species_join_number < species_info.size(); species_join_number++) {
            
            double neg_inf = -1*numeric_limits<double>::infinity();
            vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
            double cum_time = 0.0;
            int a = 0;
            
            double species_increment = species_info[species_join_number].second;
            
            tuple<string, string, string> species_joined = species_info[species_join_number].first;
            
            // update species partition with species_joined
            updateSpeciesPartition(species_joined);
            
            double species_tree_height = species_info[species_join_number].second;
            if (species_tree_height == 0.0) {
                species_tree_height = species_info[species_join_number-1].second;
            }
            
            species_tree_height = 0.0;
            for (unsigned i=0; i<species_join_number+1; i++) {
                species_tree_height += species_info[i].second;
            }
            assert (species_tree_height > 0.0);
            
            if (species_increment > 0) {
                for (unsigned i=nincrements; i < heights_and_nodes.size(); i++) {
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
                                    species = s.first;
                                    node = nd;
                                }
                            }
                        }
                        for (auto &s:_species_partition) {
                            if (s.first == species) {
                                // coalescence
                                unsigned nlineages = (int) s.second.size();
                                
                                if (nlineages == 1) {
                                    log_coalescent_likelihood = neg_inf;
                                    return log_coalescent_likelihood;
                                }
                                
                                double coalescence_rate = nlineages*(nlineages-1) / proposed_theta;
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
                                
                                double coalescence_rate = nlineages*(nlineages-1) / proposed_theta;
                                log_coalescent_likelihood -= increment * coalescence_rate;
                                
//                                cout << "pr no coalescence = " << -1 * increment * coalescence_rate << endl;
                            }
                        }
                        
                    }
                }
                double remaining_chunk_of_branch = species_increment - cum_time;
                // no coalescence
                for (auto &s:_species_partition) {
                    unsigned nlineages = (int) s.second.size();
                    double coalescence_rate = nlineages*(nlineages-1) / proposed_theta;
                    log_coalescent_likelihood -= remaining_chunk_of_branch * coalescence_rate;
//                    cout << "pr no coalescence = " << -1 * remaining_chunk_of_branch * coalescence_rate << endl;
                }
                nincrements += a;
            }
            
            else {
                // final step; no deep coalescence; one species
                assert (species_increment == 0.0);
                assert (_species_partition.size() == 1);
                
                for (unsigned i=nincrements; i < heights_and_nodes.size(); i++) {
                    // there must be coalescence at this point
                    Node* node = nullptr;
                    string species;
                    bool found = false;
                    for (auto &s:_species_partition) {
                        for (auto &nd:s.second) {
                            if (nd == heights_and_nodes[i].second->_left_child) {
                                node = nd;
                                found = true;
                            }
                        }
                        assert (found);
                        for (auto &s:_species_partition) {
                            // coalescence
                            unsigned nlineages = (int) s.second.size();
                            
                            double coalescence_rate = nlineages*(nlineages-1) / proposed_theta;
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
                assert (log_coalescent_likelihood != 0.0);
                return log_coalescent_likelihood;
            }
        }
        assert (log_coalescent_likelihood != 0.0);
        return log_coalescent_likelihood;
    }

    inline double Forest::calcCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height) {
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
                                species = s.first;
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
                            
                            double coalescence_rate = nlineages*(nlineages-1) / Forest::_theta;
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
                            
                            double coalescence_rate = nlineages*(nlineages-1) / _theta;
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
                double coalescence_rate = nlineages*(nlineages-1) / _theta;
                log_coalescent_likelihood -= remaining_chunk_of_branch * coalescence_rate;
                
//                cout << "pr no coalescence = " << -1 * remaining_chunk_of_branch * coalescence_rate << endl;
            }
            _nincrements += a;
        }
        
        // calculate coalescent likelihood for the rest of the panmictic tree
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
                double coalescence_rate = panmictic_nlineages*(panmictic_nlineages-1) / _theta;
                double nChooseTwo = panmictic_nlineages*(panmictic_nlineages-1);
                double log_prob_join = log(2/nChooseTwo);
//                log_coalescent_likelihood += log_prob_join + log(coalescence_rate) - (increment * coalescence_rate);
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
                        
                        double coalescence_rate = nlineages*(nlineages-1) / _theta;
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
        return log_coalescent_likelihood;
    }

    inline int Forest::selectPair(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u = rng.uniform();
        double cum_prob = 0.0;
        int index = 0.0;
        for (unsigned i=0; i < weight_vec.size(); i++) {
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
        for (unsigned a = 0; a < likelihood_vec.size(); a++) {
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

        double log_max_weight = *max_element(log_weight_choices.begin(), log_weight_choices.end());
        for (auto & i:log_weight_choices) {
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
        _node_choices.push_back(make_pair(subtree1, subtree2));

        return s;
    }

    inline tuple<Node*, Node*, Node*> Forest::createNewSubtreeFromPrior(pair<unsigned, unsigned> t, double increment) {
        
        Node* subtree1 = _lineages[t.first];
        Node* subtree2 = _lineages[t.second];

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
        new_nd->_partial=ps.getPartial(_npatterns*4, _index-1);
        assert(new_nd->_left_child->_right_sib);
        calcPartialArray(new_nd);

        // don't update the species list
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        _node_choices.push_back(make_pair(subtree1, subtree2));
        
        return make_tuple(subtree1, subtree2, new_nd);
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
        new_nd->_partial=ps.getPartial(_npatterns*4, _index-1);
        assert(new_nd->_left_child->_right_sib);
        calcPartialArray(new_nd);

        // don't update the species list
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        // calculate the coalescent likelihood for the first join since it's the same for all joins
        
#if defined(GENE_TREE_COALESCENT_LIKELIHOOD)
            // update increments and priors
            double log_increment_prior = 0.0;
            bool coalescence = false;
            for (auto &s:_species_partition) {
                if (s.first == species) {
                    coalescence = true;
                }
                else {
                    coalescence = false;
                }
                
                if (coalescence) {
                    // if there is coalescence, need to use number of lineages before the join
                    double coalescence_rate = (s.second.size())*(s.second.size()-1) / _theta;
                    assert (coalescence_rate > 0.0); // rate should be >0 if there is coalescence
                    double nChooseTwo = (s.second.size())*(s.second.size()-1);
                    double log_prob_join = log(2/nChooseTwo);
                    log_increment_prior += log(coalescence_rate) - (increment*coalescence_rate) + log_prob_join;
                }
                else {
                    double coalescence_rate = s.second.size()*(s.second.size() - 1) / _theta;
                    log_increment_prior -= increment*coalescence_rate;
                }
            }
            
            if (_deep_coalescent_increments.size() > 0) { // include any deep coalescence in the coalescent likelihood but do not clear yet
                for (auto &d:_deep_coalescent_increments) {
                    increment += d.first;
                    log_increment_prior += d.second;
                }
            }
            _gene_tree_log_coalescent_likelihood += log_increment_prior;
#else
            _gene_tree_log_coalescent_likelihood = 0.0;
#endif
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

    inline void Forest::createDefaultTree() {
//        clear();
        clearGeneForest();
        //create taxa
        double edge_length = rng.gamma(1.0, 1.0/_ntaxa);
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
        _count = other._count;
        _nstates = other._nstates;
        _npatterns = other._npatterns;
        _nodes.clear();
        _nodes.resize(other._nodes.size());
        _lineages.resize(other._lineages.size());
        _new_nodes.resize(other._new_nodes.size());
        _preorder.resize(other._preorder.size());
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
        _hybrid_species_joined = other._hybrid_species_joined;
        _migration_rate = other._migration_rate;
        _hybridization_rate = other._hybridization_rate;
        _last_direction = other._last_direction;
        _gamma = other._gamma;
        _gene_tree_log_weight = other._gene_tree_log_weight;
        _prev_log_likelihood = other._prev_log_likelihood;
        _log_weight_vec = other._log_weight_vec;
        _theta = other._theta;
        _increments = other._increments;
        _topology_prior = other._topology_prior;
        _log_joining_prob = other._log_joining_prob;
        _increment_choices = other._increment_choices;
        _deep_coalescent_increments = other._deep_coalescent_increments;
        _extended_increment = other._extended_increment;
        _species_join_number = other._species_join_number;
        _prev_gene_tree_log_likelihood = other._prev_gene_tree_log_likelihood;
        _ready_to_join_species = other._ready_to_join_species;
        _depths = other._depths;
        _nincrements = other._nincrements;
        _gene_tree_log_coalescent_likelihood = other._gene_tree_log_coalescent_likelihood;
        _panmictic_coalescent_likelihood = other._panmictic_coalescent_likelihood;

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
            else {
                nd->_right_sib = 0;
            }
                
            nd->_number = othernd._number;
            nd->_name = othernd._name;
            nd->_edge_length = othernd._edge_length;
            nd->_position_in_lineages = othernd._position_in_lineages;
            nd->_partial = othernd._partial;
            nd->_visited = othernd._visited;
            nd->_hybrid_newick_name = othernd._hybrid_newick_name;
            nd->_n_descendants = othernd._n_descendants;
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
//        assert (_index==0);
        assert (_nspecies == (unsigned) species_names.size());
//        clear();
        clearSpeciesForest();
        //create species
//        double edge_length = 0.0;
        for (unsigned i = 0; i < _nspecies; i++) {
            Node* nd = &*next(_nodes.begin(), i);
//            nd->_right_sib=0;
            nd->_name=species_names[i];
//            nd->_left_child=0;
//            nd->_right_sib=0;
//            nd->_parent=0;
//            nd->_number=i;
//            nd->_edge_length = edge_length;
//            nd->_position_in_lineages=i;
            }
        _nleaves=_nspecies;
        _ninternals=0;
        assert (_nodes.size() == _nspecies);
        assert (_lineages.size() == _nspecies);
//        _nodes.resize(_nspecies);
//        _lineages.resize(_nspecies);
    }

    inline void Forest::chooseSpeciesIncrement(double max_depth) {
        // conditioning on max gene tree height
        assert (max_depth >= 0.0);
        if (max_depth > 0.0) {
            // hybridization prior
            double rate = (_lambda+_hybridization_rate)*_lineages.size();
            
            double u = rng.uniform();
            double inner_term = 1-exp(-rate*max_depth);
            _last_edge_length = -log(1-u*inner_term)/rate;
            assert (_last_edge_length < max_depth);

            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            double nChooseTwo = _lineages.size()*(_lineages.size()-1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = log(rate) - (_last_edge_length*rate) - log(1 - exp(-rate*max_depth)) + log_prob_join;
            
            _increments.push_back(make_pair(_last_edge_length, increment_prior));
        }
        else {
            // hybridization prior
            double rate = (_lambda +_hybridization_rate)*_lineages.size();

            _last_edge_length = rng.gamma(1.0, 1.0/rate);

            for (auto nd:_lineages) {
                nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
            }
            
            double nChooseTwo = _lineages.size()*(_lineages.size()-1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = log(rate) - (_last_edge_length*rate) - log(1 - exp(-rate*max_depth)) + log_prob_join;
            
            _increments.push_back(make_pair(_last_edge_length, increment_prior));
        }
    }

    inline tuple<string,string, string> Forest::speciesTreeProposal() {
        // this function creates a new node and joins two species

        Node* subtree1 = nullptr;
        Node* subtree2 = nullptr;
        
        bool done = false;
        while (!done) {
            pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size());
            subtree1=_lineages[t.first];
            subtree2=_lineages[t.second];
            assert(!subtree1->_parent && !subtree2->_parent);
            
            if (_lineages.size() == 2) {
                done = true;
            }
            else if (subtree1->_name != _outgroup && subtree2->_name != _outgroup) {
                done = true;
            }
        }

        if (_lineages.size() > 2) {
            assert (subtree1->_name != _outgroup);
            assert (subtree2->_name != _outgroup);
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

        calcTopologyPrior((int) _lineages.size());

        updateNodeVector (_lineages, subtree1, subtree2, new_nd);

        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
    }

    inline void Forest::setUpGeneForest(map<string, string> &taxon_map) {
        assert (_index >0);
        _species_partition.clear();
        
        for (auto &nd:_nodes) {
            if (!nd._left_child) {
                string species_name = taxon_map[nd._name];
                _species_partition[species_name].push_back(&nd);
            }
        }
        assert (_species_partition.size() > 0);
    }

    inline void Forest::drawFromGeneTreePrior() { // TODO: don't need to join and rejoin for choice of 2, but this should rarely happen
        assert (_lineages.size() > 1);
        
        double coalescence_rate = (_lineages.size()*(_lineages.size() - 1)) / _theta;
        
        double increment = rng.gamma(1.0, 1.0/(coalescence_rate));
        
        for (auto &nd:_lineages) {
            nd->_edge_length += increment;
        }
            
        // choose taxa to join
        Node* subtree1;
        Node* subtree2;
        
        // prior-prior proposal
        if (_proposal == "prior-prior") {
            pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size());
            
            subtree1 = _lineages[t.first];
            subtree2 = _lineages[t.second];
        }
        
        else {
            if (_lineages.size() != _ntaxa) {
                _prev_gene_tree_log_likelihood = _gene_tree_log_likelihood;
            }
            else {
                _prev_gene_tree_log_likelihood = 0.0;
            }
            pair<Node*, Node*> t = chooseAllPairsFromPrior(increment);
            
            subtree1 = t.first;
            subtree2 = t.second;
        }
            
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

        //always calculating partials now
        assert (new_nd->_partial == nullptr);
        new_nd->_partial=ps.getPartial(_npatterns*4, _index-1);
        assert(new_nd->_left_child->_right_sib);
        
        //update species list
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);

        calcPartialArray(new_nd);
        calcLogLikelihood();
    }

    inline double Forest::calcLogSpeciesTreeDensity(double lambda) {
        assert (_index == 0);
        
        refreshPreorder();
        
        // Assume that this species forest is fully resolved
//        assert(_preorder.size() == 1);
        assert(_lineages.size() == 1);
                
        vector< pair<double, Node *>> heights_and_nodes = sortPreorder();
        
        // Build vector of internal node heights
        vector<double> internal_heights;
        
        for (unsigned i=0; i < heights_and_nodes.size(); i++) {
            internal_heights.push_back(heights_and_nodes[i].first);
        }
        
        // Number of internal nodes should be _nspecies - 1
        assert(internal_heights.size() == Forest::_nspecies - 1);
        
        // internal heights are already sorted
        
        double log_prob_density = 0.0;
        unsigned n = Forest::_nspecies;
        double h0 = 0.0;
        for (auto it = internal_heights.begin(); it != internal_heights.end(); ++it) {
            double h = *it;
            double r = lambda * n;
            double logr = log(r);
            double t = h - h0;
            double log_exponential_density = logr - r*t;
            log_prob_density += log_exponential_density;
            h0 = h;
            n--;
        }
        
        return log_prob_density;
    }

    inline void Forest::resetLineages() {
        assert (_lineages.size () == 1);
//        // this function rebuilds the _lineages vector, setting _position_in_lineages
        _lineages.clear();
        unsigned i=0;
        for (auto &nd:_nodes) {
            if (!nd._left_child) {
                _lineages.push_back(&nd);
                nd._position_in_lineages = i;
                i++;
            }
        }

        for (auto &nd : _nodes) {
            if (nd._left_child) {
                updateNodeVector(_lineages, nd._left_child, nd._left_child->_right_sib, &nd);
            }
        }
    }


    inline pair<double, string> Forest::chooseDelta(vector<pair<tuple<string, string, string>, double>> species_info) {
        // assert forest is not fully resolved
        
        assert (_lineages.size() > 1);
         // get species info
        double species_increment = species_info[_species_join_number].second;
        assert (species_increment >= 0.0);
        
        if (_species_join_number == 0) {
//            assert (species_increment > 0.0);
        }

         // join species if necessary
         if (_ready_to_join_species && _species_partition.size() > 1) {
             _ready_to_join_species = false;
             
             // update species partition
             _species_join_number++;
             if (_species_join_number > species_info.size() - 1) {
                 _species_join_number = (int) species_info.size() - 1;
             }

             string species1 = get<0> (species_info[_species_join_number].first);
             string species2 = get<1> (species_info[_species_join_number].first);
             string new_name = get<2> (species_info[_species_join_number].first);
             
             list<Node*> &nodes = _species_partition[new_name];
             copy(_species_partition[species1].begin(), _species_partition[species1].end(), back_inserter(nodes));
             copy(_species_partition[species2].begin(), _species_partition[species2].end(), back_inserter(nodes));
             _species_partition.erase(species1);
             _species_partition.erase(species2);

             species_increment = species_info[_species_join_number].second;
         }
         
      // calculate coalescence rate for each population
         double coalescence_rate = 0.0;
         vector<double> population_coalescent_rates;
         vector<string> eligible_species; // only possibility of coalescing if species has >1 lineage present
         
         for (auto &s:_species_partition) {
             if (s.second.size() > 1) {
                 eligible_species.push_back(s.first);
                 double population_coalescence_rate = s.second.size()*(s.second.size()-1)/_theta;
                 population_coalescent_rates.push_back(population_coalescence_rate);
                 coalescence_rate += population_coalescence_rate;
             }
         }

         assert(eligible_species.size() > 0);
         
         // use combined rate to draw an increment (delta)
         double increment = rng.gamma(1.0, 1.0/(coalescence_rate));
//        cout << increment << endl;
         
         // choose which species coalescent event occurred in
         
         for (auto &p:population_coalescent_rates) {
             p = log(p/coalescence_rate);
         }
         int index = selectPair(population_coalescent_rates);
         
         string species_for_join = eligible_species[index];
                 
         double species_tree_height = 0.0;
         for (unsigned a=0; a<_species_join_number+1; a++) {
             species_tree_height += species_info[a].second;
         }
         
         bool done = false;
//        cout << "species increment is " << species_increment << endl;
//        cout << "increment is " << increment << endl;
//        cout << "species partition size is " << _species_partition.size() << endl;
//        for (auto &s:_species_partition) {
//            cout << "s first is " << s.first << " and s second is " << s.second.size() << endl;
//        }
//        cout << "species join number is " << _species_join_number << endl;
//        for (auto &s:species_info) {
//            cout << "species info is: " << s.second << endl;
//        }
//
//        showForest();
         if (increment > species_increment && _species_partition.size() > 1 && species_increment > 0.0) {
             // deep coalescence
             double cum_time = species_increment;
             while (!done) {
                 // extend existing lineages to species barrier
                 extendGeneTreeLineages(species_tree_height, true);
                 
                 // update species partition - two species must merge now to accommodate deep coalescence
                 _species_join_number++;
                 assert (_species_join_number <= species_info.size());
                 if (_species_join_number > species_info.size()-1) {
                     _species_join_number = (int) species_info.size()-1;
                 }
                 
                 species_tree_height += species_info[_species_join_number].second;
                 
                 string species1 = get<0> (species_info[_species_join_number].first);
                 string species2 = get<1> (species_info[_species_join_number].first);
                 string new_name = get<2> (species_info[_species_join_number].first);
                 
                 list<Node*> &nodes = _species_partition[new_name];
                 copy(_species_partition[species1].begin(), _species_partition[species1].end(), back_inserter(nodes));
                 copy(_species_partition[species2].begin(), _species_partition[species2].end(), back_inserter(nodes));
                 _species_partition.erase(species1);
                 _species_partition.erase(species2);
                                  
                 species_increment = species_info[_species_join_number].second;
                 assert (species_increment >= 0.0);
                 
                 cum_time += species_increment;
                 
                 // choose a new species
                 eligible_species.clear();
                 population_coalescent_rates.clear();
                 coalescence_rate = 0.0;
                 
                 for (auto &s:_species_partition) {
                     if (s.second.size() > 1) {
                         eligible_species.push_back(s.first);
                         double population_coalescence_rate = s.second.size()*(s.second.size()-1)/_theta;
                         population_coalescent_rates.push_back(population_coalescence_rate);
                         coalescence_rate += population_coalescence_rate;
                     }
                 }
                 
                 assert (eligible_species.size() > 0);
                 
                    // choose which species coalescent event occurred in
                 for (auto &p:population_coalescent_rates) {
                     p = log(p/coalescence_rate);
                 }
                 int index = selectPair(population_coalescent_rates);
                 species_for_join = eligible_species[index];
                 
                 // draw a new increment
                 unsigned nlineages = 0;
                 for (auto &s:_species_partition) {
                     if (s.first == species_for_join) {
                         nlineages = (int) s.second.size();
                         break;
                     }
                 }
                 
                 assert (nlineages > 1);
                 
                 double deep_coalescence_rate = nlineages*(nlineages-1)/_theta;
                 double deep_coalescence_increment = rng.gamma(1.0, 1.0/(deep_coalescence_rate));
                 
                 for (auto &nd:_lineages) {
                     nd->_edge_length += deep_coalescence_increment;
                     for (int i=0; i<1; i++) {
                     }
                 }
                 
                 increment = deep_coalescence_increment;
                 
                 // check if there is another deep coalescent event
                 if (deep_coalescence_increment < species_increment || species_increment == 0.0) {
                     done = true;
                 }
             }
         }
         
         else {
             for (auto &nd:_lineages) {
                 nd->_edge_length += increment;
             }
         }
         
         return make_pair(increment, species_for_join);
     }

    inline void Forest::evolveSpeciesFor(list<Node*> &nodes, double increment, string species) {
        calcTopologyPrior((int) nodes.size());
        allowCoalescence(nodes, increment, species);
    }

    inline void Forest::allowCoalescence(list<Node*> &nodes, double increment, string species) {
        Node *subtree1 = nullptr;
        Node *subtree2 = nullptr;
        unsigned s = (unsigned) nodes.size();
        bool one_choice = false;

        // prior-prior proposal
        if (_proposal == "prior-prior") {
            pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
            auto it1 = std::next(nodes.begin(), t.first);
            subtree1 = *it1;

            auto it2 = std::next(nodes.begin(), t.second);
            subtree2 = *it2;
            assert (t.first < nodes.size() && t.second < nodes.size());
        }
        else {
            if (_lineages.size() != _ntaxa) {
                _prev_gene_tree_log_likelihood = _gene_tree_log_likelihood;
            }
            else {
                _prev_gene_tree_log_likelihood = 0.0;
            }
            
            if (nodes.size() > 3) {
                pair<Node*, Node*> t = chooseAllPairs(nodes, increment, species);
                
                subtree1 = t.first;
                subtree2 = t.second;
            }
            else {
                one_choice = true;
                
                subtree1 = nodes.front();
                subtree2 = nodes.back();
            }
                
        }
        
        // if only 1 choice, subtree has already been created in chooseAllPairs function
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

        //always calculating partials now
        assert (new_nd->_partial == nullptr);
        new_nd->_partial=ps.getPartial(_npatterns*4, _index-1);
        assert(new_nd->_left_child->_right_sib);
        
        _new_nodes.push_back(new_nd);
        calcPartialArray(new_nd);

        //update species list
        updateNodeList(nodes, subtree1, subtree2, new_nd);
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        if (_proposal == "prior-prior" || one_choice) {
            double prev_log_likelihood = _gene_tree_log_likelihood;
            _gene_tree_log_likelihood = calcLogLikelihood();
            vector<double> choice;
            choice.push_back(_gene_tree_log_likelihood);
            // don't need to recalculate this for prior-post with >1 choice
            
            vector<double> log_weight = reweightChoices(choice, prev_log_likelihood);
            double log_weight_choices_sum = getRunningSumChoices(log_weight);
            _gene_tree_log_weight = log_weight_choices_sum;
        }
        
        // update increments and priors
        double log_increment_prior = 0.0;
        for (auto &s:_species_partition) {
            bool coalescence = false;
            for (auto &nd:s.second) {
                if (nd == new_nd) {
                    coalescence = true;
                    break;
                }
                else {
                    coalescence = false;
                }
            }
            if (coalescence) {
                // if there is coalescence, need to use number of lineages before the join
                double coalescence_rate = (s.second.size()+1)*(s.second.size()) / _theta;
                assert (coalescence_rate > 0.0); // rate should be >0 if there is coalescence
                double nChooseTwo = (s.second.size()+1)*(s.second.size());
                double log_prob_join = log(2/nChooseTwo);
                log_increment_prior += log(coalescence_rate) - (increment*coalescence_rate) + log_prob_join;
            }
            else {
                double coalescence_rate = s.second.size()*(s.second.size() - 1) / _theta;
                log_increment_prior -= increment*coalescence_rate;
            }
        }
        if (_deep_coalescent_increments.size() > 0) {
            for (auto &d:_deep_coalescent_increments) {
                increment += d.first;
                log_increment_prior += d.second;
            }
        }
#if defined(GENE_TREE_COALESCENT_LIKELIHOOD)
        _gene_tree_log_coalescent_likelihood += log_increment_prior;
#else
        _gene_tree_log_coalescent_likelihood = 0.0;
#endif
        _increments.push_back(make_pair(increment, log_increment_prior));
        
        _extended_increment = 0.0;
        _deep_coalescent_increments.clear();
        
        
    }

    inline void Forest::geneTreeProposal(pair<double, string> species_info, vector<pair<tuple<string, string, string>, double>> _t) {
        string species_name = species_info.second;
        
        bool joined = false;
        double increment = species_info.first;
        
        for (auto &s:_species_partition) {
            if (s.first == species_name) {
                evolveSpeciesFor(s.second, increment, species_name);
                joined = true;
                break;
            }
        }
        assert (joined);
        
        bool extend = true;
        if (_species_partition.size() > 1) {
            for (auto &s:_species_partition) {
                if (s.second.size() != 1) {
                    extend = false;
                    break;
                }
            }
        }
        
        double species_tree_height = 0.0;
        for (unsigned a=0; a<_species_join_number+1; a++) {
            species_tree_height += _t[a].second;
        }
        
        if (extend && _species_partition.size() > 1) {
            extendGeneTreeLineages(species_tree_height, false);
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
        for (unsigned i=0; i < _lineages.size(); i++) {
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
        for (unsigned i=0; i < _lineages.size(); i++) {
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

        if (position1 > (int) node_vector.size()+1) {
            position1 = position1-1;
            if (position2 > 0) {
                position2 = position2-1;
            }
        }
        if (position2 > (int) node_vector.size()+1) {
            position2 = position2-1;
            if (position1 > 0) {
                position1 = position1-1;
            }
        }

        assert (position1 != position2);

        // lower position must be inserted first
        if (position1 < position2) {
            node_vector.insert(node_vector.begin()+position1, iter1);
            node_vector.insert(node_vector.begin()+position2, iter2);
        }
        else {
            node_vector.insert(node_vector.begin()+position2, iter2);
            node_vector.insert(node_vector.begin()+position1, iter1);
        }

        // reset _position_in_lineages
        for (unsigned i=0; i < _lineages.size(); i++) {
            _lineages[i] -> _position_in_lineages=i;
        }
    }

    inline void Forest::revertNodeList(list<Node *> &node_list, Node *addnode1, Node *addnode2, Node *delnode1) {
        // Delete delnode1 from node_vector
        auto it = find(node_list.begin(), node_list.end(), delnode1);
        assert (it != node_list.end());
        node_list.erase(it);

        // find positions of nodes to insert
        auto position1 = addnode1->_position_in_lineages;
//        auto iter1 = addnode1;
        auto iter1 = node_list.begin();
        advance(iter1, position1);

        auto position2 = addnode2->_position_in_lineages;
//        auto iter2 = addnode2;
        auto iter2 = node_list.begin();
        advance(iter2, position2);

        // lower position must be inserted first
        if (position1 < position2) {
            node_list.insert(iter1, addnode1);
            node_list.insert(iter2, addnode2);
        }
        else {
            node_list.insert(iter2, addnode1);
            node_list.insert(iter1, addnode2);
        }

        // reset _position_in_lineages
        for (unsigned i=0; i < _lineages.size(); i++) {
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

    inline void Forest::hybridizeGene(vector<string> hybridized_nodes, double species_tree_increment) {
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
//            evolveSpeciesFor(s.second, species_tree_increment);
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
        hybridGeneTreeProposal(species_tree_increment);

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
        hybridGeneTreeProposal(species_tree_increment);

        // save likelihood
        likelihood_vec.push_back(calcLogLikelihood());

        // choose major or minor path
        int index_of_choice = chooseDirectionOfHybridization(likelihood_vec);
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

            resetToMinor(minor_nodes, minor_left_children, minor_right_children, minor_left_edge_lengths, minor_right_edge_lengths);

            // revert all _lineages edge lengths to minor nodes
            for (unsigned i=0; i < _lineages.size(); i++) {
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
        for (auto &l:_lineages) {
            assert (!l->_parent);
            assert(!l->_right_sib);
        }
        }

    inline void Forest::hybridGeneTreeProposal(double species_tree_increment) {
        if (_species_partition.size() == 1) {
//            fullyCoalesceGeneTree(_species_partition.begin()->second);
        }

        else {
            for (auto &s:_species_partition) {
                assert (s.second.size()>0);
//                evolveSpeciesFor(s.second, species_tree_increment, 0.0); // TODO: fix species tree height
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

    inline int Forest::chooseDirectionOfHybridization(vector<double> likelihood_vec) {
        _node_choices.clear();
        _log_likelihood_choices.clear();

        _log_likelihood_choices.push_back(likelihood_vec[0]+log(.15)); // multiply minor likelihood by (1-gamma)
        _log_likelihood_choices.push_back(likelihood_vec[1]+log(.85)); // multiply major likelihood by (gamma)

        // reweight each choice of pairs
        vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, _prev_log_likelihood); // TODO: this is incorrect

        // sum unnormalized weights before choosing the pair
        _gene_tree_log_weight = 0.0;
        for (auto &l:log_weight_choices) {
            _gene_tree_log_weight += l;
        }

        // normalize weights
        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        for (unsigned b=0; b < log_weight_choices.size(); b++) {
            log_weight_choices[b] -= log_weight_choices_sum;
        }

        // select a direction
        _index_of_choice = selectPair(log_weight_choices);
        return _index_of_choice;
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
        for (unsigned a = (int) _new_nodes.size()-1; a>=0; a--) {
            revertNodeVector(_lineages, _new_nodes[a]->_left_child, _new_nodes[a]->_left_child->_right_sib, _new_nodes[a]);
        }

        // reset _lineages edge lengths
        assert (_lineages.size() == branch_lengths.size());
        for (unsigned i=0; i < _lineages.size(); i++) {
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
        double rate = (_lambda)*_lineages.size();

        // choose edge length but don't add it yet
        _last_edge_length = rng.gamma(1.0, 1.0/rate);

        if (_lineages.size()>1) {
            double nChooseTwo = _lineages.size()*(_lineages.size()-1);
            double log_prob_join = log(2/nChooseTwo);
            double increment_prior = (log(rate)-_last_edge_length*rate) + log_prob_join;
            
            _increments.push_back(make_pair(_last_edge_length, increment_prior));
        }

        // add the previously chosen edge length
        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
    }

    inline double Forest::calcTopologyPrior(unsigned nlineages) {
        _log_joining_prob += -log(0.5*nlineages*(nlineages-1));
        assert (!isinf(_log_joining_prob));
        return _log_joining_prob;
    }

    inline unsigned Forest::countDescendants(Node* nd, unsigned count) {
        for (Node * child = nd->_left_child; child; child=child->_right_sib) {
            if (child->_left_child) {
                count++;
                count += countDescendants(child, 0);
            }
            else {
                count += 0;
            }
        }
        return count;
    }

    inline void Forest::extendGeneTreeLineages(double species_tree_height, bool deep_coalescence) {
        // pick one node and figure out increment to add from there
//        vector<double> extended_increment_options;
//        double deep_coalescent_increment = 0.0;
        double extended_increment = 0.0;
        
        if (_lineages.size() > 1) {
            extended_increment = species_tree_height - getLineageHeight(_lineages[0]);
            
            for (auto &l:_lineages) {
                if (l->_left_child) {
                    if (getLineageHeight(l) < species_tree_height) {
                        l->_edge_length += extended_increment;
                    }
                }
                else {
                    if (l->_edge_length < species_tree_height) {
                        l->_edge_length += extended_increment;
                    }
                }
            }
        }
        else {
            extended_increment = species_tree_height - getLineageHeight(_lineages[0]->_left_child);
            _lineages[0]->_left_child->_edge_length = species_tree_height - getLineageHeight(_lineages[0]->_left_child);
            _lineages[0]->_left_child->_right_sib->_edge_length = species_tree_height - getLineageHeight(_lineages[0]->_left_child->_right_sib);
        }
        
        _ready_to_join_species = true;
        
        if (deep_coalescence) {
            
            double deep_coalescent_prior = 0.0; // must account for a deep coalescence scenario
            
            for (auto &s:_species_partition) {
                unsigned nlineages = (int) s.second.size();
                double rate = (nlineages * (nlineages-1) / _theta);
                deep_coalescent_prior -= extended_increment * rate;
            }
            assert (deep_coalescent_prior != 0.0);
            
            _deep_coalescent_increments.push_back(make_pair(extended_increment, deep_coalescent_prior));
        }
        
        _extended_increment = extended_increment;
    }

    inline void Forest::remakeGeneTree(map<string, string> &taxon_map) {
        _lineages.clear(); // TODO: can do this more efficiently without clearing and remaking the list
        _nodes.clear();
        int j=0;
        for (auto &t:taxon_map) {
            Node nd;
            nd._name = t.first;
            nd._number = j;
            nd._edge_length = 0.0;
            nd._position_in_lineages = j;
            _nodes.push_back(nd);
            j++;
        }
        _ninternals = 0;
        _log_likelihood_choices.clear();
        _new_nodes.clear();
        _species_join_number = 0;
        _gene_tree_log_weight = 0.0;
        _gene_tree_log_likelihood = 0.0;
        _node_choices.clear();
        _extended_increment = 0.0;
        _prev_gene_tree_log_likelihood = 0.0;
        _ready_to_join_species = false;
        _deep_coalescent_increments.clear();
        _preorder.clear();
        _increments.clear();
        _extended_increment = 0.0;
        
        for (auto &nd:_nodes) {
            _lineages.push_back(&nd);
        }
        
        // reset the species partition
        setUpGeneForest(taxon_map);
    }

}


