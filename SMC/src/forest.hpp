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
#include "tree_summary.hpp"

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
        void                        createDefaultTree();
        void operator=(const Forest & other);
        void                        debugForest();
        void                        debugLogLikelihood(Node* nd, double log_like);
        double                      calcTopologyPrior(int nlineages);

    private:

        typedef std::vector <double> partial_array_t;
        void                        clear();
        void                        setData(Data::SharedPtr d, int index, map<string, string> &taxon_map);
        Node *                      findNextPreorder(Node * nd);
        std::string                 makeNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>    chooseTaxaToJoin(double s);
        tuple<Node*, Node*, Node*>  createNewSubtree(pair<unsigned, unsigned> p, list<Node*> node_list, double increment);
        void                        calcPartialArray(Node* new_nd);
        void                        setUpGeneForest(map<string, string> &taxon_map);
        void                        setUpSpeciesForest(vector<string> &species_names);
        tuple<string,string, string> speciesTreeProposal();
        void                        firstGeneTreeProposal(double time_increment);
        void                        geneTreeProposal(tuple<string, string, string> &species_merge_info, double time_increment);
        void                        evolveSpeciesFor(list <Node*> &nodes, double time_increment);
        void                        fullyCoalesceGeneTree(list<Node*> &nodes);
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                        hybridizeNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node* delnode3, Node * addnode1);
        void                        revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        void                        revertNodeList(list<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                      getRunningSumChoices(vector<double> &log_weight_choices);
        double                      getRunningSumHybridChoices(vector<double> &log_weight_choices);
        vector<double>              reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        pair<Node*, Node*>          chooseAllPairs(list<Node *> &node_list, double increment);
        pair<Node*, Node*>          getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list);
        int                         selectPair(vector<double> weight_vec);
        void                        chooseSpeciesIncrement();
        void                        addSpeciesIncrement();
        string                      chooseEvent();
        void                        allowMigration(list<Node*> &nodes);
        void                        setGeneration(double g) {_generationf = g;}
        double                      chooseTaxonToMigrate(double s);
        string                      findKeyToDel(Node* taxon_to_migrate);
        void                        migrateTaxon(unsigned taxon_choice, string key_to_del, Node* taxon_to_migrate);
        string                      chooseLineage(Node* taxon_to_migrate, string key_to_del);
        void                        addMigratingTaxon(string key_to_add, string key_to_del, Node* taxon_to_migrate);
        void                        deleteTaxon(string key_to_del, unsigned taxon_choice);
        void                        allowCoalescence(list<Node*> &nodes, double increment);
        void                        allowDummyCoalescence(list<Node*> &nodes, double increment);
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
        double                      findShallowestCoalescence();
        void                        revertToShallowest(double smallest_branch);
        double                      getTreeHeight();
        double                      getLineageHeight(Node* nd);
        void                        revertNewNode(list <Node*> &nodes, Node*new_nd, Node* subtree1, Node* subtree2, double increment);
        void                        revertBranches(list <Node*> &nodes, Node* subtree1, Node* subtree2, double increment);
        void                        chooseCoalescentEvent();
        void                        mergeChosenPair(double species_tree_increment);
        void                        extendGeneTreeLineages(double species_tree_increment);

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
        map<string, list<Node*> > _species_partition;
        double                    _gene_tree_log_likelihood;
        double                      _gene_tree_log_weight;
        vector<double>              _log_weight_vec;
        vector<pair<Node*, Node*>> _node_choices;
        vector<double>              _log_likelihood_choices;
        int                         _index_of_choice;
        pair<Node*, Node*>          _species_joined;
        tuple<Node*, Node*, Node*>  _hybrid_species_joined;
        double                      _generationf = 0;
        string                      _last_direction;
        vector<double>              _rand_numbers;
        vector<pair<double, double>>         _increments;
        double                      _topology_prior;
        unsigned                    _num_coalescent_events_in_generation;
    vector<pair<double, double>> _searchable_branch_lengths; // pair is lineage height, increment
        double                      _gene_tree_cum_time; // gene tree cumulative height in this lineage only
        double                      _prev_log_likelihood;
        vector<double>               _prev_log_likelihood_choices;
        double                      _log_joining_prob;
        vector<double>              _increment_choices;

        void                        showSpeciesJoined();
        double                      calcTransitionProbability(Node* child, double s, double s_child);
        double                      calculateNewEdgeLength(string key_to_add, Node* taxon_to_migrate);
        void                        setNewEdgeLength(double difference, Node* taxon_to_migrate, string key_to_add);
        void                        hybridizeGene(vector<string> hybridized_nodes, double species_tree_increment);
        void                        resetToMinor(vector<Node*> minor_nodes, vector<Node*> minor_left_children, vector<Node*> minor_right_children, vector<double> minor_left_edge_lengths, vector<double> minor_right_edge_lengths);

    public:

        typedef std::shared_ptr<Forest> SharedPtr;
        double               _theta;
        static double                      _starting_theta;
        static double               _speciation_rate;
        static string               _proposal;
        static string               _model;
        static double               _kappa;
        static vector<double>       _base_frequencies;
        static string               _string_base_frequencies;
        static double               _migration_rate;
        static double               _hybridization_rate;
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
        _nodes.resize(_ntaxa);
        _npatterns = 0;
        _nstates = 4;
        _last_edge_length = 0.0;
        _lineages.reserve(_nodes.size());
        _rand_numbers.clear();
        _num_coalescent_events_in_generation = 0;
        _log_joining_prob = 0.0;

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
        _nleaves=_ntaxa;
        _ninternals=0;
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::setData(Data::SharedPtr d, int index, map<string, string> &taxon_map) {
        _data = d;
        _index = index;

        //don't set data for species tree
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
        double sum_height = 0.0;

        sum_height += nd->getEdgeLength();
        if (nd->_left_child) {
            for (Node* child = nd->_left_child; child; child=child->_left_child) {
                sum_height += child->getEdgeLength();
            }
        }
        return sum_height;
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
        }
        else if (_index == 0) {
            cout << " species tree: " << endl;
        }
        cout << " " << makeNewick(9, true) << "\n";
        cout << "\n";
    }

    inline string Forest::makeNewick(unsigned precision, bool use_names) {
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
            _rand_numbers.push_back(t1);
            _rand_numbers.push_back(t2);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = ::rng.randint(0, nsubtrees-1);
                _rand_numbers.push_back(t2);
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
        _rand_numbers.push_back(t1);
        _rand_numbers.push_back(t2);
        _rand_numbers.push_back(t3);

        //keep calling t2 until it doesn't equal t1 or t3
        while (t2 == t1 || t2 == t3) {
            t2 = ::rng.randint(0, nsubtrees-1);
            _rand_numbers.push_back(t2);
        }
        // keep calling t3 until it doesn't equal t1 or t2
        while (t3 == t1 || t3 == t2) {
            t3 = ::rng.randint(0, nsubtrees-1);
            _rand_numbers.push_back(t3);
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
        auto data_matrix=_data->getDataMatrix();

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

    inline pair<Node*, Node*> Forest::chooseAllPairs(list<Node*> &node_list, double increment) {
        _node_choices.clear();
        
        // TODO: not sure about this next chunk of code, check if using prior-post
        // reset _log_likelihood_choices to 0 if we are in a new lineage
        if (_log_likelihood_choices.size() == 0) {
            _prev_log_likelihood = 0.0;
        }
        else {
            _prev_log_likelihood = _log_likelihood_choices[_index_of_choice];
        }
        _log_likelihood_choices.clear();

        _log_likelihood_choices.clear();
        vector<double> topology_prior_choices;

        //choose pair of nodes to try
        for (int i = 0; i < (int) node_list.size()-1; i++) {
            for (int j = i+1; j < (int) node_list.size(); j++) {
                // createNewSubtree returns subtree1, subtree2, new_nd
                tuple<Node*, Node*, Node*> t = createNewSubtree(make_pair(i,j), node_list, increment);
                _log_likelihood_choices.push_back(calcLogLikelihood());

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
        vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, _prev_log_likelihood);

        // sum unnormalized weights before choosing the pair
        _gene_tree_log_weight = 0.0;
        for (auto &l:log_weight_choices) {
            _gene_tree_log_weight += l;
        }

        // normalize weights
        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        for (int b=0; b < (int) log_weight_choices.size(); b++) {
            log_weight_choices[b] -= log_weight_choices_sum;
        }

        // randomly select a pair
        _index_of_choice = selectPair(log_weight_choices);

        // find nodes to join in node_list
        Node *subtree1 = _node_choices[_index_of_choice].first;
        Node *subtree2 = _node_choices[_index_of_choice].second;

        // erase extra nodes created from node list
        for (int i = 0; i < _node_choices.size(); i++) {
            _nodes.pop_back();
        }

        _gene_tree_log_likelihood = _log_likelihood_choices[_index_of_choice];

        return make_pair(subtree1, subtree2);
    }

    inline void Forest::chooseCoalescentEvent() {
        _gene_tree_log_weight = 0.0;
        if (_log_likelihood_choices.size() > 0) {
            // reweight each choice of pairs
            vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, _prev_log_likelihood);
            
            // sum unnormalized weights before choosing the pair
            // must include the likelihoods of all pairs in the final particle weight
//            _gene_tree_log_weight = getRunningSumChoices(log_weight_choices);
            
            
//            for (auto &l:log_weight_choices) {
//                _gene_tree_log_weight += l;
//            }

            // normalize weights
            double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
            _gene_tree_log_weight = log_weight_choices_sum;
            for (int b=0; b < (int) log_weight_choices.size(); b++) {
                log_weight_choices[b] -= log_weight_choices_sum;
            }

//            _gene_tree_log_likelihood = _gene_tree_log_weight;
            // randomly select a pair
            _index_of_choice = selectPair(log_weight_choices);
            
            _gene_tree_log_likelihood = _log_likelihood_choices[_index_of_choice];
//            _prev_log_likelihood_choices = _log_likelihood_choices;
        }
    }

    inline void Forest::mergeChosenPair(double species_tree_increment) { // TODO: make sure be able to account for deep  coalescnece
        if (_log_likelihood_choices.size() > 0) {
            _log_likelihood_choices.clear();
            // find nodes to join in node_list
            Node *subtree1 = _node_choices[_index_of_choice].first;
            Node *subtree2 = _node_choices[_index_of_choice].second;
            
            if (_increment_choices[_index_of_choice] < species_tree_increment) { // TODO: not sure this works for 3+ individuals
                // no deep coalescence, must rejoin chosen pair
                
                _node_choices.clear();

                _gene_tree_log_likelihood = _log_likelihood_choices[_index_of_choice];
                _log_likelihood_choices.clear();
                _new_nodes.clear();
                
                // new node is always needed
                Node nd;
                _nodes.push_back(nd);
                Node* new_nd = &_nodes.back();
                new_nd->_parent=0;
                new_nd->_number=_nleaves+_ninternals;
                new_nd->_edge_length=0.0;
                new_nd->_right_sib=0;
                _ninternals++;

                new_nd->_left_child=subtree1;
                subtree1->_right_sib=subtree2;

                subtree1->_parent=new_nd;
                subtree2->_parent=new_nd;
                
                subtree1->_edge_length += _increment_choices[_index_of_choice];
                subtree2->_edge_length += _increment_choices[_index_of_choice];
                
                _num_coalescent_events_in_generation = 1;

                //always calculating partials now
                assert (new_nd->_partial == nullptr);
                new_nd->_partial=ps.getPartial(_npatterns*4);
                assert(new_nd->_left_child->_right_sib);
                calcPartialArray(new_nd);

                //update species list
                updateNodeVector(_lineages, subtree1, subtree2, new_nd);
                
                    // update species partition
                for (auto &s:_species_partition) {
                    for (auto &nd:s.second) {
                        // TODO: there must be a better way to do this
                        if (nd == subtree1 || nd == subtree2) {
    //                    if (nd->_name == subtree1->_name || nd->_name == subtree2->_name) {
                            updateNodeList(s.second, subtree1, subtree2, new_nd);
                            break;
                        }
                    }
                }
                
//                for (auto &s:_species_partition) {
//                    cout << "x";
//                }
                
                // clear new node from _nodes
//                for (auto iter = _nodes.begin(); iter != _nodes.end(); iter++) {
//                    if (iter->_name == "unused") {
//                        iter = _nodes.erase(iter);
//                        --iter;
//                        _ninternals--;
//                    }
//                }
//
//                // reset node numbers
//                int n = 0;
//                for (auto &nd:_nodes) {
//                    nd._number = n;
//                    n++;
//                }
            }
            
            else {
                // deep coalescence, don't need to join anything
                // reset branch lengths for the chosen deep coalescence event
                subtree1->_edge_length = _increment_choices[_index_of_choice];
                subtree2->_edge_length = _increment_choices[_index_of_choice];
                _num_coalescent_events_in_generation = 0;
            }
            
            // save the correct increment from _increments vector
            pair<double, double> chosen_increment;
            
            for (auto &i:_increments) {
                if (i.first == _increment_choices[_index_of_choice]) {
                    chosen_increment = i;
                }
            }
            
            // delete all new increments from _increments vector
            // don't need to do this if there was only one choice
            if (_increment_choices.size() > 1) {
                assert(_increments.size() >= _increment_choices.size());
                _increments.resize(_increments.size() - _increment_choices.size());
                
                _increments.push_back(chosen_increment);
            }
        }
        
        _increment_choices.clear();
        _node_choices.clear();
    }

    inline int Forest::selectPair(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u = rng.uniform();
        _rand_numbers.push_back(u);
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
//        assert (likelihood_vec.size() == prev_log_likelihood_choices.size());
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
        _node_choices.push_back(s);

        return s;
    }

    inline tuple<Node*, Node*, Node*> Forest::createNewSubtree(pair<unsigned, unsigned> t, list<Node*> node_list, double increment) {
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

        //update species list
        updateNodeList(node_list, subtree1, subtree2, new_nd);
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

    inline void Forest::createDefaultTree() {
        clear();
        //create taxa
        double edge_length = rng.gamma(1.0, 1.0/_ntaxa);
        _rand_numbers.push_back(edge_length);
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
        _generationf = other._generationf;
        _gamma = other._gamma;
        _rand_numbers = other._rand_numbers;
        _gene_tree_log_weight = other._gene_tree_log_weight;
        _prev_log_likelihood = other._prev_log_likelihood;
        _log_weight_vec = other._log_weight_vec;
        _theta = other._theta;
        _increments = other._increments;
        _topology_prior = other._topology_prior;
        _num_coalescent_events_in_generation = other._num_coalescent_events_in_generation;
        _searchable_branch_lengths = other._searchable_branch_lengths;
        _gene_tree_cum_time = other._gene_tree_cum_time;
        _log_joining_prob = other._log_joining_prob;
        _increment_choices = other._increment_choices;
        _prev_log_likelihood_choices = other._prev_log_likelihood_choices;

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
    }

    inline void Forest::setUpSpeciesForest(vector<string> &species_names) {
        assert (_index==0);
        assert (_nspecies = (unsigned) species_names.size());
        clear();
        //create species
        double edge_length = 0.0;
        for (unsigned i = 0; i < _nspecies; i++) {
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=species_names[i];
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length = edge_length;
            nd->_position_in_lineages=i;
            }
        _nleaves=_nspecies;
        _ninternals=0;
        _nodes.resize(_nspecies);
        _lineages.resize(_nspecies);
    }

    inline string Forest::chooseEvent() {
        string event;
        // hybridization prior
        double rate = (_speciation_rate+_hybridization_rate)*_lineages.size();

        double hybridization_prob = _hybridization_rate/(_hybridization_rate+_speciation_rate);

        double u = rng.uniform();
        _rand_numbers.push_back(u);
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

        if (_lineages.size()>2) {
            _increments.push_back(make_pair(_last_edge_length, log(rate)-_last_edge_length*rate));
        }

        return event;
    }

    inline void Forest::chooseSpeciesIncrement() {
        // hybridization prior
        double rate = (_speciation_rate+_hybridization_rate)*_lineages.size();

        _last_edge_length = rng.gamma(1.0, 1.0/rate);

        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
        _increments.push_back(make_pair(_last_edge_length, log(rate)-_last_edge_length*rate));
    }

    inline tuple<string,string, string> Forest::speciesTreeProposal() {
        // this function creates a new node and joins two species

        pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size());
        Node *subtree1=_lineages[t.first];
        Node *subtree2=_lineages[t.second];
        assert(!subtree1->_parent && !subtree2->_parent);

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

        _species_joined = make_pair(subtree1, subtree2);

        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
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
        assert (_index >0);
        _species_partition.clear();
        for (auto nd:_lineages) {
            if (!nd->_left_child) {
                string species_name = taxon_map[nd->_name];
                _species_partition[species_name].push_back(nd);
            }
        }
        assert (_species_partition.size() > 0);
    }

    inline void Forest::evolveSpeciesFor(list<Node*> &nodes, double species_tree_increment) {
        bool done = false;
        double cum_time = _gene_tree_cum_time;
//        double cum_time = _gene_tree_cum_time = getLineageHeight(nodes.front()); // TODO: make sure all nodes in the list have the same height here
        bool coalescence = false;
        
        // for prior post ish proposal, do not add increment if not attempting a coalescent event
        if (_proposal == "prior-post-ish") {
            if (nodes.size() == 1) {
                done = true;
            }
        }

        while (!done) {
            string event = "coalescence";
            double increment;

            double s = nodes.size();
            double coalescence_rate = s*(s-1)/_theta;
            increment = rng.gamma(1.0, 1.0/(coalescence_rate));

            if (_migration_rate > 0.0) {
                increment = rng.gamma(1.0, 1.0/(coalescence_rate+_migration_rate));
                double migration_prob = _migration_rate/(_migration_rate+coalescence_rate);
                double coalescence_prob = coalescence_rate/(_migration_rate+coalescence_rate);

    //             draw random number to determine coalescent or migration event
                vector<double> prob_vec;
                prob_vec.push_back(coalescence_prob);
                prob_vec.push_back(migration_prob);

                double u = rng.uniform();
                _rand_numbers.push_back(u);

                double cum_prob = 0.0;
                for (int i = 0; i<2; i++) {
                    cum_prob += prob_vec[i];
                    if (u <= cum_prob) {
                        if (i == 0) {event = "coalescence";}
                        if (i == 1) {event = "migration";}
                        break;
                    }
                }
            }

            bool time_limited = species_tree_increment > 0.0;
            bool lineages_left_to_join = s > 1;
            bool time_limit_reached = cum_time+increment > species_tree_increment;
            if ((time_limited && time_limit_reached) || !lineages_left_to_join)  {
                done = true;
                increment = species_tree_increment - cum_time;
            }

            //if increment is 0, the lineage is done and we don't consider it at all
            if (increment > 0) {
                if (_proposal == "prior-post-ish") {
                    _increment_choices.push_back(increment);
                }
                for (auto nd:nodes) {
                    nd->_edge_length += increment; //add most recently chosen branch length to each node in lineage
                }
            if (_migration_rate == 0.0 || event == "coalescence") {
                if (!done) {
                    allowCoalescence(nodes, increment);
                    coalescence = true;
                    _increments.push_back(make_pair(increment, log(coalescence_rate)-increment*coalescence_rate));
                    if (_proposal != "prior-post-ish") {
                        _searchable_branch_lengths.push_back(make_pair (getLineageHeight(_new_nodes.back()), increment));
                    }
                    done = true;
                }
                cum_time += increment;
                if (!coalescence) {
                    allowDummyCoalescence(nodes, increment);
                    _increments.push_back(make_pair(increment, log(coalescence_rate)-increment*coalescence_rate));

                }
            }
            else if (event == "migration" && nodes.size() > 1) {
                // TODO: can this happen with only one individual in a lineage?
                // if there's only one individual in a lineage, don't allow migration b/c then there will be no individuals in the lineage, which we know is not true
                assert(nodes.size()>1);
                if (!done) {
                    allowMigration(nodes);
                }
            }
            }
        }
    }

    inline void Forest::allowDummyCoalescence(list<Node*> &nodes, double increment) {
        Node *subtree1 = nullptr;
        Node *subtree2 = nullptr;
        unsigned s = (unsigned) nodes.size();
        if (nodes.size()>2) {
                pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
                auto it1 = std::next(nodes.begin(), t.first);
                subtree1 = *it1;

                auto it2 = std::next(nodes.begin(), t.second);
                subtree2 = *it2;
                assert (t.first < nodes.size() && t.second < nodes.size());
        }
        else {
            // if there are only two lineages left, there is only one choice
            subtree1 = nodes.front();
            subtree2 = nodes.back();
        }

        assert (subtree1 != subtree2);

        //new node is always needed
//        Node nd; // TODO: not sure this is necessary
//        _nodes.push_back(nd);
//        Node* new_nd = &_nodes.back();
//        new_nd->_name = "unused";
//
//        new_nd->_parent=0;
//        new_nd->_number=_nleaves+_ninternals;
//        new_nd->_edge_length=0.0;
//        _ninternals++;
//        new_nd->_right_sib=0;

//        _new_nodes.push_back(new_nd);
//        calcPartialArray(new_nd);
        
        if (_proposal == "prior-post-ish") {
            _log_likelihood_choices.push_back(calcLogLikelihood());
//            new_nd->_name = "unused";
            _node_choices.push_back(make_pair(subtree1, subtree2));
            revertBranches(nodes, subtree1, subtree2, increment);
        }
    }

    inline void Forest::allowCoalescence(list<Node*> &nodes, double increment) {
        Node *subtree1 = nullptr;
        Node *subtree2 = nullptr;
        unsigned s = (unsigned) nodes.size();
        if (nodes.size()>2) {
            // prior-post proposal
            if (_proposal == "prior-post") {
                pair<Node*, Node*> t = chooseAllPairs(nodes, increment);

                subtree1 = t.first;
                subtree2 = t.second;
            }
            
            // prior-prior proposal
            else {
                pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
                auto it1 = std::next(nodes.begin(), t.first);
                subtree1 = *it1;

                auto it2 = std::next(nodes.begin(), t.second);
                subtree2 = *it2;
                assert (t.first < nodes.size() && t.second < nodes.size());
            }
        }
        else {
            // if there are only two lineages left, there is only one choice
            // prior-prior and prior-post proposals will return the same thing
            subtree1 = nodes.front();
            subtree2 = nodes.back();
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
        new_nd->_partial=ps.getPartial(_npatterns*4);
        assert(new_nd->_left_child->_right_sib);

        _new_nodes.push_back(new_nd);
        calcPartialArray(new_nd);
        _num_coalescent_events_in_generation++;

        if (_proposal != "prior-post-ish" ) {
            calcTopologyPrior((int) nodes.size()); // TODO: this won't work with prior-post-ish
        }

        //update species list
        updateNodeList(nodes, subtree1, subtree2, new_nd);
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        if (_proposal == "prior-post-ish") {
            _log_likelihood_choices.push_back(calcLogLikelihood());
            new_nd->_name = "unused";
            _node_choices.push_back(make_pair(subtree1, subtree2));
            revertNewNode(nodes, new_nd, subtree1, subtree2, increment);
        }
    }

    inline void Forest::fullyCoalesceGeneTree(list<Node*> &nodes) { // TODO: add prior post ish here too - I think this is okay because deep coalescence is never an issue here
//        bool coalescence = false;
        assert (nodes.size()>0);
        bool done = false;
        if (nodes.size() == 1) {
            done = true;
        }
        
        // for prior post-ish propposal, do not add increment if not attempting a coalescent event
        if (_proposal == "prior-post-ish") {
            if (nodes.size() == 1) {
                done = true;
            }
        }

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
                if (_proposal == "prior-post-ish") {
                    _increment_choices.push_back(increment);
                }
                
                for (auto nd:nodes) {
                    nd->_edge_length += increment;
                }
                Node* subtree1 = nullptr;
                Node *subtree2 = nullptr;

                _increments.push_back(make_pair(increment, log(coalescence_rate)-increment*coalescence_rate));

                if (nodes.size()>2) {

// prior-post proposal
                    if (_proposal == "prior-post") {
                        pair<Node*, Node*> t = chooseAllPairs(nodes, increment);

                        subtree1 = t.first;
                        subtree2 = t.second;
                    }
                    // prior prior or prior post ish
                    else {
                        pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
                        auto it1 = std::next(nodes.begin(), t.first);
                        subtree1 = *it1;

                        auto it2 = std::next(nodes.begin(), t.second);
                        subtree2 = *it2;
                    }
                }

                else {
                    // if there are only two lineages left, there is only one choice
                    // prior-prior and prior-post proposals will return the same thing
                    subtree1 = nodes.front();
                    subtree2 = nodes.back();
                }
                
                assert (subtree1 != subtree2); // TODO: make sure this is right for deep coalescence case too
                
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

                _new_nodes.push_back(new_nd);

                assert (new_nd->_partial == nullptr);
                new_nd->_partial=ps.getPartial(_npatterns*4);
                assert(new_nd->_left_child->_right_sib);
                calcPartialArray(new_nd);
                _num_coalescent_events_in_generation++;
//                coalescence = true;
                _searchable_branch_lengths.push_back(make_pair(getLineageHeight(_new_nodes.back()), increment));
                
                if (_proposal != "prior-post-ish" ) {
                    calcTopologyPrior((int) nodes.size()); // TODO: this won't work with prior-post-ish, must do it later
                }

                //update species list
                updateNodeList(nodes, subtree1, subtree2, new_nd);
                updateNodeVector(_lineages, subtree1, subtree2, new_nd);
                
                if (_proposal == "prior-post-ish") {
                    _prev_log_likelihood_choices.push_back(_gene_tree_log_likelihood);
                    _log_likelihood_choices.push_back(calcLogLikelihood());
                    new_nd->_name = "unused";
                    _node_choices.push_back(make_pair(subtree1, subtree2));
                    revertNewNode(nodes, new_nd, subtree1, subtree2, increment);
                    done = true;
                }
            }
        }
    }

    inline void Forest::firstGeneTreeProposal(double time_increment) {
        _prev_log_likelihood = _gene_tree_log_likelihood;
        _prev_log_likelihood_choices.clear();
        
        if (_species_partition.size() == 1) {
//            _prev_log_likelihood_choices.push_back(calcLogLikelihood());
            fullyCoalesceGeneTree(_species_partition.begin()->second); // TODO: make sure this function is correct now
        }

        else {
            for (auto &s:_species_partition) {
//                _prev_log_likelihood_choices.push_back(calcLogLikelihood());
                assert (s.second.size()>0);
                _gene_tree_cum_time = getLineageHeight(s.second.front());
                evolveSpeciesFor(s.second, time_increment);
            }
        }
    }

    inline void Forest::geneTreeProposal(tuple<string, string, string> &species_merge_info, double time_increment) {
        //update species partition
        string new_name = get<2>(species_merge_info);
        string species1 = get<0>(species_merge_info);
        string species2 = get<1>(species_merge_info);
        
        _prev_log_likelihood = _gene_tree_log_likelihood;
        _prev_log_likelihood_choices.clear();
        
        bool reset_gene_tree_cum_time = true;
        _gene_tree_cum_time = 0.0;

        if (new_name != "null" ){
            // skip this for generation 0, no species have been joined yet
            list<Node*> &nodes = _species_partition[new_name];
            copy(_species_partition[species1].begin(), _species_partition[species1].end(), back_inserter(nodes));
            copy(_species_partition[species2].begin(), _species_partition[species2].end(), back_inserter(nodes));
            _species_partition.erase(species1);
            _species_partition.erase(species2);
            reset_gene_tree_cum_time = false;
        }

        if (_species_partition.size() == 1) {
//            _prev_log_likelihood_choices.push_back(calcLogLikelihood());
            fullyCoalesceGeneTree(_species_partition.begin()->second);
        }

        else {
            for (auto &s:_species_partition) {
//                _prev_log_likelihood_choices.push_back(calcLogLikelihood());
                assert (s.second.size()>0);
                if (reset_gene_tree_cum_time) {
                    _gene_tree_cum_time = getLineageHeight(s.second.front()); // TODO: needs to be within this increment only
                }
                evolveSpeciesFor(s.second, time_increment);
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

        if (position1 > node_vector.size()+1) {
            position1 = position1-1;
            if (position2 > 0) {
                position2 = position2-1;
            }
        }
        if (position2 > node_vector.size()+1) {
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
        for (int i=0; i < (int) _lineages.size(); i++) {
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

    inline void Forest::revertBranches(list <Node*> &nodes, Node* subtree1, Node* subtree2, double increment) {
        // save increment that was added to to subtrees
        vector<double> increments;
        increments.push_back(subtree1->_edge_length);
        increments.push_back(subtree2->_edge_length);

        // clear new node from _nodes
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
        
        // revert edge lengths back
        for (auto &nd:nodes) {
            nd->_edge_length -= increment;
        }
    }

    inline void Forest::revertNewNode(list <Node*> &nodes, Node* new_nd, Node* subtree1, Node* subtree2, double increment) {
        // revert the species partition
        for (auto &s:_species_partition) {
            for (auto &nd:s.second) {
                if (nd->_name == "unused") {
                    revertNodeList(s.second, nd->_left_child, nd->_left_child->_right_sib, nd);
                    break;
                }
            }
        }

        // revert _lineages vector
        revertNodeVector(_lineages, subtree1, subtree2, new_nd);

        // reset siblings and parents of new node back to 0
        new_nd->_left_child->_right_sib->resetNode(); // subtree 2
        new_nd->_left_child->resetNode(); // subtree 1

        // clear new node from _new_nodes
//        _new_nodes.erase(remove(_new_nodes.begin(), _new_nodes.end(), new_nd));

        // clear new node from _nodes
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

        // check the reversion worked
        for (auto &s:_species_partition) {
            for (auto &n:s.second) {
                assert(_lineages[n->_position_in_lineages] == n);
            }
        }
        
        // revert edge lengths back
        for (auto &nd:nodes) {
            nd->_edge_length -= increment;
        }
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
        _rand_numbers.push_back(taxon_choice);
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
        _rand_numbers.push_back(lineage_choice);
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
            _rand_numbers.push_back(u);
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
            evolveSpeciesFor(s.second, species_tree_increment);
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
        for (auto &l:_lineages) {
            assert (!l->_parent);
            assert(!l->_right_sib);
        }
        }

    inline void Forest::hybridGeneTreeProposal(double species_tree_increment) {
        if (_species_partition.size() == 1) {
            fullyCoalesceGeneTree(_species_partition.begin()->second);
        }

        else {
            for (auto &s:_species_partition) {
                assert (s.second.size()>0);
                evolveSpeciesFor(s.second, species_tree_increment);
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
        vector<double> log_weight_choices = reweightChoices(_log_likelihood_choices, _prev_log_likelihood);

        // sum unnormalized weights before choosing the pair
        _gene_tree_log_weight = 0.0;
        for (auto &l:log_weight_choices) {
            _gene_tree_log_weight += l;
        }

        // normalize weights
        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        for (int b=0; b < (int) log_weight_choices.size(); b++) {
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
        double rate = (_speciation_rate)*_lineages.size();

        double u = rng.uniform();
        _rand_numbers.push_back(u);
        // choose edge length but don't add it yet
        _last_edge_length = rng.gamma(1.0, 1.0/rate);

        if (_lineages.size()>1) {
            _increments.push_back(make_pair(_last_edge_length, log(rate)-_last_edge_length*rate));
        }

        // add the previously chosen edge length
        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
    }

    inline double Forest::calcTopologyPrior(int nlineages) {
        _log_joining_prob += -log(0.5*nlineages*(nlineages-1));
        return _log_joining_prob;
    }

    inline double Forest::findShallowestCoalescence() {
        pair<double, double> smallest_branch = *min_element(_searchable_branch_lengths.begin(), _searchable_branch_lengths.end());
        return smallest_branch.first;
    }

    inline void Forest::revertToShallowest(double smallest_branch){
        // revert joins that are not smallest join

        // save nodes that must be reverted
        vector<Node*> nodes_to_revert;
        for (auto &lineage:_new_nodes) {
            if (lineage->_left_child) {
                if (getLineageHeight(lineage->_left_child) != smallest_branch) {
                    nodes_to_revert.push_back(lineage);
                    lineage->_name = "unused";
                }
                else {
                    lineage->_edge_length = 0.0;
                }
            }
        }

        assert (nodes_to_revert.size() == _new_nodes.size()-1);
        // TODO: this will only work if there is only one unused node per lineage (which is true for now)
        for (auto &s:_species_partition) {
            list<Node*> my_list;
            for (auto &nd:s.second) {
                if (nd->_name == "unused") {
                    my_list = s.second;
                    revertNodeList(my_list, nd->_left_child, nd->_left_child->_right_sib, nd);
                    s.second = my_list;
                    break;
                }
            }
        }

        // revert nodes in reverse
        for (auto &nd : boost::adaptors::reverse(nodes_to_revert)) {
            assert (nd->_name == "unused");
            revertNodeVector(_lineages, nd->_left_child, nd->_left_child->_right_sib, nd);

            //reset siblings and parents of original nodes back to 0
            nd->_left_child->_right_sib->resetNode(); //subtree2
            nd->_left_child->resetNode(); //subtree1


            // clear new node from _new_nodes
            _new_nodes.erase(remove(_new_nodes.begin(), _new_nodes.end(), nd));
        }

        for (auto iter = _nodes.begin(); iter != _nodes.end(); iter++) {
            if (iter->_name == "unused") {
                iter = _nodes.erase(iter);
                --iter;
                _ninternals--;
            }
        }

        // reset branches to smallest coalescence
        for (auto &lineage:_lineages) {
            if (lineage->_edge_length != 0.0) {
                if (!lineage->_left_child) {
                    lineage->_edge_length = smallest_branch;
                }
                else {
                    lineage->_edge_length = smallest_branch - getLineageHeight(lineage->_left_child);
                }
            }
        }

        // delete extra increments from _branch_lengths and _branch_length_priors
        vector<double> increments_to_remove;
        for (auto &s:_searchable_branch_lengths) {
            // s.first is lineage height
            // s.second is increment
            if (s.first != smallest_branch) {
                increments_to_remove.push_back(s.second);
            }
        }

        for (auto &i:_increments) {
            for (int a = 0; a<increments_to_remove.size(); a++) {
                if (i.first == increments_to_remove[a]) {
                    _increments.erase(remove(_increments.begin(), _increments.end(), i), _increments.end());
                }
            }
        }

        // reset _searchable_branch_lengths
        _searchable_branch_lengths.clear();

        // reset node numbers
        int n = 0;
        for (auto &nd:_nodes) {
            nd._number = n;
            n++;
        }

        // clear new nodes to reset for next generation
        _new_nodes.clear();

        for (auto &s:_species_partition) {
            for (auto &n:s.second) {
                assert(_lineages[n->_position_in_lineages] == n);
            }
        }

        // calculate gene tree height
        _gene_tree_cum_time = getLineageHeight(_lineages[0]);
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

    inline void Forest::extendGeneTreeLineages(double species_tree_height) {
        if (_lineages.size() > 1) {
            for (auto &l:_lineages) {
                if (l->_left_child) {
                    if (getLineageHeight(l) < species_tree_height) {
                        l->_edge_length += species_tree_height - getLineageHeight(l);
                    }
                }
                else {
                    if (l->_edge_length < species_tree_height) {
                        l->_edge_length += species_tree_height - getLineageHeight(l);
                    }
                }
            }
        }
        else {
            _lineages[0]->_left_child->_edge_length = species_tree_height - getLineageHeight(_lineages[0]->_left_child);
            _lineages[0]->_left_child->_right_sib->_edge_length = species_tree_height - getLineageHeight(_lineages[0]->_left_child->_right_sib);
        }
    }
}


