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
        void                        geneTreeProposal(tuple<string, string, string> &species_merge_info, double time_increment);
        void                        evolveSpeciesFor(list <Node*> &nodes, double time_increment);
        void                        fullyCoalesceGeneTree(list<Node*> &nodes);
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                        revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                      getRunningSumChoices(vector<double> &log_weight_choices);
        vector<double>              reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        pair<Node*, Node*>          chooseAllPairs(list<Node *> &node_list, double increment);
        pair<Node*, Node*>          getSubtreeAt(pair<unsigned, unsigned> t, list<Node*> node_list);
        int                         selectPair(vector<double> weight_vec);
        void                        chooseSpeciesIncrement();
        void                        allowMigration(list<Node*> &nodes);
        void                        setGeneration(double g) {_generation = g;}
        double                      chooseTaxonToMigrate(double s);
        string                      findKeyToDel(Node* taxon_to_migrate);
        void                        migrateTaxon(unsigned taxon_choice, string key_to_del, Node* taxon_to_migrate);
        void                        allowHybridization(list<Node*> &nodes);
        string                      chooseLineage(Node* taxon_to_migrate, string key_to_del);
        void                        addMigratingTaxon(string key_to_add, string key_to_del, Node* taxon_to_migrate);
        void                        deleteTaxon(string key_to_del, unsigned taxon_choice);

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
        vector<pair<Node*, Node*>> _node_choices;
        vector<double>              _log_likelihood_choices;
        double                      _prev_log_likelihood;
        int                         _index_of_choice;
        pair<Node*, Node*>          _species_joined;
        double                      _generation = 0;
        void                        showSpeciesJoined();
        double                      calcTransitionProbability(Node* child, double s, double s_child);
        void                        calculateNewEdgeLength(string key_to_add, Node* taxon_to_migrate);


    public:

        typedef std::shared_ptr<Forest> SharedPtr;
        static double               _theta;
        static double               _speciation_rate;
        static string               _proposal;
        static string               _model;
        static double               _kappa;
        static vector<double>       _base_frequencies;
        static string               _string_base_frequencies;
        static double               _migration_rate;
        double                      _hybridization_rate;
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
        _nodes.resize(2*_ntaxa);
        _npatterns = 0;
        _nstates = 4;

        //create taxa
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd=&_nodes[i];
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_position_in_lineages=i;
            }
        _nleaves=_ntaxa;
        _ninternals=0;

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

    inline void Forest::calcPartialArray(Node* new_nd) {
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

        // reset _log_likelihood_choices to 0 if we are in a new lineage
        if (_log_likelihood_choices.size() == 0) {
            _prev_log_likelihood = 0.0;
        }
        else {
            _prev_log_likelihood = _log_likelihood_choices[_index_of_choice];
        }
        _log_likelihood_choices.clear();

        //choose pair of nodes to try
        for (int i = 0; i < node_list.size()-1; i++) {
            for (int j = i+1; j < node_list.size(); j++) {
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

        // normalize weights
        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        for (int b=0; b<log_weight_choices.size(); b++) {
            log_weight_choices[b] -= log_weight_choices_sum;
        }

        // randomly select a pair
        _index_of_choice = selectPair(log_weight_choices);

        // find nodes to join in node_list
        Node *subtree1 = _node_choices[_index_of_choice].first;
        Node *subtree2 = _node_choices[_index_of_choice].second;
        return make_pair(subtree1, subtree2);
    }

    inline int Forest::selectPair(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u = rng.uniform();
        double cum_prob = 0.0;
        int index = 0.0;
        for (int i=0; i<weight_vec.size(); i++) {
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
        for (int a = 0; a<likelihood_vec.size(); a++) {
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
        Node* new_nd=&_nodes[_nleaves+_ninternals];

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
        for (unsigned i = 0; i < _ntaxa; i++) {
            Node* nd=&_nodes[i];
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
    }

    inline void Forest::operator=(const Forest & other) {
        _nstates = other._nstates;
        _npatterns = other._npatterns;
        _nodes.clear();
        _nodes.resize(other._nodes.size());
        _lineages.resize(other._lineages.size());
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
        _prev_log_likelihood = other._prev_log_likelihood;
        _species_joined = other._species_joined;

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
                if (othernd._parent) {
                    unsigned parent_number = othernd._parent->_number;
//                    _nodes[k]._partial = othernd._partial;

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
            _nodes[k]._position_in_lineages = othernd._position_in_lineages;
            _nodes[k]._partial = othernd._partial;
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
            nd->_right_sib=0;
            nd->_name=species_names[i];
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length = edge_length;
            _lineages.push_back(nd);
            nd->_position_in_lineages=i;
            }
        _nleaves=_nspecies;
        _ninternals=0;
        _nodes.resize(_nspecies*2);
        _lineages.resize(_nspecies);
    }

    inline void Forest::chooseSpeciesIncrement() {
        double rate = _speciation_rate*_lineages.size();
        _last_edge_length = rng.gamma(1.0, 1.0/rate);

        for (auto nd:_lineages) {
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each species node
        }
    }

    inline tuple<string,string, string> Forest::speciesTreeProposal() {
        chooseSpeciesIncrement();

        pair<unsigned, unsigned> t = chooseTaxaToJoin(_lineages.size());
        Node *subtree1=_lineages[t.first];
        Node *subtree2=_lineages[t.second];
        assert(!subtree1->_parent && !subtree2->_parent);
        assert(!subtree1->_right_sib && !subtree2->_right_sib);

        Node* new_nd = &_nodes[_nleaves+_ninternals];
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
//
//        cout << "species tree increment: " << "\t" << "\t" << "\t" << _last_edge_length << endl;
//        cout << "\t" << "species joined are " << subtree1->_name << " and " << subtree2->_name << endl;

        _species_joined = make_pair(subtree1, subtree2);

        return make_tuple(subtree1->_name, subtree2->_name, new_nd->_name);
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
        assert (_index >0);
        _species_partition.clear();
        for (auto nd:_lineages) {
            if (!nd->_left_child) {
                string species_name = taxon_map[nd->_name];
                _species_partition[species_name].push_back(nd);
            }
        }
    }

    inline void Forest::evolveSpeciesFor(list<Node*> &nodes, double species_tree_increment) {
        // reset _log_likelihood_choices and _prev_log_likelihood each time a new lineage is entered
        _log_likelihood_choices.clear();
        _prev_log_likelihood = 0.0;
        bool done = false;
        double cum_time = 0.0;
        bool migration;
        if (_migration_rate > 0) {migration = true;}

        while (!done) {
            string event;
            double increment;

            double s = nodes.size();
            double coalescence_rate = s*(s-1)/_theta;
            _hybridization_rate = 0.1;
            increment = rng.gamma(1.0, 1.0/(coalescence_rate + _hybridization_rate));

            // _hybridization_rate chance of of staying in lineage
            // (1-_hybridization_rate) chance of moving to another lineage
            double r = rng.uniform();
            if (r < _hybridization_rate) {event = "hybridization";}
            else if (r > _hybridization_rate) {event = "coalescence";}

# if migration == true
            increment = rng.gamma(1.0, 1.0/(coalescence_rate+_migration_rate));
            double migration_prob = _migration_rate/(_migration_rate+coalescence_rate);
            double coalescence_prob = coalescence_rate/(_migration_rate+coalescence_rate);

//             draw random number to determine coalescent or migration event
            vector<double> prob_vec;
            prob_vec.push_back(coalescence_prob);
            prob_vec.push_back(migration_prob);

            double u = rng.uniform();

            double cum_prob = 0.0;
            for (int i = 0; i<2; i++) {
                cum_prob += prob_vec[i];
                if (u <= cum_prob) {
                    if (i == 0) {event = "coalescence";}
                    if (i == 1) {event = "migration";}
                    break;
                }
            }
# endif
            if (event == "coalescence" || nodes.size() == 1) {
                bool time_limited = species_tree_increment > 0.0;
                bool lineages_left_to_join = s > 1;
                bool time_limit_reached = cum_time+increment > species_tree_increment;
                if ((time_limited && time_limit_reached) || !lineages_left_to_join)  {
                    done = true;
                    increment = species_tree_increment - cum_time;
                }

                //add increment to each lineage
                for (auto nd:nodes) {
                    nd->_edge_length += increment; //add most recently chosen branch length to each node in lineage
                }

                if (!done) {
                    Node *subtree1;
                    Node *subtree2;

                    if (nodes.size()>2) {
                        // prior-prior proposal
                        if (_proposal == "prior-prior") {
                            pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
                            auto it1 = std::next(nodes.begin(), t.first);
                            subtree1 = *it1;

                            auto it2 = std::next(nodes.begin(), t.second);
                            subtree2 = *it2;
                            assert (t.first < nodes.size() && t.second < nodes.size());
                        }

                        // prior-post proposal
                        if (_proposal == "prior-post") {
                            pair<Node*, Node*> t = chooseAllPairs(nodes, increment);

                            subtree1 = t.first;
                            subtree2 = t.second;
                        }
                    }
                    else {
                        // if there are only two lineages left, there is only one choice
                        // prior-prior and prior-post proposals will return the same thing
                        subtree1 = nodes.front();
                        subtree2 = nodes.back();
                    }

                    //new node is always needed
                    Node* new_nd=&_nodes[_nleaves+_ninternals];

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
                    calcPartialArray(new_nd);

                    //update species list
                    updateNodeList(nodes, subtree1, subtree2, new_nd);
                    updateNodeVector(_lineages, subtree1, subtree2, new_nd);
                }
                cum_time += increment;
            }
            else if (event == "migration" && nodes.size()>1) {
                // TODO: not sure about this
                // if there's only one individual in a lineage, don't allow migration b/c then there will be no individuals in the lineage, which we know is not true
                assert(nodes.size()>1);
                allowMigration(nodes);
            }
            else if (event == "hybridization") {
                allowHybridization(nodes);
            }
        }
    }

    inline void Forest::fullyCoalesceGeneTree(list<Node*> &nodes) {
        _prev_log_likelihood = 0.0;
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
                    nd->_edge_length += increment;
                }
                Node* subtree1;
                Node *subtree2;

                if (nodes.size()>2) {

// prior-prior proposal
                    if (_proposal == "prior-prior") {
                        pair<unsigned, unsigned> t = chooseTaxaToJoin(s);
                        auto it1 = std::next(nodes.begin(), t.first);
                        subtree1 = *it1;

                        auto it2 = std::next(nodes.begin(), t.second);
                        subtree2 = *it2;
                    }

// prior-post proposal
                    if (_proposal == "prior-post") {
                        pair<Node*, Node*> t = chooseAllPairs(nodes, increment);

                        subtree1 = t.first;
                        subtree2 = t.second;
                    }
                }

                else {
                    // if there are only two lineages left, there is only one choice
                    // prior-prior and prior-post proposals will return the same thing
                    subtree1 = nodes.front();
                    subtree2 = nodes.back();
                }


                Node* new_nd=&_nodes[_nleaves+_ninternals];

                new_nd->_parent=0;
                new_nd->_number=_nleaves+_ninternals;
                new_nd->_edge_length=0.0;
                _ninternals++;
                new_nd->_right_sib=0;

                new_nd->_left_child=subtree1;
                subtree1->_right_sib=subtree2;

                subtree1->_parent=new_nd;
                subtree2->_parent=new_nd;

                assert (new_nd->_partial == nullptr);
                new_nd->_partial=ps.getPartial(_npatterns*4);
                assert(new_nd->_left_child->_right_sib);
                calcPartialArray(new_nd);

                //update species list
                updateNodeList(nodes, subtree1, subtree2, new_nd);
                updateNodeVector(_lineages, subtree1, subtree2, new_nd);
            }
        }
    }

    inline void Forest::geneTreeProposal(tuple<string, string, string> &species_merge_info, double time_increment) {
        //update species partition
        string new_name = get<2>(species_merge_info);
        string species1 = get<0>(species_merge_info);
        string species2 = get<1>(species_merge_info);

        if (new_name != "null" ){
            // skip this for generation 0, no species have been joined yet
            list<Node*> &nodes = _species_partition[new_name];
            copy(_species_partition[species1].begin(), _species_partition[species1].end(), back_inserter(nodes));
            copy(_species_partition[species2].begin(), _species_partition[species2].end(), back_inserter(nodes));
            _species_partition.erase(species1);
            _species_partition.erase(species2);
        }

        if (_species_partition.size() == 1) {
            fullyCoalesceGeneTree(_species_partition.begin()->second);
        }
        else {
            for (auto & s:_species_partition) {
                assert (s.second.size()>0);
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
        for (int i=0; i<_lineages.size(); i++) {
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
        for (int i=0; i<_lineages.size(); i++) {
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

        // reset _position_in_lineages
        for (int i=0; i<_lineages.size(); i++) {
            _lineages[i] -> _position_in_lineages=i;
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
        mtx.unlock();
        return taxon_choice;
    }

    inline string Forest::findKeyToDel(Node* taxon_to_migrate) {
        string key_to_del;

        // find taxon to migrate in species partition
        // then key_to_del is s.first

        bool done = false;
        while (!done) {
        for (auto &s:_species_partition) {
            for (auto &t:s.second) {
                if (t == taxon_to_migrate) {
                    key_to_del = s.first;
                    done = true;
                    }
                }
            }
        }
        return key_to_del;
    }

    inline void Forest::migrateTaxon(unsigned taxon_choice, string key_to_del, Node* taxon_to_migrate) {
        // create vector of species names
        vector<string> species_names;
//        for (auto & s:_species_partition) {
//            species_names.push_back(s.first);
//        }
        
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
        
        // set edge length of migrating taxon or target edge length (lineage taxon is migrating into)
        calculateNewEdgeLength(key_to_add, taxon_to_migrate);
        
        // set edge length of migrating taxon or target edge length (lineage taxon is migrating into)
        
//        // find target lineage
//        Node* taxon_in_target_lineage;
//        double target_edge_length = 0.0;
//        for (auto &s:_species_partition) {
//            if (s.first == key_to_add) {
//                // TODO: check that this works if migrating taxon is the only taxon in the lineage
//                taxon_in_target_lineage = s.second.front();
//            }
//        }
        
//        // get target lineage edge length
//        for (Node* child = taxon_in_target_lineage; child; child=child->_left_child) {
//            target_edge_length += child->_edge_length;
//        }
        
        // get edge length of migrating lineage
//        double migrating_edge_length = 0.0;
//        for (Node* child = taxon_to_migrate; child; child = child->_left_child) {
//            migrating_edge_length += child->_edge_length;
//        }
//
//        // compare target edge length to migrating edge length
//        double difference = target_edge_length - migrating_edge_length;
//
//        // if migrating lineage is shorter than target lineage, extend migrating taxon edge length
//        if (difference > 0.0) {
//            taxon_to_migrate->_edge_length = difference;
//        }
//
//        // if migrating lineage is longer than target lineage, extend target edge length
//        // TODO: can this be simplified?
//        else if (difference < 0.0) {
//            for (auto &s:_species_partition) {
//                if (s.first == key_to_add) {
//                    for (auto iter = s.second.begin(); iter != s.second.end(); iter++) {
//                    Node* taxon = *iter;
//                    taxon->_edge_length += -1.0*difference;
//                    }
//                }
//            }
//        }
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
        double lineage_choice = ::rng.randint(0, _species_partition.size()-1);
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

    inline void Forest::calculateNewEdgeLength(string key_to_add, Node* taxon_to_migrate) {
        // find target lineage
        Node* taxon_in_target_lineage;
        double target_edge_length = 0.0;
        for (auto &s:_species_partition) {
            if (s.first == key_to_add) {
                // TODO: check that this works if migrating taxon is the only taxon in the lineage
                taxon_in_target_lineage = s.second.front();
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
        
        // if migrating lineage is shorter than target lineage, extend migrating taxon edge length
        if (difference > 0.0) {
            taxon_to_migrate->_edge_length = difference;
        }
        
        // if migrating lineage is longer than target lineage, extend target edge length
        // TODO: can this be simplified?
        else if (difference < 0.0) {
            for (auto &s:_species_partition) {
                if (s.first == key_to_add) {
                    for (auto iter = s.second.begin(); iter != s.second.end(); iter++) {
                    Node* taxon = *iter;
                    taxon->_edge_length += -1.0*difference;
                    }
                }
            }
        }
    }

    inline void Forest::allowHybridization(list<Node*> &nodes) {
        unsigned taxon_choice;
        Node* taxon_to_hybridize;
        taxon_choice = chooseTaxonToMigrate(nodes.size());
        auto iter = std::next(nodes.begin(), taxon_choice);
        taxon_to_hybridize = *iter;

        // find taxa to delete and lineage to add to
        string key_to_del = findKeyToDel(taxon_to_hybridize);

        // delete hybridizing taxon from its original lineage and add to new lineage
        migrateTaxon(taxon_choice, key_to_del, taxon_to_hybridize);
    }
}
