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
std::mutex mtx;

extern proj::StopWatch stopwatch;

#include "partial_store.hpp"
extern proj::PartialStore ps;

#include "node.hpp"

namespace proj {

//#if defined(LAZY_COPYING)
//    class ForestExtension;
//#endif

using namespace std;

//class Likelihood;
class Particle;

class Forest {

//        friend class Likelihood;
        friend class Particle;
    
#if defined(LAZY_COPYING)
        friend class ForestExtension;
#endif

    public:
            typedef std::shared_ptr<const Forest> ConstSharedPtr;

                                    Forest();
                                    ~Forest();
                                    Forest(const Forest & other);
    
    typedef tuple<double, unsigned, vector<G::species_t> >  coalinfo_t;

        unsigned                        numInternals() const;
        unsigned                        numNodes() const;
        void                            showForest();
        double                          calcLogLikelihood();
        void                            createDefaultTree(Lot::SharedPtr lot);
        void operator=(const Forest & other);
        void                            debugForest();
        void                            debugLogLikelihood(Node* nd, double log_like);
        double                          getForestHeight() {return _forest_height;}

    private:

        void                            clear();
        void                            setData(Data::SharedPtr d, int index, map<string, string> &taxon_map, bool partials);
        void                            setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa);
        Node *                          findNextPreorder(Node * nd);
        string                          makeNewick(unsigned precision, bool use_names);
        string                          makeAltNewick(unsigned precision, bool use_names);
        string                          makePartialNewick(unsigned precision, bool use_names);
        pair<unsigned, unsigned>        chooseTaxaToJoin(double s, Lot::SharedPtr lot);
#if defined (LAZY_COPYING)
        double                          calcPartialArrayJC(Node * new_nd, const Node * lchild, const Node * rchild) const;
#else
        void                            calcPartialArrayJC(Node* new_nd);
        void                            calcPartialArrayHKY(Node* new_nd);
#endif

        void                            setUpGeneForest(map<string, string> &taxon_map);
        void                            updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                            updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void                            revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        double                          getRunningSumChoices(vector<double> &log_weight_choices);
        vector<double>                  reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        int                             selectPair(vector<double> weight_vec, Lot::SharedPtr lot);
        void                            allowCoalescence(string species_name, double increment, Lot::SharedPtr lot);
    
#if defined (LAZY_COPYING)
        vector<pair<double, unsigned long>>    calcForestRate(Lot::SharedPtr lot, unordered_map<G::species_t, double> theta_map);
        void                                   coalesce(G::species_t species_name, Lot::SharedPtr lot, double prev_log_likelihood);
        void                                   coalescePreviousJoin(pair<Node*, Node*> subtrees);
//        G::species_t                           getPrevSpecies();
        vector<G::species_t>                   getAllPrevSpecies();

#else
        vector<pair<double, string>>    calcForestRate(Lot::SharedPtr lot, unordered_map<string, double> theta_map);
#endif
    
#if defined (LAZY_COPYING)
        vector<pair<double, string>>    calcForestRateSim(Lot::SharedPtr lot, unordered_map<G::species_t, double> theta_map);
#else
        vector<pair<double, string>>    calcForestRateSim(Lot::SharedPtr lot, unordered_map<string, double> theta_map);
#endif
        void                            updateSpeciesPartition(tuple<string, string, string> species_info);
    
#if defined (LAZY_COPYING)
        void                            updateSpeciesPartitionSim(tuple<string, string, string> species_info, G::species_t new_species);
        void                            mergeSpecies(G::species_t left_spp, G::species_t right_spp);
#endif
    
        double                          calcTopologyPrior(unsigned nlineages);
    
        void                            clearPartials();
        double                          getLogLikelihood() const;
        unsigned                        getDeepCoal(tuple <string, string, string> species_joined);
        unsigned                        getMaxDeepCoal(tuple <string, string, string> species_joined);
        void                            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
        vector< pair<double, Node *>>   sortPreorder();
        void                            refreshPreorder();
        void                            buildFromNewick(const string newick, bool rooted, bool allow_polytomies);
        void                            extractNodeNumberFromName(Node * nd, std::set<unsigned> & used);
        void                            stripOutNexusComments(std::string & newick);
        unsigned                        countNewickLeaves(const std::string newick);
        unsigned                        countNewickInternals(const std::string newick);
        void                            extractEdgeLen(Node * nd, std::string edge_length_string);
        void                            renumberInternals();
        bool                            canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);
        vector<tuple<string, string, string>>              buildFromNewickTopology(const string newick);
        pair<Node*, Node*>              chooseAllPairs(vector<Node*> &nodes, double increment, string species, Lot::SharedPtr lot);
        tuple<Node*, Node*, Node*>      createNewSubtree(pair<unsigned, unsigned> t, vector<Node*> node_list, double increment, string species);
        pair<Node*, Node*>              getSubtreeAt(pair<unsigned, unsigned> t, vector<Node*> node_list);
        void                            debugShowDistanceMatrix(const vector<double> & d) const;
    
#if defined (UPGMA)
        void                            constructUPGMA();
        void                            destroyUPGMA();
        void                            mergeDMatrixPair(vector<Split> & dmatrows, vector<double> & dmatrix, Split & s1, Split & s2);
        void                            buildStartingUPGMAMatrix();
#endif
    
#if defined (USING_MPI)
        vector<pair<tuple<string, string, string>, double>>                            resetLineages(Lot::SharedPtr lot);
#endif
        vector<pair<tuple<string, string, string>, double>> resetT();
    
#if defined (REUSE_PARTIALS)
        void                            stowPartial(Node *nd);
#endif
    
        void                            saveCoalInfoInitial();
        void                            saveCoalInfoGeneForest(vector<Forest::coalinfo_t> & coalinfo_vect) const;
        void                            addCoalInfoElem(const Node *, vector<coalinfo_t> & recipient);
        void                            fixupCoalInfo(vector<coalinfo_t> & coalinfo_vect, vector<coalinfo_t> & sppinfo_vect) const;
        static bool                     subsumed(G::species_t test_species, G::species_t subtending_species);
        void                            refreshPreorderNew(vector<Node*> & preorder) const;
        Node *                          findNextPreorderNew(Node * nd) const;
        void                            setTreeHeight();
    

        vector<coalinfo_t>                  _coalinfo;
        mutable unsigned                    _next_node_number;

        void                            setNodeHeights();
#if defined (USING_MPI)
        void                            resetSpeciesPartition(string species_partition);
        map<string, vector<string>>     saveSpeciesPartition();
    vector<pair<tuple<string, string, string>, double>> buildFromNewickMPI(const string newick, bool rooted, bool allow_polytomies, Lot::SharedPtr lot);
#endif
    
#if defined (UPGMA)
        void                            setGeneUPGMAMatrices();
#endif
    
        static bool                     compareNodeHeights(const Node& a, const Node& b) ;
    
#if defined(LAZY_COPYING)
        PartialStore::partial_t         pullPartial();
        void                            buildCoalInfoVect();
        virtual pair<double, double>    getBoundaryExtension() const;
        void                            joinLineagePair(Node * new_nd, Node * subtree1, Node * subtree2);
        Node *                          getMostRecentAncestor();
        void                            refreshAllPreorders() const;
        void                            advanceAllLineagesBy(double increment);
        void                            addIncrAndJoin(double incr, const Split & lsplit, const Split & rsplit, ForestExtension & gfx);
        void                            copyLineageSpecies(vector<G::species_t> & species_of_lineages) const;
        void                            debugCheckBleedingEdge(string msg, double anc_height) const;
        void                            mergeSpecies(G::species_t left_species, G::species_t right_species, G::species_t anc_species);
        void                            setLogLikelihood(double log_likelihood) {_gene_tree_log_likelihood = log_likelihood;}
        unsigned                        checkNumberOfUniqueSpeciesInExistence();
        pair<Node*, Node*>              getPrevJoin();
        void                            revertGeneForest(double prev_gene_increment, vector<G::species_t> prev_species_assignments);
        unsigned                        getNRemainingSpecies();
    
        mutable vector<Node::ptr_vect_t> _preorders;
#endif
    
        std::vector<Node *>             _lineages;
        vector<Node>                    _nodes;

        unsigned                        _ninternals;
        unsigned                        _npatterns;
        double                          _last_edge_length;

        Data::SharedPtr                 _data;
    
        unsigned                        _first_pattern = 0;
        unsigned                        _index; // gene indices start at 1, not 0
        map<string, vector<Node*> >     _species_partition;
        double                          _gene_tree_log_likelihood;
        double                          _log_joining_prob;
        vector<pair<double, double>>    _increments_and_priors;
    
#if defined (DEBUG_MODE)
        pair<Node*, Node*>                  _species_joined;
#endif

        vector<Node*>                   _preorder;
    
        double                          _log_weight;
        vector<pair<Node*, Node*>>      _node_choices;
        vector<double>                  _log_likelihood_choices;
    
#if defined (UPGMA)
        stack<Node *>                   _upgma_additions;
        map<Node *, double>             _upgma_starting_edgelen;
        
        // Local copy of the precalculated global JC pairwise distance matrix
        vector<double>                  _dmatrix;
        vector<Split>                   _dmatrix_rows;
    
        vector<double>                  _starting_dij;
#endif
    
        vector<pair<string, unsigned>>  _lineages_per_species;
        double                          _forest_height;
        double                          _forest_length;
    
#if defined (DEBUG_MODE)
        void                            showSpeciesJoined();
#endif
        double                          calcTransitionProbabilityJC(double s, double s_child, double edge_length) const;
        double                          calcTransitionProbabilityHKY(double s, double s_child, double edge_length);
        double                          calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length);
        double                          getTreeLength();
        void                            calcTreeLength();
        double                          getLineageHeight(Node* nd);
        void                            addIncrement(double increment);
        void                            simulateData(Lot::SharedPtr lot, unsigned starting_site, unsigned nsites);
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
        _ninternals=0;
        _preorder.clear();
        _forest_length = 0.0;
      
#if defined (UPGMA)
        _upgma_additions = stack<Node*>();
        _upgma_starting_edgelen.clear();

        _starting_dij.clear();
#endif
        
        _lineages_per_species.clear();
        _forest_height = 0.0;
        _coalinfo.clear();
        _next_node_number = 0;
        
#if defined (LAZY_COPYING)
        _preorders.clear();
#endif
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

#if defined (UPGMA)
    inline void Forest::setGeneUPGMAMatrices() {
        assert (_index > 0);
        _dmatrix      = G::_dmatrix[_index-1];
        _dmatrix_rows = G::_dmatrix_rows;
    }
#endif

    inline void Forest::setSimData(Data::SharedPtr d, int index, map<string, string> &taxon_map, unsigned ntaxa) {
        _index = index;
        assert (index > 0);         //don't set data for species tree, though it doesn't really matter for simulations
        G::_ntaxa = ntaxa;
        
        _data = d;
        
        _nodes.resize(2*G::_ntaxa - 1);
//        _nodes.reserve(2*G::_ntaxa - 1);
        _lineages.reserve(_nodes.size());
        unsigned i= 0;
        
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
//#if defined (LAZY_COPYING)
//            string taxon_name = G::_taxon_names[i];
//#endif
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
            
#if defined (LAZY_COPYING)
//            _nodes[i]._name = taxon_name;
            _nodes[i]._split.resize(G::_ntaxa);
            _nodes[i]._split.setBitAt(i);
//            if (G::_taxon_to_species.count(taxon_name) == 0)
//                throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % taxon_name));
//            else {
//                Node::setSpeciesBit(_nodes[i]._species, G::_taxon_to_species.at(taxon_name), /*init_to_zero_first*/true);
//            }
#endif
        }
        
        vector<string> taxon_names;
        for (auto &t:taxon_map) {
            taxon_names.push_back(t.first);
        }
        
//#if defined (LAZY_COPYING)
//        G::species_t curr_species = 1;
//#endif
        
        for (auto &nd:_lineages) {
            if (!nd->_left_child) {
                // replace all spaces with underscores so that other programs do not have
                  // trouble parsing your tree descriptions
                  std::string name = taxon_names[i++];
                  boost::replace_all(name, " ", "_");
                nd->_name = name;
//#if defined (LAZY_COPYING)
//                nd->_species = curr_species;
//#endif
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
        
#if !defined (LAZY_COPYING)
        const Data::taxon_names_t & taxon_names = _data->getTaxonNames();
        unsigned i = 0;
#endif
        auto &data_matrix=_data->getDataMatrix();
        
        _nodes.resize(2*G::_ntaxa - 1);
        _lineages.reserve(_nodes.size());
        
//        _nodes.reserve(2*G::_ntaxa - 1);
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            string taxon_name = G::_taxon_names[i];
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
            _lineages.push_back(nd);
#if defined (LAZY_COPYING)
            _nodes[i]._name = taxon_name;
            _nodes[i]._split.resize(G::_ntaxa);
            _nodes[i]._split.setBitAt(i);
            if (G::_taxon_to_species.count(taxon_name) == 0)
                throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % taxon_name));
            else {
                Node::setSpeciesBit(_nodes[i]._species, G::_taxon_to_species.at(taxon_name), /*init_to_zero_first*/true);
            }
#endif
        }
        
        if (G::_upgma) {
            for (unsigned i = 0; i < G::_ntaxa; i++) {
                Node * nd = &(_nodes[i]);
                // set splits
                nd->_split.resize(G::_ntaxa);
                nd->_split.setBitAt(i);
            }
        }

        for (auto &nd:_lineages) {
            if (!nd->_left_child) {
#if !defined (LAZY_COPYING)
                // replace all spaces with underscores so that other programs do not have
                  // trouble parsing your tree descriptions
                  std::string name = taxon_names[i++];
                  boost::replace_all(name, " ", "_");
                nd->_name = name;
#endif

                if (!G::_save_memory || (G::_save_memory && partials)) { // if save memory setting, don't set tip partials yet
#if defined (REUSE_PARTIALS)
                    mtx.lock();
                    nd->_partial=ps.getPartial(_npatterns*4, _index);
                    mtx.unlock();
#else
                    nd->_partial=ps.getPartial(_npatterns*4);
#endif
                    for (unsigned p=0; p<_npatterns; p++) {
                        
                        unsigned pp = _first_pattern+p;
                        unsigned pxnstates = p*G::_nstates;
#if defined (UNROLL_LOOPS)
                        
                        // loop 1
                        unsigned s = 0;
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[nd->_number][pp];
                        double result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 2
                        s = 1;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 3
                        s = 2;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 4
                        s = 3;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
#else
                        
#if defined (UNROLL_LOOPS)
                        // loop 1
                        unsigned s = 0;
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[nd->_number][pp];
                        double result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 2
                        s = 1;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 3
                        s = 2;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 4
                        s = 3;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
#else
#if defined (UNROLL_LOOPS)
                        assert (1 == 2);
                        // loop 1
                        
                        unsigned s = 0;
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[nd->_number][pp];
                        double result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 2
                        s = 1;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 3
                        s = 2;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        
                        // loop 4
                        s = 3;
                        state = (Data::state_t)1 << s;
                        d = data_matrix[nd->_number][pp];
                        result = state & d;
#if defined (REUSE_PARTIALS)
                        (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        

#else
                        for (unsigned s=0; s<G::_nstates; s++) {
                            Data::state_t state = (Data::state_t)1 << s;
                            Data::state_t d = data_matrix[nd->_number][pp];
                            double result = state & d;
#if defined (REUSE_PARTIALS)
                            (nd->_partial->_v)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                            (*nd->_partial)[pxnstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                        }
#endif
#endif
#endif
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
        new_nd->_number=G::_ntaxa+_ninternals;
        new_nd->_name+=boost::str(boost::format("node-%d")%new_nd->_number);
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
        
        calcPartialArray(new_nd);
        
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
        new_nd2->_number=G::_ntaxa+_ninternals;
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
        new_nd3->_number=G::_ntaxa+_ninternals;
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
        new_nd4->_number=G::_ntaxa+_ninternals;
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
        
        calcPartialArray(new_nd4);
        
        calcLogLikelihood();
#endif
        
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

#if defined (LAZY_COPYING)
    inline double Forest::calcPartialArrayJC(Node * new_nd, const Node * lchild, const Node * rchild) const {
        // Computes the partial array for new_nd and returns the difference in
        // log likelihood due to the addition of new_nd
        //char base[] = {'A','C','G','T'};
        
        // Get pattern counts
        auto counts = _data->getPatternCounts();

        // Get the first and last pattern index for this gene's data
        Data::begin_end_pair_t be = _data->getSubsetBeginEnd(_index - 1);
        unsigned first_pattern = be.first;
                
        auto & parent_partial_array = new_nd->_partial->_v;
        unsigned npatterns = _data->getNumPatternsInSubset(_index - 1);
    #if 1
        // Determine if there is an edge length extension (this would be the
        // case if new_nd comes from a gene forest extension)
        double lchild_stem_height = lchild->_height + lchild->_edge_length;
        double rchild_stem_height = rchild->_height + rchild->_edge_length;
        assert(fabs(lchild_stem_height - rchild_stem_height) < G::_small_enough);
        
        // Calculate the edge length extension
        double edgelen_extension = new_nd->_height - lchild_stem_height;
        
        // Edge length extension may be slightly negative due to roundoff
        assert(edgelen_extension >= -G::_small_enough);
        if (edgelen_extension < 0.0)
            edgelen_extension = 0.0;
            
        for (const Node * child : {lchild, rchild})  {
            assert(child->_partial);
            auto & child_partial_array = child->_partial->_v;
                
            double pr_same = calcTransitionProbabilityJC(0, 0, child->_edge_length + edgelen_extension);
            double pr_diff = calcTransitionProbabilityJC(0, 1, child->_edge_length + edgelen_extension);
            for (unsigned p = 0; p < npatterns; p++) {
                unsigned pxnstates = p*G::_nstates;
                //unsigned pp = first_pattern + p;

#if defined (UNROLL_LOOPS)
                // unroll parent loop
                assert (G::_nstates == 4);
                unsigned s = 0;
                double sum_over_child_states = 0.0;
                    
                    // child state subloop 0
                unsigned s_child = 0;
                double child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_same * child_partial;
                    
                    // child state subloop 1
                s_child = 1;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;

                    
                    // child state subloop 2
                s_child = 2;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 3
                s_child = 3;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;

                parent_partial_array[pxnstates + s] *= sum_over_child_states;
                
                s = 1;
                sum_over_child_states = 0.0;
                // child state subloop 0
                s_child = 0;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 1
                s_child = 1;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_same * child_partial;
                    
                    // child state subloop 2
                s_child = 2;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 3
                s_child = 3;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;

                parent_partial_array[pxnstates + s] *= sum_over_child_states;
                
                s = 2;
                sum_over_child_states = 0.0;
                // child state subloop 0
                s_child = 0;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 1
                s_child = 1;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 2
                s_child = 2;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_same * child_partial;
                
                    // child state subloop 3
                s_child = 3;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;

                parent_partial_array[pxnstates + s] *= sum_over_child_states;
                
                s = 3;
                sum_over_child_states = 0.0;
                // child state subloop 0
                s_child = 0;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 1
                s_child = 1;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 2
                s_child = 2;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_diff * child_partial;
                    
                    // child state subloop 3
                s_child = 3;
                child_partial = child_partial_array[pxnstates + s_child];

                sum_over_child_states += pr_same * child_partial;

                parent_partial_array[pxnstates + s] *= sum_over_child_states;
#else
                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? pr_same : pr_diff);
                        double child_partial = child_partial_array[pxnstates + s_child];
                                                
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    
                    parent_partial_array[pxnstates + s] *= sum_over_child_states;
                }   // parent state loop
#endif
            }   // pattern loop
        }
    #else
        for (Node * child = new_nd->_left_child; child; child = child->_right_sib) {
            assert(child->_partial);
            auto & child_partial_array = child->_partial->_v;

            // If this gene forest is an extension, check to see if edge length
            // needs to be extended to account for both the edge length in the
            // parent forest as well as the delta accumulated in the extension
            double edgelen_extension = 0.0;
            bool straddler = (new_nd->_height > _starting_height) && (child->_height < _starting_height);
            if (_is_extension && straddler) {
                edgelen_extension = _proposed_delta;
            }
                
            double pr_same = calcTransitionProbability(0, 0, child->_edge_length + edgelen_extension);
            double pr_diff = calcTransitionProbability(0, 1, child->_edge_length + edgelen_extension);
            for (unsigned p = 0; p < npatterns; p++) {
                unsigned pxnstates = p*G::_nstates;
                //unsigned pp = first_pattern + p;

                for (unsigned s = 0; s < G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? pr_same : pr_diff);
                        double child_partial = child_partial_array[pxnstates + s_child];
                                                
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
                    
                    parent_partial_array[pxnstates + s] *= sum_over_child_states;
                }   // parent state loop
            }   // pattern loop
        }   // child loop
    #endif

        // Compute the ratio of after to before likelihoods
        //TODO: make more efficient
        double prev_loglike = 0.0;
        double curr_loglike = 0.0;
        auto & newnd_partial_array = new_nd->_partial->_v;
        auto & lchild_partial_array = lchild->_partial->_v;
        auto & rchild_partial_array = rchild->_partial->_v;
        for (unsigned p = 0; p < npatterns; p++) {
            unsigned pxnstates = p*G::_nstates;
            
            unsigned pp = first_pattern + p;
            //unsigned count = counts[pp];
            double left_sitelike = 0.0;
            double right_sitelike = 0.0;
            double newnd_sitelike = 0.0;
#if defined (UNROLL_LOOPS)
            // loop 0
            unsigned s = 0;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];
            
            // loop 1
            s = 1;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];
            
            // loop 2
            s = 2;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];
            
            // loop 3
            s = 3;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];

#else
            for (unsigned s = 0; s < G::_nstates; s++) {
                left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
                right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
                newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];
            }
#endif
            prev_loglike += log(left_sitelike)*counts[pp];
            prev_loglike += log(right_sitelike)*counts[pp];
            curr_loglike += log(newnd_sitelike)*counts[pp];
        }
        
        return curr_loglike - prev_loglike;
    }
#else
    inline void Forest::calcPartialArrayJC(Node * new_nd) {
        assert (_index > 0);
    
        if (!new_nd->_left_child) {
            auto &data_matrix=_data->getDataMatrix();
            assert (G::_save_memory || G::_start_mode == "sim");
            if (!new_nd->_left_child) {
#if defined (REUSE_PARTIALS)
                mtx.lock();
                new_nd->_partial=ps.getPartial(_npatterns*4, _index);
                mtx.unlock();
#else
                new_nd->_partial=ps.getPartial(_npatterns*4);
#endif
                for (unsigned p=0; p<_npatterns; p++) {
                    unsigned pp = _first_pattern+p;
                    for (unsigned s=0; s<G::_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[new_nd->_number][pp];
                        double result = state & d;
#if defined (REUSE_PARTIALS)
                        (new_nd->_partial->_v)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*new_nd->_partial)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                    }
                }
            }
        }
    
#if defined (REUSE_PARTIALS)
        auto & parent_partial_array = new_nd->_partial->_v;
#else
        auto & parent_partial_array = *(new_nd->_partial);
#endif
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {

            if (child->_partial == nullptr) {
#if defined (REUSE_PARTIALS)
                mtx.lock();
                child->_partial = ps.getPartial(_npatterns*4, _index);
                mtx.unlock();
#else
                child->_partial = ps.getPartial(_npatterns*4);
#endif
                calcPartialArrayJC(child);
            }
            assert (child->_partial != nullptr);
#if defined (REUSE_PARTIALS)
            auto & child_partial_array = child->_partial->_v;
#else
            auto & child_partial_array = *(child->_partial);
#endif
            
            // calculate both transition probabilities before loop, then choose 1 based on s == s_child
            double transition_prob_same = calcTransitionProbabilityJC(0, 0, child->_edge_length);
            double transition_prob_dif = calcTransitionProbabilityJC(0, 1, child->_edge_length);
            
#if defined (UNROLL_LOOPS)
            for (unsigned p = 0; p < _npatterns; p++) {
                unsigned pxnstates = p*G::_nstates;
                
                // loop 1
                unsigned s = 0;
                double sum_over_child_states = 0.0;
                unsigned index = pxnstates+s;
                
                // subloop 1
                unsigned s_child = 0;
                double child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                double child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 2
                s_child = 1;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 3
                s_child = 2;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 4
                s_child = 3;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;

                parent_partial_array[index] *= sum_over_child_states;
            
                // loop 2
                s = 1;
                sum_over_child_states = 0.0;
                index = pxnstates+s;
                
                // subloop 1
                s_child = 0;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 2
                s_child = 1;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 3
                s_child = 2;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 4
                s_child = 3;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                parent_partial_array[index] *= sum_over_child_states;
            
                // loop 3
                s = 2;
                sum_over_child_states = 0.0;
                index = pxnstates+s;
                
                // subloop 1
                s_child = 0;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 2
                s_child = 1;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 3
                s_child = 2;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 4
                s_child = 3;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;

                parent_partial_array[index] *= sum_over_child_states;
            
                // loop 4
                s = 3;
                sum_over_child_states = 0.0;
                index = pxnstates+s;
                
                // subloop 1
                s_child = 0;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 2
                s_child = 1;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 3
                s_child = 2;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                // subloop 4
                s_child = 3;
                child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                child_partial = child_partial_array[pxnstates + s_child];
                sum_over_child_states += child_transition_prob * child_partial;
                
                parent_partial_array[index] *= sum_over_child_states;
            }
        } // child loop
#else
            for (unsigned p = 0; p < _npatterns; p++) {
                for (unsigned s = 0; s <G::_nstates; s++) {
                    double sum_over_child_states = 0.0;
                    for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                        double child_transition_prob = (s == s_child ? transition_prob_same : transition_prob_dif);
                        double child_partial = child_partial_array[p*G::_nstates + s_child];
                        sum_over_child_states += child_transition_prob * child_partial;
                    }   // child state loop
//                    if (child == new_nd->_left_child) {
//                        parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
//                    }
//                    else {
                        parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
//                    }
                }   // parent state loop
            }   // pattern loop
        }   // child loop
#endif
    }
#endif


#if !defined (LAZY_COPYING)
    inline void Forest::calcPartialArrayHKY(Node * new_nd) {
        assert (_index > 0);

        if (!new_nd->_left_child) {
            auto &data_matrix=_data->getDataMatrix();
            assert (G::_save_memory || G::_start_mode_type == G::StartModeType::START_MODE_SIM);
            if (!new_nd->_left_child) {
#if defined (REUSE_PARTIALS)
                mtx.lock();
                new_nd->_partial=ps.getPartial(_npatterns*4, _index);
                mtx.unlock();
#else
                new_nd->_partial=ps.getPartial(_npatterns*4);
#endif
                for (unsigned p=0; p<_npatterns; p++) {
                    unsigned pp = _first_pattern+p;
                    
                    // TODO: unroll this loop
                    for (unsigned s=0; s<G::_nstates; s++) {
                        Data::state_t state = (Data::state_t)1 << s;
                        Data::state_t d = data_matrix[new_nd->_number][pp];
                        double result = state & d;
#if defined (REUSE_PARTIALS)
                        (new_nd->_partial->_v)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
#else
                        (*new_nd->_partial)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
#endif
                    }
                }
            }
        }

#if defined (REUSE_PARTIALS)
        auto & parent_partial_array = new_nd->_partial->_v;
#else
        auto & parent_partial_array = *(new_nd->_partial);
#endif
        for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {

            if (child->_partial == nullptr) {
#if defined (REUSE_PARTIALS)
                mtx.lock();
                child->_partial = ps.getPartial(_npatterns*4, _index);
                mtx.unlock();
#else
                child->_partial = ps.getPartial(_npatterns*4);
#endif
                calcPartialArrayHKY(child);
            }
            assert (child->_partial != nullptr);
#if defined (REUSE_PARTIALS)
            auto & child_partial_array = child->_partial->_v;
#else
            auto & child_partial_array = *(child->_partial);
#endif

            for (unsigned p = 0; p < _npatterns; p++) {
                // TODO: unroll this loop
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
#endif

    inline double Forest::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) {
        double relative_rate = G::_double_relative_rates[_index-1];
        assert(pi.size() == 4);
        assert(fabs(accumulate(pi.begin(), pi.end(), 0.0) - 1.0) < G::_small_enough);
        assert(relative_rate > 0.0);
        
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
        double betat = 0.5*relative_rate*edge_length/((pi[0] + pi[2])*(pi[1] + pi[3]) + kappa*(pi[0]*pi[2] + pi[1]*pi[3]));
        
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

    inline double Forest::calcTransitionProbabilityJC(double s, double s_child, double edge_length) const {
        double relative_rate = G::_double_relative_rates[_index-1];
        assert (relative_rate > 0.0);
        
        double child_transition_prob = 0.0;

            if (s == s_child) {
                child_transition_prob = 0.25 + 0.75*exp(-4.0*relative_rate*edge_length/3.0);
            }
            
            else {
                child_transition_prob = 0.25 - 0.25*exp(-4.0*relative_rate*edge_length/3.0);
            }
            return child_transition_prob;
    }

    inline double Forest::calcTransitionProbabilityHKY(double s, double s_child, double edge_length) {
        double relative_rate = G::_double_relative_rates[_index-1];
        assert (relative_rate > 0.0);
        
        double child_transition_prob = 0.0;
    
        double pi_A = G::_base_frequencies[0];
        double pi_C = G::_base_frequencies[1];
        double pi_G = G::_base_frequencies[2];
        double pi_T = G::_base_frequencies[3];

        double pi_j = 0.0;
        double PI_J = 0.0;

        double phi = (pi_A+pi_G)*(pi_C+pi_T)+_kappa*(pi_A*pi_G+pi_C*pi_T);
        double beta_t = 0.5*(edge_length * relative_rate )/phi;

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
        
#if defined (UPGMA)
        if (G::_upgma) {
            if (_forest_height > 0) {
                constructUPGMA();
            }
        }
#endif

#if defined (UNROLL_LOOPS)
        for (auto &nd:_lineages) {
            double log_like = 0.0;
            for (unsigned p=0; p<_npatterns; p++) {
                double site_like = 0.0;
                unsigned pxnstates = p*G::_nstates;
                
                // loop 1
                unsigned s = 0;
#if defined (REUSE_PARTIALS)
                double partial = (nd->_partial->_v)[pxnstates+s];
#else
                double partial = (*nd->_partial)[pxnstates+s];
#endif
                site_like += 0.25*partial;
                
                // loop 2
                s = 1;
#if defined (REUSE_PARTIALS)
                partial = (nd->_partial->_v)[pxnstates+s];
#else
                partial = (*nd->_partial)[pxnstates+s];
#endif
                site_like += 0.25*partial;
                
                // loop 3
                s = 2;
#if defined (REUSE_PARTIALS)
                partial = (nd->_partial->_v)[pxnstates+s];
#else
                partial = (*nd->_partial)[pxnstates+s];
#endif
                site_like += 0.25*partial;
                
                // loop 4
                s = 3;
#if defined (REUSE_PARTIALS)
                partial = (nd->_partial->_v)[pxnstates+s];
#else
                partial = (*nd->_partial)[pxnstates+s];
#endif
                site_like += 0.25*partial;
            
                
                assert(site_like>0);
                log_like += log(site_like)*counts[_first_pattern+p];
            }

            _gene_tree_log_likelihood += log_like;
            
//            debugLogLikelihood(nd, log_like);
        }
#else
        for (auto &nd:_lineages) {
            double log_like = 0.0;
            for (unsigned p=0; p<_npatterns; p++) {
                double site_like = 0.0;
                for (unsigned s=0; s<G::_nstates; s++) {
#if defined (REUSE_PARTIALS)
                    double partial = (nd->_partial->_v)[p*G::_nstates+s];
#else
                    double partial = (*nd->_partial)[p*G::_nstates+s];
#endif
                    site_like += 0.25*partial;
                }
                assert(site_like>0);
                log_like += log(site_like)*counts[_first_pattern+p];
            }

            _gene_tree_log_likelihood += log_like;
            
//            debugLogLikelihood(nd, log_like);
        }
#endif
        
#if defined (UPGMA)
        if (G::_upgma) {
            if (_forest_height) {
                destroyUPGMA();
            }
        }
#endif
        
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

    inline void Forest::operator=(const Forest & other) {
        if (other._nodes.size() > 0) {
            _nodes.resize(other._nodes.size()); // don't need to clear these members because they will get overwritten
            _lineages.resize(other._lineages.size());
        }
        
        _gene_tree_log_likelihood = other._gene_tree_log_likelihood;
        
        // the following data members apply only to the first round
        if (!G::_in_second_level) {
            _ninternals         = other._ninternals;
            _last_edge_length   = other._last_edge_length;
            
            _preorder.resize(other._preorder.size());
            
            _log_joining_prob = other._log_joining_prob;
            _npatterns = other._npatterns;
            _edge_rate_variance = other._edge_rate_variance;
            _asrv_shape = other._asrv_shape;
            _comphet = other._comphet;
            _first_pattern      = other._first_pattern;
            _data               = other._data;
            
#if defined (UPGMA)
            if (G::_upgma) {
                _upgma_additions = other._upgma_additions;
                _upgma_starting_edgelen = other._upgma_starting_edgelen;

                _dmatrix = other._dmatrix;
                _dmatrix_rows = other._dmatrix_rows;
                
                _starting_dij = other._starting_dij;
            }
#endif
            
#if defined (LAZY_COPYING)
            _preorders.clear();
#endif
        }
        
        _index              = other._index; // TODO: things like this will be the same as the copied particle - need to copy them in the constructor still?
        _increments_and_priors = other._increments_and_priors;
        _forest_length = other._forest_length;
            
        _forest_height = other._forest_height;
        
        // the following data members apply only when simulating and do not need to be copied because simulating data only deals with one particle at a time
//        _lineages_per_species = other._lineages_per_species;
        
        // the following data members get reset if using prior-post and do not need to be copied
//       _node_choices = other._node_choices;
//       _log_likelihood_choices = other._log_likelihood_choices;
        
        
        // the following data members apply only to the second level
        if (G::_in_second_level) {
            _coalinfo.reserve(other._coalinfo.size());
            _coalinfo = other._coalinfo;

            _next_node_number = other._next_node_number;
        }
        
#if defined (DEBUG_MODE)
        _species_joined = other._species_joined;
#endif

        // copy tree itself
        
        if (other._nodes.size() > 0) { // otherwise, there is no forest and nothing needs to be copied
            
#if !defined (LAZY_COPYING) // no species partition to copy
            if (!G::_in_second_level) {
                
                bool faster_spp_partition = true;
                
                if (!faster_spp_partition) {
                        
                    _species_partition.clear();

                    for (auto spiter : other._species_partition) {
                        for (auto s : spiter.second) {
                            unsigned number = s->_number;
                            Node* nd = &_nodes[number];
                            _species_partition[spiter.first].push_back(nd);
                        }
                    }
                }
                    
                else {
                    vector<string> names_to_erase;
                    names_to_erase.reserve(_species_partition.size());
                    for (auto &s:_species_partition) {
                        names_to_erase.push_back(s.first);
                    }

                    unsigned count = 0;
                    for (auto &spiter : other._species_partition) {
                        // if key is found in names_to_erase, delete it from names_to_erase because it shouldn't be erased
                        auto position = std::find(names_to_erase.begin(), names_to_erase.end(), spiter.first);
                        if (position != names_to_erase.end()) {
                            names_to_erase.erase(position);
                        }
                        count = 0;
                        for (auto &s : spiter.second) {
                            unsigned number = s->_number;
                            Node* nd = &_nodes[number];
                            if (count == 0) {
                                _species_partition[spiter.first].clear();
                                _species_partition[spiter.first].push_back(nd);
                            }
                            else {
                                _species_partition[spiter.first].push_back(nd);
                            }
                            count++;
                        }
                    }
                    
                    for (unsigned n=0; n<names_to_erase.size(); n++) {
                        _species_partition.erase(names_to_erase[n]);
                    }
                }
        }
#endif

            for (auto &othernd : other._nodes) {
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
                    nd->_partial = othernd._partial;
                    nd->_height = othernd._height;
                    nd->_split = othernd._split;
                    nd->_species = othernd._species;
                }
            }

            unsigned j = 0;
            for (auto &othernd : other._lineages) {
                unsigned k = othernd->_number;
                Node * nd = &_nodes[k];
                _lineages[j] = nd;
                j++;
            }
            
            if (other._preorder.size() > 0) {
                unsigned m = 0;
                for (auto &othernd : other._preorder) {
                    unsigned n = othernd->_number;
                    Node * nd = &_nodes[n];
                    _preorder[m] = nd;
                    m++;
                }
            }
        }
    }

    inline void Forest::updateSpeciesPartition(tuple<string, string, string> species_info) {
        string spp1 = get<0>(species_info);
        string spp2 = get<1>(species_info);
        string new_spp = get<2>(species_info);
        
        unsigned before = (int) _species_partition.size();

        vector<Node*> &nodes = _species_partition[new_spp];
        
        copy(_species_partition[spp1].begin(), _species_partition[spp1].end(), back_inserter(nodes));
        copy(_species_partition[spp2].begin(), _species_partition[spp2].end(), back_inserter(nodes));
        _species_partition.erase(spp1);
        _species_partition.erase(spp2);
        
        if (spp1 != "null") {
            assert (_species_partition.size() == before - 1);
        }
    }

#if defined (LAZY_COPYING)
    inline void Forest::updateSpeciesPartitionSim(tuple<string, string, string> species_info, G::species_t new_species_number) {
        string spp1 = get<0>(species_info);
        string spp2 = get<1>(species_info);
        string new_spp = get<2>(species_info);
        
        unsigned before = (int) _species_partition.size();

        vector<Node*> &nodes = _species_partition[new_spp];
        
        copy(_species_partition[spp1].begin(), _species_partition[spp1].end(), back_inserter(nodes));
        copy(_species_partition[spp2].begin(), _species_partition[spp2].end(), back_inserter(nodes));
        _species_partition.erase(spp1);
        _species_partition.erase(spp2);
        
        if (spp1 != "null") {
            assert (_species_partition.size() == before - 1);
        }
        
        for (auto &nd:nodes) {
            nd->_species = new_species_number;
        }
    }
#endif

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
        assert (_index >0);
        assert (_species_partition.size() == 0);
//        _species_partition.clear();
        
        unsigned count = 0;
        
#if defined (LAZY_COPYING)
            G::species_t curr_species = 0;
            string prev_species_name = "";
#endif
        for (auto &nd:_nodes) {
            count++;
            assert (!nd._left_child);
            string species_name = taxon_map[nd._name];
            
#if defined (LAZY_COPYING)
            if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
            if (prev_species_name != species_name) {
                if (curr_species == 0) {
                    curr_species = 1;
                }
                else {
                    curr_species += curr_species;
                }
            }
        }
#endif
            
            if (!G::_gene_newicks_specified) {
                _species_partition[species_name].push_back(&nd); // TODO: may need this for old second level even if starting from gene newicks?
#if defined (LAZY_COPYING)
                if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
                    nd._species = curr_species;
                    prev_species_name = species_name;
                }
#endif
            }
            
            if (count == G::_ntaxa) {
                break;
            }
            if (G::_start_mode_type != G::StartModeType::START_MODE_SIM) {
                if (G::_taxon_to_species.count(nd._name) == 0) {
                    throw XProj(str(format("Could not find an index for the taxon name \"%s\"") % nd._name));
                }
                else {
                    Node::setSpeciesBit(nd._species, G::_taxon_to_species.at(nd._name), /*init_to_zero_first*/true);
                }
            }
        }
        
        if (!G::_gene_newicks_specified) {
            assert (_species_partition.size() == G::_nspecies);
        }
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

#if defined (UPGMA)
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
            
#if defined (REUSE_PARTIALS)
            mtx.lock();
            stowPartial(parent);
            mtx.unlock();
#endif

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
//            _nodes.pop_back(); // remove unused node from node list
        
            _ninternals--;
            
        }
        
        // Restore starting edge lengths
        for (auto &nd : _lineages) {
            nd->_edge_length = _upgma_starting_edgelen.at(nd);
        }
        
        _upgma_starting_edgelen.clear();
        
        if (_lineages.size() == 1) {
            _lineages.back()->_partial = nullptr; // last step, only node with partials should be the new node, and partials are no longer needed if likelihood has been calculated
        }
    }
#endif

#if defined (REUSE_PARTIALS)
    inline void Forest::stowPartial(Node * nd) {
        if (nd && nd->_left_child && nd->_partial) {
            // Nothing to do if nd or nd->_partial is null
            // Only internal partials are ever reset/stowed;
            // leaf partials are calculated at the beginning
            // and never changed.
            if (nd->_partial.use_count() == 1) {
                // Partial is not being used by any other node, so it is
                // safe to stow it in PartialStore for reuse later
                ps.putPartial(_index, nd->_partial);
            }
        
            nd->_partial.reset();
        }
    }
#endif

#if defined (UPGMA)
    inline void Forest::constructUPGMA() {
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
            for (auto &nd : _lineages) {
                nd->_edge_length += dh;
            }
            
            //debugShowLineages();

            // Join lineages i and j
//            Node nd;
//            _nodes.push_back(nd);
//            Node* new_nd = &_nodes.back();
            Node * new_nd = &_nodes[G::_ntaxa + _ninternals];
            
            assert (new_nd->_parent==0);
           assert (new_nd->_number == -1);
           assert (new_nd->_right_sib == 0);

            Node * subtree1 = node_for_row[i];
            Node * subtree2 = node_for_row[j];
            
//            new_nd->_parent=0;
            new_nd->_number=G::_ntaxa+_ninternals;
//            new_nd->_right_sib=0;

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
#if defined (REUSE_PARTIALS)
            mtx.lock();
            new_nd->_partial=ps.getPartial(_npatterns*4, _index);
            mtx.unlock();
#else
            new_nd->_partial=ps.getPartial(_npatterns*4);
#endif
            assert(new_nd->_left_child->_right_sib);
            
#if defined (LAZY_COPYING)
            calcPartialArrayJC(new_nd, new_nd->_left_child, new_nd->_left_child->_right_sib); // use JC model for all UPGMA
#else
            calcPartialArrayJC(new_nd);
#endif
            
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
    }
#endif

#if defined (UPGMA)
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

    inline pair<Node*, Node*> Forest::chooseAllPairs(vector<Node*> &node_list, double increment, string species, Lot::SharedPtr lot) {
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

    inline pair<Node*, Node*> Forest::getSubtreeAt(pair<unsigned, unsigned> t, vector<Node*> node_list) {
          Node *subtree1 = nullptr;
          Node *subtree2 = nullptr;

          unsigned a = 0; // TODO: can make this fater because node_list is actually a vector now
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


    inline tuple<Node*, Node*, Node*> Forest::createNewSubtree(pair<unsigned, unsigned> t, vector<Node*> node_list, double increment, string species) {
         pair<Node*, Node*> p = getSubtreeAt(t, node_list);

         Node* subtree1 = p.first;
         Node* subtree2 = p.second;

    //        new node is always needed
//         Node nd;
//         _nodes.push_back(nd);
//         Node* new_nd = &_nodes.back();
        
        Node * new_nd = &_nodes[G::_ntaxa + _ninternals];
        
        assert (new_nd->_parent==0);
       assert (new_nd->_number == -1);
       assert (new_nd->_right_sib == 0);
        
//         new_nd->_parent=0;
         new_nd->_number=G::_ntaxa+_ninternals;
         new_nd->_edge_length=0.0;
//         new_nd->_right_sib=0;

         new_nd->_left_child=subtree1;
         subtree1->_right_sib=subtree2;

         subtree1->_parent=new_nd;
         subtree2->_parent=new_nd;

         //always calculating partials now
         assert (new_nd->_partial == nullptr);
#if defined (REUSE_PARTIALS)
        mtx.lock();
        new_nd->_partial=ps.getPartial(_npatterns*4, _index);
        mtx.unlock();
#else
         new_nd->_partial=ps.getPartial(_npatterns*4);
#endif
         assert(new_nd->_left_child->_right_sib);
#if defined (LAZY_COPYING)
        calcPartialArrayJC(new_nd, new_nd->_left_child, new_nd->_left_child->_right_sib);
#else
        if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
            calcPartialArrayJC(new_nd);
        }
        else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
            calcPartialArrayHKY(new_nd);
        }
        else {
            throw XProj("model must be either JC or HKY");
        }
#endif
        new_nd->_height = _forest_height;

         // don't update the species list
         updateNodeVector(_lineages, subtree1, subtree2, new_nd);
                  
         return make_tuple(subtree1, subtree2, new_nd);
     }

    inline void Forest::allowCoalescence(string species_name, double increment, Lot::SharedPtr lot) {
         double prev_log_likelihood = _gene_tree_log_likelihood;
        
         Node *subtree1 = nullptr;
         Node *subtree2 = nullptr;
        vector<Node*> nodes = _species_partition[species_name];
        
        assert (nodes.size() > 0);

         unsigned s = (unsigned) nodes.size();
         calcTopologyPrior(s);

         assert (s > 1);
         bool one_choice = false;
         if (nodes.size() == 2) {
             one_choice = true;
         }

         if (!G::_prior_prior && (!one_choice)) {
             if (G::_save_memory) {
                 for (auto &nd:_lineages) {
                     if (nd->_partial == nullptr) {
#if defined (REUSE_PARTIALS)
                         mtx.lock();
                         nd->_partial = ps.getPartial(_npatterns*4, _index);
                         mtx.unlock();
#else
                         nd->_partial = ps.getPartial(_npatterns*4);
#endif
#if defined (LAZY_COPYING)
                         calcPartialArrayJC(nd, nd->_left_child, nd->_left_child->_right_sib);
#else
                         if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                             calcPartialArrayJC(nd);
                         }
                         else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
                             calcPartialArrayHKY(nd);
                         }
                         else {
                             throw XProj("model must be either HKY or JC");
                         }
#endif
                     }
                 }
             }
             
             pair<Node*, Node*> t = chooseAllPairs(nodes, increment, species_name, lot);
             
             subtree1 = t.first;
             subtree2 = t.second;
         }
         
         else {
             assert (G::_prior_prior || one_choice);
             // prior-prior proposal
             pair<unsigned, unsigned> t = chooseTaxaToJoin(s, lot);
             subtree1 = nodes[t.first];
             subtree2 = nodes[t.second];

             assert (t.first < nodes.size());
             assert (t.second < nodes.size());

             assert (subtree1 != subtree2);
         }

         // access next unused node
        
        Node * new_nd = &_nodes[G::_ntaxa + _ninternals];

        assert (new_nd->_parent==0);
        assert (new_nd->_number == -1);
        assert (new_nd->_right_sib == 0);
        
         new_nd->_number=G::_ntaxa+_ninternals;
         new_nd->_edge_length=0.0;
         _ninternals++;
         new_nd->_name += to_string(new_nd->_number);

         new_nd->_left_child=subtree1;
         subtree1->_right_sib=subtree2;

         subtree1->_parent=new_nd;
         subtree2->_parent=new_nd;

         if (!G::_run_on_empty) {
             //always calculating partials now
             assert (new_nd->_partial == nullptr);
#if defined (REUSE_PARTIALS)
             mtx.lock();
             new_nd->_partial=ps.getPartial(_npatterns*4, _index);
             mtx.unlock();
#else
             new_nd->_partial=ps.getPartial(_npatterns*4);
#endif
             assert(new_nd->_left_child->_right_sib);

             if (G::_save_memory) {
                 for (auto &nd:_lineages) {
                     if (nd->_partial == nullptr) {
#if defined (REUSE_PARTIALS)
                         mtx.lock();
                         nd->_partial = ps.getPartial(_npatterns*4, _index);
                         mtx.unlock();
#else
                         nd->_partial = ps.getPartial(_npatterns*4);
#endif
#if defined (LAZY_COPYING)
                         calcPartialArrayJC(nd, nd->_left_child, nd->_left_child->_right_sib);
#else
                         if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                             calcPartialArrayJC(nd);
                         }
                         else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
                             calcPartialArrayHKY(nd);
                         }
                         else {
                             throw XProj("model specified must be JC or HKY");
                         }
#endif
                     }
                 }
             }
             
#if defined (LAZY_COPYING)
             calcPartialArrayJC(new_nd, new_nd->_left_child, new_nd->_left_child->_right_sib);
#else
             if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                 calcPartialArrayJC(new_nd);
             }
             else if (G::_model_type == G::ModelType::MODEL_TYPE_HKY) {
                 calcPartialArrayHKY(new_nd);
             }
             else {
                 throw XProj("model specified must be JC or HKY");
             }
#endif
             new_nd->_height = _forest_height;
             
             new_nd->_species = 0;

             subtree1->_partial=nullptr; // throw away subtree partials now, no longer needed
             subtree2->_partial=nullptr;
             
#if defined (UPGMA)
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
#endif

         }
        
        //update species list
        updateNodeVector(nodes, subtree1, subtree2, new_nd);
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        if (G::_upgma) {
            new_nd->_split.resize(G::_ntaxa);
            new_nd->_split += subtree1->_split;
            new_nd->_split += subtree2->_split;
        }
        
#if defined (LAZY_COPYING)
        if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
            new_nd->_species = new_nd->_left_child->_species;
        }
#endif

        _species_partition[species_name] = nodes;

        if ((G::_prior_prior || one_choice) && (!G::_run_on_empty) ) {
            _gene_tree_log_likelihood = calcLogLikelihood();
            _log_weight = _gene_tree_log_likelihood - prev_log_likelihood;
        }

        if (G::_save_memory) {
            for (auto &nd:_nodes) {
                nd._partial=nullptr;
            }
        }
     }

#if defined (LAZY_COPYING)
    void Forest::mergeSpecies(G::species_t left_spp, G::species_t right_spp) {
//        showForest();
        // Merge species in _species_vect
        G::species_t anc_spp = (left_spp | right_spp);
        for (auto &nd:_lineages) {
            if (nd->_species == left_spp || nd->_species == right_spp) {
                nd->_species = anc_spp;
            }
        }
    }
#endif

//#if defined (LAZY_COPYING)
//    inline G::species_t Forest::getPrevSpecies() {
//        return _lineages.back()->_species;
//    }
//#endif

#if defined (LAZY_COPYING)
    inline vector<G::species_t> Forest::getAllPrevSpecies() {
        vector<G::species_t> species;
        for (auto &nd:_nodes) {
            species.push_back(nd._species);
        }
        return species;
    }
#endif

#if defined (LAZY_COPYING)
    inline void Forest::coalescePreviousJoin(pair<Node*, Node*> subtrees) {
        double prev_log_likelihood = calcLogLikelihood(); // TODO: save this instead of calculating it
        Node *subtree1 = subtrees.first;
        Node *subtree2 = subtrees.second;

        // access next unused node
        Node * new_nd = &_nodes[G::_ntaxa + _ninternals];

        assert (new_nd->_parent==0);
        assert (new_nd->_number == -1);
        assert (new_nd->_right_sib == 0);
        
        new_nd->_number=G::_ntaxa+_ninternals;
        new_nd->_edge_length=0.0;
        _ninternals++;
        new_nd->_name += to_string(new_nd->_number);

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
//        // set new node species
//        new_nd->setSpecies(species_name);

        // Set new node split to union of the two child splits
        new_nd->_split.resize(G::_ntaxa);
        
        new_nd->_split += subtree1->_split;
        new_nd->_split += subtree2->_split;
        
        new_nd->_height = _forest_height;
        
        new_nd->_partial=ps.getPartial(_npatterns*4, _index);
        
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        calcPartialArrayJC(new_nd, subtree1, subtree2); // TODO: add HKY
        
        _gene_tree_log_likelihood = calcLogLikelihood();
        _log_weight = _gene_tree_log_likelihood - prev_log_likelihood;
    }
#endif

#if defined (LAZY_COPYING)
    inline void Forest::coalesce(G::species_t species_name, Lot::SharedPtr lot, double prev_log_likelihood) {
        // create a species map
        
       vector<unsigned> eligible_lineages;
        
        for (auto &nd:_lineages) {
            if (nd->_species == species_name) {
                eligible_lineages.push_back(nd->_position_in_lineages);
            }
        }
        
        // choose two nodes to join
        unsigned s = (unsigned) eligible_lineages.size();
        assert (s > 1);
        pair<unsigned, unsigned> t = chooseTaxaToJoin(s, lot);
        
        unsigned child1_num = eligible_lineages[t.first];
        unsigned child2_num = eligible_lineages[t.second];
        
        Node *subtree1 = _lineages[child1_num];
        Node *subtree2 = _lineages[child2_num];

        // access next unused node
        Node * new_nd = &_nodes[G::_ntaxa + _ninternals];

        assert (new_nd->_parent==0);
        assert (new_nd->_number == -1);
        assert (new_nd->_right_sib == 0);
        
        new_nd->_number=G::_ntaxa+_ninternals;
        new_nd->_edge_length=0.0;
        _ninternals++;
        new_nd->_name += to_string(new_nd->_number);

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
        // set new node species
        new_nd->setSpecies(species_name);

        // Set new node split to union of the two child splits
        new_nd->_split.resize(G::_ntaxa);
        
        new_nd->_split += subtree1->_split;
        new_nd->_split += subtree2->_split;
        
        new_nd->_height = _forest_height;
        
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        new_nd->_partial=ps.getPartial(_npatterns*4, _index);
        calcPartialArrayJC(new_nd, subtree1, subtree2); // TODO: add HKY
        
        _gene_tree_log_likelihood = calcLogLikelihood();
        _log_weight = _gene_tree_log_likelihood - prev_log_likelihood;
    }
#endif

    inline pair<Node*, Node*> Forest::getPrevJoin() {
        return make_pair(_lineages.back()->_left_child, _lineages.back()->_left_child->_right_sib);
    }

    inline unsigned Forest::getNRemainingSpecies() {
        vector<G::species_t> existing_species;
        for (auto &nd:_lineages) {
            existing_species.push_back(nd->_species);
        }
        std::sort(existing_species.begin(), existing_species.end());
        unsigned n_unique_species = unique(existing_species.begin(), existing_species.end()) - existing_species.begin();
        return n_unique_species;
    }

    inline void Forest::revertGeneForest(double prev_gene_increment, vector<G::species_t> prev_species_assignments) {
        _forest_height -= prev_gene_increment;
        Node* child1 = _lineages.back()->_left_child;
        Node* child2 = _lineages.back()->_left_child->_right_sib;
        
        Node* ancestor = _lineages.back();
        revertNodeVector(_lineages, child1, child2, ancestor);
        
        child1->_parent = 0;
        child1->_right_sib = 0;
        child2->_parent = 0;
        child2->_right_sib = 0;
        ancestor->_species = 0;
        ancestor->_number = -1;
        ancestor->_name = "";
        
        for (unsigned i=0; i<_nodes.size(); i++) {
            _nodes[i]._species = prev_species_assignments[i];
        }
                
        for (auto &nd:_lineages) {
            nd->_edge_length -= prev_gene_increment;
            if (nd->_height > 0) { // if height is 0, node is a tip and height should remain 0
                nd->_height = getLineageHeight(nd->_left_child);
                assert (nd->_height > 0.0);
            }
        }
        
        stowPartial(ancestor);
        
        _ninternals--;
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

    inline void Forest::calcTreeLength() {
        // sum of all edge lengths in tree
        double sum_height = 0.0;
        
        for (auto &nd:_nodes) {
            // sum edge lengths from all nodes
            sum_height += nd._edge_length;
        }
        _forest_length = sum_height;
    }


    inline double Forest::getTreeLength() {
        // sum of all edge lengths in tree
        return _forest_length;
    }

#if defined (LAZY_COPYING)
    inline  vector<pair<double, unsigned long>> Forest::calcForestRate(Lot::SharedPtr lot, unordered_map<G::species_t, double> theta_map) {
        vector<pair<double, unsigned long>> rates;
        
        map<G::species_t, unsigned> species_counts;
        
        for (auto &nd:_lineages) {
            if (species_counts.count(nd->_species)) {
                species_counts[nd->_species]++;
            }
            else {
                species_counts[nd->_species];
            }
        }
        
        // counts start at 0, so 0 means 1 entry
        pair<double, unsigned long> rate_and_name;
        
        for (auto &s:species_counts) {
            if (s.second > 0) { // if size == 1, no possibility of coalescence and rate is 0
                double population_coalescence_rate = 0.0;
    #if defined (DRAW_NEW_THETA)
                double population_theta = theta_map[s.first];
                population_coalescence_rate = (s.second+1)*(s.second)/population_theta;
    #else
                population_coalescence_rate = (s.second+1)*(s.second)/G::_theta;
    #endif
                G::species_t name = s.first;
                rate_and_name = make_pair(population_coalescence_rate, name);
                rates.push_back(rate_and_name);
            }
        }
        return rates;
    }
#else
    inline vector<pair<double, string>> Forest::calcForestRate(Lot::SharedPtr lot, unordered_map<string, double> theta_map) {
        vector<pair<double, string>> rates;
        pair<double, string> rate_and_name;

        for (auto &s:_species_partition) {
            if (s.second.size() > 1) { // if size == 0, no possibility of coalescence and rate is 0
                double population_coalescence_rate = 0.0;
    #if defined (DRAW_NEW_THETA)
                double population_theta = theta_map[s.first];
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
#endif

#if defined (LAZY_COPYING)
    inline vector<pair<double, string>> Forest::calcForestRateSim(Lot::SharedPtr lot, unordered_map<G::species_t, double> theta_map) {
#else
        inline vector<pair<double, string>> Forest::calcForestRateSim(Lot::SharedPtr lot, unordered_map<string, double> theta_map) {
#endif
        vector<pair<double, string>> rates;
        pair<double, string> rate_and_name;

        for (auto &s:_species_partition) {
            if (s.second.size() > 1) { // if size == 0, no possibility of coalescence and rate is 0
                double population_coalescence_rate = 0.0;
    #if defined (DRAW_NEW_THETA)
                G::species_t species = s.second[0]->_species;
#if defined (LAZY_COPYING)
                double population_theta = theta_map[species];
#else
                double population_theta = theta_map[s.first];
#endif
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
        return (unsigned) (count - 1);
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
        unsigned curr_internal = G::_ntaxa;
        for (auto & nd : boost::adaptors::reverse(_preorder)) {
            if (nd->_left_child) {
                // nd is an internal node
                nd->_number = curr_internal++;
            }
        }

        _ninternals = curr_internal - G::_ntaxa;

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
//        _nleaves = countNewickLeaves(commentless_newick);
        
    //        if (_nleaves < 4) {
    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //        }
        unsigned max_nodes = 2*G::_ntaxa - (rooted ? 0 : 2);
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
            for (auto &ch : commentless_newick) {
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
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % G::_ntaxa));

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

    inline void Forest::buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies) {
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
        G::_ntaxa = countNewickLeaves(commentless_newick);
    //        if (_nleaves < 4) {
    //            throw XProj("Expecting newick tree description to have at least 4 leaves");
    //        }
        unsigned max_nodes = 2*G::_ntaxa - (rooted ? 0 : 2);
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
            for (auto &ch : commentless_newick) {
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
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % G::_ntaxa));

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

#if defined (USING_MPI)
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
        G::_ntaxa = countNewickLeaves(commentless_newick);
        unsigned max_nodes = 2*G::_ntaxa - (rooted ? 0 : 2);
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
            for (auto &ch : commentless_newick) {
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
                            throw XProj(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _nodes.size() % G::_ntaxa));
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

    //            unsigned max_nodes = countNewickInternals(newick) + G::_ntaxa;
            unsigned ninternals = countNewickInternals(newick);
            unsigned max_internals = G::_ntaxa-1;
            unsigned max_nodes = ninternals + G::_ntaxa + 1;
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
            
            unsigned num = G::_ntaxa;
            
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
#endif

#if defined (USING_MPI)
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
                    
//                    Node nd;
//                    _nodes.push_back(nd);
//                    Node* new_nd = &_nodes.back();
                    assert (1 == 2); // TODO: fix for _nodes as vector
                    Node * new_nd = &_nodes[G::_ntaxa + _ninternals];
                    
                    new_nd->_parent=0;
                    new_nd->_number=G::_ntaxa+_ninternals;
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
                    Node * new_nd = &_nodes[G::_ntaxa + _ninternals];
                    new_nd->_parent=0;
                    new_nd->_number=G::_ntaxa+_ninternals;
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
#endif

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

    inline void Forest::saveCoalInfoGeneForest(vector<Forest::coalinfo_t> & coalinfo_vect) const {
        // Appends to coalinfo_vect; clear before calling if desired
        // GeneForest version ignores cap argument.
        // Assumes heights and preorders are up-to-date
        
        // coalinfo_t is a tuple with these elements:
        // - height of node
        // - 1-offset gene index (0 means speciation)
        // - vector of child species
        
        // Should only be called for complete gene trees
    //        assert(_lineages.size() == 1);

        // Copy tuples stored in _coalinfo to end of coalinfo_vect
        coalinfo_vect.insert(coalinfo_vect.end(), _coalinfo.begin(), _coalinfo.end());
    }

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
            for (auto&  nd : boost::adaptors::reverse(_preorder)) {
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

    inline bool Forest::subsumed(G::species_t test_species, G::species_t subtending_species) {
        bool not_equal = (test_species != subtending_species);
        bool is_subset = ((test_species & subtending_species) == test_species);
        if (not_equal && is_subset)
            return true;
        else
            return false;
    }

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

    inline void Forest::refreshPreorderNew(vector<Node*> & preorder) const {
        // Assumes preorder just contains the root node when this function is called
        // Also assumes that _next_node_number was initialized prior to calling this function
        assert(preorder.size() == 1);
        
        Node * nd = preorder[0];
        while (true) {
            nd = findNextPreorderNew(nd);
            if (nd) {
                preorder.push_back(nd);
//                if (nd->_left_child) {
//                    nd->_number = _next_node_number++;
//                }
            }
            else
                break;
        }
    }

    inline void Forest::setTreeHeight() {
        _forest_height = _lineages.back()->_height;
    }

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

#if defined (USING_MPI)
    inline void Forest::resetSpeciesPartition(string species_partition_string) {
        _species_partition.clear();
        cout << species_partition_string << endl;
        for (auto &s:species_partition_string) {
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
#endif

#if defined (USING_MPI)
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
#endif

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

#if defined(LAZY_COPYING)
    inline Node * Forest::getMostRecentAncestor() {
        assert(_lineages.size() == 1); // assumes forest is a complete tree
        assert(_preorders[0].size() > 0); // TODO: need to make sure _preorders has been created
        Node * anc = nullptr;
        for (auto nd : boost::adaptors::reverse(_preorders[0])) {
            if (nd->_left_child) {
                anc = nd;
                break;
            }
        }
        assert(anc);
        return anc;
    }
#endif

#if defined (LAZY_COPYING)
    inline void Forest::refreshAllPreorders() const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        _preorders.clear();
        if (_lineages.size() == 0) {
            return;
        }

        for (auto & nd : _lineages) {
//            if (nd->_left_child) {
////                nd->_number = _next_node_number++;
//            }

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

#if defined(LAZY_COPYING)
    inline pair<double, double> Forest::getBoundaryExtension() const {
        return make_pair(-1.0, 0.0);
    }
#endif

#if defined (LAZY_COPYING)
    inline void Forest::advanceAllLineagesBy(double increment) {
        // Add dt to the edge length of all lineage root nodes, unless
        // 1. there is just one lineage or
        // 2. this is a gene forest extension and the
        //    lineage root node belongs to the parent,
        // in which case do nothing
        unsigned n = (unsigned)_lineages.size();
        if (n > 1) {
            for (auto &nd : _lineages) {
                double edge_len = nd->_edge_length + increment;
                assert(edge_len >= 0.0 || fabs(edge_len) < Node::_smallest_edge_length);
                nd->_edge_length = edge_len;
            }
        
            // Add to to the current forest height
            _forest_height += increment;
        }
    }
#endif

#if defined(LAZY_COPYING)
    inline void Forest::buildCoalInfoVect() {
        // Assumes heights of all nodes are accurate
        
        // Assumes this is not a gene forest extension
        assert(getBoundaryExtension().first < 0);
        
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

#if defined(LAZY_COPYING)
    inline void Forest::mergeSpecies(G::species_t left_species, G::species_t right_species, G::species_t anc_species) {
        // Every node previously assigned to left_species
        // or right_species should be reassigned to anc_species
        
        // Create a functor that assigns anc_species to the
        // supplied nd if it is currently in either left_species
        // or right_species
        for (auto &nd:_lineages) {
            if (nd->_species == left_species || nd->_species == right_species) {
                nd->setSpecies(anc_species);
            }
        }
//        auto reassign = [left_species, right_species, anc_species, this](Node * nd) {
//            G::species_t ndspp = nd->getSpecies();
//            if (ndspp == left_species || ndspp == right_species) {
//                nd->setSpecies(anc_species);
//            }
//        };
//
//        // Apply functor reassign to each node in _lineages
//        for_each(_lineages.begin(), _lineages.end(), reassign);
    }
#endif

#if defined(LAZY_COPYING)
    inline void Forest::copyLineageSpecies(vector<G::species_t> & species_of_lineages) const {
        species_of_lineages.resize(_lineages.size());
        unsigned i = 0;
        for (auto nd : _lineages) {
            species_of_lineages[i++] = nd->_species;
        }
    }
#endif

#if defined (LAZY_COPYING)
    inline unsigned Forest::checkNumberOfUniqueSpeciesInExistence() {
        vector<G::species_t> species_in_existence;
        for (auto &nd:_lineages)
            if (std::find(species_in_existence.begin(), species_in_existence.end(), nd->_species) == species_in_existence.end()) {
                species_in_existence.push_back(nd->_species);
            }
        return (unsigned) species_in_existence.size();
        }
#endif

#if defined(LAZY_COPYING)
    inline void Forest::debugCheckBleedingEdge(string msg, double anc_height) const {
        // Find maximum height of all nodes in _lineages vector
        double maxh = 0.0;
        for (auto nd : _lineages) {
            if (nd->_height > maxh)
                maxh = nd->_height;
        }

        // Find maximum height + edgelen of all nodes in _lineages vector
        double maxhplus = 0.0;
        for (auto nd : _lineages) {
            if (nd->_height + nd->_edge_length > maxhplus)
                maxhplus = nd->_height + nd->_edge_length;
        }
        
        double diff = fabs(maxhplus - _forest_height);
        //if (diff > G::_small_enough) {
        //    cerr << endl;
        //}
        assert(diff <= G::_small_enough);
        
        double diff2 = fabs(_forest_height - anc_height);
        //if (diff2 > G::_small_enough) {
        //    cerr << endl;
        //}
        assert(diff2 <= G::_small_enough);
    }
#endif

    inline double Forest::getLogLikelihood() const {
        return _gene_tree_log_likelihood;
    }

//#if defined(LAZY_COPYING)
//    inline void Forest::computeAllPartials() {
//        // Assumes _leaf_partials have been computed but that every node in the tree
//        // has _partial equal to nullptr.
//        assert(_data);
//        assert(_gene_index >= 0);
//
//        if (_preorders.size() == 0) {
//            refreshAllPreorders();
//        }
//
//        for (auto & preorder : _preorders) {
//            // Visit nodes in the subtree rooted at preorder in post-order sequence
//            for (auto nd : boost::adaptors::reverse(preorder)) {
//                assert(nd->_partial == nullptr);
//                if (nd->_left_child) {
//                    assert(nd->_left_child->_partial != nullptr);
//                    assert(nd->_left_child->_right_sib->_partial != nullptr);
//                    nd->_partial = pullPartial();
//                    calcPartialArrayJC(nd, nd->_left_child, nd->_left_child->_right_sib);
//                }
//                else {
//                    nd->_partial = _leaf_partials[_index][nd->_number];
//                }
//            }
//        }
//    }
//#endif

#if defined(LAZY_COPYING)
    inline PartialStore::partial_t Forest::pullPartial() {
    #if defined(USING_MULTITHREADING)
        lock_guard<mutex> guard(mutex);
    #endif
        assert(_index >= 0);
        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
        ptr = ps.getPartial(_npatterns*4, _index);
        return ptr;
    }
#endif

#if defined (LAZY_COPYING)
    inline void Forest::joinLineagePair(Node * new_nd, Node * subtree1, Node * subtree2) {
        // Assumes pullNode has been called to obtain new_nd before calling this function
        // and that new_nd already has correct _height, _edge_length, and _species.
        
        // Perform sanity checks
        assert(new_nd);
        assert(subtree1);
        assert(subtree2);
        assert(new_nd->_number > -1);
        assert(new_nd->_edge_length == 0.0);
        assert(new_nd->_species == (subtree1->_species | subtree2->_species));

        // Check whether the left and right subtrees imply the same ancestral node height
        double h1 = subtree1->_height + subtree1->_edge_length;
        double h2 = subtree2->_height + subtree2->_edge_length;
        assert(fabs(h1 - h2) < G::_small_enough);

        // Check whether height of the new node is consistent with average of h1 and h2
        double havg = (h1 + h2)/2.0;
        assert(fabs(havg - new_nd->_height) < G::_small_enough);

        new_nd->_name = "anc-" + to_string(new_nd->_number);
        new_nd->_left_child  = subtree1;
                
        subtree1->_right_sib = subtree2;
        subtree1->_parent    = new_nd;
        subtree2->_parent    = new_nd;
    }
#endif

}

