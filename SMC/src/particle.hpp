#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace boost;

#include "lot.hpp"

extern proj::Lot rng;

namespace proj {

class Particle {
    public:

        Particle();
        Particle(const Particle & other);
        ~Particle();
        typedef std::shared_ptr<Particle>               SharedPtr;


        void                                    debugParticle(std::string name);
        void                                    showParticle();
        double                                  proposal(bool gene_trees_only, bool deconstruct, vector<pair<tuple<string, string, string>, double>> species_joined);
        void                                    setData(Data::SharedPtr d, map<string, string> &taxon_map, int index) {
                                                    _nsubsets = d->getNumSubsets();
                                                    _data = d;
                                                    _forest.setData(d, index, taxon_map);
                                                }
        void                                    resetGeneTreePartials(Data::SharedPtr d, map<string, string> taxon_map, int i);
        void                                    mapSpecies(map<string, string> &taxon_map, vector<string> &species_names, int n);
        void                                    mapGeneTrees(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        void                                    setLogLikelihood(double log_likelihood);
        void                                    setLogCoalescentLikelihood(double log_coal_likelihood);
        void                                    setParticleGeneration(int n);
        double                                  calcHeight();
        void                                    buildEntireSpeciesTree();
        string                                  getSpeciesNewick() {return _forest.makeNewick(9, true);}
        double                                  getLogWeight(string a);
        void                                    setLogWeight(double w, string a);
        void                                    operator=(const Particle & other);
        std::string                             saveForestNewick() {
            return _forest.makeNewick(25, true);
        }
        bool operator<(const Particle::SharedPtr & other) const {
            return _log_weight<other->_log_weight;
        }

        bool operator>(const Particle::SharedPtr & other) const {
            return _log_weight>other->_log_weight;
        }

        static void                                     setNumSubsets(unsigned n);
        void                                            showSpeciesIncrement();
        void                                            showSpeciesTree();
        void                                            showHybridNodes();
        string                                          saveHybridNodes();
        double                                          saveParticleWeights();
        double                                          saveParticleLikelihoods();
        void                                            showGamma();
        string                                          saveGamma();
        void                                            calculateGamma();
        void                                            setRunOnEmpty(bool a) {_running_on_empty = a;}
        vector<double>                                  getBranchLengths();
        vector<double>                                  getBranchLengthPriors();
        vector<double>                                  getGeneTreeLogLikelihoods();
        vector<double>                                  getGeneTreeLogCoalescentLikelihood();
        vector<double>                                  getTopologyPriors();
        void                                            hybridizationProposal();
        void                                            calcParticleWeight();
        void                                            speciesJoinedProposal();
        vector<string>                                  getGeneTreeNewicks();
        vector<pair<tuple<string, string, string>, double>>                    getSpeciesJoined(){return _t;}
        double                                          getCoalescentLikelihood(){return _log_coalescent_likelihood;}
        void                                            geneTreeProposal(bool deconstruct, vector<pair<tuple<string, string, string>, double>> species_joined);
        void                                            speciesProposal(vector<double> max_depths, tuple<string, string, string> species_joined);
        void                                            resetSpeciesInfo(){_t.clear();}
        void                                            resetSpeciesTreeHeight(){ _species_tree_height = 0.0;}
        void                                            resetSpecies();
        void                                            resetGeneIncrements();
        void                                            processSpeciesNewick(vector<string> newicks, bool topology_only);
        void                                            processGeneNewicks(vector<string> newicks, int gene_number);
        void                                            remakeGeneTrees(map<string, string> &taxon_map) ;
        vector<pair<double, pair<string, string>>>      getMinDepth(){return _forest.getMinDepths();}
        void                                            calcGeneTreeMinDepth();
        void                                            refreshGeneTreePreorder();
        tuple<string, string, string>                   speciesTopologyProposal();
        double                                          calcConstrainedProposal(tuple<string, string, string> species_joined);
        double                                          calcGeneCoalescentLikelihood(double last_edge_len, tuple<string, string, string> species_joined, double species_tree_height);
        double                                          getSpeciesTreeHeight();
        double                                          getLastEdgeLen();
        void                                            calcSpeciesParticleWeight(double log_coalescent_likelihood);
        void                                            drawHeightsFromPrior();
        void                                            resetLogTopologyPrior(){_forest._log_joining_prob = 0.0;}
        double                                          calcCoalLikeForNewTheta(double proposed_theta, vector<pair<tuple<string, string, string>, double>> species_info, bool both);
        void                                            buildEntireGeneTree();
        void                                            sampleGeneTreePrior();
        double                                          calcLogSpeciesTreeDensityGivenLambda(double lambda);
    
    private:

        static unsigned                         _nsubsets;
        Forest                                  _forest;
        double                                  _log_weight;
        double                                  _species_log_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        double                                  _log_coalescent_likelihood;
        int                                     _generation = 0;
        bool                                    _running_on_empty = false;
        double                                  _species_tree_height;
        vector<pair<tuple<string, string, string>, double>> _t;
        bool                                    _inf = false;
        string                                  _name;
};

    inline Particle::Particle() {
        //log weight and log likelihood are 0 for first generation
        _log_weight = 0.0;
        _log_likelihood = 0.0;
        _log_coalescent_likelihood = 0.0;
        _species_tree_height = 0.0;
        _species_tree_height = 0.0;
        _species_log_weight = 0.0;
    };

    inline Particle::~Particle() {
//        cout << "destroying a particle" << endl;
    }

    inline void Particle::showParticle() {
        //print out weight of each particle
        cout << "\nParticle:\n";
        double log_weight = _log_weight;
        if (_name == "species") {
            log_weight = _species_log_weight;
        }
        cout << "  _log_weight: " << log_weight << "\n" ;
        if (_name != "species") {
            cout << " _log_likelihood: " << _log_likelihood << "\n";
        }
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forest.showForest();
    }

    inline void Particle::showSpeciesTree() {
        assert (_name == "species");
        //print out weight of each particle
        cout << "\nParticle:\n";
        cout << "  _log_weight: " << _species_log_weight << "\n" ;
        cout << " _log_coalescent_likelihood: " << _log_coalescent_likelihood << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forest.showForest();
    }

    //more detailed version of showParticle
    inline void Particle::debugParticle(std::string name) {
        cout << "debugging particle" << endl;
        //print out weight of each particle
        cout << "\nParticle " << name << ":\n";
        cout << "  _log_weight:               " << _log_weight                 << "\n" ;
        cout << "  _log_likelihood:           " << _log_likelihood             << "\n";
        cout << "  _forest._nleaves:          " << _forest._nleaves            << "\n";
        cout << "  _forest._ninternals:       " << _forest._ninternals         << "\n";
        cout << "  _forest._npatterns:        " << _forest._npatterns          << "\n";
        cout << "  _forest._nstates:          " << _forest._nstates            << "\n";
        cout << "  _forest._last_edge_length: " << _forest._last_edge_length   << "\n";
        cout << "  newick description:        " << _forest.makeNewick(5,false) << "\n";
    }

    inline double Particle::calcLogLikelihood() {
        //calculate likelihood for each gene tree
        assert (_name != "species");
        
        double log_likelihood = 0.0;
        double gene_tree_log_likelihood = 0.0;
        
        gene_tree_log_likelihood = _forest.calcLogLikelihood();
        assert(!isnan (log_likelihood));
        assert(!isnan (gene_tree_log_likelihood));
        //total log likelihood is sum of gene tree log likelihoods
        log_likelihood += gene_tree_log_likelihood;
        log_likelihood += _forest._gene_tree_log_coalescent_likelihood;

        _generation++;
        return log_likelihood;
    }

    inline void Particle::remakeGeneTrees(map<string, string> &taxon_map) {
        assert (_name != "species");
        _forest.clear();
        _forest.remakeGeneTree(taxon_map);
    }

    inline double Particle::proposal(bool gene_trees_only, bool deconstruct, vector<pair<tuple<string, string, string>, double>> species_joined) {
        // this function proposes gene trees, not species trees
        
        string event;
        
        if (_generation == 0 || _generation == Forest::_ntaxa - 1) {
            _forest._nincrements = 0;
            _forest._gene_tree_log_coalescent_likelihood = 0.0;
        }
        if (gene_trees_only) {
            geneTreeProposal(deconstruct, species_joined);
        }
        
        if (_running_on_empty) {
            _generation++;
            _log_weight = 0.0;
        }
        
        return _log_weight;
    }

    inline void Particle::geneTreeProposal(bool deconstruct, vector<pair<tuple<string, string, string>, double>> species_joined) {
        assert(_name != "species");
        
        if (deconstruct) {
            _forest.deconstructGeneTree();
        }
        
        pair<double, string> species_info = _forest.chooseDelta(species_joined);
        _forest.geneTreeProposal(species_info, species_joined);
        if (_forest._lineages.size() == 1) {
            _forest.refreshPreorder();
        }

        if (!_running_on_empty) {
            double prev_log_likelihood = _log_likelihood;
            _log_likelihood = _forest._gene_tree_log_likelihood + _forest._gene_tree_log_coalescent_likelihood;

            if (Forest::_proposal == "prior-prior") {
                _log_weight = _log_likelihood - prev_log_likelihood;
            }
            else {
                calcParticleWeight();
                _generation++;
            }
        }
    }

    inline double Particle::getLogWeight(string a) {
        if (a == "s") {
            return _species_log_weight;
        }
        else {
            return _log_weight;
        }
    }

    inline void Particle::setLogWeight(double w, string a) {
        if (a == "s") {
            _species_log_weight = w;
        }
        else {
            _log_weight = w;
        }
    }

    inline tuple<string, string, string> Particle::speciesTopologyProposal() {
        assert (_name == "species");
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        
        if (_forest._last_edge_length > 0.0) {
        // choose species to join if past the first species generation for each forest vector
            species_joined = _forest.speciesTreeProposal();
        }
        
        _t.push_back(make_pair(species_joined, 0.0));
        
        return species_joined;
    }

    inline double Particle::calcConstrainedProposal(tuple<string, string, string> species_joined) {
        assert (!_inf);
        
        _forest.updateSpeciesPartition(species_joined);
        
        string species1 = get<0>(species_joined);
        string species2 = get<1>(species_joined);
        
        double max_depth = 0.0;
        if (_forest._species_partition.size() > 1) {
            _forest.resetDepthVector(species_joined);
            max_depth = _forest.getMinDepths()[0].first;
            assert (max_depth > 0.0);
        }
        
        return max_depth;
    }

    inline double Particle::getSpeciesTreeHeight() {
        assert (_name == "species");
        return _species_tree_height;
    }

    inline double Particle::getLastEdgeLen() {
        assert (_name == "species");
        return _forest._last_edge_length;
    }

    inline double Particle::calcGeneCoalescentLikelihood(double last_edge_len, tuple<string, string, string> species_joined, double species_tree_height) {
        // reset panmictic part of coalescent likelihood to 0
        _forest._panmictic_coalescent_likelihood = 0.0;
        
        double coal_like_increment = _forest.calcCoalescentLikelihood(last_edge_len, species_joined, species_tree_height);
        _log_coalescent_likelihood += coal_like_increment;
        
        return _log_coalescent_likelihood + _forest._panmictic_coalescent_likelihood;
    }


    inline void Particle::speciesProposal(vector<double> max_depth_vector, tuple<string, string, string> species_joined ) {
        assert (!_inf);
        
        double max_depth = 0.0;

        if (_forest._lineages.size() > 1) {
            max_depth = *min_element(max_depth_vector.begin(), max_depth_vector.end());
            max_depth -= _forest.getTreeHeight();
            // choose a species tree increment
        }
        
        if (_forest._lineages.size() > 1) {
            assert (max_depth > 0.0);
            _forest.chooseSpeciesIncrement(max_depth);
            _species_tree_height += _forest._last_edge_length;
        }
        
        if (_forest._lineages.size() == 1) {
            _forest._last_edge_length = 0.0;
        }
        
        _t[_forest._species_join_number].second = _forest._last_edge_length;
        _forest._species_join_number++;
        
        }

    inline void Particle::calcSpeciesParticleWeight(double log_coalescent_likelihood) {
        double prev_log_coalescent_likelihood = _log_coalescent_likelihood;
        _log_coalescent_likelihood = log_coalescent_likelihood;
        
        if (!_running_on_empty) {
            _species_log_weight = log_coalescent_likelihood - prev_log_coalescent_likelihood;
            double test = 1 / _species_log_weight;
            if (test == -0) {
                _inf = true;
            }
            assert (!_inf);
        }
    }

    inline void Particle::calcParticleWeight() {
        // particle weight is the forest weight
        _log_weight = _forest._gene_tree_log_weight;
    }

    inline Particle::Particle(const Particle & other) {
        *this = other;
    }

    inline void Particle::calculateGamma() {
        double major = 0.0;
        if (_forest._last_direction == "major") {
            major++;
        }
        double gamma = major;
        _forest._gamma.push_back(gamma);
    }

    inline void Particle::saveForest(std::string treefilename)  {
        ofstream treef(treefilename);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree test = [&R] " << _forest.makeNewick(8, true)  << ";\n";
        treef << "end;\n";
        treef.close();
    }

    inline void Particle::savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const {
        ofstream paupf(paupfilename);
        paupf << "#nexus\n\n";
        paupf << "begin paup;\n";
        paupf << "exe " << datafilename << ";\n";
        paupf << "set crit=like forcepoly;\n";
        paupf << "lset nst=1 basefreq=equal rates=equal pinvar=0;\n";
        paupf << "gettrees file=" << treefilename << " storebrlen;\n";
        paupf << "lscores 1 / userbrlen;\n";
        paupf << "[!expected lnL = " << expected_lnL << "]\n";
        paupf << "end;\n";
        paupf.close();
    }

    inline void Particle::processSpeciesNewick(vector<string> newicks, bool topology_only) {
        assert (_name == "species");
        
        if (newicks.size() > 1) {
            throw XProj(boost::str(boost::format("only one species tree newick may be specified")));
        }
        
        if (newicks.size() > 0) {
            string newick = newicks[0];
                _t = _forest.buildFromNewickTopology(newick, topology_only);
            if (topology_only) {
                drawHeightsFromPrior();
            }
        }
        
        else {
            buildEntireSpeciesTree();
        }
    }

    inline void Particle::drawHeightsFromPrior() {
        assert (_name == "species");
        
        _forest.resetIncrements();
        
       vector<string> existing_lineages = _forest.setUpExistingLineagesVector();
        
        for (int a=0; a<Forest::_nspecies-1; a++) {
            existing_lineages = _forest.updateExistingLineagesVector(existing_lineages, _t[a].first);
            _forest.chooseSpeciesIncrementFromNewick(existing_lineages);
            _t[a].second = _forest._last_edge_length;
        }
        
        assert (_t[Forest::_nspecies-1].second == 0); // last element of _t should not draw a branch length
    }

    inline void Particle::processGeneNewicks(vector<string> newicks, int gene_number) {
        assert (_name != "species");
        _forest.buildFromNewick(newicks[gene_number], true, false);
        _forest.refreshPreorder();
    }

    inline void Particle::buildEntireGeneTree() {
        _forest.drawFromGeneTreePrior();
    }

    inline void Particle::sampleGeneTreePrior() {
        double prev_log_likelihood = _forest._gene_tree_log_likelihood;
        
        _forest.drawFromGeneTreePrior();
        if (Forest::_proposal == "prior-prior") {
            _log_weight = _forest._gene_tree_log_likelihood - prev_log_likelihood;
        }
        else {
            _log_weight = _forest._gene_tree_log_weight;
        }
    }

    inline double Particle::calcHeight() {
        assert (_name == "species");
        //species tree
        double sum_height = 0.0;

        // calculate height of lineage
        Node* base_node = _forest._lineages[0];
        sum_height += base_node->getEdgeLength();
        for (Node* child=base_node->_left_child; child; child=child->_left_child) {
            sum_height += child->getEdgeLength();
        }
        return sum_height;
    }

    inline void Particle::hybridizationProposal() {
        vector<string> hybridized_nodes = _forest.hybridizeSpecies();
        if (_forest._lineages.size()>1) {
            _forest.addSpeciesIncrement();
        }
        _forest.hybridizeGene(hybridized_nodes, _forest._last_edge_length);
        calculateGamma();
    }

    inline void Particle::setNumSubsets(unsigned n) {
        _nsubsets = n;
    }

    inline void Particle::mapSpecies(map<string, string> &taxon_map, vector<string> &species_names, int n) {
        _forest.setIndex(n);
        
        //species tree
        if (n == 0) {
            _name = "species";
            _forest.setUpSpeciesForest(species_names);
        }
        else {
            //gene tree
            string gene_index = to_string(n);
            _name = "gene" + gene_index;
            _forest.setUpGeneForest(taxon_map);
        }
    }

    inline void Particle::mapGeneTrees(map<string, string> &taxon_map, vector<string> &species_names) {
        //gene trees
        assert (_name != "species");
        _forest.setUpGeneForest(taxon_map);
    }

    inline void Particle::showSpeciesIncrement(){
        assert (_name == "species");
        cout << "species tree increment: " << "     " << _forest._last_edge_length << endl;
    }

    inline void Particle::showHybridNodes() {
        cout << "particle" << endl;
        showGamma();
        for (auto &nd:_forest._nodes) {
            if (nd._major_parent) {
                cout << "       " << "hybridized node is: " << nd._name << " with minor parent " << nd._minor_parent->_name << " and major parent " << nd._major_parent->_name << endl;
            }
        }
    }

    inline string Particle::saveHybridNodes() {
        string nodes = "";
        int i = 0;
        for (auto &nd:_forest._nodes) {
            if (nd._major_parent) {
                string gammastr = to_string(_forest._gamma[i]);
                nodes +=  "hybridized node is: " + nd._name + " with minor parent " + nd._minor_parent->_name + " and major parent " + nd._major_parent->_name + "\n" + "gamma is: " + gammastr + "\n";
                i++;
            }
        }
        return nodes;
    }

    inline double Particle::saveParticleWeights() {
        return _log_weight;
    }

    inline double Particle::saveParticleLikelihoods() {
        return _log_likelihood;
    }

    inline void Particle::showGamma() {
        if (_forest._gamma.size() > 0) {
            cout << "   " << "gamma is: " << endl;
            for (auto &g:_forest._gamma) {
                cout << g << "   ";
            }
            cout << "\n";
        }
    }

    inline void Particle::setLogLikelihood(double log_likelihood) {
        _log_likelihood = log_likelihood;
    }

    inline void Particle::setLogCoalescentLikelihood(double log_coal_likelihood) {
        _log_coalescent_likelihood = log_coal_likelihood;
    }

    inline void Particle::setParticleGeneration(int n) {
        _generation = n;
    }

    inline vector<double> Particle::getBranchLengths() {
        vector<double> divergence_times;
        for (auto &b:_forest._increments) {
            divergence_times.push_back(b.first);
        }
        return divergence_times;
    }

    inline vector<double> Particle::getBranchLengthPriors() {
        vector<double> priors;
        for (auto &b:_forest._increments) {
            priors.push_back(b.second);
        }
        return priors;
    }

    inline vector<double> Particle::getGeneTreeLogLikelihoods() {
        vector<double> gene_tree_log_likelihoods;
        gene_tree_log_likelihoods.push_back(_forest._gene_tree_log_likelihood);
        return gene_tree_log_likelihoods;
    }

    inline vector<double> Particle::getGeneTreeLogCoalescentLikelihood() {
        vector<double> gene_tree_log_coalescent_likelihoods;
        gene_tree_log_coalescent_likelihoods.push_back(_forest._gene_tree_log_coalescent_likelihood);
        return gene_tree_log_coalescent_likelihoods;
    }

    inline vector<double> Particle::getTopologyPriors() {
        // calculate species tree topology probability
        vector<double> topology_priors;
        topology_priors.push_back(_forest._log_joining_prob);
        return topology_priors;
    }

    inline vector<string> Particle::getGeneTreeNewicks() {
        vector<string> newicks;
        string newick = _forest.makeNewick(9, true);
        newicks.push_back(newick);
        return newicks;
    }

    inline void Particle::resetSpecies() {
        setLogLikelihood(0.0);
        setLogWeight(0.0, "s");
        resetSpeciesInfo();
        resetSpeciesTreeHeight();
        if (_name == "species") {
            _forest._increments.clear();
        }
        _forest._log_joining_prob = 0.0;
        _log_coalescent_likelihood = 0.0;
    }

    inline void Particle::resetGeneIncrements() {
        assert (_name != "species");
        _forest._increments.clear();
    }

    inline double Particle::calcCoalLikeForNewTheta(double proposed_theta,  vector<pair<tuple<string, string, string>, double>> species_info, bool both) {
        return _forest.calcLogCoalLikeGivenTheta(proposed_theta, species_info, both);
    }

    inline void Particle::buildEntireSpeciesTree() {
        assert (_name == "species");
        
        double max_depth = 0.0;
        _forest.chooseSpeciesIncrement(max_depth);
        double edge_len = _forest._last_edge_length;
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _t.push_back(make_pair(species_joined, edge_len));

        for (int i=0; i < _forest._nspecies-1; i++) {
            if (_forest._lineages.size() > 1) {
                species_joined = _forest.speciesTreeProposal();
                
                // if the species tree is not finished, add another species increment
                if (_forest._lineages.size()>1) {
                    _forest.addSpeciesIncrement();
                }
                
                double edge_len = 0.0;
                if (_forest._lineages.size() > 1) {
                    edge_len = _forest._last_edge_length;
                }
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::resetGeneTreePartials(Data::SharedPtr d, map<string, string> taxon_map, int i) {
        assert (_name != "species");
        _nsubsets = d->getNumSubsets();
        _data = d;
        _forest.setData(d, i, taxon_map);
    }

    inline void Particle::calcGeneTreeMinDepth() {
        assert (_name != "species");
        _forest.calcMinDepth();
    }

    inline void Particle::refreshGeneTreePreorder() {
        assert (_name != "species");
        _forest.refreshPreorder();
    }

    inline double Particle::calcLogSpeciesTreeDensityGivenLambda(double lambda) {
        assert (_name == "species");
        
         double prev_lambda = Forest::_lambda;
         Forest::_lambda = lambda;
         
         double log_density = _forest.calcLogSpeciesTreeDensity(lambda);
         
         Forest::_lambda = prev_lambda;
        
         return log_density;
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forest         = other._forest;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _t              = other._t;
        _running_on_empty = other._running_on_empty;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _species_tree_height = other._species_tree_height;
        _inf = other._inf;
        _species_log_weight = other._species_log_weight;
        _name = other._name;
    };
}
