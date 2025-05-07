#pragma once

#if defined(LAZY_COPYING)

namespace proj {

    class ForestExtension {
        public:
            
                            ForestExtension();
                                
            void            dock(const Forest::SharedPtr gf, PartialStore::partial_t partial, Lot::SharedPtr lot);
            void            dockSim(const Forest::SharedPtr gf, Lot::SharedPtr lot);
            void            undock();

            G::uint_pair_t  chooseNodesToJoin(const vector<unsigned> & node_indices) const;
            G::species_t    chooseSpecies(double total_rate, double theta) const;
            unsigned        speciesCount(G::species_t spp) const;
            double          getProposedDelta() const;
            double          getHeight() const;
            double          getLogWeight() const;
            const Node *    getProposedAnc() const;
            const Node *    getProposedLChild() const;
            const Node *    getProposedRChild() const;
            
            const G::merge_vect_t & getMergers() const;
            
            double          calcTotalRate(double theta);
            void            mergeSpecies(G::species_t left_spp, G::species_t right_spp);
            void            addIncrement(double increment);
            void            coalesce(double total_rate, G::species_t species_name);
            void            debugCheckSpeciesVect() const;
            unsigned        getDeepCoal(tuple <G::species_t, G::species_t, G::species_t> species_joined);
            unsigned        getMaxDeepCoal(tuple <G::species_t, G::species_t, G::species_t> species_joined);
            void            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
        
            inline vector<pair<double, unsigned long>> calcForestRate(Lot::SharedPtr lot, unordered_map<G::species_t, double> theta_map);
                    
            PartialStore::partial_t getExtensionPartial();
        
            unsigned        getSpeciesPartitionSize(){return (unsigned) _species_partition.size();}

        private:
        
            // key: species ==> value: vector of indices into docked forest _lineages
            typedef unordered_map<G::species_t, vector<unsigned> > sppmap_t;

            Forest::ConstSharedPtr        _docked_gene_forest;
            double                        _log_weight;
            double                        _proposed_delta;
            Node                          _proposed_anc;
            const Node *                  _proposed_lchild;
            const Node *                  _proposed_rchild;
            vector<G::species_t>          _species_vect;
            G::merge_vect_t               _mergers;
            sppmap_t                      _species_partition;
            Lot::SharedPtr                _lot;
            vector<pair<G::species_t, unsigned>>  _lineages_per_species;
    };
    
    inline ForestExtension::ForestExtension() {
        undock();
    }
    
    inline void ForestExtension::dock(const Forest::SharedPtr gf, PartialStore::partial_t partial, Lot::SharedPtr lot) {
        // Check to make sure this extension was previously undocked
        assert(gf);
        assert(_docked_gene_forest == nullptr);
        assert(_species_vect.empty());
        
        // Reset the random number generator
        _lot = lot;

        // Attach Forest
        _docked_gene_forest = gf;
        
        // Create vector of species corresponding to lineages in _docked_gene_forest
        _docked_gene_forest->copyLineageSpecies(_species_vect);
        
        _proposed_delta = 0.0;
        _proposed_anc._height = gf->getForestHeight();
        _proposed_anc._edge_length = 0.0;
        _proposed_anc._partial = partial;
    }

    inline void ForestExtension::dockSim(const Forest::SharedPtr gf, Lot::SharedPtr lot) {
        // Check to make sure this extension was previously undocked
        assert(gf);
        assert(_docked_gene_forest == nullptr);
        assert(_species_vect.empty());
        
        // Reset the random number generator
        _lot = lot;

        // Attach Forest
        _docked_gene_forest = gf;
        
        // Create vector of species corresponding to lineages in _docked_gene_forest
        _docked_gene_forest->copyLineageSpecies(_species_vect);
        
        _proposed_delta = 0.0;
        _proposed_anc._height = gf->getForestHeight();
        _proposed_anc._edge_length = 0.0;
    }
    
    inline void ForestExtension::undock() {
        _lot.reset();
        _docked_gene_forest.reset();
        _species_vect.clear();
        _mergers.clear();
        _species_partition.clear();
        _log_weight = 0.0;
        _proposed_delta = 0.0;
        _proposed_anc.clearPointers();
        _proposed_anc._number = -2;
        _proposed_anc._name = "fake";
        _proposed_anc._height = 0.0;
        _proposed_anc._partial.reset();
        _proposed_anc._edge_length = 0.0;
        _proposed_anc._flags = 0;
        _proposed_anc._species = 0;
        _proposed_anc._split.clear();
        _proposed_lchild = nullptr;
        _proposed_rchild = nullptr;
    }

    inline unsigned ForestExtension::speciesCount(G::species_t spp) const {
        return (unsigned)count(_species_vect.begin(), _species_vect.end(), spp);
    }

    inline const Node * ForestExtension::getProposedAnc() const {
        return &_proposed_anc;
    }
    
    inline const Node * ForestExtension::getProposedLChild() const {
        return _proposed_lchild;
    }
    
    inline const Node * ForestExtension::getProposedRChild() const {
        return _proposed_rchild;
    }
    
    inline double ForestExtension::getProposedDelta() const {
        return _proposed_delta;
    }
    
    inline double ForestExtension::getHeight() const {
        return _proposed_anc._height;
    }
    
    inline double ForestExtension::getLogWeight() const {
        return _log_weight;
    }
    
    inline const G::merge_vect_t & ForestExtension::getMergers() const {
        return _mergers;
    }

    inline void ForestExtension::addIncrement(double increment) {
        assert (increment > 0.0);
        _proposed_delta += increment;
        _proposed_anc._height += increment;
    }
    
    inline void ForestExtension::mergeSpecies(G::species_t left_spp, G::species_t right_spp) {
        // Assumes height of _proposed_anc is the level at which the merger occurs
        
        // Add element to _mergers vector
        _mergers.push_back(make_tuple(_proposed_anc._height, left_spp, right_spp));
            
        // Merge species in _species_vect
        G::species_t anc_spp = (left_spp | right_spp);
        for (auto & s : _species_vect) {
            if (s == left_spp || s == right_spp) {
                s = anc_spp;
            }
        }
    }
    
    inline double ForestExtension::calcTotalRate(double theta) {
        // Enumerate lineages belonging to each current species
        _species_partition.clear();
        unsigned i = 0;
        for (auto s : _species_vect) {
            _species_partition[s].push_back(i++);
        }
        
        // Calculate total rate and return it
        double total_rate = 0.0;
        for (auto & sv : _species_partition) {
            unsigned n = (unsigned)sv.second.size();
            double r = 1.0*n*(n-1)/theta;
            total_rate += r;
        }
        return total_rate;
    }
    
    inline G::species_t ForestExtension::chooseSpecies(double total_rate, double theta) const {
        unsigned nspecies = (unsigned)_species_partition.size();
        assert(nspecies > 0);
        
        vector<double> probs(nspecies, 0.0);
        vector<G::species_t> species_vect;
        species_vect.reserve(nspecies);
        
        // If n_i = number of lineages in species i, the probability of coalescence is
        //   r_i / total_rate = (n_i*(n_i-1)/theta) / (sum_i r_i)
        unsigned i = 0;
        for (auto & sv : _species_partition) {
            species_vect.push_back(sv.first);
            double n = (double)sv.second.size();
            probs[i++] = n*(n - 1.0)/(theta*total_rate);
        }
        
        // Choose a species using a multinomial draw from probs
        unsigned which = G::multinomialDraw(_lot, probs);
        assert(which < probs.size());
        G::species_t spp = species_vect[which];
        return spp;
    }
    
    inline G::uint_pair_t ForestExtension::chooseNodesToJoin(const vector<unsigned>  & node_indices) const {
        double nsubtrees = node_indices.size();
        assert (nsubtrees>1);
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        
        if (nsubtrees > 2) {
            t1 = _lot->randint(0, nsubtrees-1);
            t2 = _lot->randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = _lot->randint(0, nsubtrees-1);
            }
        }
        assert(t1 < nsubtrees);
        assert (t2 < nsubtrees);

        return make_pair(t1, t2);
    }

    inline void ForestExtension::coalesce(double total_rate, G::species_t species_name) {
        _proposed_anc.setSpecies(species_name);
                    
        // Choose the two nodes to join
         G::uint_pair_t chosen = chooseNodesToJoin(_species_partition[species_name]);
        
        // Get pointers to the two nodes to join
        _proposed_lchild = _docked_gene_forest->_lineages[_species_partition[species_name][chosen.first]];
        _proposed_rchild = _docked_gene_forest->_lineages[_species_partition[species_name][chosen.second]];
        
        // Set _proposed_anc's split to union of the two child splits
        _proposed_anc._split.resize(G::_ntaxa);
        _proposed_anc._split += _proposed_lchild->_split;
        _proposed_anc._split += _proposed_rchild->_split;

        // Compute partial likelihood array of ancestral node
        _log_weight = _docked_gene_forest->calcPartialArrayJC(&_proposed_anc, _proposed_lchild, _proposed_rchild);

        assert(!isnan(_log_weight));
        assert(!isinf(_log_weight));
    }
    
    inline PartialStore::partial_t ForestExtension::getExtensionPartial() {
        return _proposed_anc._partial;
    }

    inline void ForestExtension::debugCheckSpeciesVect() const {
        unsigned n = (unsigned)_species_vect.size();
        assert(_docked_gene_forest->_lineages.size() == n);
        for (unsigned i = 0; i < n; i++) {
            assert(_docked_gene_forest->_lineages[i]->_species == _species_vect[i]);
        }
    }

    inline unsigned ForestExtension::getDeepCoal(tuple <G::species_t, G::species_t, G::species_t> species_joined) {
        unsigned num_deep_coal = 0;
        
        G::species_t spp1 = get<0> (species_joined);
        G::species_t spp2 = get<1> (species_joined);
        
        unsigned nlineages1 = (unsigned) _species_partition[spp1].size();
        unsigned nlineages2 = (unsigned) _species_partition[spp2].size();
        
        num_deep_coal += nlineages1 - 1;

        num_deep_coal += nlineages2 - 1;
        
        return num_deep_coal;
    }

    inline unsigned ForestExtension::getMaxDeepCoal(tuple <G::species_t, G::species_t, G::species_t> species_joined) {
        G::species_t species1 = get<0>(species_joined);
        G::species_t species2 = get<1>(species_joined);
        G::species_t species3 = get<2>(species_joined);
        
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

    inline void ForestExtension::setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species) {
        unsigned i = 0;
        for (auto &s:G::_species_names_typed) {
            _lineages_per_species.push_back(make_pair(s, ntaxa_per_species[i]));
            i++;
        }
    }


#if defined (LAZY_COPYING)
    inline vector<pair<double, unsigned long>> ForestExtension::calcForestRate(Lot::SharedPtr lot, unordered_map<G::species_t, double> theta_map) {
#else
    inline vector<pair<double, unsigned long>> ForestExtension::calcForestRate(Lot::SharedPtr lot, unordered_map<string, double> theta_map) {
#endif
        vector<pair<double, unsigned long>> rates;
        pair<double, unsigned long> rate_and_name;
        
        _species_partition.clear();
        unsigned i = 0;
        for (auto &s : _species_vect) {
            _species_partition[s].push_back(i++);
        }

        for (auto &s:_species_partition) {
            if (s.second.size() > 1) { // if size == 0, no possibility of coalescence and rate is 0
                double population_coalescence_rate = 0.0;
    #if defined (DRAW_NEW_THETA)
                    double population_theta = theta_map[s.first];
                population_coalescence_rate = s.second.size()*(s.second.size()-1)/population_theta;
    #else
                population_coalescence_rate = s.second.size()*(s.second.size()-1)/G::_theta;
    #endif
                G::species_t name = s.first;
                rate_and_name = make_pair(population_coalescence_rate, name);
                rates.push_back(rate_and_name);
            }
        }
        return rates;
    }
 }
#endif


