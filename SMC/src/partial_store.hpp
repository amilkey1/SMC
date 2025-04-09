#pragma once

namespace proj {

#if defined (REUSE_PARTIALS)
    struct Partial {
        Partial(unsigned g, unsigned n);
        ~Partial();
        unsigned        _g; // the locus
        vector<double>  _v; // the partial array: length = _nstates*<no. patterns>
    };
#endif

#if defined (REUSE_PARTIALS)
    inline Partial::Partial(unsigned g, unsigned n) {
        _g = g;
        _v.resize(n);
    }
#endif

#if defined (REUSE_PARTIALS)
    inline Partial::~Partial() {
    }
#endif

    class PartialStore {
        public:
            PartialStore();
            ~PartialStore();
#if defined (REUSE_PARTIALS)
        typedef std::shared_ptr<Partial>          partial_t;
#else
        typedef std::shared_ptr <vector<double>>   partial_t;
#endif
        
#if defined (REUSE_PARTIALS)
        void            setNLoci(unsigned nloci);
        partial_t       getPartial(unsigned locus);
        void            putPartial(unsigned locus, partial_t partial);
        void            setNElements(unsigned nelements, unsigned locus);
        unsigned        getNElements(unsigned locus) const {return _nelements[locus-1];}

        
        typedef vector<partial_t>            vect_partial_t;
        typedef vector<vect_partial_t>       storage_t;
        typedef vector<vect_partial_t>       leaf_partials_t;
        
        unsigned         _total_partials_created;
        unsigned         _total_elements_created;
        unsigned         _total_partials_reused;
        unsigned         _total_elements_reused;
#endif
        
#if defined (REUSE_PARTIALS)
        PartialStore::partial_t getPartial(unsigned nelements, unsigned locus);
#else
        partial_t getPartial(unsigned nelements);
#endif
        
#if !defined (REUSE_PARTIALS)
        void setnelements(unsigned nelements) {_nelements = nelements;}
#endif
        
        void resetPartial();
        private:
#if defined (REUSE_PARTIALS)
            storage_t        _storage;
            vector<unsigned> _nelements;
#else
        unsigned _nelements;
#endif
    };

    inline PartialStore::PartialStore() {
        //std::cout << "Constructing a partial store" << std::endl;
#if defined (REUSE_PARTIALS)
        _total_partials_created = 0;
        _total_elements_created = 0;
        _total_partials_reused = 0;
        _total_elements_reused = 0;
#else
        _nelements=0;
#endif
    }

    inline PartialStore::~PartialStore() {
        //std::cout << "Destroying a partial store" << std::endl;
#if defined (REUSE_PARTIALS)
        _nelements.clear();
        _storage.clear();
#endif
    }

#if defined (REUSE_PARTIALS)
    inline void PartialStore::setNLoci(unsigned nloci) {
        // Should be called before any partials are stored
        assert(_nelements.empty());
        assert(_storage.empty());
        
        // Resize both containers
        _nelements.resize(nloci);
        _storage.resize(nloci);
    }
#endif

#if defined (REUSE_PARTIALS)
    inline void PartialStore::setNElements(unsigned nelements, unsigned locus) {
        assert(_nelements.size() > locus-1);
        _nelements[locus-1] = nelements;
    }
#endif

#if defined (REUSE_PARTIALS)
    inline PartialStore::partial_t PartialStore::getPartial(unsigned nelements, unsigned locus) {
#else
    inline PartialStore::partial_t PartialStore::getPartial(unsigned nelements) {
#endif
#if defined(REUSE_PARTIALS)
        // Check to make sure supplied value of locus is valid
        assert(_nelements.size() > locus-1);
        assert(_nelements[locus-1] > 0);
        
        partial_t partial;
        if (_storage[locus-1].empty()) {
            // No stored partials for this locus, so allocate one
            partial = partial_t(new Partial(locus, _nelements[locus-1]));
            _total_partials_created++;
            _total_elements_created += _nelements[locus-1];
        }
        else {
            size_t n = _storage[locus-1].size();
            partial = _storage[locus-1].at(n-1);
            _storage[locus-1].pop_back();
            _total_partials_reused++;
            _total_elements_reused += _nelements[locus-1];
        }
        
        return partial;
#else
        // allocate a new partial
        assert(nelements>0);
        return partial_t(new std::vector<double> (nelements));
#endif
    }
        
#if defined(REUSE_PARTIALS)
    inline void PartialStore::putPartial(unsigned locus, partial_t partial) {
        // Check to make sure supplied value of locus is valid
        assert(_nelements.size() > locus-1);
        assert(_nelements[locus-1] > 0);

        // Store the partial for later
        assert(partial->_v.size() == _nelements[locus-1]);
        partial->_v.assign(_nelements[locus-1], 0.0);
        _storage[locus-1].push_back(partial);
    }
#endif
}
