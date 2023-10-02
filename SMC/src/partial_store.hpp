#pragma once

using namespace std;

namespace proj {
    class PartialStore {
        public:
            PartialStore();
            ~PartialStore();
            typedef std::shared_ptr < std::vector<double> >   partial_t;
            partial_t getPartial(unsigned nelements, unsigned gene);
//            void setnelements(unsigned nelements) {_nelements = nelements;}
            typedef vector<partial_t>            vect_partial_t;
            typedef vector<vect_partial_t>       storage_t;

            void           memoryReport(ofstream & memf) const;
            void resetPartial();
            void        setNGenes(unsigned ngenes);
            void       setNElements(unsigned nelements, unsigned gene);
            unsigned long   getTotalBytesAllocated() const;
            void       stowPartial(partial_t & p, unsigned gene);

        private:
//            unsigned _nelements;
            vector<unsigned> _nelements;
            vector<unsigned> _nallocated;
            storage_t        _storage;
    };

    inline PartialStore::PartialStore() {
        //std::cout << "Constructing a partial store" << std::endl;
//        _nelements=0;
    }

    inline PartialStore::~PartialStore() {
        //std::cout << "Destroying a partial store" << std::endl;
        _nelements.clear();
        _nallocated.clear();
        _storage.clear();
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned nelements, unsigned gene) {
        assert(nelements>0);
//        _nelements[gene] = nelements;
        assert(_nelements.size() > gene);
        assert(_nallocated.size() > gene);
        assert(_storage.size() > gene);
        assert(_nelements[gene] > 0);
        
        if (_storage[gene].size() > 0) {
            partial_t last = _storage[gene].back();
            _storage[gene].pop_back();
            return last;
        }
        
        else {
            partial_t ptr = partial_t(new vector<double>(nelements));
            _nallocated[gene]++;
            return ptr;
        }
        
//        _storage[gene].push_back(partial_t(new std::vector<double> (nelements)));
//        _nallocated[gene]++;
//        return partial_t(new std::vector<double> (nelements));
    }

    inline unsigned long PartialStore::getTotalBytesAllocated() const {
        unsigned long total_bytes = 0L;
        for (unsigned g = 0; g < _nallocated.size(); ++g) {
            total_bytes += _nelements[g]*_nallocated[g];
        }
        return total_bytes;
    }

    inline void PartialStore::setNGenes(unsigned ngenes) {
        assert(_nelements.empty());
        assert(_nallocated.empty());
        assert(_storage.empty());

        // Resize both containers
        _nelements.resize(ngenes);
        _nallocated.resize(ngenes, 0);
        _storage.resize(ngenes);
        }

    inline void PartialStore::setNElements(unsigned nelements, unsigned gene) {
        assert(_nelements.size() > gene);
        _nelements[gene] = nelements*4;
    }

    inline void PartialStore::stowPartial(partial_t & p, unsigned gene) {
        assert(_storage.size() > gene);
        fill(p->begin(), p->end(), 0.0);

        // //temporary!
        //cerr << str(format("~~~> stowing partial for gene %d\n") % gene);

        _storage[gene].push_back(p);
    }

    inline void PartialStore::memoryReport(ofstream & memf) const {
        getTotalBytesAllocated();
        memf << "\nPartialStore memory report:\n\n";
        memf << str(boost::format("  %12s %12s %12s %12s %12s %12s\n") %         "gene" %         "size" %       "in use" %       "stored" %        "total" %        "bytes");
        memf << str(boost::format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
        unsigned long total_bytes = 0L;
        unsigned total_in_use = 0;
        unsigned total_stored = 0;
        unsigned total_allocated = 0;
        for (unsigned g = 0; g < _nallocated.size(); ++g) {
            unsigned long nbytes = _nelements[g]*_nallocated[g];
            unsigned nallocated = _nallocated[g];
            total_allocated += nallocated;
            unsigned nstored = (unsigned)_storage[g].size();
            total_stored += nstored;
            unsigned in_use = nallocated - nstored;
            total_in_use += in_use;
            memf << str(boost::format("  %12d %12d %12d %12d %12d %12d\n") % g % _nelements[g] % in_use % nstored % nallocated % nbytes);
            total_bytes += nbytes;
        }
        memf << str(boost::format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
        memf << str(boost::format("  %12s %12s %12s %12d %12d %12d\n") % " " % " " % total_in_use % total_stored % total_allocated % total_bytes);
        memf << str(boost::format("  Total megabytes: %.5f\n") % (1.0*total_bytes/1048576));
    }
}
