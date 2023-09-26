#pragma once

namespace proj {
    class PartialStore {
        public:
            PartialStore();
            ~PartialStore();
            typedef std::shared_ptr < std::vector<double> >   partial_t;
            partial_t getPartial(unsigned nelements);
//            void setnelements(unsigned nelements) {_nelements = nelements;}
//        void           memoryReport(ofstream & memf) const;
        void resetPartial();
//        private:
//            unsigned _nelements;
    };

    inline PartialStore::PartialStore() {
        //std::cout << "Constructing a partial store" << std::endl;
//        _nelements=0;
    }

    inline PartialStore::~PartialStore() {
        //std::cout << "Destroying a partial store" << std::endl;
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned nelements) {
        assert(nelements>0);
        return partial_t(new std::vector<double> (nelements));
    }

//    inline void PartialStore::memoryReport(ofstream & memf) const {
//        memf << "\nPartialStore memory report:\n\n";
//        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") %         "gene" %         "size" %       "in use" %       "stored" %        "total" %        "bytes");
//        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
//        unsigned long total_bytes = 0L;
//        unsigned total_in_use = 0;
//        unsigned total_stored = 0;
//        unsigned total_allocated = 0;
//        for (unsigned g = 0; g < _nallocated.size(); ++g) {
//            unsigned long nbytes = _nelements[g]*_nallocated[g];
//            unsigned nallocated = _nallocated[g];
//            total_allocated += nallocated;
//            unsigned nstored = (unsigned)_storage[g].size();
//            total_stored += nstored;
//            unsigned in_use = nallocated - nstored;
//            total_in_use += in_use;
//            memf << str(format("  %12d %12d %12d %12d %12d %12d\n") % g % _nelements[g] % in_use % nstored % nallocated % nbytes);
//            total_bytes += nbytes;
//        }
//        memf << str(format("  %12s %12s %12s %12s %12s %12s\n") % " -----------" % " -----------" % " -----------" % " -----------" % " -----------" % " -----------");
//        memf << str(format("  %12s %12s %12s %12d %12d %12d\n") % " " % " " % total_in_use % total_stored % total_allocated % total_bytes);
//        memf << str(format("  Total megabytes: %.5f\n") % (1.0*total_bytes/1048576));
//    }
}
