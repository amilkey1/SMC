#pragma once

namespace proj {
    struct Partial {
        Partial(unsigned g, unsigned n);
        ~Partial();
        unsigned        _g; // the gene
        vector<double>  _v; // the partial array: length = _nstates*<no. patterns>
    };

    inline Partial::Partial(unsigned g, unsigned n) {
        _g = g;
        _v.resize(n);
    }

    inline Partial::~Partial() {
    }

    class PartialStore {
        public:
            PartialStore();
            ~PartialStore();
            typedef std::shared_ptr<Partial>          partial_t;
//            typedef std::shared_ptr <vector<double>>   partial_t;
//            partial_t getPartial(unsigned nelements);
            partial_t getPartial(unsigned index, unsigned nelements);
//            void setnelements(unsigned nelements) {_nelements = nelements;}
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

    inline PartialStore::partial_t PartialStore::getPartial(unsigned index, unsigned nelements) {
        assert(nelements>0);
        return partial_t(new Partial(index, nelements));
//        return partial_t(new std::vector<double> (nelements));
    }
}
