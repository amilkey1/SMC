#pragma once

#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;
using namespace boost;

namespace proj {

    class Lot {
        public:
                                            Lot();
                                            ~Lot();
            
            void                            setSeed(unsigned seed);
            double                          uniform();
            int                             randint(int low, int high);
            pair<unsigned, unsigned>        nchoose2(unsigned n);
            double                          normal();
            double                          gamma(double shape, double scale);
            double                          logUniform();
            
            typedef std::shared_ptr<Lot>    SharedPtr;

        private:
        
            typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> >                uniform_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::random::normal_distribution<> >       normal_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::random::gamma_distribution<> >        gamma_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> >  uniform_int_generator_t;

            unsigned                                        _seed;
            boost::mt19937                                  _generator;
            std::shared_ptr<uniform_variate_generator_t>    _uniform_variate_generator;
            std::shared_ptr<normal_variate_generator_t>     _normal_variate_generator;
            std::shared_ptr<gamma_variate_generator_t>      _gamma_variate_generator;
            std::shared_ptr<uniform_int_generator_t>        _uniform_int_generator;

            double                                          _gamma_shape;
            int                                             _low;
            int                                             _high;
    };
    
    inline Lot::Lot() : _seed(0), _gamma_shape(1.0), _low(0), _high(100) {
        //std::cout << "Constructing a Lot" << std::endl;
        _generator.seed(static_cast<unsigned int>(std::time(0)));
        _uniform_variate_generator = std::shared_ptr<uniform_variate_generator_t>(new uniform_variate_generator_t(_generator, boost::random::uniform_01<>()));
        _normal_variate_generator = std::shared_ptr<normal_variate_generator_t>(new normal_variate_generator_t(_generator, boost::random::normal_distribution<>()));
        _gamma_variate_generator = std::shared_ptr<gamma_variate_generator_t>(new gamma_variate_generator_t(_generator, boost::random::gamma_distribution<>(_gamma_shape)));
        _uniform_int_generator = std::shared_ptr<uniform_int_generator_t>(new uniform_int_generator_t(_generator, boost::random::uniform_int_distribution<>(_low, _high)));
    }
        
    inline Lot::~Lot() {
        //std::cout << "Destroying a Lot" << std::endl;
        _uniform_variate_generator.reset();
        _normal_variate_generator.reset();
        _gamma_variate_generator.reset();
        _uniform_int_generator.reset();
    }
        
    inline void Lot::setSeed(unsigned seed) {
        _seed = seed;
        _generator.seed(_seed > 0 ? _seed : static_cast<unsigned int>(std::time(0)));
    }
        
    inline double Lot::uniform() {
        double u = (*_uniform_variate_generator)();
        while (u <= 0.0)
            u = (*_uniform_variate_generator)();
        return u;
    }

    inline double Lot::logUniform() {
        double u = (*_uniform_variate_generator)();
        while (u <= 0.0)
            u = (*_uniform_variate_generator)();
        return std::log(u);
    }
    
    inline double Lot::normal() {
        return (*_normal_variate_generator)();
    }

    inline double Lot::gamma(double shape, double scale) {
        assert(shape > 0.0);
        assert(scale > 0.0);
        if (shape != _gamma_shape) {
            _gamma_shape = shape;
            _gamma_variate_generator.reset(new gamma_variate_generator_t(_generator, boost::random::gamma_distribution<>(_gamma_shape,1)));
        }
        double deviate = (*_gamma_variate_generator)();
        return scale*deviate;
    }

    inline int Lot::randint(int low, int high) {
        if (low != _low || high != _high) {
            _low  = low;
            _high = high;
            _uniform_int_generator.reset(new uniform_int_generator_t(_generator, boost::random::uniform_int_distribution<>(_low,_high)));
        }
        return (*_uniform_int_generator)();
    }

    inline pair<unsigned, unsigned> Lot::nchoose2(unsigned n) {
        assert (n>1);
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        if (n > 2) {
            t1 = randint(0, n-1);
            t2 = randint(0, n-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = randint(0, n-1);
            }
        }
        assert (t1 < n);
        assert (t2 < n);

        return make_pair(t1, t2);
    }
    }
