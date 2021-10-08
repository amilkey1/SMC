#pragma once

#include <boost/format.hpp>

namespace proj {

    class XProj : public std::exception {
        public:
                                XProj() throw() {}
                                XProj(const std::string s) throw() : _msg() {_msg = s;}
                                XProj(const boost::format & f) throw() : _msg() {_msg = boost::str(f);}
            virtual             ~XProj() throw() {}
            const char *        what() const throw() {return _msg.c_str();}

        private:

            std::string         _msg;
    };

}
