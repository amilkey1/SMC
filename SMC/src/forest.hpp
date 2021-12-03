#pragma once    

#include <stack>    
#include <memory>
#include <iostream>
#include <boost/format.hpp>
#include <vector>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/algorithm/string/replace.hpp>


#include "lot.hpp"
extern proj::Lot rng;

#include "node.hpp"

namespace proj {

using namespace std;

class Likelihood;
class Updater;
class TreeUpdater;
class PolytomyUpdater;
class Particle;

class Forest {

        friend class Likelihood;
        friend class Updater;
        friend class TreeUpdater;
        friend class PolytomyUpdater;
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
        double                      calcLogLikelihood();
        void                        createDefaultTree();
        void operator=(const Forest & other);
        
    
    private:

        //POL placed all private methods together below here
        void                        clear();
        void                        setData(Data::SharedPtr d);
        void                        refreshPreorder();
        Node *                      findNextPreorder(Node * nd);
        std::string                 makeNewick(unsigned precision, bool use_names) const;
        pair<unsigned, unsigned>    chooseTaxaToJoin();
        void                        detachSubtree(Node * s);
        void                        insertSubtreeOnLeft(Node * s, Node * u);
        Node *                      findLeftSib(Node * nd);
        Node *                      getSubtreeAt(unsigned i);
        unsigned                    getNumSubtrees();
        void                        setEdgeLength(Node * nd);
        void                        createNewSubtree(unsigned t1, unsigned t2);
        double                      calcTransitionProbability(unsigned from, unsigned to, double edge_length);
        pair<double, double>        proposeBasalHeight();
        unsigned                    countNumberTrees() const;

        //POL placed all private data members together below here
        Node *                      _root;
        std::vector<Node *>         _preorder;
        std::vector<Node *>         _levelorder;
        std::vector<Node>           _nodes;
        std::vector< std::vector <double> > _partials;
        std::vector<Node *>         _unused_nodes;

        unsigned                    _nleaves;
        unsigned                    _ninternals;
        unsigned                    _npatterns;
        unsigned                    _nstates;
        unsigned                    _nsubtrees;
        double                      _last_edge_length;
        double                      _speciation_rate;
        pair <double, double>       _new_basal_height;
        pair <double, double>       _old_basal_height;
        
        Data::SharedPtr             _data;
        static unsigned             _nspecies;
        
    public:
        typedef std::shared_ptr<Forest> SharedPtr;
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
    _preorder.clear();
    _levelorder.clear();
    _partials.clear();
    _nodes.resize(2*_nspecies);
    _partials.resize(2*_nspecies);
    _npatterns = 0;
    _nstates = 4;
    _nsubtrees = _nspecies;
    _speciation_rate = 10.9; //temporary
    //creating root node
    _root = &_nodes[_nspecies];
    _root->_name="root";
    _root->_left_child=0;
    _root->_right_sib=0;
    _root->_parent=0;
    _root->_number=_nspecies;
    _root->_edge_length=double(0.0);
    
    //creating subroot node
    Node* subroot=&_nodes[_nspecies+1];
    subroot->_name="subroot";
    subroot->_left_child=0;
    subroot->_right_sib=0;
    subroot->_parent=_root;
    subroot->_number=_nspecies+1;
    subroot->_edge_length=double(0.0);
    _root->_left_child=subroot;
    
    //create species
    for (unsigned i = 0; i < _nspecies; i++) {
        Node* nd=&_nodes[i];
        if (i==0) {
            subroot->_left_child=nd;
        }
        else {
            _nodes[i-1]._right_sib=nd;
        }
        nd->_name=" ";
        nd->_left_child=0;
        nd->_right_sib=0;
        nd->_parent=subroot;
        nd->_number=i;
        nd->_edge_length=0.0;
        }
    _nleaves=_nspecies;
    _ninternals=2;
    refreshPreorder();
    
    //initialize basal edge lengths
    _old_basal_height = _new_basal_height;
    _new_basal_height = proposeBasalHeight();
//    _old_basal_height = _new_basal_height;
    
    //initialize species nodes with basal height
    for (auto nd:_preorder) {
        if (nd->_parent==_root->_left_child) {
        nd->_edge_length=_new_basal_height.first;
//          nd->_edge_length=_old_basal_height.first;
        }
    }
    }

inline Forest::Forest(const Forest & other) {
    clear();
    *this = other;
}

inline void Forest::setData(Data::SharedPtr d) {
    _data = d;
    _npatterns = d->getNumPatterns();
    const Data::taxon_names_t & taxon_names = _data->getTaxonNames();
    unsigned i = 0;
    for (unsigned j=0; j<_nodes.size(); j++) {
        _partials[j].resize(_npatterns * _nstates);
    }
    for (auto nd:_preorder) {
//        _partials[nd->_number].resize(_npatterns * _nstates);
        if (!nd->_left_child) {
            // replace all spaces with underscores so that other programs do not have
              // trouble parsing your tree descriptions
              std::string name = taxon_names[i++];
              boost::replace_all(name, " ", "_");
              nd->_name = name;
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

inline void Forest::refreshPreorder() {
    // Create vector of node pointers in preorder sequence
    _preorder.clear();
    _preorder.reserve(_nodes.size() - 1); // _preorder does not include root node

    if (!_root)
        return;

    Node * first_preorder = _root->_left_child;

    // sanity check: first preorder node should be the only child of the root node
    assert(first_preorder->_right_sib == 0);

    Node * nd = first_preorder;
    _preorder.push_back(nd);

    while (true) {
        nd = findNextPreorder(nd);
        if (nd)
            _preorder.push_back(nd);
        else
            break;
    }   // end while loop
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
    cout << " " << makeNewick(9, true) << "\n"; //POL: changed precision from 3 to 9
}

inline std::string Forest::makeNewick(unsigned precision, bool use_names) const {
    std::string newick;
    const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
    const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
    const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
    std::stack<Node *> node_stack;


    for (auto nd : _preorder) {
        if (nd->_left_child) {
            newick += "(";
            node_stack.push(nd);
        }
        else {
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
                        newick += ")";
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
        }
    }

    return newick;
}

inline void Forest::setNumSpecies(unsigned n){
    //TODO be careful, this function should not be called after forests have already been created
    _nspecies=n;
}

inline pair<unsigned, unsigned> Forest::chooseTaxaToJoin(){
    unsigned t1=0;
    unsigned t2=1;
    unsigned nsubtrees = getNumSubtrees();
    //don't use this when there's only one choice
        //i.e. # of trees = 2
    if (nsubtrees > 2) {
        t1 = ::rng.randint(0, nsubtrees-1);
        t2 = ::rng.randint(0, nsubtrees-1);

        //keep calling t2 until it doesn't equal t1
        while (t2 == t1) {
            t2 = ::rng.randint(0, nsubtrees-1);
        }
    }
    return make_pair(t1, t2);
}

inline unsigned Forest::countNumberTrees() const{
    Node* subroot = _root->_left_child;
    assert (subroot);
    unsigned n = 0;
    for (Node* child=subroot->_left_child; child; child=child->_right_sib) {
        n++;
    }
    return n;
}

inline pair<double, double> Forest::proposeBasalHeight() {
    //draw random Uniform(0,1)
    double new_basal_height = 0.0;
    double log_prob_basal_height = 0.0;
    unsigned s = countNumberTrees();
    for (unsigned k=2; k<=s; k++) {
        double u = rng.uniform();
        //transform uniform deviate to exponential, number of basal lineages = nspecies - ninternals + 2 (from root & subroot)
        double t = -log(1-u)/(_speciation_rate*k);
        new_basal_height += t;
        log_prob_basal_height += log(k*_speciation_rate)-(k*_speciation_rate*t);
    }
    
    return make_pair(new_basal_height, log_prob_basal_height);
}

inline void Forest::createNewSubtree(unsigned t1, unsigned t2) {
    _last_edge_length = rng.gamma(1.0, 1.0/(_nsubtrees*_speciation_rate));
    for (auto nd:_preorder) {
        //if node's parent is subroot, subtract basal height and add new branch length
        if (nd->_parent==_root->_left_child){
            nd->_edge_length -= _new_basal_height.first; //subtract old basal height from each node whose parent is subroot
            nd->_edge_length += _last_edge_length; //add most recently chosen branch length to each node whose parent is subroot
        }
    }
    
    assert(_nsubtrees>1);
    _nsubtrees--;
    Node * subtree1=getSubtreeAt(t1);
    Node * subtree2 = getSubtreeAt(t2);

    detachSubtree(subtree1);
    detachSubtree(subtree2);
    //creating new node
    Node* new_nd=&_nodes[_nleaves+_ninternals];
    new_nd->_name=" ";
    new_nd->_left_child=subtree1;
    new_nd->_right_sib=0;
    new_nd->_parent=_root->_left_child;
    new_nd->_number=_nleaves+_ninternals;
    
    new_nd->_edge_length=0.0;

//    cout << "New node branch length is: " << new_nd->_edge_length << endl;

    _ninternals++;
    subtree1 -> _right_sib=subtree2;
    subtree1->_parent=new_nd;
    subtree2->_parent=new_nd;

    insertSubtreeOnLeft(new_nd, _root->_left_child);

    refreshPreorder();
    
    //POL cout << makeNewick(3, true) << endl;
    
    //once new taxa have been joined, add new basal height to finish off tree
    _old_basal_height = _new_basal_height;
    _new_basal_height = proposeBasalHeight();
    
    for (auto nd:_preorder) {
        if (nd->_parent==_root->_left_child){
            nd->_edge_length += _new_basal_height.first;
            }
    }
    
//    _old_basal_height = _new_basal_height;
}

inline void Forest::detachSubtree(Node * s) {
    assert(s);
    assert(s->_parent);
    
    // Save pointers to relevant nodes
    Node * s_leftsib  = findLeftSib(s);
    Node * s_rightsib = s->_right_sib;
    Node * s_parent   = s->_parent;

    // Completely detach s and seal up the wound
    s->_parent = 0;
    s->_right_sib = 0;
    if (s_leftsib)
        s_leftsib->_right_sib = s_rightsib;
    else
        s_parent->_left_child = s_rightsib;
}

inline void Forest::insertSubtreeOnLeft(Node * s, Node * u) {
    assert(u);
    assert(s);
    s->_right_sib  = u->_left_child;
    s->_parent     = u;
    u->_left_child = s;
}

inline Node * Forest::findLeftSib(Node * nd) {
    assert(nd);
    assert(nd->_parent);
    Node * child = nd->_parent->_left_child;
    while (child && child->_right_sib != nd)
        child = child->_right_sib;
    return child;
}

inline Node * Forest::getSubtreeAt(unsigned i) {
    unsigned nsubtrees = 0;
    Node * the_nd = 0;
    for (auto nd : _preorder) {
        if (nd->_parent == _root->_left_child) {
            if (nsubtrees == i) {
                the_nd = nd;
                break;
            }
            nsubtrees++;
        }
    }
    assert(the_nd);
    return the_nd;
}

inline unsigned Forest::getNumSubtrees() {
    unsigned nsubtrees = 0;
    for (auto nd : _preorder) {
        if (nd->_parent == _root->_left_child) {
            nsubtrees++;
        }
    }
//    cout << "Number of subtrees is: " << nsubtrees << endl;
    return nsubtrees;
}

inline double Forest::calcTransitionProbability(unsigned from, unsigned to, double edge_length){
    double transition_prob = 0.0;
    if (from == to) {
        transition_prob = 0.25 + 0.75*exp(-4.0*edge_length/3.0);
    }
    
    else {
        transition_prob = 0.25 - 0.25*exp(-4.0*edge_length/3.0);
    }
    return transition_prob;
}

inline double Forest::calcLogLikelihood() {
//    return 0.0;
    auto data_matrix=_data->getDataMatrix();
    for (auto nd : boost::adaptors::reverse(_preorder)) {
        if (nd ==_root) {
            
            //if root is reached
            break;
        }
        if (nd->_left_child){
            //internal node
            assert(nd->_left_child->_right_sib);
            for (unsigned p=0; p<_npatterns; p++) {
                for (unsigned s=0; s<_nstates; s++) {
                    double partial=1.0;
                    for (Node * child=nd->_left_child; child; child=child->_right_sib){
                        double child_sum=0.0;
                        for (unsigned s_child=0; s_child<_nstates; s_child++) {
                            double child_transition_prob = calcTransitionProbability(s, s_child, child->_edge_length);
                            double child_partial = _partials[child->_number][p*_nstates + s_child];
                            child_sum += child_transition_prob * child_partial;
                        }
                        partial *= child_sum;
                    }
//                    assert(_partials[nd->_number][0]);
                    _partials[nd->_number][p*_nstates+s]= partial;
                }
            }
        }
        else {
            //leaf node
            for (unsigned p=0; p<_npatterns; p++) {
                for (unsigned s=0; s<_nstates; s++) {
                    Data::state_t state = (Data::state_t)1 << s;
                    Data::state_t d = data_matrix[nd->_number][p];
                    double result = state & d;
                    _partials[nd->_number][p*_nstates+s]= (result == 0.0 ? 0.0:1.0);
                }
            }
        }
    }
//    compute log likelihood of every subtree whose parent is subroot
    auto counts = _data->getPatternCounts();
//    Node* subroot = _root->_left_child;
    double composite_log_likelihood = 0.0;
    for (auto nd=_root->_left_child; nd; nd=nd->_right_sib) {
        double log_like = 0.0;
        for (unsigned p=0; p<_npatterns; p++) {
            double site_like = 0.0;
            for (unsigned s=0; s<_nstates; s++) {
                for (unsigned ss=0; ss<_nstates; ss++) {
                    double child_transition_prob = calcTransitionProbability(s, ss, nd->_edge_length);
                    site_like += 0.25*child_transition_prob*_partials[nd->_number][p*_nstates+ss];
                }
            }
            log_like += log(site_like)*counts[p];
        }
        composite_log_likelihood += log_like;
    }
    
    return composite_log_likelihood;
}

inline void Forest::createDefaultTree() {
    clear();
    
    //creating root node
    _root = &_nodes[_nspecies];
    _root->_name="root";
    _root->_left_child=0;
    _root->_right_sib=0;
    _root->_parent=0;
    _root->_number=_nspecies;
    _root->_edge_length=double(0.0);
    
    //creating subroot node
    Node* subroot=&_nodes[_nspecies+1];
    subroot->_name="subroot";
    subroot->_left_child=0;
    subroot->_right_sib=0;
    subroot->_parent=_root;
    subroot->_number=_nspecies + 1;
    subroot->_edge_length=double(0.0);
    _root->_left_child=subroot;
    
    //create species
    double edge_length = rng.gamma(1.0, 1.0/_nspecies);
    for (unsigned i = 0; i < _nspecies; i++) {
        Node* nd=&_nodes[i];
        if (i==0) {
            subroot->_left_child=nd;
        }
        else {
            _nodes[i-1]._right_sib=nd;
        }
        nd->_name="";
        nd->_left_child=0;
        nd->_right_sib=0;
        nd->_parent=subroot;
        nd->_number=i;
        nd->_edge_length = edge_length;
        }
    _nleaves=_nspecies;
    _ninternals=2;
    refreshPreorder();
}

inline void Forest::operator=(const Forest & other) {
    _nstates = other._nstates;
    _npatterns = other._npatterns;
    _nodes.resize(other._nodes.size());
    _preorder.resize(other._preorder.size());
    _partials.resize(other._partials.size());
    copy(other._partials.begin(), other._partials.end(), _partials.begin());
    
    _speciation_rate  = other._speciation_rate;     //POL added
    _nsubtrees        = other._nsubtrees;           //POL added
    _nleaves          = other._nleaves;             //POL added
    _ninternals       = other._ninternals;          //POL added
    _last_edge_length = other._last_edge_length;    //POL added
    _old_basal_height = other._old_basal_height;
    _new_basal_height = other._new_basal_height;
    // copy tree itself
    assert(other._root->_number == _nspecies);
    assert(other._root->_left_child->_number == _nspecies + 1);
    _root               = &_nodes[_nspecies];
    _root->_number      = _nspecies;
    _root->_name        = "root";
    _root->_edge_length = 0.0;
    _root->_parent      = 0;
    _root->_left_child  = &_nodes[_nspecies + 1];
    _root->_right_sib   = 0;

    unsigned i = 0;
    for (auto othernd : other._preorder) {
        // get number of next node in preorder sequence (serves as index of node in _nodes vector)
        unsigned k = othernd->_number;

        // update preorder vector
        _preorder[i] = &_nodes[k];
        
        // copy parent
        assert(othernd->_parent);
        unsigned parent_number = othernd->_parent->_number;
        _nodes[k]._parent = &_nodes[parent_number];

        // copy left child
        if (othernd->_left_child) {
            unsigned left_child_number = othernd->_left_child->_number;
            _nodes[k]._left_child = &_nodes[left_child_number];
        }
        else
            _nodes[k]._left_child = 0;
        
        // copy right sibling
        if (othernd->_right_sib) {
            unsigned right_sib_number = othernd->_right_sib->_number;
            _nodes[k]._right_sib = &_nodes[right_sib_number];
        }
        else
            _nodes[k]._right_sib = 0;
            
        _nodes[k]._number      = othernd->_number;
        _nodes[k]._name        = othernd->_name;
        _nodes[k]._edge_length = othernd->_edge_length;
        
        i++;
    }
}

}
