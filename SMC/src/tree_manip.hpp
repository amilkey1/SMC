//#pragma once
//
////#include "conditionals.hpp"
//
//#include <cassert>
//#include <memory>
//#include <stack>
//#include <queue>
//#include <set>
//#if defined(USE_BOOST_REGEX)
//#   include <boost/regex.hpp>
//#else
//#   include <regex>
//#endif
//#include <boost/range/adaptor/reversed.hpp>
//#include <boost/format.hpp>
//#if defined(HPD_PWK_METHOD)
//#   include <Eigen/Dense>
//#endif
//#include "particle.hpp"
//#include "lot.hpp"
//#include "xproj.hpp"
//
//namespace proj {
//
//    class TreeManip {
//        public:
//                                        TreeManip();
//                                        TreeManip(Particle::SharedPtr p);
//                                        ~TreeManip();
//
//            void                        setTree(Particle::SharedPtr p);
//            Particle::SharedPtr             getTree();
//
//            double                      calcTreeLength() const;
//            unsigned                    calcResolutionClass() const;
//            unsigned                    countEdges() const;
//            unsigned                    countInternals() const;
//            void                        scaleAllEdgeLengths(double scaler);
//
//            void                        createTestTree();
//            std::string                 makeNewick(unsigned precision, bool use_names = false) const;
//
//            void                        buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies);
//            void                        storeSplits(std::set<Split> & splitset);
//            void                        rerootAtNodeNumber(int node_number);
//
//            Node *                      randomEdge(Lot::SharedPtr lot);
//            Node *                      randomInternalEdge(Lot::SharedPtr lot);
//            Node *                      randomChild(Lot::SharedPtr lot, Node * x, Node * avoid, bool parent_included);
//            void                        LargetSimonSwap(Node * a, Node * b);
//
//            bool                        isPolytomy(Node * nd) const;
//            void                        nniNodeSwap(Node * a, Node * b);
//            unsigned                    countChildren(Node * nd) const;
//            Node *                      findLeftSib(Node * nd);
//            Node *                      findRightmostChild(Node * nd);
//            Node *                      findLastPreorderInClade(Node * start);
//            void                        insertSubtreeOnLeft(Node * s, Node * u);
//            void                        insertSubtreeOnRight(Node * s, Node * u);
//            void                        detachSubtree(Node * s);
//            void                        rectifyNumInternals(int incr);
//            void                        refreshNavigationPointers();
//            Node *                      getUnusedNode(Node * sought = 0);
//            void                        putUnusedNode(Node * nd);
//
//            void                        selectAll();
//            void                        deselectAll();
//            void                        selectAllPartials();
//            void                        deselectAllPartials();
//            void                        selectAllTMatrices();
//            void                        deselectAllTMatrices();
//            void                        selectPartialsHereToRoot(Node * a);
//            void                        flipPartialsAndTMatrices();
//
//            void                        clear();
//
//            double                      copyEdgeProportionsTo(std::vector<double> & receptacle) const;
//            void                        copyEdgeProportionsFrom(double TL, const std::vector<double> & new_props);
//
//#if defined(HPD_PWK_METHOD)
//            void                        saveParamNames(std::vector<std::string> & param_name_vect) const;
//            double                      logTransformEdgeLengths(std::vector<double> & param_vect) const;
//            double                      setEdgeLengthsFromLogTransformed(Eigen::VectorXd & param_vect, double TL, unsigned first, unsigned nedges);
//#endif
//
//        private:
//
//            Node *                      findNextPreorder(Node * nd);
//            void                        refreshPreorder();
//            void                        refreshLevelorder();
//            void                        renumberInternals();
//            void                        rerootAtNode(Node * prospective_root);
//            void                        extractNodeNumberFromName(Node * nd, std::set<unsigned> & used);
//            void                        extractEdgeLen(Node * nd, std::string edge_length_string);
//            unsigned                    countNewickLeaves(const std::string newick);
//            void                        stripOutNexusComments(std::string & newick);
//            bool                        canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);
//
//            Particle::SharedPtr             _particle;
//
//        public:
//
//            typedef std::shared_ptr< TreeManip > SharedPtr;
//    };
//
//    inline TreeManip::TreeManip() {
//        //std::cout << "Constructing a TreeManip" << std::endl;
//        clear();
//    }
//
//    inline TreeManip::TreeManip(Particle::SharedPtr t) {
//        //std::cout << "Constructing a TreeManip with a supplied tree" << std::endl;
//        clear();
//        setTree(t);
//    }
//
//    inline TreeManip::~TreeManip() {
//        //std::cout << "Destroying a TreeManip" << std::endl;
//    }
//
//    inline void TreeManip::clear() {
//        _particle.reset();
//    }
//
//    inline void TreeManip::setTree(Particle::SharedPtr t) {
//        assert(t);
//        _particle = t;
//    }
//
//    inline Particle::SharedPtr TreeManip::getTree() {
//        return _particle;
//    }
//
////    inline double TreeManip::calcTreeLength() const {
////        double TL = 0.0;
////        for (auto nd : _tree->_preorder) {
////            TL += nd->_edge_length;
////        }
////        return TL;
////    }
//
////    inline unsigned TreeManip::calcResolutionClass() const {
////        return _particle->_ninternals;
////    }
//
////    inline unsigned TreeManip::countEdges() const {
////        return (unsigned)_particle->_preorder.size();
////    }
//
////    inline unsigned TreeManip::countInternals() const {
////        unsigned m = 0;
////        for (auto nd : _tree->_preorder) {
////            if (nd->_left_child)
////                m++;
////        }
////        return m;
////    }
//
////    inline unsigned TreeManip::countChildren(Node * nd) const {
////        assert(nd);
////        unsigned nchildren = 0;
////        Node * child = nd->getLeftChild();
////        while (child) {
////            nchildren++;
////            child = child->getRightSib();
////        }
////        return nchildren;
////    }
//
////    inline void TreeManip::scaleAllEdgeLengths(double scaler) {
////        for (auto nd : _particle->_preorder) {
////            nd->setEdgeLength(scaler*nd->_edge_length);
////        }
////    }
//
////    inline void TreeManip::createTestTree() {
////        clear();
////        _particle = Particle::SharedPtr(new Particle());
////        _particle->_nodes.resize(6);
////
////        Node * root_node       = &_particle->_nodes[0];
////        Node * first_internal  = &_particle->_nodes[1];
////        Node * second_internal = &_particle->_nodes[2];
////        Node * first_leaf      = &_particle->_nodes[3];
////        Node * second_leaf     = &_particle->_nodes[4];
////        Node * third_leaf      = &_particle->_nodes[5];
////
////        // Here is the structure of the tree (numbers in
////        // parentheses are node numbers, other numbers
////        // are edge lengths):
////        //
////        // first_leaf (0)   second_leaf (1)   third_leaf (2)
////        //      \              /                  /
////        //       \ 0.1        / 0.1              /
////        //        \          /                  /
////        //     second_internal (3)             / 0.2
////        //             \                      /
////        //              \ 0.1                /
////        //               \                  /
////        //                first_internal (4)
////        //                        |
////        //                        | 0.1
////        //                        |
////        //                    root_node (5)
////        //
////        root_node->_parent = 0;
////        root_node->_left_child = first_internal;
////        root_node->_right_sib = 0;
////        root_node->_number = 5;
////        root_node->_name = "root_node";
////        root_node->_edge_length = 0.0;
////
////        first_internal->_parent = root_node;
////        first_internal->_left_child = second_internal;
////        first_internal->_right_sib = 0;
////        first_internal->_number = 4;
////        first_internal->_name = "first_internal_node";
////        first_internal->_edge_length = 0.1;
////
////        second_internal->_parent = first_internal;
////        second_internal->_left_child = first_leaf;
////        second_internal->_right_sib = third_leaf;
////        second_internal->_number = 3;
////        second_internal->_name = "second_internal_node";
////        second_internal->_edge_length = 0.1;
////
////        first_leaf->_parent = second_internal;
////        first_leaf->_left_child = 0;
////        first_leaf->_right_sib = second_leaf;
////        first_leaf->_number = 0;
////        first_leaf->_name = "first_leaf";
////        first_leaf->_edge_length = 0.1;
////
////        second_leaf->_parent = second_internal;
////        second_leaf->_left_child = 0;
////        second_leaf->_right_sib = 0;
////        second_leaf->_number = 1;
////        second_leaf->_name = "second_leaf";
////        second_leaf->_edge_length = 0.1;
////
////        third_leaf->_parent = first_internal;
////        third_leaf->_left_child = 0;
////        third_leaf->_right_sib = 0;
////        third_leaf->_number = 2;
////        third_leaf->_name = "third_leaf";
////        third_leaf->_edge_length = 0.2;
////
////        _tree->_is_rooted = true;
////        _tree->_root = root_node;
////        _tree->_nleaves = 3;
////
////        // Note that root node is not included in _preorder
////        _tree->_preorder.push_back(first_internal);
////        _tree->_preorder.push_back(second_internal);
////        _tree->_preorder.push_back(first_leaf);
////        _tree->_preorder.push_back(second_leaf);
////        _tree->_preorder.push_back(third_leaf);
////
////        _tree->_levelorder.push_back(first_internal);
////        _tree->_levelorder.push_back(second_internal);
////        _tree->_levelorder.push_back(third_leaf);
////        _tree->_levelorder.push_back(first_leaf);
////        _tree->_levelorder.push_back(second_leaf);
////    }
//
////    inline std::string TreeManip::makeNewick(unsigned precision, bool use_names) const {
////        std::string newick;
////        const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
////        const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
////        const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
////        std::stack<Node *> node_stack;
////
////        Node * root_tip = (_tree->_is_rooted ? 0 : _tree->_root);
////        for (auto nd : _tree->_preorder) {
////            if (nd->_left_child) {
////                newick += "(";
////                node_stack.push(nd);
////                if (root_tip) {
////                    if (use_names) {
////                        newick += boost::str(boost::format(tip_node_name_format)
////                            % root_tip->_name
////                            % nd->_edge_length);
////                    } else {
////                        newick += boost::str(boost::format(tip_node_number_format)
////                            % (root_tip->_number + 1)
////                            % nd->_edge_length);
////                    }
////                    newick += ",";
////                    root_tip = 0;
////                }
////            }
////            else {
////                if (use_names) {
////                    newick += boost::str(boost::format(tip_node_name_format)
////                        % nd->_name
////                        % nd->_edge_length);
////                } else {
////                    newick += boost::str(boost::format(tip_node_number_format)
////                        % (nd->_number + 1)
////                        % nd->_edge_length);
////                }
////                if (nd->_right_sib)
////                    newick += ",";
////                else {
////                    Node * popped = (node_stack.empty() ? 0 : node_stack.top());
////                    while (popped && !popped->_right_sib) {
////                        node_stack.pop();
////                        if (node_stack.empty()) {
////                            newick += ")";
////                            popped = 0;
////                        }
////                        else {
////                            newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
////                            popped = node_stack.top();
////                        }
////                    }
////                    if (popped && popped->_right_sib) {
////                        node_stack.pop();
////                        newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
////                        newick += ",";
////                    }
////                }
////            }
////        }
////
////        return newick;
////    }
//
////    inline void TreeManip::extractNodeNumberFromName(Node * nd, std::set<unsigned> & used) {
////        assert(nd);
////        bool success = true;
////        unsigned x = 0;
////        try {
////            x = std::stoi(nd->_name);
////        }
////        catch(std::invalid_argument &) {
////            // node name could not be converted to an integer value
////            success = false;
////        }
////
////        if (success) {
////            // conversion succeeded
////            // attempt to insert x into the set of node numbers already used
////            std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
////            if (insert_result.second) {
////                // insertion was made, so x has NOT already been used
////                nd->_number = x - 1;
////            }
////            else {
////                // insertion was not made, so set already contained x
////                throw XStrom(boost::str(boost::format("leaf number %d used more than once") % x));
////            }
////        }
////        else
////            throw XStrom(boost::str(boost::format("node name (%s) not interpretable as a positive integer") % nd->_name));
////    }
//
////    inline void TreeManip::extractEdgeLen(Node * nd, std::string edge_length_string) {
////        assert(nd);
////        bool success = true;
////        double d = 0.0;
////        try {
////            d = std::stod(edge_length_string);
////        }
////        catch(std::invalid_argument &) {
////            // edge_length_string could not be converted to a double value
////            success = false;
////        }
////
////        if (success) {
////            // conversion succeeded
////            nd->setEdgeLength(d);
////        }
////        else
////            throw XProj(boost::str(boost::format("%s is not interpretable as an edge length") % edge_length_string));
////
////    }
//
////    inline unsigned TreeManip::countNewickLeaves(const std::string newick) {
////#if defined(USE_BOOST_REGEX)
////        boost::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
////        boost::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
////        boost::sregex_iterator m2;
////#else
////        std::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
////        std::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
////        std::sregex_iterator m2;
////#endif
////        return (unsigned)std::distance(m1, m2);
////    }
//
//    inline void TreeManip::stripOutNexusComments(std::string & newick) {
//#if defined(USE_BOOST_REGEX)
//        boost::regex commentexpr("\\[.*?\\]");
//        newick = boost::regex_replace(newick, commentexpr, std::string(""));
//#else
//        std::regex commentexpr("\\[.*?\\]");
//        newick = std::regex_replace(newick, commentexpr, std::string(""));
//#endif
//    }
//
//    inline void TreeManip::refreshPreorder() {
//        // Create vector of node pointers in preorder sequence
//        _tree->_preorder.clear();
//        _tree->_preorder.reserve(_tree->_nodes.size() - 1); // _preorder does not include root node
//
//        if (!_tree->_root)
//            return;
//
//        Node * first_preorder = _tree->_root->_left_child;
//
//        // sanity check: first preorder node should be the only child of the root node
//        assert(first_preorder->_right_sib == 0);
//
//        Node * nd = first_preorder;
//        _tree->_preorder.push_back(nd);
//
//        while (true) {
//            nd = findNextPreorder(nd);
//            if (nd)
//                _tree->_preorder.push_back(nd);
//            else
//                break;
//        }   // end while loop
//    }
//
//    //                            1. start by adding only descendant of root node to buffer queue
//    //                               queue = [1], stack = []
//    //                            2. move node at front of buffer queue to back of stack vector
//    //                               queue = [], stack = [1]
//    //                 8    9     3. add this node's immediate children to back of buffer queue
//    //                  \  /         queue = [2,3], stack = [1]
//    //                   \/       4. repeat 2 and 3 until all nodes are processed
//    //      4  5    6    7           (2) queue = [3], stack = [1,2]
//    //       \ |     \  /            (3) queue = [3,4,5], stack = [1,2]
//    //        \|      \/             (2) queue = [4,5], stack = [1,2,3]
//    //         2      3              (3) queue = [4,5,6,7], stack = [1,2,3]
//    //          \    /               (2) queue = [5,6,7], stack = [1,2,3,4]
//    //           \  /                (3) no-op: 4 has no children
//    //            \/                 (2) queue = [6,7], stack = [1,2,3,4,5]
//    //            1                  (3) no-op: 5 has no children
//    //            |                  (2) queue = [7], stack = [1,2,3,4,5,6]
//    //            0                  (3) no-op: 6 has no children
//    //                               (2) queue = [], stack = [1,2,3,4,5,6,7]
//    //                               (3) queue = [8,9], stack = [1,2,3,4,5,6,7]
//    //                               (2) queue = [9], stack = [1,2,3,4,5,6,7,8]
//    //                               (3) no-op: 8 has no children
//    //                               (2) queue = [], stack = [1,2,3,4,5,6,7,8,9]
//    //                               (3) no-op: 9 has no children
//    //                            5. stack vector is now in level order
//    inline void TreeManip::refreshLevelorder() {
//        if (!_tree->_root)
//            return;
//
//        // q is the buffer queue
//        std::queue<Node *> q;
//
//        // _tree->_levelorder is the stack vector
//        _tree->_levelorder.clear();
//        _tree->_levelorder.reserve(_tree->_nodes.size() - 1);
//
//        Node * nd = _tree->_root->_left_child;
//
//        // sanity check: first node should be the only child of the root node
//        assert(nd->_right_sib == 0);
//
//        // Push nd onto back of queue
//        q.push(nd);
//
//        while (!q.empty()) {
//            // pop nd off front of queue
//            nd = q.front(); q.pop();
//
//            // and push it onto the stack
//            _tree->_levelorder.push_back(nd);
//
//            // add all children of nd to back of queue
//            Node * child = nd->_left_child;
//            if (child) {
//                q.push(child);
//                child = child->_right_sib;
//                while (child) {
//                    q.push(child);
//                    child = child->_right_sib;
//                }
//            }
//        }   // end while loop
//    }
//
//    inline void TreeManip::renumberInternals() {
//        assert(_tree->_preorder.size() > 0);
//
//        // Renumber internal nodes in postorder sequence
//        unsigned curr_internal = _tree->_nleaves;
//        for (auto nd : boost::adaptors::reverse(_tree->_preorder)) {
//            if (nd->_left_child) {
//                // nd is an internal node
//                nd->_number = curr_internal++;
//            }
//        }
//
//        // Root node is not included in _tree->_preorder, so if the root node
//        // is an internal node we need to number it here
//        if (_tree->_is_rooted)
//            _tree->_root->_number = curr_internal++;
//
//        _tree->_ninternals = curr_internal - _tree->_nleaves;
//
//        _tree->_unused_nodes.clear();
//        for (; curr_internal < (unsigned)_tree->_nodes.size(); curr_internal++) {
//            Node * curr = &_tree->_nodes[curr_internal];
//            putUnusedNode(curr);
//            assert(curr->_number == -1);
//            curr->_number = curr_internal;
//        }
//
//    }
//
//    inline bool TreeManip::canHaveSibling(Node * nd, bool rooted, bool allow_polytomies) {
//        assert(nd);
//        if (!nd->_parent) {
//            // trying to give root node a sibling
//            return false;
//        }
//
//        if (allow_polytomies)
//            return true;
//
//        bool nd_can_have_sibling = true;
//        if (nd != nd->_parent->_left_child) {
//            if (nd->_parent->_parent) {
//                // trying to give a sibling to a sibling of nd, and nd's parent is not the root
//                nd_can_have_sibling = false;
//            }
//            else {
//                if (rooted) {
//                    // root node has exactly 2 children in rooted trees
//                    nd_can_have_sibling = false;
//                }
//                else if (nd != nd->_parent->_left_child->_right_sib) {
//                    // trying to give root node more than 3 children
//                    nd_can_have_sibling = false;
//                }
//            }
//        }
//
//        return nd_can_have_sibling;
//    }
//
//    inline void TreeManip::rerootAtNodeNumber(int node_number) {
//        // Locate node having _number equal to node_number
//        Node * nd = 0;
//        for (auto & curr : _tree->_nodes) {
//            if (curr._number == node_number) {
//                nd = &curr;
//                break;
//            }
//        }
//
//        if (!nd)
//            throw XStrom(boost::str(boost::format("no node found with number equal to %d") % node_number));
//
//        if (nd != _tree->_root) {
//            if (nd->_left_child)
//                throw XStrom(boost::str(boost::format("cannot currently root trees at internal nodes (e.g. node %d)") % nd->_number));
//            rerootAtNode(nd);
//        }
//    }
//
//    inline void TreeManip::rerootAtNode(Node * prospective_root) {
//        Node * a = prospective_root;
//        Node * b = prospective_root->_parent;
//        Node * c = 0;
//        Node * d = 0;
//        Node * p = 0;
//        a->_parent = 0;
//        double tmp_edgelen  = 0.0;
//        double prev_edgelen = a->getEdgeLength();
//
//        while (b) {
//            // Prune node a from b
//            if (a == b->_left_child) {
//                if (a->_right_sib) {
//                    b->_left_child = a->_right_sib;
//                    a->_right_sib = 0;
//                }
//                else {
//                    b->_left_child = 0;
//                }
//            }
//            else {
//                c = b->_left_child;
//                while (c->_right_sib != a)
//                    c = c->_right_sib;
//                d = a->_right_sib;
//                c->_right_sib = d;
//                a->_right_sib = 0;
//            }
//
//            // Graft node b onto node a (but don't unhook node b from its parent just yet)
//            if (a->_left_child) {
//                c = a->_left_child;
//                while (c->_right_sib)
//                    c = c->_right_sib;
//                c->_right_sib = b;
//            }
//            else {
//                a->_left_child = b;
//            }
//
//            // Rotate
//            p = a;
//            a = b;
//            b = b->_parent;
//            a->_parent = p;
//
//            // Swap nd's edge length with its new parent's edge length
//            tmp_edgelen = a->getEdgeLength();
//            a->setEdgeLength(prev_edgelen);
//            prev_edgelen = tmp_edgelen;
//        }
//        prospective_root->setEdgeLength(0.0);
//        _tree->_root = prospective_root;
//        refreshPreorder();
//        refreshLevelorder();
//    }
//
//    inline void TreeManip::buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies) {
//        _tree.reset(new Tree());
//        _tree->_is_rooted = rooted;
//
//        std::set<unsigned> used; // used to ensure that no two leaf nodes have the same number
//        unsigned curr_leaf = 0;
//        unsigned num_edge_lengths = 0;
//        unsigned curr_node_index = 0;
//
//        // Remove comments from the supplied newick string
//        std::string commentless_newick = newick;
//        stripOutNexusComments(commentless_newick);
//
//        // Resize the _nodes vector
//        _tree->_nleaves = countNewickLeaves(commentless_newick);
//        if (_tree->_nleaves < 4)
//            throw XStrom("Expecting newick tree description to have at least 4 leaves");
//        unsigned max_nodes = 2*_tree->_nleaves - (rooted ? 0 : 2);
//        _tree->_nodes.resize(max_nodes);
//        for (auto & nd : _tree->_nodes )
//            nd._number = -1;
//
//        try {
//            // Root node
//            Node * nd = &_tree->_nodes[curr_node_index];
//            _tree->_root = nd;
//
//            if (_tree->_is_rooted) {
//                nd = &_tree->_nodes[++curr_node_index];
//                nd->_parent = &_tree->_nodes[curr_node_index - 1];
//                nd->_parent->_left_child = nd;
//            }
//
//            // Some flags to keep track of what we did last
//            enum {
//                Prev_Tok_LParen		= 0x01,	// previous token was a left parenthesis ('(')
//                Prev_Tok_RParen		= 0x02,	// previous token was a right parenthesis (')')
//                Prev_Tok_Colon		= 0x04,	// previous token was a colon (':')
//                Prev_Tok_Comma		= 0x08,	// previous token was a comma (',')
//                Prev_Tok_Name		= 0x10,	// previous token was a node name (e.g. '2', 'P._articulata')
//                Prev_Tok_EdgeLen	= 0x20	// previous token was an edge length (e.g. '0.1', '1.7e-3')
//            };
//            unsigned previous = Prev_Tok_LParen;
//
//            // Some useful flag combinations
//            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
//            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
//            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
//            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
//            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);
//
//            // Set to true while reading an edge length
//            bool inside_edge_length = false;
//            std::string edge_length_str;
//            unsigned edge_length_position = 0;
//
//            // Set to true while reading a node name surrounded by (single) quotes
//            bool inside_quoted_name = false;
//
//            // Set to true while reading a node name not surrounded by (single) quotes
//            bool inside_unquoted_name = false;
//
//            // Set to start of each node name and used in case of error
//            unsigned node_name_position = 0;
//
//            // loop through the characters in newick, building up tree as we go
//            unsigned position_in_string = 0;
//            for (auto ch : commentless_newick) {
//                position_in_string++;
//
//                if (inside_quoted_name) {
//                    if (ch == '\'') {
//                        inside_quoted_name = false;
//                        node_name_position = 0;
//                        if (!nd->_left_child) {
//                            extractNodeNumberFromName(nd, used);
//                            curr_leaf++;
//                        }
//                        previous = Prev_Tok_Name;
//                    }
//                    else if (iswspace(ch))
//                        nd->_name += ' ';
//                    else
//                        nd->_name += ch;
//
//                    continue;
//                }
//                else if (inside_unquoted_name) {
//                    if (ch == '(')
//                        throw XStrom(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));
//
//                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
//                        inside_unquoted_name = false;
//
//                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
//                        if (!(previous & Name_Valid))
//                            throw XStrom(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));
//
//                        if (!nd->_left_child) {
//                            extractNodeNumberFromName(nd, used);
//                            curr_leaf++;
//                        }
//
//                        previous = Prev_Tok_Name;
//                    }
//                    else {
//                        nd->_name += ch;
//                        continue;
//                    }
//                }
//                else if (inside_edge_length) {
//                    if (ch == ',' || ch == ')' || iswspace(ch)) {
//                        inside_edge_length = false;
//                        edge_length_position = 0;
//                        extractEdgeLen(nd, edge_length_str);
//                        ++num_edge_lengths;
//                        previous = Prev_Tok_EdgeLen;
//                    }
//                    else {
//                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
//                        if (!valid)
//                            throw XStrom(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
//                        edge_length_str += ch;
//                        continue;
//                    }
//                }
//
//                if (iswspace(ch))
//                    continue;
//
//                switch(ch) {
//                    case ';':
//                        break;
//
//                    case ')':
//                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
//                        if (!nd->_parent)
//                            throw XStrom(boost::str(boost::format("Too many right parentheses at position %d in tree description") % position_in_string));
//
//                        // Expect right paren only after an edge length, a node name, or another right paren
//                        if (!(previous & RParen_Valid))
//                            throw XStrom(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));
//
//                        // Go down a level
//                        nd = nd->_parent;
//                        if (!nd->_left_child->_right_sib)
//                            throw XStrom(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
//                        previous = Prev_Tok_RParen;
//                        break;
//
//                    case ':':
//                        // Expect colon only after a node name or another right paren
//                        if (!(previous & Colon_Valid))
//                            throw XStrom(boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));
//                        previous = Prev_Tok_Colon;
//                        break;
//
//                    case ',':
//                        // Expect comma only after an edge length, a node name, or a right paren
//                        if (!nd->_parent || !(previous & Comma_Valid))
//                            throw XStrom(boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));
//
//                        // Check for polytomies
//                        if (!canHaveSibling(nd, rooted, allow_polytomies)) {
//                            throw XStrom(boost::str(boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % newick));
//                        }
//
//                        // Create the sibling
//                        curr_node_index++;
//                        if (curr_node_index == _tree->_nodes.size())
//                            throw XStrom(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _tree->_nodes.size() % _tree->_nleaves));
//                        nd->_right_sib = &_tree->_nodes[curr_node_index];
//                        nd->_right_sib->_parent = nd->_parent;
//                        nd = nd->_right_sib;
//                        previous = Prev_Tok_Comma;
//                        break;
//
//                    case '(':
//                        // Expect left paren only after a comma or another left paren
//                        if (!(previous & LParen_Valid))
//                            throw XStrom(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % position_in_string));
//
//                        // Create new node above and to the left of the current node
//                        assert(!nd->_left_child);
//                        curr_node_index++;
//                        if (curr_node_index == _tree->_nodes.size())
//                            throw XStrom(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _tree->_nodes.size()));
//                        nd->_left_child = &_tree->_nodes[curr_node_index];
//                        nd->_left_child->_parent = nd;
//                        nd = nd->_left_child;
//                        previous = Prev_Tok_LParen;
//                        break;
//
//                    case '\'':
//                        // Encountered an apostrophe, which always indicates the start of a
//                        // node name (but note that node names do not have to be quoted)
//
//                        // Expect node name only after a left paren (child's name), a comma (sib's name)
//                        // or a right paren (parent's name)
//                        if (!(previous & Name_Valid))
//                            throw XStrom(boost::str(boost::format("Not expecting node name at position %d in tree description") % position_in_string));
//
//                        // Get the rest of the name
//                        nd->_name.clear();
//
//                        inside_quoted_name = true;
//                        node_name_position = position_in_string;
//
//                        break;
//
//                    default:
//                        // Get here if ch is not one of ();:,'
//
//                        // Expecting either an edge length or an unquoted node name
//                        if (previous == Prev_Tok_Colon) {
//                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
//                            inside_edge_length = true;
//                            edge_length_position = position_in_string;
//                            edge_length_str = ch;
//                        }
//                        else {
//                            // Get the node name
//                            nd->_name = ch;
//
//                            inside_unquoted_name = true;
//                            node_name_position = position_in_string;
//                        }
//                }   // end of switch statement
//            }   // loop over characters in newick string
//
//            if (inside_unquoted_name)
//                throw XStrom(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
//            if (inside_edge_length)
//                throw XStrom(boost::str(boost::format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
//            if (inside_quoted_name)
//                throw XStrom(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));
//
//            if (_tree->_is_rooted) {
//                refreshPreorder();
//                refreshLevelorder();
//            }
//            else {
//                // Root at leaf whose _number = 0
//                // refreshPreorder() and refreshLevelorder() called after rerooting
//                rerootAtNodeNumber(0);
//            }
//            renumberInternals();
//        }
//        catch(XStrom x) {
//            clear();
//            throw x;
//        }
//    }
//
//    inline void TreeManip::storeSplits(std::set<Split> & splitset) {
//        // Start by clearing and resizing all splits
//        for (auto & nd : _tree->_nodes) {
//            nd._split.resize(_tree->_nleaves);
//        }
//
//        // Now do a postorder traversal and add the bit corresponding
//        // to the current node in its parent node's split
//        for (auto nd : boost::adaptors::reverse(_tree->_preorder)) {
//            if (nd->_left_child) {
//                // add this internal node's split to splitset
//                splitset.insert(nd->_split);
//            }
//            else {
//                // set bit corresponding to this leaf node's number
//                nd->_split.setBitAt(nd->_number);
//            }
//
//            if (nd->_parent) {
//                // parent's bits are the union of the bits set in all its children
//                nd->_parent->_split.addSplit(nd->_split);
//            }
//        }
//    }
//
//    inline Node * TreeManip::randomEdge(Lot::SharedPtr lot) {
//        // Unrooted case:                        Rooted case:
//        //
//        // 2     3     4     5                   1     2     3     4
//        //  \   /     /     /                     \   /     /     /
//        //   \ /     /     /                       \ /     /     /
//        //    8     /     /                         7     /     /
//        //     \   /     /                           \   /     /
//        //      \ /     /                             \ /     /
//        //       7     /                               6     /
//        //        \   /                                 \   /
//        //         \ /                                   \ /
//        //          6   nleaves = 5                       5    nleaves = 4
//        //          |   preorder length = 7               |    preorder length = 7
//        //          |   num_edges = 7 - 0 = 7             |    num_edges = 7 - 1 = 6
//        //          1   choices: 2,3,4,5,6,7,8           root  choices: 1,2,3,4,6,7
//        //
//        // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
//        //
//        // Note: _preorder is actually a vector of Node *, but is shown here as a
//        // vector of integers solely to illustrate the algorithm below.
//
//        // If tree is rooted, add one to skip first node in _preorder vector (and to subtract
//        // one from number of valid edges) because, for rooted tree case, this first noe is
//        // an internal node whose edge is not a valid choice.
//        int rooted_offset = (_tree->_is_rooted ? 1 : 0);
//
//        int num_edges = (unsigned)_tree->_preorder.size() - rooted_offset;
//        double uniform_deviate = lot->uniform();
//        unsigned index_of_chosen = rooted_offset + (unsigned)std::floor(uniform_deviate*num_edges);
//
//        unsigned nodes_visited = 0;
//        Node * chosen_node = 0;
//        for (auto nd : _tree->_preorder) {
//            if (nodes_visited == index_of_chosen) {
//                chosen_node = nd;
//                break;
//            }
//            else
//                ++nodes_visited;
//        }
//        assert(chosen_node);
//        return chosen_node;
//    }
//
//    inline Node * TreeManip::randomInternalEdge(Lot::SharedPtr lot) {
//        // Unrooted case:                        Rooted case:
//        //
//        // 2     3     4     5                   1     2     3     4
//        //  \   /     /     /                     \   /     /     /
//        //   \ /     /     /                       \ /     /     /
//        //    8     /     /                         7     /     /
//        //     \   /     /                           \   /     /
//        //      \ /     /                             \ /     /
//        //       7     /                               6     /
//        //        \   /                                 \   /
//        //         \ /                                   \ /
//        //          6   nleaves = 5                       5    nleaves = 4
//        //          |   preorder length = 7               |    preorder length = 7
//        //          |   num_internal_edges = 7 - 5 = 2    |    num_internal_edges = 7 - 4 - 1 = 2
//        //          1   choose node 7 or node 8          root  choose node 6 or node 7
//        //
//        // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
//        //
//        // Note: _preorder is actually a vector of T *, but is shown here as a
//        // vector of integers solely to illustrate the algorithm below
//
//        int num_internal_edges = (unsigned)_tree->_preorder.size() - _tree->_nleaves - (_tree->_is_rooted ? 1 : 0);
//        if (num_internal_edges == 0) {
//            // Star tree: return hub node, which is the first node in the preorder sequence
//            return _tree->_preorder[0];
//        }
//
//        // Add one to skip first node in _preorder vector, which is an internal node whose edge
//        // is either a terminal edge (if tree is unrooted) or invalid (if tree is rooted)
//        double uniform_deviate = lot->uniform();
//        unsigned index_of_chosen = 1 + (unsigned)std::floor(uniform_deviate*num_internal_edges);
//
//        unsigned internal_nodes_visited = 0;
//        Node * chosen_node = 0;
//        for (auto nd : _tree->_preorder) {
//            if (nd->_left_child) {
//                if (internal_nodes_visited == index_of_chosen) {
//                    chosen_node = nd;
//                    break;
//                }
//                else
//                    ++internal_nodes_visited;
//            }
//        }
//        assert(chosen_node);
//        return chosen_node;
//    }
//
//    inline Node * TreeManip::randomChild(Lot::SharedPtr lot, Node * x, Node * avoid, bool parent_included) {
//        // Count number of children of x
//        unsigned n = 0;
//        Node * child = x->getLeftChild();
//        while (child) {
//            if (child != avoid)
//                n++;
//            child = child->getRightSib();
//    }
//
//        // Choose random child index
//        unsigned upper = n + (parent_included ? 1 : 0);
//        unsigned chosen = lot->randint(0,upper - 1);
//
//        // If chosen < n, then find and return that particular child
//        if (chosen < n) {
//            child = x->getLeftChild();
//            unsigned i = 0;
//            while (child) {
//                if (child != avoid && i == chosen)
//                    return child;
//                else if (child != avoid)
//                    i++;
//                child = child->getRightSib();
//            }
//        }
//
//        // If chosen equals n, then the parent was chosen, indicated by returning NULL
//        return NULL;
//    }
//
//    inline void TreeManip::LargetSimonSwap(Node * a, Node * b) {
//        // a and b are the ends of the selected 3-edge path in a Larget-Simon move
//        // The 3-edge path is indicated by parentheses around the nodes involved.
//        // x is always the parent of a
//        // y can be the parent of b (case 1) or the child of b (case 2)
//
//        Node * x = a->_parent;
//        assert(x);
//
//        Node * y = x->_parent;
//        assert(y);
//
//        if (y == b->_parent) {
//            // Case 1: y is the parent of b
//            //
//            //    (a) d  e             (b) d  e
//            //      \ | /                \ | /
//            //       \|/                  \|/
//            //       (x) f (b)            (x) f (a)    Swap a and b, leaving everything
//            //         \ | /                \ | /      else as is
//            //          \|/     ==>          \|/
//            //          (y)                  (y)
//            //           |                    |
//            //           |                    |
//            //           c                    c
//            //
//
//            // Detach a from tree
//            if (a == x->_left_child) {
//                x->_left_child = a->_right_sib;
//            } else {
//                Node * child = x->_left_child;
//                while (child->_right_sib != a)
//                    child = child->_right_sib;
//                child->_right_sib = a->_right_sib;
//            }
//            a->_parent = 0;
//            a->_right_sib = 0;
//
//            // Detach b from tree
//            if (b == y->_left_child) {
//                y->_left_child = b->_right_sib;
//            } else {
//                Node * child = y->_left_child;
//                while (child->_right_sib != b)
//                    child = child->_right_sib;
//                child->_right_sib = b->_right_sib;
//            }
//            b->_parent = 0;
//            b->_right_sib = 0;
//
//            // Reattach a to y
//            a->_right_sib = y->_left_child;
//            y->_left_child = a;
//            a->_parent = y;
//
//            // Reattach b to x
//            b->_right_sib = x->_left_child;
//            x->_left_child = b;
//            b->_parent = x;
//        }
//        else {
//            // Case 2: y is the child of b
//            //
//            //    (a) d  e             (a) f  c
//            //      \ | /                \ | /
//            //       \|/                  \|/
//            //       (x) f  c            (x) d  e    swap x's children (except a)
//            //         \ | /               \ | /     with y's children (except x)
//            //          \|/     ==>         \|/
//            //          (y)                 (y)
//            //           |                   |
//            //           |                   |
//            //          (b)                 (b)
//            assert(b == y->_parent);
//
//            // Remove x's children from tree and store in xchildren stack
//            std::stack<Node *> xchildren;
//            Node * child = x->_left_child;
//            Node * prevchild = 0;
//            while (child) {
//                if (child == a) {
//                    prevchild = child;
//                    child = child->_right_sib;
//                } else {
//                    if (child == x->_left_child) {
//                        x->_left_child = child->_right_sib;
//                        child->_right_sib = 0;
//                        child->_parent = 0;
//                        xchildren.push(child);
//                        child = x->_left_child;
//                    } else if (child->_right_sib) {
//                        prevchild->_right_sib = child->_right_sib;
//                        child->_right_sib = 0;
//                        child->_parent = 0;
//                        xchildren.push(child);
//                        child = prevchild->_right_sib;
//                    } else {
//                        assert(prevchild == a);
//                        a->_right_sib = 0;
//                        child->_parent = 0;
//                        xchildren.push(child);
//                        child = 0;
//                        prevchild = 0;
//                    }
//                }
//            }
//
//            // Remove y's children from tree and store in ychildren stack
//            std::stack<Node *> ychildren;
//            child = y->_left_child;
//            prevchild = 0;
//            while (child) {
//                if (child == x) {
//                    prevchild = child;
//                    child = child->_right_sib;
//                } else {
//                    if (child == y->_left_child) {
//                        y->_left_child = child->_right_sib;
//                        child->_right_sib = 0;
//                        child->_parent = 0;
//                        ychildren.push(child);
//                        child = y->_left_child;
//                    } else if (child->_right_sib) {
//                        prevchild->_right_sib = child->_right_sib;
//                        child->_right_sib = 0;
//                        child->_parent = 0;
//                        ychildren.push(child);
//                        child = prevchild->_right_sib;
//                    } else {
//                        assert(prevchild == x);
//                        x->_right_sib = 0;
//                        child->_parent = 0;
//                        ychildren.push(child);
//                        child = 0;
//                        prevchild = 0;
//                    }
//                }
//            }
//
//            // Reattach xchildren to y
//            while (!xchildren.empty()) {
//                Node * popped = xchildren.top();
//                xchildren.pop();
//                popped->_right_sib = y->_left_child;
//                y->_left_child = popped;
//                popped->_parent = y;
//            }
//
//            // Reattach ychildren to x
//            while (!ychildren.empty()) {
//                Node * popped = ychildren.top();
//                ychildren.pop();
//                popped->_right_sib = x->_left_child;
//                x->_left_child = popped;
//                popped->_parent = x;
//            }
//        }
//
//        refreshPreorder();
//        refreshLevelorder();
//    }
//
//    inline void TreeManip::selectAll() {
//        for (auto & nd : _tree->_nodes) {
//            nd.select();
//        }
//    }
//
//    inline void TreeManip::deselectAll() {
//        for (auto & nd : _tree->_nodes) {
//            nd.deselect();
//        }
//    }
//
//    inline void TreeManip::selectAllPartials() {
//        for (auto & nd : _tree->_nodes)
//            nd.selectPartial();
//    }
//
//    inline void TreeManip::deselectAllPartials() {
//        for (auto & nd : _tree->_nodes) {
//            nd.deselectPartial();
//        }
//    }
//
//    inline void TreeManip::selectAllTMatrices() {
//        for (auto & nd : _tree->_nodes)
//            nd.selectTMatrix();
//    }
//
//    inline void TreeManip::deselectAllTMatrices() {
//        for (auto & nd : _tree->_nodes) {
//            nd.deselectTMatrix();
//        }
//    }
//
//    inline void TreeManip::selectPartialsHereToRoot(Node * a) {
//        a->selectPartial();
//        while (a->_parent) {
//            a = a->_parent;
//            a->selectPartial();
//        }
//    }
//
//    inline void TreeManip::flipPartialsAndTMatrices() {
//        for (auto & nd : _tree->_nodes) {
//            if (nd.isSelPartial())
//                nd.flipPartial();
//
//            if (nd.isSelTMatrix())
//                nd.flipTMatrix();
//        }
//    }
//
//    inline Node * TreeManip::findNextPreorder(Node * nd) {
//        assert(nd);
//        Node * next = 0;
//        if (!nd->_left_child && !nd->_right_sib) {
//            // nd has no children and no siblings, so next preorder is the right sibling of
//            // the first ancestral node that has a right sibling.
//            Node * anc = nd->_parent;
//            while (anc && !anc->_right_sib)
//                anc = anc->_parent;
//            if (anc) {
//                // We found an ancestor with a right sibling
//                next = anc->_right_sib;
//            }
//            else {
//                // nd is last preorder node in the tree
//                next = 0;
//            }
//        }
//        else if (nd->_right_sib && !nd->_left_child) {
//            // nd has no children (it is a tip), but does have a sibling on its right
//            next = nd->_right_sib;
//        }
//        else if (nd->_left_child && !nd->_right_sib) {
//            // nd has children (it is an internal node) but no siblings on its right
//            next = nd->_left_child;
//        }
//        else {
//            // nd has both children and siblings on its right
//            next = nd->_left_child;
//        }
//        return next;
//    }
//
//    inline Node * TreeManip::findLeftSib(Node * nd) {
//        assert(nd);
//        assert(nd->_parent);
//        Node * child = nd->_parent->_left_child;
//        while (child && child->_right_sib != nd)
//            child = child->_right_sib;
//        return child;
//    }
//
//    inline Node * TreeManip::findRightmostChild(Node * nd) {
//        assert(nd);
//        Node * child = nd->getLeftChild();
//        while (child->getRightSib())
//            child = child->getRightSib();
//        return child;
//    }
//
//    inline Node * TreeManip::findLastPreorderInClade(Node * start) {
//        assert(start);
//        Node * curr = start;
//        Node * rchild = findRightmostChild(curr);
//        while (rchild) {
//            curr = rchild;
//            rchild = findRightmostChild(curr);
//        }
//        return curr;
//    }
//
//    inline void TreeManip::insertSubtreeOnLeft(Node * s, Node * u) {
//        assert(u);
//        assert(s);
//        s->_right_sib  = u->_left_child;
//        s->_parent     = u;
//        u->_left_child = s;
//    }
//
//    inline void TreeManip::insertSubtreeOnRight(Node * s, Node * u) {
//        assert(u);
//        assert(s);
//
//        s->_right_sib = 0;
//        s->_parent    = u;
//        if (u->_left_child) {
//            Node * u_rchild = findRightmostChild(u);
//            u_rchild->_right_sib = s;
//        }
//        else
//            u->_left_child = s;
//    }
//
//    inline void TreeManip::detachSubtree(Node * s) {
//        assert(s);
//        assert(s->_parent);
//
//        // Save pointers to relevant nodes
//        Node * s_leftsib  = findLeftSib(s);
//        Node * s_rightsib = s->_right_sib;
//        Node * s_parent   = s->_parent;
//
//        // Completely detach s and seal up the wound
//        s->_parent = 0;
//        s->_right_sib = 0;
//        if (s_leftsib)
//            s_leftsib->_right_sib = s_rightsib;
//        else
//            s_parent->_left_child = s_rightsib;
//    }
//
//    inline void TreeManip::rectifyNumInternals(int incr) {
//        assert(_tree->_nodes.size() == _tree->_unused_nodes.size() + _tree->_nleaves + _tree->_ninternals + incr);
//        _tree->_ninternals += incr;
//    }
//
//    inline void TreeManip::refreshNavigationPointers() {
//        refreshPreorder();
//        refreshLevelorder();
//    }
//
//    inline Node * TreeManip::getUnusedNode(Node * sought) {
//        assert(!_tree->_unused_nodes.empty());
//        Node * nd = 0;
//        if (sought) {
//            unsigned i = 0;
//            for (Node * und : _tree->_unused_nodes) {
//                if (und == sought) {
//                    nd = und;
//                    _tree->_unused_nodes.erase(_tree->_unused_nodes.begin()+i);
//                    break;
//                }
//                ++i;
//            }
//            assert(nd);
//        }
//        else {
//            nd = _tree->_unused_nodes.back();
//            _tree->_unused_nodes.pop_back();
//        }
//        nd->clearPointers();
//        return nd;
//    }
//
//    inline void TreeManip::putUnusedNode(Node * nd) {
//        nd->clearPointers();
//        _tree->_unused_nodes.push_back(nd);
//    }
//
//    inline bool TreeManip::isPolytomy(Node * nd) const {
//        Node * lchild = nd->_left_child;
//        assert(lchild);    // should only call this function for internal nodes
//
//        Node * rchild = lchild->_right_sib;
//        if (rchild && rchild->_right_sib)
//            return true;
//        return false;
//    }
//
//    inline void TreeManip::copyEdgeProportionsFrom(double TL, const std::vector<double> & new_props) {
//        assert(new_props.size() == _tree->_preorder.size());
//        unsigned i = 0;
//        for (auto nd : _tree->_preorder) {
//            nd->setEdgeLength(TL*new_props[i++]);
//        }
//    }
//
//    inline double TreeManip::copyEdgeProportionsTo(std::vector<double> & receptacle) const {
//        receptacle.resize(_tree->_preorder.size());
//        unsigned i = 0;
//        double TL = 0.0;
//        for (auto nd : _tree->_preorder) {
//            assert(nd->_edge_length > 0.0);
//            receptacle[i++] = nd->_edge_length;
//            TL += nd->_edge_length;
//        }
//
//        // Convert edge lengths into proportions by dividing each by TL
//        std::transform(receptacle.begin(), receptacle.end(), receptacle.begin(), [TL](double v) {return v/TL;});
//
//        return TL;
//    }
//
//#if defined(HPD_PWK_METHOD)
//
//    inline void TreeManip::saveParamNames(std::vector<std::string> & param_name_vect) const {
//        param_name_vect.push_back("TL");
//
//        // the number of transformed edge length proportions is one less than the number of edges
//        unsigned num_edgelen_proportions = (unsigned)(_tree->_preorder.size() - 1);
//        for (unsigned i = 0; i < num_edgelen_proportions; ++i) {
//            param_name_vect.push_back(boost::str(boost::format("edgeprop-%d") % (i+1)));
//        }
//    }
//
//    inline double TreeManip::logTransformEdgeLengths(std::vector<double> & param_vect) const {
//        double TL = 0.0;
//        std::vector<double> edge_length_proportions(_tree->_preorder.size());
//        unsigned i = 0;
//        for (auto nd : _tree->_preorder) {
//            assert(nd->_edge_length > 0.0);
//            edge_length_proportions[i++] = nd->_edge_length;
//            TL += nd->_edge_length;
//        }
//
//        // Convert edge lengths into proportions by dividing each by TL
//        std::transform(edge_length_proportions.begin(), edge_length_proportions.end(), edge_length_proportions.begin(), [TL](double v) {return v/TL;});
//
//        // Record everything in param_vect and compute the log of the Jacobian
//        double log_TL = log(TL);
//        double log_jacobian = log_TL;
//        param_vect.push_back(log_TL);
//
//        double first = edge_length_proportions[0];
//        double log_first = log(first);
//        log_jacobian += log_first;
//        for (unsigned i = 1; i < edge_length_proportions.size(); ++i) {
//            double p = edge_length_proportions[i];
//            double logp = log(p);
//            double transformed = logp - log_first;
//            log_jacobian += logp;
//            param_vect.push_back(transformed);
//        }
//
//        return log_jacobian;
//    }
//
//    inline double TreeManip::setEdgeLengthsFromLogTransformed(Eigen::VectorXd & param_vect, double TL, unsigned start_at, unsigned nedges) {
//        double phi = 1.0;
//        Node * nd = _tree->_preorder[0];
//        nd->setEdgeLength(1.0);
//        for (unsigned i = 0; i < nedges; ++i) {
//            double t = param_vect[start_at + i];
//            double r = exp(t);
//            nd = _tree->_preorder[i + 1];
//            nd->setEdgeLength(r);
//            phi += r;
//        }
//        double first = 1.0/phi;
//        double log_jacobian = log(TL);
//        for (auto nd : _tree->_preorder) {
//            double proportion = first*nd->getEdgeLength();
//            double edgelen = TL*proportion;
//            nd->setEdgeLength(edgelen);
//            log_jacobian += log(proportion);
//        }
//        return log_jacobian;
//    }
//
//#endif  //HPD
//
//}
