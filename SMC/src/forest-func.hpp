#pragma once

namespace proj {

// put some forest functions in here to avoid circular include issues
    
#if defined(LAZY_COPYING)
    inline void Forest::addIncrAndJoin(double incr, const Split & lsplit, const Split & rsplit, ForestExtension & gfx) {
        // Identify the two nodes to join
        Node * first_node = nullptr;
        Node * second_node = nullptr;
        for (auto nd : _lineages) {
            const Split & s = nd->_split;
            if (s == lsplit) {
                first_node = nd;
            }
            else if (s == rsplit) {
                second_node = nd;
            }
            if (first_node && second_node)
                break;
        }
        
        assert(first_node);
        assert(second_node);
        
        // Increment all lineages
        const G::merge_vect_t mergers = gfx.getMergers();
        if (mergers.size() > 0) {
            double starting_forest_height = _forest_height;
            double partial_incr = 0.0;
            for (auto m : mergers) {
                // m is a tuple comprising:
                // 0: height
                // 1: first species to merge
                // 2: second species to merge
                double merge_height = get<0>(m);
                G::species_t left_spp = get<1>(m);
                G::species_t right_spp = get<2>(m);
                G::species_t anc_spp = (left_spp | right_spp);
                
                // Bring all _lineages (and _forest_height) up to merge_height
                partial_incr = merge_height - _forest_height;
                advanceAllLineagesBy(partial_incr);
                                
                // Merge left_spp and right_spp
                mergeSpecies(left_spp, right_spp, anc_spp);
            }
            
            // Bring all _lineages (and _forest_height) up to final height
            partial_incr = starting_forest_height + incr - _forest_height;
            assert(partial_incr > 0.0);
            advanceAllLineagesBy(partial_incr);
        }
        else {
            advanceAllLineagesBy(incr);
        }
                
        assert(first_node->getSpecies() == second_node->getSpecies());
        assert(first_node->getSpecies() > 0);

        // get next available node to serve as ancestral node
        Node * anc_node = &_nodes[G::_ntaxa + _ninternals];

        assert (anc_node->_parent==0);
        assert (anc_node->_number == -1);
        assert (anc_node->_right_sib == 0);
        
        anc_node->_number=G::_ntaxa+_ninternals;
        anc_node->_edge_length=0.0;
         _ninternals++;
        anc_node->_name += to_string(anc_node->_number);

        anc_node->_left_child=first_node;
        first_node->_right_sib=second_node;
        
        first_node->_parent = anc_node;
        second_node->_parent = anc_node;
        
        anc_node->_height = _forest_height;
        anc_node->_species = first_node->getSpecies();

        // Set anc_node split to union of the two child splits
        anc_node->_split.resize(G::_ntaxa);
        anc_node->_split += first_node->_split;
        anc_node->_split += second_node->_split;
        
        // Set partial to supplied (already-calculated) partial
        anc_node->_partial = gfx.getExtensionPartial();
        
        // Set species of anc node
        assert(first_node->_species == second_node->_species);
        assert(first_node->_species > 0);
        assert(second_node->_species > 0);
        anc_node->_species = first_node->_species;

        // Fix up _lineages
        updateNodeVector(_lineages, first_node, second_node, anc_node);
        refreshAllPreorders();
    }
#endif

}

