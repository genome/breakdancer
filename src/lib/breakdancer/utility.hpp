#pragma once

#include <boost/container/flat_map.hpp>

#include <utility>

template<typename K, typename V>
struct ConfigMap {
    typedef boost::container::flat_map<K, V> type;
};

// This function is similar to a.insert(b.begin(), b.end()), but it
// resolves conflicts by applying the supplied Combiner operation.
// The combiner is expected to be a binary function taking two
// MapType::mapped_type objects and returning one.
template<typename MapType, typename Combiner>
void merge_maps(MapType& a, MapType const& b, Combiner const& combiner) {
    typedef typename MapType::iterator IterType;
    typedef typename MapType::const_iterator ConstIterType;
    for (ConstIterType b_iter = b.begin(); b_iter != b.end(); ++b_iter) {
        std::pair<IterType, bool> inserted = a.insert(*b_iter);
        if (!inserted.second)
            inserted.first->second = combiner(inserted.first->second, b_iter->second);
    }
}
