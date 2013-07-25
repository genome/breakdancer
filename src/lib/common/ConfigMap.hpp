#pragma once

#include "common/namespace.hpp"

#include <boost/container/flat_map.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>

template<typename K, typename V>
struct ConfigMap {
    typedef boost::container::flat_map<K, V> type;
};

BEGIN_NAMESPACE(boost)
BEGIN_NAMESPACE(serialization)

template<class Archive, class Key, class Value, class Compare, class Allocator>
inline void save(
        Archive& ar,
        boost::container::flat_map<Key, Value, Compare, Allocator> const& m,
        const unsigned int /* version */)
{
    boost::serialization::stl::save_collection(ar, m);
}

template<class Archive, class Key, class Value, class Compare, class Allocator>
inline void load(
        Archive& ar,
        boost::container::flat_map<Key, Value, Compare, Allocator>& m,
        const unsigned int /* version */)
{
    typedef boost::container::flat_map<Key, Value, Compare, Allocator> MapType;

    boost::serialization::stl::load_collection<
        Archive,
        MapType,
        boost::serialization::stl::archive_input_map<Archive, MapType>,
        boost::serialization::stl::no_reserve_imp<MapType>
        >
    (ar, m);

}


template<class Archive, class Key, class Value, class Compare, class Allocator>
inline void serialize(
        Archive& ar,
        boost::container::flat_map<Key, Value, Compare, Allocator>& m,
        const unsigned int version
        )
{
    boost::serialization::split_free(ar, m, version);
}

END_NAMESPACE(serialization)
END_NAMESPACE(boost)

