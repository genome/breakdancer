#pragma once

#include <map>

template<typename VertexType, typename WeightType>
class UndirectedWeightedGraph {
public:
    typedef std::map<VertexType, WeightType> EdgeMap;
    typedef std::map<VertexType, EdgeMap> VertexMap;

    typedef typename VertexMap::size_type size_type;
    typedef typename VertexMap::iterator iterator;
    typedef typename VertexMap::const_iterator const_iterator;

    size_type size() const {
        return _graph.size();
    }

    EdgeMap& operator[](VertexType const& v) {
        return _graph[v];
    }

    iterator begin() { return _graph.begin(); }
    iterator end() { return _graph.end(); }

    iterator find(VertexType const& v) { return _graph.find(v); }
    const_iterator find(VertexType const& v) const { return _graph.find(v); }
    const_iterator begin() const { return _graph.begin(); }
    const_iterator end() const { return _graph.end(); }

    void erase(iterator const& iter) {
        _graph.erase(iter);
    }

    void erase(VertexType const& v) {
        _graph.erase(v);
    }

private:
    VertexMap _graph;
};
