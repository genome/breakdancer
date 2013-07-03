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

    void erase_edge(VertexType const& src, VertexType const& dst) {
        iterator iter = find(src);
        if (iter != end())
            iter->second.erase(dst);
    }

    void increment_edge_weight(VertexType const& v1, VertexType const& v2) {
        ++_graph[v1][v2];
        if (v1 != v2) {
            ++_graph[v2][v1];
        }
    }

    void clear() {
        _graph.clear();
    }

private:
    VertexMap _graph;
};
