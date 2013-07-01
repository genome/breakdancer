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
        //track the number of links between two nodes
        //
        // This doesn't make a lot of sense to me. When v1 == v2 and v1 is not
        // in the map, both are set to one. If v1 is in the map, then we increment
        // twice. We should either double count or not. Doing a mixture of both is
        // silly. -ta
        //
        // We'll carry this forward until we get the green light to fix the bug.
        // -ta

        iterator v1_iter = find(v1);
        typename EdgeMap::iterator v2_iter;
        if (v1_iter != end() && (v2_iter = v1_iter->second.find(v2)) != v1_iter->second.end()) {
            ++v2_iter->second;
            ++_graph[v2][v1];
        }
        else {
            _graph[v1][v2] = 1;
            _graph[v2][v1] = 1;
        }
    }

private:
    VertexMap _graph;
};
