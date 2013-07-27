#pragma once

#include <map>
#include <ostream>

template<typename VertexType, typename WeightType>
class UndirectedWeightedGraph {
public:
    typedef std::map<VertexType, WeightType> EdgeMap;
    typedef std::map<VertexType, EdgeMap> VertexMap;

    typedef typename VertexMap::size_type size_type;
    typedef typename VertexMap::iterator iterator;
    typedef typename VertexMap::const_iterator const_iterator;

    size_type num_vertices() const {
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
        erase_edge_impl(src, dst);
        erase_edge_impl(dst, src);
    }

    void increment_edge_weight(VertexType const& v1, VertexType const& v2) {
        ++_graph[v1][v2];
        if (v1 != v2) {
            ++_graph[v2][v1];
        }
    }

    WeightType const& get_edge_weight_default(
            VertexType const& v1,
            VertexType const& v2,
            WeightType const& dflt)
    {
        const_iterator x = find(v1);
        if (x == end())
            return dflt;

        typename EdgeMap::const_iterator y = x->second.find(v2);
        if (y == x->second.end())
            return dflt;

        return y->second;
    }

    void clear() {
        _graph.clear();
    }

    template<typename VT, typename WT>
    friend std::ostream& operator<<(std::ostream& s, UndirectedWeightedGraph<VT, WT> const& g) {
        typedef UndirectedWeightedGraph<VT, WT> Graph;
        typedef typename Graph::const_iterator EI;
        typedef typename Graph::EdgeMap::const_iterator EJ;
        for (EI i = g.begin(); i != g.end(); ++i) {
            for (EJ j = i->second.begin(); j != i->second.end(); ++j) {
                s << "(" << i->first << ", " << j->first << "): " << j->second << "\n";
            }
        }

        return s;
    }

private:
    void erase_edge_impl(VertexType const& src, VertexType const& dst) {
        iterator iter = find(src);
        if (iter != end()) {
            iter->second.erase(dst);
            //if (iter->second.empty())
                //_graph.erase(iter);
        }
    }
private:
    VertexMap _graph;
};
