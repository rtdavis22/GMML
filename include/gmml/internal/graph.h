// Author: Robert Davis

#ifndef GMML_INTERNAL_GRAPH_H_
#define GMML_INTERNAL_GRAPH_H_

#include <algorithm>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

// This class represents a graph: a set of vertices and edges. Vertices are
// represented by the numbers 0...n-1. Each vertex has a random access list of
// other vertices that is guaranteed to be in increasing order.
class Graph {
  public:
    // AdjList is a random access list representing an adjacency list.
    typedef std::vector<size_t> AdjList;

    struct BFSResults {
        //~BFSResults();

        std::vector<size_t> *found;
        std::vector<int> *previous;
    };

    // Create an empty graph.
    Graph() : edges_(0) {}

    // Create a graph with the given number of vertices and no edges.
    explicit Graph(size_t size) : edges_(size) {}

    // TODO: Get rid of this if possible. If it stays, we must check that the
    // argument is valid.
    explicit Graph(const std::vector<std::vector<size_t> >& edges)
            : edges_(edges) {}

    virtual Graph *clone() const { return new Graph(edges_); }

    // The number of vertices of the graph.
    size_t size() const { return edges_.size(); }

    // Returns true if there is an edge between the given vertices.
    // Invariant: is_adjacent(v1, v2) <==> is_adjacent(v2, v1).
    bool is_adjacent(size_t v1, size_t v2) const {
        return v1 < size() && std::binary_search(edges_[v1].begin(),
                                                 edges_[v1].end(), v2);
    }

    // The function appends the vertices of another graph, and updates their
    // adjacency lists to reflect their new vertex indices.
    void append(const Graph& graph);

    // Remove the vertex from the graph, along with all incident edges.
    // It is more efficient to remove a vertex with a large index.
    void remove_vertex(size_t index);

    // Remove a range of vertices from the graph.
    void remove_vertex_range(size_t start_index, size_t count);

    // Remove a list of vertices. This function is more efficient when removing
    // many vertices.
    void remove_vertex_list(const std::vector<size_t>& vertices);

    // The function performs a breadth-first search, starting at the given
    // index. It returns a list of all the vertex indices it finds.
    BFSResults *bfs(size_t start_index) const {
        return bfs(start_index, edges_.size());
    }

    // The function performs a breadth-first search, but only goes out to a
    // specified distance.
    BFSResults *bfs(size_t start_index, size_t distance) const;

    // This is another form of breadth-first search. It only finds vertices
    // on the end_index side of the edge from start_index to end_index.
    std::vector<size_t> *edge_bfs(size_t start_index, size_t end_index) const;

    // This function returns a list of all simple cycles of the graph. It uses
    // an iterative depth-first search algorithm to avoid what could be a
    // large recursion stack.
    std::vector<std::vector<size_t> > *get_cycles() const;

    // Add an edge. If the indices are invalid or if the edge already exists,
    // false is returned.
    bool add_edge(size_t start_index, size_t end_index);

    // Remove an edge. If it doesn't exist, false is returned.
    bool remove_edge(size_t start_index, size_t end_index);

    // This function applies a map to all vertices, so that vertex v's label
    // is changed to permutation[v]. The argument must be a permutation of
    // 0...n-1. This doesn't affect the topology of the graph.
    void apply_map(const std::vector<size_t>& permutation);

    // Returns the adjacency list of a particular vertex.
    const AdjList& edges(size_t v) const { return edges_[v]; }

    const std::vector<std::vector<size_t> >& edges() const { return edges_; }

    void print() const;

  private:
    std::vector<AdjList> edges_;
};

inline void Graph::append(const Graph& graph) {
    size_t current_size = edges_.size();
    edges_.insert(edges_.end(), graph.edges_.begin(), graph.edges_.end());
    for (size_t i = current_size; i < edges_.size(); i++)
        for (size_t j = 0; j < edges_[i].size(); j++)
            edges_[i][j] += current_size;
}

inline bool Graph::add_edge(size_t start_index, size_t end_index) {
    size_t size = this->size();
    if (start_index >= size || end_index >= size)
        return false;
    std::vector<size_t>::iterator it;
    it = std::lower_bound(edges_[start_index].begin(),
                          edges_[start_index].end(),
                          end_index);
    if (it != edges_[start_index].end() && *it == end_index)
        return false;
    edges_[start_index].insert(it, end_index);
    it = std::lower_bound(edges_[end_index].begin(),
                          edges_[end_index].end(),
                          start_index);
    edges_[end_index].insert(it, start_index);
    return true;
}


inline bool Graph::remove_edge(size_t start_index, size_t end_index) {
    if (!is_adjacent(start_index, end_index))
        return false;
    std::vector<size_t>::iterator it;
    it = std::lower_bound(edges_[start_index].begin(),
                          edges_[start_index].end(),
                          end_index);
    edges_[start_index].erase(it);
    it = std::lower_bound(edges_[end_index].begin(),
                          edges_[end_index].end(),
                          start_index);
    edges_[end_index].erase(it);
    return true;
}

}  // namespace gmml

#endif  // GMML_INTERNAL_GRAPH_H_
