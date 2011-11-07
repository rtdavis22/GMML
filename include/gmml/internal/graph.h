// Author: Robert Davis

#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <vector>

#include "utilities.h"

namespace gmml {

// This class represents a graph: a set of vertices and edges. Vertices are
// represented by the numbers 0...n-1. Each vertex has a random access list of
// other vertices that is guaranteed to be in increasing order.
class Graph {
  public:
    // AdjList is a random access list representing an adjacency list.
    // It may be appropriate for it to be a template argument.
    typedef std::vector<size_t> AdjList;

    Graph() : edges_(0) {}
    explicit Graph(size_t size) : edges_(size) {}
    explicit Graph(const std::vector<std::vector<size_t> >& edges)
            : edges_(edges) {}

    Graph *clone() const { return new Graph(edges_); }

    size_t size() const { return edges_.size(); }

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
    std::vector<size_t> *bfs(size_t start_index) const {
        return bfs(start_index, edges_.size());
    }

    // The function performs a breadth-first search, but only goes out to a
    // specified distance.
    std::vector<size_t> *bfs(size_t start_index, size_t distance) const;

    // This is another form of breadth-first search. It only finds vertices
    // on the end_index side of the edge from start_index to end_index.
    std::vector<size_t> *edge_bfs(size_t start_index, size_t end_index) const;

    // This function returns a list of all simple cycles of the graph. It uses
    // an iterative depth-first search algorithm to avoid what could be a
    // large recursion stack.
    std::vector<std::vector<size_t> > *get_cycles() const;

    // Add an edge, and maintain sorted order
    void add_edge(size_t start_index, size_t end_index) {
        std::vector<size_t>::iterator it;
        it = std::lower_bound(edges_[start_index].begin(),
                              edges_[start_index].end(),
                              end_index);
        edges_[start_index].insert(it, end_index);
        it = std::lower_bound(edges_[end_index].begin(),
                              edges_[end_index].end(),
                              start_index);
        edges_[end_index].insert(it, start_index);
    }

    // This function applies a map to all vertices, so that vertex v's label
    // is changed to permutation[v]. The argument must be a permutation of
    // 0...n-1. This doesn't affect the topology of the graph.
    void apply_map(const std::vector<size_t>& permutation);

    const std::vector<size_t>& edges(size_t v) const { return edges_[v]; }
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

}  // namespace gmml

#endif  // GRAPH_H
