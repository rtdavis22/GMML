// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Robert Davis

#include "gmml/internal/graph.h"

#include <algorithm>
#include <deque>
#include <iostream>  // for print()
#include <queue>
#include <vector>

namespace gmml {

using std::deque;
using std::queue;
using std::vector;

void Graph::remove_vertex_list(const vector<size_t>& vertices) {
    vector<int> new_indices(edges_.size());
    for (int i = 0; i < vertices.size(); i++)
        new_indices[vertices[i]] = -1;

    int count = 0;
    for (int i = 0; i < new_indices.size(); i++) {
        if (new_indices[i] == -1)
            count++;
        else
            new_indices[i] = i - count;
    }

    vector<AdjList> new_edges(edges_.size() - count);
    int insert_index = 0;
    for (int i = 0; i < edges_.size(); i++) {
        if (new_indices[i] == -1)
            continue;
        for (int j = 0; j < edges_[i].size(); j++) {
            int new_dest = new_indices[edges_[i][j]];
            if (new_dest != -1)
                new_edges[insert_index].push_back(new_dest);
        }
        insert_index++;
    }
    edges_ = new_edges;
}

// This could probable made to be more efficient.
void Graph::remove_vertex(size_t vertex) {
    vector<size_t> vertices(1, vertex);
    remove_vertex_list(vertices);
}

vector<size_t> *Graph::bfs(size_t start_index, size_t distance) const {
    deque<bool> marked(edges_.size());
    vector<int> distances(edges_.size(), edges_.size());
    vector<int> *previous = new vector<int>(edges_.size(), -1);
    marked[start_index] = true;
    distances[start_index] = 0;
    queue<size_t> q;
    vector<size_t> *vertices_found = new vector<size_t>();
    vertices_found->push_back(start_index);
    q.push(start_index);
    while (!q.empty()) {
        size_t current_vertex = q.front();
        q.pop();
        for (size_t i = 0; i < edges_[current_vertex].size(); i++) {
            if (!marked[edges_[current_vertex][i]]) {
                marked[edges_[current_vertex][i]] = true;
                (*previous)[edges_[current_vertex][i]] = current_vertex;
                size_t current_distance = distances[current_vertex] + 1;
                distances[edges_[current_vertex][i]] = current_distance;
                if (current_distance <= distance)
                    vertices_found->push_back(edges_[current_vertex][i]);
                q.push(edges_[current_vertex][i]);
            }
        }
    }
    delete previous;
    return vertices_found;
}

vector<size_t> *Graph::edge_bfs(size_t start_index, size_t end_index) const {
    deque<bool> marked(edges_.size());
    marked[end_index] = true;
    queue<size_t> q;
    vector<size_t> *vertices_found = new vector<size_t>();
    vertices_found->push_back(end_index);
    for (size_t i = 0; i < edges_[end_index].size(); i++) {
        if (edges_[end_index][i] != start_index) {
            size_t current_vertex = edges_[end_index][i];
            vertices_found->push_back(current_vertex);
            marked[current_vertex] = true;
            q.push(current_vertex);
        }
    }
    while (!q.empty()) {
        size_t current_vertex = q.front();
        q.pop();
        for (size_t i = 0; i < edges_[current_vertex].size(); i++) {
            if (!marked[edges_[current_vertex][i]]) {
                vertices_found->push_back(edges_[current_vertex][i]);
                marked[edges_[current_vertex][i]] = true;
                q.push(edges_[current_vertex][i]);
            }
        }
    }
    return vertices_found;
}

// Iterative depth-first search by double-threading all vertices
vector<vector<size_t> > *Graph::get_cycles() const {
    vector<vector<size_t> > *cycles = new vector<vector<size_t> >();
    vector<int> parents(edges_.size());
    deque<bool> is_visited(edges_.size(), false);
    vector<deque<bool> > is_threaded;
    is_threaded.reserve(edges_.size());
    for (size_t i = 0; i < edges_.size(); i++) {
        is_threaded.push_back(deque<bool>(edges_[i].size(), false));
    }
    parents[0] = -1;

    int cur_vertex = 0;
    is_visited[cur_vertex] = true;
    while (true) {
        bool done = false;
        // Look for nodes with no threads
        for (size_t i = 0; i < edges_[cur_vertex].size(); i++) {
            vector<size_t>::const_iterator it;
            it = std::find(edges_[edges_[cur_vertex][i]].begin(),
                           edges_[edges_[cur_vertex][i]].end(),
                           cur_vertex);
            // The index of cur_vertex in edges_[cur_vertex][i]
            size_t index = std::distance(
                    edges_[edges_[cur_vertex][i]].begin(), it);
            if (!is_threaded[cur_vertex][i] &&
                    !is_threaded[edges_[cur_vertex][i]][index]) {
                is_threaded[cur_vertex][i] = true;
                if (is_visited[edges_[cur_vertex][i]]) {
                    is_threaded[edges_[cur_vertex][i]][index] = true;
                    vector<size_t> cycle;
                    cycle.push_back(cur_vertex);
                    while (cycle.back() != edges_[cur_vertex][i])
                        cycle.push_back(parents[cycle.back()]);
                    cycles->push_back(cycle);
                } else {
                    parents[edges_[cur_vertex][i]] = cur_vertex;
                    cur_vertex = edges_[cur_vertex][i];
                    is_visited[cur_vertex] = true;
                }
                done = true;
                break;
            }
        }
        if (done) {
            continue;
        }
        // Look for nodes with an incoming but no outgoing thread
        for (size_t i = 0; i < edges_[cur_vertex].size(); i++) {
            if (!is_threaded[cur_vertex][i]) {
                is_threaded[cur_vertex][i] = true;
                cur_vertex = edges_[cur_vertex][i];
                is_visited[cur_vertex] = true;
                done = true;
                break;
            }
        }
        // Everything is already double-threaded
        if (!done) {
            break;
        }
    }
    return cycles;
}

void Graph::print() const {
    for (int i = 0; i < edges_.size(); i++) {
        std::cout << i << ": ";
        for (int j = 0; j < edges_[i].size(); j++)
            std::cout << edges_[i][j] << " ";
        std::cout << std::endl;
    }
}

}  // namespace gmml
