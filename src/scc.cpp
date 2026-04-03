// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/scc.hpp"

#include <unordered_set>
#include <vector>

namespace SparseLib::SCC {


// Iterative DFS with hooks for pre/post visit and edge events
template <typename PreVisit, typename PreEdge, typename PostEdge,
          typename PostVisit>
void genericDFS(int start, const Sparse::CSCMatrix &matrix, PreVisit &pre_visit,
                PreEdge &pre_edge, PostEdge &post_edge, PostVisit &post_visit) {

  // tracks vertex, adjacency position, and parent for stack-based DFS
  struct Frame {
    int u;
    int p;
    int parent;
  };

  std::vector<Frame> stack;

  stack.reserve(1024);

  stack.push_back({start, matrix.col_ptrs[start], -1});
  pre_visit(start);

  while (!stack.empty()) {
    Frame &f = stack.back();
    int u = f.u;
    int end = matrix.col_ptrs[u + 1];

    if (f.p < end) {
      int v = matrix.row_indices[f.p];
      f.p++;

      if (pre_edge(u, v)) {
        pre_visit(v);
        stack.push_back({v, matrix.col_ptrs[v], u});
      }
    } else {  
      int parent = f.parent;
      stack.pop_back();

      post_visit(u);

      if (parent != -1) {
        post_edge(parent, u);
      }
    }
  }
}

std::vector<std::vector<int>> tarjanSCC(const Sparse::CSCMatrix &matrix) {
  const int n = matrix.cols;
  std::vector<int> index(n, -1);
  std::vector<int> lowLink(n, 0);
  std::vector<bool> onstack(n, false);
  std::vector<int> stack;
  std::vector<std::vector<int>> sccs;

  int currentIndex = 0;

  // node discovered: assign index and push to stack
  auto pre_visit = [&](int v) {
    index[v] = lowLink[v] = currentIndex++;
    stack.push_back(v);
    onstack[v] = true;
  };

  // edge discovered: recurse if unvisited, otherwise update lowLink
  auto pre_edge = [&](int v, int w) -> bool {
    if (index[w] == -1) {
      return true;
    } else if (onstack[w]) {
      if (index[w] < lowLink[v])
        lowLink[v] = index[w];
    }
    return false;
  };

  // returning from child: propagate lowLink
  auto post_edge = [&](int v, int w) {
    if (lowLink[w] < lowLink[v])
      lowLink[v] = lowLink[w];
  };

  // node finished: if root, pop stack to form an SCC
  auto post_visit = [&](int v) {
    if (lowLink[v] == index[v]) {
      std::vector<int> component;
      while (true) {
        int w = stack.back();
        stack.pop_back();
        onstack[w] = false;
        component.push_back(w);
        if (v == w)
          break;
      }
      sccs.push_back(std::move(component));
    }
  };

  for (int v = 0; v < n; v++) {
    if (index[v] == -1) {
      genericDFS(v, matrix, pre_visit, pre_edge, post_edge, post_visit);
    }
  }
  return sccs;
}

std::vector<std::vector<int>> kosarajuSCC(const Sparse::CSCMatrix &matrix) {
  const int n = matrix.cols;
  std::vector<int> finish_stack;
  std::vector<bool> visited(n, false);

  // Phase 1
  auto p1_pre_visit = [&](int v) { visited[v] = true; };
  auto p1_pre_edge  = [&](int /*u*/, int v) { return !visited[v]; };
  auto p1_post_edge = [&](int /*parent*/, int /*u*/) {};
  auto p1_post_visit = [&](int v) { finish_stack.push_back(v); };

  for (int i = 0; i < n; ++i) {
    if (!visited[i]) {
      genericDFS(i, matrix, p1_pre_visit, p1_pre_edge, p1_post_edge, p1_post_visit);
    }
  }

  // reverse the graph
  Sparse::CSCMatrix reversed = transpose(Sparse::CSCMatrix(matrix));

  // Phase 2
  std::fill(visited.begin(), visited.end(), false);
  std::vector<std::vector<int>> sccs;
  std::vector<int> current_scc;

  auto p2_pre_visit = [&](int v) { visited[v] = true; current_scc.push_back(v); };
  auto p2_post_visit = [&](int /*v*/) {};

  while (!finish_stack.empty()) {
    int v = finish_stack.back();
    finish_stack.pop_back();

    if (!visited[v]) {
      current_scc.clear();
      genericDFS(v, reversed, p2_pre_visit, p1_pre_edge, p1_post_edge, p2_post_visit);
      sccs.push_back(current_scc);
    }
  }

  return sccs;
}

Sparse::CSCMatrix condensationGraph(const Sparse::CSCMatrix &matrix,
                                    const std::vector<std::vector<int>> &sccs) {
  int n = matrix.cols;
  int m = static_cast<int>(sccs.size());
  std::vector<int> node_to_scc(n, -1);

  // compute a node -> component_id mapping
  for (int cid = 0; cid < m; cid++) {
    for (int v : sccs[cid])
      node_to_scc[v] = cid;
  }

  std::vector<std::unordered_set<int>> outs(m);

  for (int u = 0; u < n; u++) {
    int su = node_to_scc[u];
    if (su < 0)
      continue;
    for (int p = matrix.col_ptrs[u]; p < matrix.col_ptrs[u + 1]; p++) {
      int v = matrix.row_indices[p];
      int sv = node_to_scc[v];
      if (sv >= 0 && sv != su)
        outs[su].insert(sv);
    }
  }

  // build the CSC matrix of the condensation graph
  Sparse::CSCMatrix condense(m, m, 0);
  condense.col_ptrs.assign(m + 1, 0);
  for (int i = 0; i < m; i++) {
    condense.col_ptrs[i + 1] =
        condense.col_ptrs[i] + static_cast<int>(outs[i].size());
  }
  condense.nnz = condense.col_ptrs[m];
  condense.row_indices.assign(condense.nnz, -1);

  int i = 0;
  for (int col = 0; col < m; col++) {
    // Sorting for deterministic output
    std::vector<int> tmp(outs[col].begin(), outs[col].end());
    std::sort(tmp.begin(), tmp.end());
    for (int sv : tmp) {
      condense.row_indices[i] = sv;
      i++;
    }
  }

  return condense;
}

std::vector<int> topologicalSort(const Sparse::CSCMatrix &matrix) {
  // Kahn's algorithm: repeatedly remove zero-indegree nodes
  std::vector<int> result;

  int num_nodes = matrix.col_ptrs.size() - 1;
  std::vector<int> incoming_count(num_nodes, 0);

  for (int row_index : matrix.row_indices) {
    incoming_count[row_index] += 1;
  }

  std::unordered_set<int> s;

  for (int i = 0; i < num_nodes; i++) {
    if (!incoming_count[i]) {
      s.insert(i);
    }
  }

  while (!s.empty()) {
    auto it = s.begin();
    int n = *it;
    s.erase(n);
    result.push_back(n);

    for (int i = matrix.col_ptrs[n]; i < matrix.col_ptrs[n + 1]; i++) {
      int m = matrix.row_indices[i];
      incoming_count[m] -= 1;
      if (incoming_count[m] == 0) {
        s.insert(m);
      }
    }
  }

  if (result.size() != static_cast<size_t>(num_nodes)) {
    throw std::runtime_error("Graph is cyclic; topological sort impossible.");
  }

  return result;
}

} // namespace SparseLib::SCC
