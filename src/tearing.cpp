// CPSC 482 : DSM Project
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/tearing.hpp"
#include <queue>
#include <unordered_map>

namespace SparseLib::Tearing {

std::vector<int> elsOrder(const std::vector<int> &scc_vertices,
                          const Sparse::CSCMatrix &matrix);

// Apply Algorithm GR to each non-trivial SCC.
std::vector<std::vector<int>>
tearAllSCCs(const std::vector<std::vector<int>> &sccs,
            const Sparse::CSCMatrix &matrix) {
  std::vector<std::vector<int>> result;
  result.reserve(sccs.size());

  for (const auto &scc : sccs) {
    if (scc.size() <= 1)
      result.push_back(scc);
    else
      result.push_back(elsOrder(scc, matrix));
  }

  return result;
}

// Algorithm GR -- greedy FAS heuristic (Eades, Lin & Smyth, 1993)
std::vector<int> elsOrder(const std::vector<int> &scc_vertices,
                          const Sparse::CSCMatrix &matrix) {
  int k = static_cast<int>(scc_vertices.size());
  if (k <= 1)
    return scc_vertices;

  // Build local subgraph for the SCC 
  std::unordered_map<int, int> g2l;
  g2l.reserve(k);
  for (int i = 0; i < k; ++i)
    g2l[scc_vertices[i]] = i;

  std::vector<std::vector<int>> succ(k);
  std::vector<std::vector<int>> pred(k);
  std::vector<int> out_deg(k, 0);
  std::vector<int> in_deg(k, 0);

  for (int i = 0; i < k; ++i) {
    int u = scc_vertices[i];
    for (int p = matrix.col_ptrs[u]; p < matrix.col_ptrs[u + 1]; ++p) {
      int v = matrix.row_indices[p]; 
      if (v == u)
        continue; 
      auto it = g2l.find(v);
      if (it == g2l.end())
        continue;
      int li = i, lj = it->second;
      succ[li].push_back(lj); 
      pred[lj].push_back(li);
      ++out_deg[li];           
      ++in_deg[lj];          
    }
  }

  // Initialise s_1, s_2, and the sink/source queues
  std::vector<bool> alive(k, true);
  int remaining = k;

  std::queue<int> sink_q;
  std::queue<int> source_q;

  for (int i = 0; i < k; ++i) {
    if (out_deg[i] == 0)
      sink_q.push(i);
    else if (in_deg[i] == 0)
      source_q.push(i);
  }

  std::vector<int> s_left; 
  std::vector<int> s_right; 
  s_left.reserve(k);
  s_right.reserve(k);

  // Remove vertex u from G and update neighbour degrees
  auto remove_vertex = [&](int u) {
    alive[u] = false;
    --remaining;
    for (int v : succ[u]) {
      if (!alive[v])
        continue;
      --in_deg[v];
      if (in_deg[v] == 0)
        source_q.push(v);
    }
    for (int v : pred[u]) {
      if (!alive[v])
        continue;
      --out_deg[v];
      if (out_deg[v] == 0)
        sink_q.push(v);
    }
  };

  // Main loop
  while (remaining > 0) {

    while (!sink_q.empty()) {
      int u = sink_q.front();
      sink_q.pop();
      if (!alive[u])
        continue; 
      remove_vertex(u);
      s_right.push_back(u);
    }

    while (!source_q.empty()) {
      int u = source_q.front();
      source_q.pop();
      if (!alive[u])
        continue;
      remove_vertex(u);
      s_left.push_back(u);
    }

    if (remaining > 0) {
      int best = -1;
      int best_delta = INT_MIN;
      for (int i = 0; i < k; ++i) {
        if (!alive[i])
          continue;
        int delta = out_deg[i] - in_deg[i];
        if (delta > best_delta) {
          best_delta = delta;
          best = i;
        }
      }
      remove_vertex(best);
      s_left.push_back(best);
    }
  }

  // Build final sequence s = s_1 s_2 
  std::vector<int> result;
  result.reserve(k);
  for (int u : s_left)
    result.push_back(scc_vertices[u]);
  for (int i = static_cast<int>(s_right.size()) - 1; i >= 0; --i)
    result.push_back(scc_vertices[s_right[i]]);

  return result;
}

} // namespace SparseLib::Tearing
