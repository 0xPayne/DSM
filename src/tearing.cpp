#include "../include/tearing.hpp"
#include <queue>
#include <unordered_map>

namespace SparseLib::Tearing {

std::vector<int> elsOrder(const std::vector<int> &scc_vertices,
                          const Sparse::CSCMatrix &matrix);

/// Apply Algorithm GR (Eades, Lin & Smyth, 1993) to each non-trivial SCC.
/// Trivial SCCs (single vertices) are passed through unchanged.
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

/// Algorithm GR -- greedy heuristic for the Feedback Arc Set (FAS) problem.
///
/// From: Eades, Lin & Smyth, "A fast and effective heuristic for the feedback
///       arc set problem", Information Processing Letters 47(6), 1993, pp. 319-323.
///       DOI: 10.1016/0020-0190(93)90079-O
///       https://scispace.com/pdf/a-fast-and-effective-heuristic-for-the-feedback-arc-set-100uza4fhy.pdf
///
/// Given a simple connected digraph G = (V, A), a vertex sequence
/// s = v_1, v_2, ..., v_n induces a feedback arc set R(s) consisting of all
/// leftward arcs v_j v_i where j > i (Section 2, p. 3). The FAS problem --
/// computing a minimum feedback arc set R*(G) -- is NP-hard [Karp, 1972].
///
/// Algorithm GR computes a vertex sequence s whose induced FAS satisfies
/// |R(s)| <= m/2 - n/6  (Theorem 2.1) in O(m) time and O(m) space
/// (Theorem 2.3), where m = |A| and n = |V|.
///
/// Procedure GR (Section 2, p. 4):
///   s_1 <- {}; s_2 <- {};
///   while G != {} do
///     while G contains a sink do
///       choose a sink u; s_2 <- u s_2; G <- G - u;
///     while G contains a source do
///       choose a source u; s_1 <- s_1 u; G <- G - u;
///     choose a vertex u for which delta(u) is a maximum;
///       s_1 <- s_1 u; G <- G - u;
///   s <- s_1 s_2.
///
/// Notation (per the paper):
///   d+(u) -- outdegree of u (number of arcs leaving u in current G)
///   d-(u) -- indegree  of u (number of arcs entering u in current G)
///   delta(u) = d+(u) - d-(u)
///   sink   -- vertex with d+(u) = 0
///   source -- vertex with d-(u) = 0

std::vector<int> elsOrder(const std::vector<int> &scc_vertices,
                          const Sparse::CSCMatrix &matrix) {
  int k = static_cast<int>(scc_vertices.size());
  if (k <= 1)
    return scc_vertices;

  // -- Build local subgraph G' for the SCC ------------------------------------
  // Map global vertex IDs -> local indices 0..k-1 so we can work with compact
  // arrays instead of sparse global IDs.
  std::unordered_map<int, int> g2l;
  g2l.reserve(k);
  for (int i = 0; i < k; ++i)
    g2l[scc_vertices[i]] = i;

  // succ[u] -- local successors of u  (arcs u -> v within the SCC)
  // pred[u] -- local predecessors of u (arcs v -> u within the SCC)
  // out_deg[u] = d+(u),  in_deg[u] = d-(u)
  // Self-loops and arcs leaving the SCC are excluded since they do not
  // affect the feedback arc set within this component.
  std::vector<std::vector<int>> succ(k);
  std::vector<std::vector<int>> pred(k);
  std::vector<int> out_deg(k, 0);
  std::vector<int> in_deg(k, 0);

  for (int i = 0; i < k; ++i) {
    int u = scc_vertices[i]; // global ID of vertex i
    for (int p = matrix.col_ptrs[u]; p < matrix.col_ptrs[u + 1]; ++p) {
      int v = matrix.row_indices[p]; // global ID of neighbour
      if (v == u)
        continue; // self-loop -- skip
      auto it = g2l.find(v);
      if (it == g2l.end())
        continue; // arc leaves the SCC -- skip
      int li = i, lj = it->second;
      succ[li].push_back(lj); // u -> v
      pred[lj].push_back(li); // v <- u
      ++out_deg[li];           // d+(u)++
      ++in_deg[lj];            // d-(v)++
    }
  }

  // -- Initialise s_1, s_2, and the sink/source queues
  std::vector<bool> alive(k, true); // tracks vertices still in G
  int remaining = k;

  // Seed queues with initial sinks (d+(u)=0) and sources (d-(u)=0).
  std::queue<int> sink_q;   // vertices with d+(u) = 0
  std::queue<int> source_q; // vertices with d-(u) = 0

  for (int i = 0; i < k; ++i) {
    if (out_deg[i] == 0)
      sink_q.push(i);
    else if (in_deg[i] == 0)
      source_q.push(i);
  }

  std::vector<int> s_left;  // s_1 -- built by appending sources / max-delta(u) vertices
  std::vector<int> s_right; // s_2 -- sinks appended here, reversed at the end
  s_left.reserve(k);
  s_right.reserve(k);

  // -- Vertex removal: G <- G - u (per the paper's notation) 
  // When vertex u is removed from G, update degrees of adjacent vertices:
  //   - For each successor v: d-(v) decreases. If d-(v) = 0, v is now a source.
  //   - For each predecessor v: d+(v) decreases. If d+(v) = 0, v is now a sink.
  // This is the mechanism by which removing a max-delta(u) vertex can cascade
  // into further sink/source removals in the next iteration of the main loop.
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

  // -- Main loop: "while G != {} do"
  while (remaining > 0) {

    // "while G contains a sink do { choose a sink u; s_2 <- u s_2; G <- G - u }"
    // A sink has d+(u) = 0: all its arcs point inward, so placing it at the
    // end of the vertex sequence ensures those arcs are not leftward (not in R(s)).
    // (We append here and reverse s_2 at the end to achieve prepend semantics.)
    while (!sink_q.empty()) {
      int u = sink_q.front();
      sink_q.pop();
      if (!alive[u])
        continue; // already removed (stale queue entry)
      remove_vertex(u);
      s_right.push_back(u);
    }

    // "while G contains a source do { choose a source u; s_1 <- s_1 u; G <- G - u }"
    // A source has d-(u) = 0: all its arcs point outward, so placing it at
    // the beginning of the vertex sequence ensures those arcs are not leftward.
    while (!source_q.empty()) {
      int u = source_q.front();
      source_q.pop();
      if (!alive[u])
        continue;
      remove_vertex(u);
      s_left.push_back(u);
    }

    // "choose a vertex u for which delta(u) is a maximum; s_1 <- s_1 u; G <- G - u"
    // This vertex has the greatest excess of outgoing over incoming arcs,
    // so appending it to s_1 minimises the number of leftward arcs in R(s).
    if (remaining > 0) {
      int best = -1;
      int best_delta = INT_MIN;
      for (int i = 0; i < k; ++i) {
        if (!alive[i])
          continue;
        int delta = out_deg[i] - in_deg[i]; // delta(u) = d+(u) - d-(u)
        if (delta > best_delta) {
          best_delta = delta;
          best = i;
        }
      }
      remove_vertex(best);
      s_left.push_back(best);
    }
  }

  // -- Build final vertex sequence s = s_1 s_2
  // s_left already holds s_1 in order.
  // s_right holds s_2 in reverse,
  // so we iterate it backwards to reconstruct the correct s_2 order.
  // Map local indices back to global vertex IDs.
  std::vector<int> result;
  result.reserve(k);
  for (int u : s_left)
    result.push_back(scc_vertices[u]);
  for (int i = static_cast<int>(s_right.size()) - 1; i >= 0; --i)
    result.push_back(scc_vertices[s_right[i]]);

  return result;
}

} // namespace SparseLib::Tearing
