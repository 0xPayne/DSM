// CPSC 482 : DSM Project 
// Authors: Simon Kraft, Joshua Payne, El Sall

#include "../include/banding.hpp"
#include <algorithm>
#include <cassert>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace SparseLib::Banding {

namespace {

struct LocalGraph {
  int k;                                  
  std::unordered_map<int, int> g2l;    
  std::vector<int> globals;                 
  std::vector<std::vector<int>> adj;      
  std::vector<int> deg;                     

  struct Edge { int src, dst; };
  std::vector<Edge> directed_edges;
};

// Build local undirected adjacency and directed edge list for an SCC.
LocalGraph buildLocalGraph(const std::vector<int> &scc_vertices,
                           const Sparse::CSCMatrix &matrix) {
  LocalGraph G;
  G.k = static_cast<int>(scc_vertices.size());
  G.globals = scc_vertices;
  G.g2l.reserve(G.k);
  for (int i = 0; i < G.k; ++i)
    G.g2l[scc_vertices[i]] = i;

  std::vector<std::unordered_set<int>> adj_set(G.k);

  for (int i = 0; i < G.k; ++i) {
    int u = scc_vertices[i];
    for (int p = matrix.col_ptrs[u]; p < matrix.col_ptrs[u + 1]; ++p) {
      int v = matrix.row_indices[p];
      if (v == u)
        continue;
      auto it = G.g2l.find(v);
      if (it == G.g2l.end())
        continue;
      int li = i, lj = it->second;
      G.directed_edges.push_back({li, lj});
      adj_set[li].insert(lj);
      adj_set[lj].insert(li);
    }
  }

  // convert to sorted vectors
  G.adj.resize(G.k);
  G.deg.resize(G.k);
  for (int i = 0; i < G.k; ++i) {
    G.adj[i].assign(adj_set[i].begin(), adj_set[i].end());
    G.deg[i] = static_cast<int>(G.adj[i].size());
  }

  // sort adj lists by degree ascending 
  for (int i = 0; i < G.k; ++i) {
    std::sort(G.adj[i].begin(), G.adj[i].end(),
              [&](int a, int b) { return G.deg[a] < G.deg[b]; });
  }

  return G;
}

int findPseudoPeripheral(const LocalGraph &G) {
  auto bfs = [&](int start) -> std::pair<std::vector<int>, int> {
    std::vector<int> dist(G.k, -1);
    std::queue<int> q;
    dist[start] = 0;
    q.push(start);
    std::vector<int> last_level;
    int max_dist = 0;
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      if (dist[u] > max_dist) {
        max_dist = dist[u];
        last_level.clear();
      }
      if (dist[u] == max_dist)
        last_level.push_back(u);
      for (int v : G.adj[u]) {
        if (dist[v] == -1) {
          dist[v] = dist[u] + 1;
          q.push(v);
        }
      }
    }
    return {last_level, max_dist};
  };

  int current = 0; 
  int prev_ecc = -1;

  for (int round = 0; round < 5; ++round) {
    auto [last_level, ecc] = bfs(current);
    if (ecc <= prev_ecc)
      break; 
    prev_ecc = ecc;
    int best = last_level[0];
    for (int v : last_level)
      if (G.deg[v] < G.deg[best])
        best = v;
    current = best;
  }

  return current;
}

// Reverse Cuthill-McKee (RCM)
std::vector<int> rcmOrder(const LocalGraph &G) {
  int start = findPseudoPeripheral(G);

  std::vector<int> order;
  order.reserve(G.k);
  std::vector<bool> visited(G.k, false);

  std::queue<int> q;
  visited[start] = true;
  q.push(start);

  while (!q.empty()) {
    int u = q.front();
    q.pop();
    order.push_back(u);

    for (int v : G.adj[u]) {
      if (!visited[v]) {
        visited[v] = true;
        q.push(v);
      }
    }
  }

  for (int i = 0; i < G.k; ++i) {
    if (!visited[i])
      order.push_back(i);
  }

  std::reverse(order.begin(), order.end());

  return order;
}

// Local FBM and TFBD computation for an SCC ordering
struct Metrics {
  int fbm;
  long long tfbd;
};

Metrics localMetrics(const std::vector<int> &order,
                     const std::vector<LocalGraph::Edge> &edges, int k) {
  std::vector<int> pos(k);
  for (int i = 0; i < k; ++i)
    pos[order[i]] = i;

  Metrics m{0, 0};
  for (const auto &e : edges) {
    int p_src = pos[e.src];
    int p_dst = pos[e.dst];
    if (p_dst < p_src) { 
      ++m.fbm;
      m.tfbd += p_src - p_dst;
    }
  }
  return m;
}

// Adjacent-swap hill climbing 
std::vector<int>
localRefine(const std::vector<int> &initial_order,
            const std::vector<LocalGraph::Edge> &edges, int k,
            int max_passes = 20) {
  std::vector<int> order = initial_order;

  std::vector<int> pos(k);
  for (int i = 0; i < k; ++i)
    pos[order[i]] = i;

  // incidence lists per vertex
  std::vector<std::vector<int>> incident(k);
  for (int ei = 0; ei < static_cast<int>(edges.size()); ++ei) {
    incident[edges[ei].src].push_back(ei);
    incident[edges[ei].dst].push_back(ei);
  }

  for (int pass = 0; pass < max_passes; ++pass) {
    bool improved = false;

    for (int i = 0; i < k - 1; ++i) {
      int a = order[i];
      int b = order[i + 1];

      int delta_fbm = 0;
      long long delta_tfbd = 0;

      auto process_edges = [&](int vertex) {
        for (int ei : incident[vertex]) {
          int s = edges[ei].src;
          int d = edges[ei].dst;

          if ((s == a && d == b) || (s == b && d == a))
            continue;

          int ps_old = pos[s];
          int pd_old = pos[d];

          int ps_new = (s == a) ? i + 1 : (s == b) ? i : ps_old;
          int pd_new = (d == a) ? i + 1 : (d == b) ? i : pd_old;

          bool fbm_old = (pd_old < ps_old);
          long long dist_old = fbm_old ? (ps_old - pd_old) : 0;

          bool fbm_new = (pd_new < ps_new);
          long long dist_new = fbm_new ? (ps_new - pd_new) : 0;

          if (fbm_new && !fbm_old)
            ++delta_fbm;
          else if (!fbm_new && fbm_old)
            --delta_fbm;

          delta_tfbd += dist_new - dist_old;
        }
      };

      process_edges(a);
      process_edges(b);

      for (int ei : incident[a]) {
        int s = edges[ei].src;
        int d = edges[ei].dst;
        if (!((s == a && d == b) || (s == b && d == a)))
          continue;

        int ps_old = pos[s], pd_old = pos[d];
        int ps_new = (s == a) ? i + 1 : i;
        int pd_new = (d == a) ? i + 1 : i;

        bool fbm_old = (pd_old < ps_old);
        long long dist_old = fbm_old ? (ps_old - pd_old) : 0;
        bool fbm_new = (pd_new < ps_new);
        long long dist_new = fbm_new ? (ps_new - pd_new) : 0;

        if (fbm_new && !fbm_old)
          ++delta_fbm;
        else if (!fbm_new && fbm_old)
          --delta_fbm;

        delta_tfbd += dist_new - dist_old;
      }

      if (delta_tfbd < 0 && delta_fbm <= 0) {
        std::swap(order[i], order[i + 1]);
        pos[a] = i + 1;
        pos[b] = i;
        improved = true;
      }
    }

    if (!improved)
      break;
  }

  return order;
}

} // anonymous namespace

std::vector<std::vector<int>>
bandAllSCCs(const std::vector<std::vector<int>> &torn_sccs,
            const Sparse::CSCMatrix &matrix) {
  std::vector<std::vector<int>> result;
  result.reserve(torn_sccs.size());

  for (const auto &scc : torn_sccs) {
    int k = static_cast<int>(scc.size());

    // Trivial SCCs:
    if (k <= 2) {
      result.push_back(scc);
      continue;
    }

    LocalGraph G = buildLocalGraph(scc, matrix);

    std::vector<int> torn_local(k);
    for (int i = 0; i < k; ++i)
      torn_local[i] = i;

    Metrics m_torn = localMetrics(torn_local, G.directed_edges, k);

    std::vector<int> refined_torn =
        localRefine(torn_local, G.directed_edges, k);
    Metrics m_refined_torn =
        localMetrics(refined_torn, G.directed_edges, k);

    // Phase 1: Reverse Cuthill-McKee.
    std::vector<int> rcm = rcmOrder(G);
    Metrics m_rcm = localMetrics(rcm, G.directed_edges, k);

    // Phase 2b: Refine RCM ordering if it did not increase FBM.
    // Then pick whichever refined result is better (lower FBM first,
    // then lower TFBD as tiebreaker).  This guarantees the output is
    // never worse than the refined torn ordering.
    std::vector<int> refined;
    if (m_rcm.fbm <= m_torn.fbm) {
      std::vector<int> refined_rcm =
          localRefine(rcm, G.directed_edges, k);
      Metrics m_refined_rcm =
          localMetrics(refined_rcm, G.directed_edges, k);

      if (m_refined_rcm.fbm < m_refined_torn.fbm ||
          (m_refined_rcm.fbm == m_refined_torn.fbm &&
           m_refined_rcm.tfbd < m_refined_torn.tfbd)) {
        refined = std::move(refined_rcm);
      } else {
        refined = std::move(refined_torn);
      }
    } else {
      refined = std::move(refined_torn);
    }

    std::vector<int> global_order;
    global_order.reserve(k);
    for (int li : refined)
      global_order.push_back(scc[li]);

    assert(static_cast<int>(global_order.size()) == k);

    result.push_back(std::move(global_order));
  }

  return result;
}

} // namespace SparseLib::Banding
