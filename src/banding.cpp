#include "../include/banding.hpp"
#include <algorithm>
#include <cassert>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace SparseLib::Banding {

// ============================================================================
// References
// ============================================================================
//
// [1] Cuthill, E. & McKee, J., "Reducing the bandwidth of sparse symmetric
//     matrices", Proceedings of the 1969 24th National Conference, ACM, 1969,
//     pp. 157-172.  DOI: 10.1145/800195.805928
//
// [2] George, A., "Computer implementation of the finite element method",
//     PhD thesis, Stanford University, 1971.
//     (Introduced the Reverse Cuthill-McKee (RCM) variant.)
//
// [3] George, A. & Liu, J. W. H., "An implementation of a pseudoperipheral
//     node finder", ACM Transactions on Mathematical Software, 5(3), 1979,
//     pp. 284-295.  DOI: 10.1145/355841.355845
//     (Pseudo-peripheral node heuristic used for RCM starting vertex.)
//
// [4] Gibbs, N. E., Poole, W. G. & Stockmeyer, P. K., "An algorithm for
//     reducing the bandwidth and profile of a sparse matrix", SIAM Journal
//     on Numerical Analysis, 13(2), 1976, pp. 236-250.
//     DOI: 10.1137/0713023
//     (GPS algorithm; establishes the pseudo-peripheral starting strategy
//      and proves RCM reduces profile compared to CM.)

// ============================================================================
// Internal helpers (anonymous namespace -- internal linkage)
// ============================================================================
namespace {

// ---------------------------------------------------------------------------
// Local SCC subgraph: undirected adjacency + directed edge list
// ---------------------------------------------------------------------------
struct LocalGraph {
  int k;                                     // number of vertices in the SCC
  std::unordered_map<int, int> g2l;          // global vertex ID -> local index
  std::vector<int> globals;                  // local index -> global vertex ID
  std::vector<std::vector<int>> adj;         // undirected adjacency (sorted)
  std::vector<int> deg;                      // undirected degree

  // Directed edges within the SCC: (source_local, target_local).
  // Needed for FBM / TFBD computation.
  struct Edge { int src, dst; };
  std::vector<Edge> directed_edges;
};

/// Build the local undirected adjacency list and directed edge list for an SCC.
/// Follows the same global-to-local mapping pattern as elsOrder() in tearing.cpp.
///
/// The undirected graph is the symmetric closure of the directed subgraph
/// induced by scc_vertices.  Self-loops and inter-SCC arcs are excluded.
/// Adjacency lists are sorted by degree (ascending) for Cuthill-McKee [1].
LocalGraph buildLocalGraph(const std::vector<int> &scc_vertices,
                           const Sparse::CSCMatrix &matrix) {
  LocalGraph G;
  G.k = static_cast<int>(scc_vertices.size());
  G.globals = scc_vertices;
  G.g2l.reserve(G.k);
  for (int i = 0; i < G.k; ++i)
    G.g2l[scc_vertices[i]] = i;

  // Temporary set-based adjacency to deduplicate undirected edges.
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

  // Convert sets to sorted vectors and compute degrees.
  G.adj.resize(G.k);
  G.deg.resize(G.k);
  for (int i = 0; i < G.k; ++i) {
    G.adj[i].assign(adj_set[i].begin(), adj_set[i].end());
    G.deg[i] = static_cast<int>(G.adj[i].size());
  }

  // Sort each adjacency list by degree ascending (Cuthill-McKee rule [1]).
  for (int i = 0; i < G.k; ++i) {
    std::sort(G.adj[i].begin(), G.adj[i].end(),
              [&](int a, int b) { return G.deg[a] < G.deg[b]; });
  }

  return G;
}

// ---------------------------------------------------------------------------
// Pseudo-peripheral vertex finder  [3]
// ---------------------------------------------------------------------------
// Uses iterated BFS to locate a vertex of high eccentricity, which serves
// as a good starting node for the Cuthill-McKee algorithm [1].
//
// George & Liu [3] showed that starting BFS from an arbitrary vertex,
// selecting a vertex from the last level, and repeating until the
// eccentricity (number of BFS levels) stops increasing, yields a
// pseudo-peripheral vertex in practice.  Two to three rounds suffice.

int findPseudoPeripheral(const LocalGraph &G) {
  // BFS returning (last-level vertices, eccentricity).
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

  int current = 0; // start from local vertex 0
  int prev_ecc = -1;

  for (int round = 0; round < 5; ++round) {
    auto [last_level, ecc] = bfs(current);
    if (ecc <= prev_ecc)
      break; // eccentricity did not increase -- converged
    prev_ecc = ecc;
    // Pick the vertex from the last level with minimum degree (heuristic [3]).
    int best = last_level[0];
    for (int v : last_level)
      if (G.deg[v] < G.deg[best])
        best = v;
    current = best;
  }

  return current;
}

// ---------------------------------------------------------------------------
// Reverse Cuthill-McKee (RCM) ordering  [1, 2, 4]
// ---------------------------------------------------------------------------
// The Cuthill-McKee algorithm [1] performs a BFS from a starting vertex,
// enqueuing unvisited neighbours in order of ascending degree.  George [2]
// observed that reversing this ordering consistently produces lower
// bandwidth and profile.  Gibbs, Poole & Stockmeyer [4] provided further
// theoretical justification for the RCM variant.
//
// Complexity: O(k + m_scc) where k = |SCC| and m_scc = edges within the SCC.

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

    // Enqueue unvisited neighbours in ascending degree order.
    // G.adj[u] is already sorted by degree ascending (see buildLocalGraph).
    for (int v : G.adj[u]) {
      if (!visited[v]) {
        visited[v] = true;
        q.push(v);
      }
    }
  }

  // Handle disconnected components within the SCC subgraph's undirected view.
  // (Rare: an SCC is strongly connected, but the undirected symmetrisation
  //  can still leave isolated components in degenerate cases.)
  for (int i = 0; i < G.k; ++i) {
    if (!visited[i])
      order.push_back(i);
  }

  // Reverse to obtain RCM [2].
  std::reverse(order.begin(), order.end());

  return order;
}

// ---------------------------------------------------------------------------
// Local FBM and TFBD computation for an SCC
// ---------------------------------------------------------------------------
// Given a vertex ordering (local indices in the desired sequence) and the
// directed edge list, compute FBM and TFBD restricted to this SCC.
//
// In the CSC convention of this project, entry (row, col) with row < col is
// a feedback mark.  After permutation, col = pos(source), row = pos(target)
// for a directed edge source -> target.  A feedback mark occurs when
// pos(target) < pos(source), i.e. the target appears before the source.

struct Metrics {
  int fbm;
  long long tfbd;
};

Metrics localMetrics(const std::vector<int> &order,
                     const std::vector<LocalGraph::Edge> &edges, int k) {
  // pos[local_vertex] = position in the ordering
  std::vector<int> pos(k);
  for (int i = 0; i < k; ++i)
    pos[order[i]] = i;

  Metrics m{0, 0};
  for (const auto &e : edges) {
    int p_src = pos[e.src];
    int p_dst = pos[e.dst];
    if (p_dst < p_src) { // feedback mark: target before source
      ++m.fbm;
      m.tfbd += p_src - p_dst;
    }
  }
  return m;
}

// ---------------------------------------------------------------------------
// Adjacent-swap hill climbing on TFBD
// ---------------------------------------------------------------------------
// After the initial RCM (or tearing) ordering, iteratively scan all adjacent
// pairs and swap them if doing so strictly reduces TFBD without increasing
// FBM.  Each swap delta is computed in O(deg(a) + deg(b)) by examining only
// the edges incident to the two vertices involved, rather than recomputing
// the full metric.
//
// Convergence: each pass is O(k * d_avg).  Passes repeat until no improving
// swap exists, capped at max_passes to bound runtime on large SCCs.
//
// This local-search strategy is a standard neighbourhood-based improvement
// heuristic for combinatorial optimisation on permutations; see e.g.:
//
// [5] Kernighan, B. W. & Lin, S., "An efficient heuristic procedure for
//     partitioning graphs", Bell System Technical Journal, 49(2), 1970,
//     pp. 291-307.  DOI: 10.1002/j.1538-7305.1970.tb01770.x
//     (Foundational pairwise-swap improvement for graph partitioning.)

std::vector<int>
localRefine(const std::vector<int> &initial_order,
            const std::vector<LocalGraph::Edge> &edges, int k,
            int max_passes = 20) {
  std::vector<int> order = initial_order;

  // pos[local_vertex] = position in the ordering
  std::vector<int> pos(k);
  for (int i = 0; i < k; ++i)
    pos[order[i]] = i;

  // Build incidence lists: for each local vertex, the indices into `edges`
  // where that vertex appears as src or dst.
  std::vector<std::vector<int>> incident(k);
  for (int ei = 0; ei < static_cast<int>(edges.size()); ++ei) {
    incident[edges[ei].src].push_back(ei);
    incident[edges[ei].dst].push_back(ei);
  }

  for (int pass = 0; pass < max_passes; ++pass) {
    bool improved = false;

    for (int i = 0; i < k - 1; ++i) {
      int a = order[i];     // vertex at position i
      int b = order[i + 1]; // vertex at position i+1

      // Compute delta_tfbd and delta_fbm if we swap a and b.
      // After swap: a moves to i+1, b moves to i.
      int delta_fbm = 0;
      long long delta_tfbd = 0;

      // Process edges incident to a or b.  Use a merged scan over both
      // incidence lists.  For each edge, compute contribution change.
      //
      // For an edge (src, dst): its TFBD contribution is
      //   max(0, pos[src] - pos[dst])   if it's a feedback mark
      //
      // We need to check: does the edge's FBM status or distance change
      // when a moves from i to i+1 and b moves from i+1 to i?

      auto process_edges = [&](int vertex) {
        for (int ei : incident[vertex]) {
          int s = edges[ei].src;
          int d = edges[ei].dst;

          // Skip the edge between a and b -- handle it separately below.
          if ((s == a && d == b) || (s == b && d == a))
            continue;

          int ps_old = pos[s];
          int pd_old = pos[d];

          // New positions after swapping a (at i) and b (at i+1).
          int ps_new = (s == a) ? i + 1 : (s == b) ? i : ps_old;
          int pd_new = (d == a) ? i + 1 : (d == b) ? i : pd_old;

          // Old FBM/TFBD contribution of this edge.
          bool fbm_old = (pd_old < ps_old);
          long long dist_old = fbm_old ? (ps_old - pd_old) : 0;

          // New contribution.
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

      // Handle the direct edge(s) between a and b.
      for (int ei : incident[a]) {
        int s = edges[ei].src;
        int d = edges[ei].dst;
        if (!((s == a && d == b) || (s == b && d == a)))
          continue;

        int ps_old = pos[s], pd_old = pos[d];
        // After swap: a goes to i+1, b goes to i.
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

      // Accept swap only if TFBD strictly improves and FBM does not worsen.
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

// ============================================================================
// Public API
// ============================================================================

/// Apply bandwidth-reduction banding to each SCC.
///
/// Two-phase refinement per non-trivial SCC:
///   Phase 1 -- Reverse Cuthill-McKee (RCM) [1,2,4] reordering of the
///              undirected symmetrisation of the SCC subgraph.
///   Phase 2 -- Adjacent-swap hill climbing [5] on TFBD, constrained to
///              never increase FBM.
///
/// If RCM increases FBM relative to the input (tearing) ordering, Phase 1
/// is skipped and Phase 2 refines the tearing ordering directly.  This
/// guarantees FBM is never worse than the input.
std::vector<std::vector<int>>
bandAllSCCs(const std::vector<std::vector<int>> &torn_sccs,
            const Sparse::CSCMatrix &matrix) {
  std::vector<std::vector<int>> result;
  result.reserve(torn_sccs.size());

  for (const auto &scc : torn_sccs) {
    int k = static_cast<int>(scc.size());

    // Trivial SCCs: nothing to reorder.
    if (k <= 2) {
      result.push_back(scc);
      continue;
    }

    // Build local subgraph (undirected adjacency + directed edge list).
    LocalGraph G = buildLocalGraph(scc, matrix);

    // Identity ordering (local index i is at position i) represents the
    // incoming tearing ordering, since scc[i] = global vertex at position i.
    std::vector<int> torn_local(k);
    for (int i = 0; i < k; ++i)
      torn_local[i] = i;

    Metrics m_torn = localMetrics(torn_local, G.directed_edges, k);

    // Phase 2a: Always refine the torn ordering as a baseline.
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

    // Map local indices back to global vertex IDs.
    std::vector<int> global_order;
    global_order.reserve(k);
    for (int li : refined)
      global_order.push_back(scc[li]);

    // Sanity: output must be a permutation of the input SCC.
    assert(static_cast<int>(global_order.size()) == k);

    result.push_back(std::move(global_order));
  }

  return result;
}

} // namespace SparseLib::Banding
