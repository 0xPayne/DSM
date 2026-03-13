#include "../include/scc.hpp"

namespace SparseLib {
namespace SCC {

std::vector<std::vector<int>> tarjanSCC(const Sparse::CSCMatrix& matrix) {
    const int n = matrix.cols;
    std::vector<int> index(n, -1);
    std::vector<int> lowLink(n, 0);
    std::vector<bool> onstack(n, false);
    std::vector<int> stack;
    std::vector<std::vector<int>> sccs;

    int currentIndex = 0;

    std::function<void(int)> dfs = [&](int v) {
        // Set the depth index for v to the smallest unused index
        index[v] = lowLink[v] = currentIndex++;
        stack.push_back(v);
        onstack[v] = true;

        // Consider successors w of v
        for (int p = matrix.col_ptrs[v]; p < matrix.col_ptrs[v + 1]; p++) {
            int w = matrix.row_indices[p]; // edge v -> w
            
            if (index[w] == -1) {
                // if successor w was not yet visited, do dfs on it
                dfs(w);
                if (lowLink[w] < lowLink[v]) lowLink[v] = lowLink[w];
            } else if (onstack[w]) {
                // Successor w is in stack S and hence in the current SCC
                // If w is not on stack, then (v, w) is an edge pointing to an SCC already found and must be ignored
                    if (index[w] < lowLink[v]) lowLink[v] = index[w];
            }
        }

        // If v is a root node, pop the stack and generate an SCC
        if (lowLink[v] == index[v]) {
            std::vector<int> component;
            while (true) {
                int w = stack.back();
                stack.pop_back();
                onstack[w] = false;
                component.push_back(w);
                if (v == w) break;
            }
            // appends comp to sccs by moving its contents into the new element instead of copying them
            sccs.push_back(std::move(component));
        }
    };

    for (int v = 0; v < n; v++) {
        if (index[v] == -1) dfs(v);
    }
    return sccs;
}
}
}