#ifndef UTILS_H
#define UTILS_H

#include <list>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>
#include "data.h"
#include "hungarian.h"

#define DEBUGGING_LEVEL 0

typedef unsigned int uint;

template <class T>
using Vec2D = std::vector<std::vector<T>>;

const double eps = 0.00001;

struct Node {
    std::vector<std::pair<int, int>> forbidden_arcs; // arcos já proibidos do no
    std::vector<std::pair<int, int>> to_forbid; // arcos que devem ser proibidos nos nós filhos
    std::vector<std::pair<int, int>> one_tree;
    std::vector<double> lambda; // penalizadores lagrangianos
    double lower_bound; // custo total da solucao
    bool feasible; // indica se a solucao eh viavel
};

int choose(const Vec2D<int>& subtours);

/**
 * @brief DFS find cycles/subtours.
*/
enum Status {is_visited, not_visited, in_stack};
Vec2D<int> findSubtours(const std::shared_ptr<Data>& pData,
                        const hungarian_problem_t& rHung);
void processDFSTree(Vec2D<int>& rSubtours, 
                    const std::vector<std::list<int>>& rAdjacent, 
                    std::list<int>& rStack, 
                    std::vector<Status>& rVisited);
void recoverSubtour(Vec2D<int>& rSubtours, 
                    std::list<int>& rStack, int v);

/**
* @brief Any vector of v2 is identical to v1?
*/
bool contains(Vec2D<int> v2, std::vector<int> v1);

inline double gap(double uAst, double sAst) {
    return 100.0 * (sAst - uAst) / uAst;
}

/**
* @brief Initial greedy solution to TSP.
*/
double solveGreedyTSP(const std::shared_ptr<Data>& pData);

/**
 * @brief Return the solution feasibility status and the node with the maximum
 * degree.
*/
std::pair<bool, int> isFeasible(const std::vector<int>& rDegree);

inline bool isl(double a, double b) {
    return a < b - eps;
}

inline bool isg(double a, double b) {
    return a > b + eps;
}

inline bool iseq(double a, double b) {
    return fabs(a - b) < eps;
}

#endif