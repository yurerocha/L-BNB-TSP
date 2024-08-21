#ifndef BRANCH_AND_BOUND_H
#define BRANCH_AND_BOUND_H

#include <iomanip>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>
#include "data.h"
#include "hungarian.h"
#include "Kruskal.h"
#include "Lagrange.h"
#include "Timer.h"
#include "Utils.h"

class Lagrange;

class BNB {
public:
    BNB(const std::shared_ptr<Data>& pData, const Vec2D<double>& pCost);

    /**
     * @brief Run branch-and-bound algorithm with DFS or BFS branching strategy.
    */
    void run(bool isDFS, double upperBound);
    /**
     * @brief Run branch-and-bound algorithm with lower bound branching 
     * strategy.
    */
    void runLB(double upperBound);

    auto branchingStrategy(const std::list<Node>& tree, bool isDFS) const {
        if(isDFS) {
            // tree.end() points to last pos + 1
            return std::make_pair(tree.back(), --tree.end());
        } else {
            return std::make_pair(tree.front(), tree.begin());
        }
    };

    /**
     * @brief populate the node with the lagrangian solution.
    */
    void computeSolutionLagrange(Node& rNode, double upperBound);

    /**
     * @brief Compute the forbidden arcs, i.e., those that are incident to node
     * iMaxDegree.
    */
    std::vector<std::pair<int, int>> 
             computeArcsToForbid(const std::vector<std::pair<int, int>>& rEdges,
                                 int iMaxDegree) const;
    
    /**
     * @brief 
    */
   void addBranches(std::list<Node>& rTree, const Node& rNode) const;

    void log(std::ostream& rOs, double lb, double up, double time) const;

    // ~BNB();
private:
    std::shared_ptr<Data> mpData;
    Vec2D<double> mpCost;
    Vec2D<double> mpCostCopy;
};

#endif