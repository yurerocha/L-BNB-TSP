#ifndef LAGRANGE_H
#define LAGRANGE_H

#include <tuple>
#include <vector>
#include "data.h"
#include "Kruskal.h"
#include "Utils.h"

class Lagrange {
public:
    Lagrange(const std::shared_ptr<Data>& pData, 
             const Vec2D<double>& pCost,
             double upperBound);

    /**
     * @brief Compute the number of incoming and oucoming arcs of each node.
    */
    std::vector<int> computeDegrees(
                          const std::vector<std::pair<int, int>>& rEdges) const;

    /**
     * @brief Update the lagrangean subgradient vector.
    */
    std::tuple<std::vector<std::pair<int, int>>,
               double,
               std::vector<double>> subgradient(std::vector<double>& r_lambda);

private:
    /**
     * @brief Update the costs of the edges.
    */
    void updateCosts(Vec2D<double>& rNewCost, 
                     const std::vector<double>& r_lambda);

    /**
     * @brief Compute the cost of the lagrangean relaxation.
    */
    double computeCost(double cost, const std::vector<double>& r_lambda) const;

    /**
     * @brief Update the lagrangean penalizers.
    */
    void updatePenalizers(double eps, double omega, 
                          std::vector<double>& r_lambda,
                          const std::vector<std::pair<int, int>>& rEdges);

    std::shared_ptr<Data> mpData;
    Vec2D<double> mpCost;
    double mUpperBound;
    double m_epsilonMin;
    uint m_kMax;
};

#endif