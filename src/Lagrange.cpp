#include "Lagrange.h"

Lagrange::Lagrange(const std::shared_ptr<Data>& pData, 
                   const Vec2D<double>& pCost,
                   double upperBound) 
    : mpData(pData), mpCost(pCost), mUpperBound(upperBound),
      m_epsilonMin(0.0005), m_kMax(30) {

}

std::vector<int> Lagrange::computeDegrees(
                         const std::vector<std::pair<int, int>>& rEdges) const {
    std::vector<int> degree(mpData->getDimension(), 0);
    for(const auto& e : rEdges) {
        degree[e.first]++;
        degree[e.second]++;
    }
    return degree;
}

std::tuple<std::vector<std::pair<int, int>>,
           double, std::vector<double>> 
                          Lagrange::subgradient(std::vector<double>& r_lambda) {
    double eps = 1.0;
    uint k = 0;

    auto best_lambda = r_lambda;
    auto kruskal = std::make_unique<Kruskal>(mpCost);
    double bestCost = kruskal->oneTree(mpData->getDimension(), mpCost);
    auto bestOneTree = kruskal->getEdges();
    double bestLowerBound = computeCost(bestCost, r_lambda);

    updatePenalizers(eps, bestLowerBound, r_lambda, kruskal->getEdges());
    auto degree = computeDegrees(kruskal->getEdges());

    #if DEBUGGING_LEVEL == 1
        std::cout << "Subgradient method" << std::endl;
    #endif
    auto newCost = mpCost;
    while(true) {
        updateCosts(newCost, r_lambda);
        // optimal solution of PR_lambda
        kruskal = std::make_unique<Kruskal>(newCost);
	    double cost = kruskal->oneTree(mpData->getDimension(), newCost);
        // update omega
        double lowerBound = computeCost(cost, r_lambda);
        // lambda
        updatePenalizers(eps, lowerBound, r_lambda, kruskal->getEdges());
        degree = computeDegrees(kruskal->getEdges());

        if(isg(lowerBound, bestLowerBound)) {
            bestLowerBound = lowerBound;
            bestCost = cost;
            best_lambda = r_lambda;
            bestOneTree = kruskal->getEdges();
            #if DEBUGGING_LEVEL == 1
                std::cout << "it: " << k 
                          << " cost:" << bestLowerBound
                          << " size:" << bestOneTree.size() << std::endl;
            #endif
            k = 0;
        } else {
            k++;
            if(k >= m_kMax) {
                k = 0;
                eps /= 2.0;
            }
        }
        // implement stop criterion
        if(isl(eps, m_epsilonMin) || isFeasible(degree).first) {
            break;
        }
    }

    return std::make_tuple(bestOneTree, bestCost, best_lambda);
}

void Lagrange::updateCosts(Vec2D<double>& rNewCost,
                           const std::vector<double>& r_lambda) {
    for(int i = 0; i < mpData->getDimension(); ++i) {
        for(int j = 0; j < mpData->getDimension(); ++j) {
            rNewCost[i][j] = mpCost[i][j] - r_lambda[i] - r_lambda[j];
        }
    }
}

double Lagrange::computeCost(double cost, 
                             const std::vector<double>& r_lambda) const {
    double omega = 0.0;
    for(auto l : r_lambda) {
        omega += l;
    }
    return cost + 2.0 * omega;
}

void Lagrange::updatePenalizers(double eps, double omega, 
                               std::vector<double>& r_lambda,
                               const std::vector<std::pair<int, int>>& rEdges) {
    auto degree = computeDegrees(rEdges);
    double sum = 0.0;
    for(auto d : degree) {
        sum += std::pow(2.0 - d, 2);
    }
    
    double mu = eps * (mUpperBound - omega) / sum;
    for(uint i = 0; i < r_lambda.size(); ++i) {
        r_lambda[i] += mu * (2.0 - degree[i]);
    }
}