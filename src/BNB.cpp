#include "BNB.h"

BNB::BNB(const std::shared_ptr<Data>& pData, const Vec2D<double>& pCost) 
    : mpData(pData), mpCost(pCost) {
    // hungarian_init(&mHungarian, mpCost, mpData->getDimension(), 
    //                mpData->getDimension(), mode);
    mpCostCopy = Vec2D<double>(mpData->getDimension());
    for(int i = 0; i < mpData->getDimension(); ++i){
        mpCostCopy[i] = std::vector<double>(mpData->getDimension());
        for(int j = 0; j < mpData->getDimension(); ++j){
            mpCostCopy[i][j] = mpCost[i][j];
        }
    }
}

void BNB::run(bool isDFS, double upperBound) {
    Timer timer;
    timer.start();

    std::cout << "Start execution of BNB algorithm (lb ub gap time)" 
              << std::endl;
    Node root; // no raiz
    root.lambda = std::vector<double>(mpData->getDimension(), 0.0);
    double lowerBound = 0.0;

    /* criacao da arvore */
    auto tree = std::list<Node>();
    tree.push_back(root);

    // gerar solução inicial
    while(!tree.empty()) {
        // escolher um dos nos da arvore
        auto [node, itNode] = branchingStrategy(tree, isDFS);

        computeSolutionLagrange(node, upperBound);
        if(iseq(lowerBound, 0.0)) {
            lowerBound = node.lower_bound;
        }
        if(isg(node.lower_bound, upperBound)) {
            tree.erase(itNode);
            continue;
        }
        if(node.feasible) {
            if(isl(node.lower_bound, upperBound)) {
                upperBound = node.lower_bound;
                log(std::cout, lowerBound, upperBound, timer.count());
            }
        } else {
            /* Adicionando os filhos */
            // iterar por todos os arcos do subtour escolhido
            addBranches(tree, node);
        }
        tree.erase(itNode);
    }
}

void BNB::runLB(double upperBound) {
    Timer timer;
    timer.start();

    std::cout << "Start execution of BNB algorithm (lb ub gap time)" 
              << std::endl;
    Node root; // no raiz
    root.lambda = std::vector<double>(mpData->getDimension(), 0.0);
    double lowerBound = 0.0;

    /* criacao da arvore */
    auto comp = [](const Node& n1, const Node& n2) {
                    return isg(n1.lower_bound, n2.lower_bound);
                };
    std::priority_queue<Node, std::vector<Node>, decltype(comp)> tree(comp);

    tree.push(root);
    // gerar solução inicial
    while(!tree.empty()) {
        // escolher um dos nos da arvore
        auto node = tree.top();

        computeSolutionLagrange(node, upperBound);
        if(iseq(lowerBound, 0.0)) {
            lowerBound = node.lower_bound;
        }
        if(isg(node.lower_bound, upperBound)) {
            tree.pop();
            continue;
        }
        if(node.feasible) {
            // upper_bound = std::min(upper_bound, node.lower_bound);
            if(isl(node.lower_bound, upperBound)) {
                upperBound = node.lower_bound;
                log(std::cout, lowerBound, upperBound, timer.count());
            }
        } else {
            /* Adicionando os filhos */
            // iterar por todos os arcos do subtour escolhido
            for(uint i = 0; i < node.to_forbid.size(); ++i) {
                Node n;
                n.lambda = node.lambda;
                n.forbidden_arcs = node.forbidden_arcs;
                n.forbidden_arcs.push_back(node.to_forbid[i]);
                tree.push(n);
            }
        }
        tree.pop();
    }
}

void BNB::computeSolutionLagrange(Node& rNode, double upperBound) {
    #if DEBUGGING_LEVEL == 2
        std::cout << "Forbidden arcs: " << rNode.forbidden_arcs.size() 
                  << std::endl;
    #endif
    for(const auto& forb_arc : rNode.forbidden_arcs) {
        #if DEBUGGING_LEVEL == 2
            std::cout << "(" << forb_arc.first << ", " << forb_arc.second 
                      << ") ";
        #endif
        mpCost[forb_arc.first][forb_arc.second] = INFINITE;
    }
    #if DEBUGGING_LEVEL == 2
        std::cout << std::endl;
    #endif

    Lagrange lagrange(mpData, mpCost, upperBound);
    #if DEBUGGING_LEVEL == 2
        std::cout << "Lambda: ";
        for(uint i = 0; i < rNode.lambda.size(); ++i) {
            std::cout << rNode.lambda[i] << " ";
        }
        std::cout << std::endl;
    #endif
    auto [oneTree, cost, lambda] = lagrange.subgradient(rNode.lambda);
    auto [isFeas, iMaxDegree] = isFeasible(lagrange.computeDegrees(oneTree));
    rNode.one_tree = oneTree;
    rNode.lower_bound = cost;
    rNode.lambda = lambda;
    rNode.feasible = isFeas;
    rNode.to_forbid = computeArcsToForbid(oneTree, iMaxDegree);
    #if DEBUGGING_LEVEL == 2
        if(isFeas) {
            std::cout << "Feasible solution: " << cost << std::endl;
            for(uint i = 0; i < oneTree.size(); ++i){
                cout << oneTree[i].first << " " << oneTree[i].second << endl;
            }
            // getchar();
        } else {
            std::cout << "Infeasible solution: " << cost << std::endl;
            std::cout << "iMaxDegree: " << iMaxDegree << std::endl;
            std::cout << "oneTree.size: " << oneTree.size() << std::endl;
            for(uint i = 0; i < oneTree.size(); ++i){
                cout << oneTree[i].first << " " << oneTree[i].second << endl;
            }
            std::cout << "toForbid.size: " << rNode.to_forbid.size() 
                      << std::endl;
            for(uint i = 0; i < rNode.to_forbid.size(); ++i){
                cout << rNode.to_forbid[i].first << " " 
                     << rNode.to_forbid[i].second << endl;
            }
        }
    #endif

    for(const auto& forb_arc : rNode.forbidden_arcs) {
        mpCost[forb_arc.first][forb_arc.second] 
                                  = mpCostCopy[forb_arc.first][forb_arc.second];
    }
}

std::vector<std::pair<int, int>> 
        BNB::computeArcsToForbid(const std::vector<std::pair<int, int>>& rEdges,
                                 int iMaxDegree) const {
    std::vector<std::pair<int, int>> toForbid;
    for(auto& e : rEdges) {
        if(e.first == iMaxDegree || e.second == iMaxDegree) {
            toForbid.push_back(e);
        }
    }
    return toForbid;
}

void BNB::addBranches(std::list<Node>& rTree, const Node& rNode) const {
    for(uint i = 0; i < rNode.to_forbid.size(); ++i) {
        Node n;
        n.lambda = rNode.lambda;
        n.forbidden_arcs = rNode.forbidden_arcs;
        n.forbidden_arcs.push_back(rNode.to_forbid[i]);
        rTree.push_back(n);
    }
}

void BNB::log(std::ostream& rOs, double lb, double up, double time) const {
    rOs << std::fixed << std::setprecision(1) << lb << "\t"
        << up << "\t"
        << gap(lb, up) << "\t"
        << time << std::endl;
}

// BNB::~BNB() {
//   hungarian_free(&hungarian);
// }
