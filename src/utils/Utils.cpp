#include "Utils.h"

int choose(const Vec2D<int>& subtours) {
    auto m = INFINITE;
    auto i = 0;
    auto chosen = 0;
    for(auto& st : subtours) {
        if(st.size() < m) {
            m = st.size();
            chosen = i;
        }
        ++i;
    }

    return chosen;
}

bool contains(Vec2D<int> v2, std::vector<int> v1) {
    for(auto vec : v2) {
        if(v1.size() == vec.size()) {
            auto c = true;
            for(uint i = 0; i < v1.size(); ++i) {
                if(v1[i] != vec[i]) {
                    c = false;
                    break;
                }
            }
            if(c) {
                return true;
            }
        }
    }

    return false;
}

Vec2D<int> findSubtours(const std::shared_ptr<Data>& pData,
                        const hungarian_problem_t& rHungarian) {
    // std::cout << "Find subtours" << std::endl;
    std::vector<std::list<int>> adjacent(pData->getDimension());
    for(int i = 0; i < pData->getDimension(); ++i) {
        for(int j = 0; j < pData->getDimension(); ++j) {
            if(rHungarian.assignment[i][j] == 1) {
                adjacent[i].push_back(j);
            }
        }
    }

    std::vector<Status> visited(pData->getDimension(), not_visited);
    Vec2D<int> subtours;
    for(int v = 0; v < pData->getDimension(); ++v) {
        if(visited[v] == not_visited) {
            std::list<int> stack;
            stack.push_back(v);
            visited[v] = in_stack;
            processDFSTree(subtours, adjacent, stack, visited);
        }
    }

    return subtours;
}

void processDFSTree(Vec2D<int>& rSubtours,
                    const std::vector<std::list<int>>& crAdjacent, 
                    std::list<int>& rStack, 
                    std::vector<Status>& rVisited) {
    for(auto v : crAdjacent[rStack.back()]) {
        if(rVisited[v] == in_stack) {
            recoverSubtour(rSubtours, rStack, v);
        } else if(rVisited[v] == not_visited) {
            rStack.push_back(v);
            rVisited[v] = in_stack;
            processDFSTree(rSubtours, crAdjacent, rStack, rVisited);
        }
    }
    rVisited[rStack.back()] = is_visited;
    rStack.pop_back();
}

void recoverSubtour(Vec2D<int>& rSubtours, 
                    std::list<int>& rStack, int v) {
    // std::cout << "Recover subtour" << std::endl;
    // Theoretically, the following is not a behavior of a stack,
    // but it suffices for the sake of implementation.
    auto it = rStack.end();
    --it;
    while(true) {
        if(*it == v) {
            break;
        }
        --it;
    }
    std::vector<int> subtour(it, rStack.end());
    subtour.push_back(v);
    // Check if the new subtour is already in the subtours.
    if(!contains(rSubtours, subtour)) {
        rSubtours.push_back(subtour);
    }
}

double solveGreedyTSP(const std::shared_ptr<Data>& pData) {
    auto c = 0.0;
    auto i = 0;
    for(; i < pData->getDimension() - 1; ++i) {
        c += pData->getDistance(i, i+1);
    }
    return c + pData->getDistance(i, 0);
}

std::pair<bool, int> isFeasible(const std::vector<int>& rDegree) {
    bool isFeas = true;
    int sel = 0;
    int maxDegree = 0;
    for(uint i = 0; i < rDegree.size(); ++i) {
        if(rDegree[i] != 2) {
            isFeas = false;
        }
        if(rDegree[i] > maxDegree) {
            maxDegree = rDegree[i];
            sel = i;
        }
    }
    return std::make_pair(isFeas, sel);
}