#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "BNB.h"
#include "data.h"
#include "hungarian.h"
#include "Kruskal.h"
#include "Utils.h"

using namespace std;

int main(int argc, char** argv) {
	auto pData = std::make_shared<Data>(argc, argv[1]);
	pData->readData();
	#if DEBUGGING_LEVEL == 2
		Vec2D<double> dist = {{ 0, 30, 26, 50, 40},
							  {30,  0, 24, 40, 50},
							  {26, 24,  0, 24, 26},
							  {50, 40, 24,  0, 30},
							  {40, 50, 26, 30,  0}};
		pData->dimension = 5;
	#endif

	vector<vector<double>> cost(pData->getDimension(), vector<double>(pData->getDimension()));
	for(int i = 0; i < pData->getDimension(); i++) {
		for(int j = 0; j < pData->getDimension(); j++) {
			#if DEBUGGING_LEVEL == 2
				cost[i][j] = dist[i][j];
			#else
				cost[i][j] = pData->getDistance(i,j);
			#endif
		}
	}

	// auto kruskal = std::make_shared<Kruskal>(cost);

	// kruskal->oneTree(pData->getDimension(), cost);
	// kruskal->printEdges();

	auto bnb = std::make_shared<BNB>(pData, cost);
	std::cout << "Branching strategy: " << argv[2] << std::endl;
	double ub = atof(argv[3]);
	if(string(argv[2]) == "dfs") {
		bnb->run(true, ub);
	} else if(string(argv[2]) == "bfs") {
		bnb->run(false, ub);
	} else if(string(argv[2]) == "lb") {
		bnb->runLB(ub);
	}

	return 0;
}

// int main(int argc, char** argv) {
// 	auto files = std::vector<string>({"bayg29", "bays29", "burma14", "fri26", 
// 									  "gr17", "gr21", "gr24", "ulysses16", 
// 									  "ulysses22"});
// 	for(auto f : files) {
// 		auto inst = "input/" + f + ".tsp";
// 		std::cout << "Read instance: " << inst << std::endl;
// 		// auto pData = std::make_shared<Data>(argc, argv[1]);
// 		auto pData = std::make_shared<Data>(argc, inst.c_str());
// 		pData->readData();

// 		double **cost = new double*[pData->getDimension()];
// 		for(int i = 0; i < pData->getDimension(); i++) {
// 			cost[i] = new double[pData->getDimension()];
// 			for(int j = 0; j < pData->getDimension(); j++) {
// 				cost[i][j] = pData->getDistance(i,j);
// 			}
// 		}

// 		auto bnb = std::make_shared<BNB>(pData, cost);
// 		std::cout << "Branching strategy: " << argv[2] << std::endl;
// 		if(string(argv[2]) == "dfs") {
// 			bnb->run(true);
// 		} else if(string(argv[2]) == "bfs") {
// 			bnb->run(false);
// 		} else if(string(argv[2]) == "lb") {
// 			bnb->runLB();
// 		}

// 		// hungarian_problem_t p;
// 		// int mode = HUNGARIAN_MODE_MINIMIZE_COST;
// 		// double obj_value = hungarian_solve(&p);
// 		// cout << "Obj. value: " << obj_value << endl;

// 		// cout << "Assignment" << endl;
// 		// hungarian_print_assignment(&p);

// 		// hungarian_free(&p);
// 		for(int i = 0; i < pData->getDimension(); i++) delete [] cost[i];
// 		delete [] cost;
// 		// delete pData.get();
// 	}

// 	return 0;
// }
