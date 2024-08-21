#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include "data.h"
#include "Utils.h"

#define FIRST_NODE 1

using namespace std;

typedef pair<int, int> ii;
typedef vector <vector<double> >vvi;
typedef vector<ii> vii;

class Kruskal{
public:
	Kruskal(const vvi& dist);

	double MST(int nodes);
	double oneTree(int nodes, const vvi& dist);
	vii getEdges();
	void printEdges();

private:
	priority_queue <pair<double,ii> > graph;
	vector <int> pset;
	vii edges;

	void initDisjoint(int n);
	int findSet(int i);
	void unionSet(int i, int j);
	bool isSameSet(int i, int j);
};

#endif