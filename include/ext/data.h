#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

#define NB_PARAMS 4
#define INFINITE 99999999.0

using namespace std;

class Data{
public:
	Data( int, const char * );
	~Data();

	void readData();
	void printMatrixDist();
	inline int getDimension(){ return dimension; };
	inline double getDistance(int i, int j){return distMatrix[i][j]; };
	inline double **getMatrixCost(){return distMatrix; }
	inline double getXCoord(int i){return xCoord[i];}
	inline double getYCoord(int i){return yCoord[i];}
	inline bool getExplicitCoord(){return explicitCoord; };

	string getInstanceName();                    //Get instance's name

	int dimension;
private:
	string instaceName;
	int nbOfPar;


	double *xCoord, *yCoord;
	double **distMatrix;

	//Computing Distances
	static double CalcDistEuc ( double *, double *, int , int );
	static double CalcDistAtt ( double *, double *, int , int );
	static double CalcDistGeo ( double *, double *, int , int );
	static void CalcLatLong ( double *, double *, int , double *, double* );

	bool explicitCoord;
};

#endif
