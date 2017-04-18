#include <vector>
#include <iostream>
#include <algorithm>

#include "ppa.h"

using namespace std;

int main()
{
	vector<CoordBound> coordinateBounds{
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0)
	};

	Member knownOptimum = { 1.0, 1.0 };

	PPA ppa;
	
	ppa.coordinateBounds = coordinateBounds;
	ppa.initialSize = 30;
	ppa.maxGenerations = 30;
	ppa.nMax = 5;
	ppa.objectiveFunctionName = "rosenbrock2d";
	ppa.objectiveFunction = rosenbrock2d;
	ppa.knownOptimum = &knownOptimum;

	return runExperiments(10, &ppa);
}