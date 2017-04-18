#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

#include "framework.h"

using namespace std;

int main()
{
	vector<CoordBound> coordinateBounds{
		CoordBound(-500.0, 500.0),
		CoordBound(-500.0, 500.0)
	};

	Member knownOptimumMember = { 420.9687, 420.9687 };
	
	PPA ppa;
	ppa.coordinateBounds = coordinateBounds;
	ppa.initialSize = 30;
	ppa.maxGenerations = 30;
	ppa.nMax = 10;
	ppa.objectiveFunction = schwefel2d;
	ppa.objectiveFunctionName = "schwefel2d";
	ppa.knownOptimum = &knownOptimumMember;

	return runExperiments(10, &ppa);
}