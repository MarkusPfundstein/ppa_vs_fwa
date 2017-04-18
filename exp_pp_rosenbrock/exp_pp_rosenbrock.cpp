#include <vector>
#include <iostream>
#include <algorithm>

#include "ppa.h"

using namespace std;

int main()
{
	auto coordinateBounds = createUniformCoordinateBounds(2, -5, 10);

	Member knownOptimum = { 1.0, 1.0 };

	PPA ppa;
	ppa.coordinateBounds = coordinateBounds;
	ppa.initialSize = 30;
	ppa.maxGenerations = 30;
	ppa.nMax = 5;
	ppa.objectiveFunctionName = "rosenbrock2d";
	ppa.knownOptimum = &knownOptimum;

	return runExperiments(10, &ppa);
}