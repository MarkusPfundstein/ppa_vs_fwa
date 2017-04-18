#include <vector>
#include <iostream>
#include <algorithm>

#include "framework.h"

int main()
{
	vector<CoordBound> coordinateBounds{
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0)
	};
	Member knownOptimum = { 1.0, 1.0 };

	FWA params;
	params.coordinateBounds = coordinateBounds;
	params.initialSize = 30;
	params.maxGenerations = 30;
	params.objectiveFunction = rosenbrock2d;
	params.objectiveFunctionName = "rosenbrock";
	params.knownOptimum = &knownOptimum;

	return runExperiments(10, &params);
}