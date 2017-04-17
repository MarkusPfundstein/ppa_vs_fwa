#include <vector>
#include <iostream>
#include <algorithm>

#include "fwa.h"


// https://en.wikipedia.org/wiki/Rosenbrock_function
double rosenbrock(const Member& member) {
	double a = 1;
	double b = 100.0;

	double x = member[0];
	double y = member[1];
	return pow((a - x), 2) + b * pow((y - pow(x, 2)), 2);
};

// from paper
double minMaxFitness(double v, double min, double max)
{
	return (max - v) / (max - min);
}

int main()
{
	vector<CoordBound> coordinateBounds{
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0)/*,
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0),
		CoordBound(-5.0, 10.0)*/
	};

	// Generate a population P = {p_i, i = 1, ...., m}
	const size_t n = 30;

	// iteration g and max number of iterations
	const size_t gMax = 30;

	const size_t nRuns = 10;

	vector<MemberWithValue> results;
	results.reserve(nRuns);

	auto knownOptimumMember = Member{ 1.0, 1.0 };
	double knownOptimum = rosenbrock(knownOptimumMember);

	string funName = "rosenbrock";

	cout << "fwa(" << funName << ") with parameters: " << endl;
	cout << "\t" << "n:\t\t " << n << endl;
	cout << "\t" << "g_max:\t\t " << gMax << endl;
	cout << "\t" << "dims:\t\t " << coordinateBounds.size() << endl;
	cout << "\t" << "bounds:\t\t " << printBounds(coordinateBounds) << endl;
	cout << "\t" << "knownOptimum:\t " << printMember(knownOptimumMember) << endl;
	cout << "\t" << "globalMinimum:\t " << knownOptimum << endl;

	for (size_t i = 0; i < nRuns; ++i) {

		auto finalPopulation = fwa(
			coordinateBounds,
			n,
			gMax,
			minMaxFitness,
			rosenbrock
		);

		auto bestSolution = finalPopulation[0];

		double objectiveValue = rosenbrock(bestSolution);

		cout << "run: " << i << ", objectiveValue=" << objectiveValue << endl;

		results.push_back(MemberWithValue(bestSolution, objectiveValue));
	}

	auto bestSolution = min_element(results.begin(), results.end(), compareMemberWithValueLower);

	cout << "--------" << endl;
	cout << "bestSolution: " << printMember(get<0>(*bestSolution)) << endl;
	cout << "bestObjective: " << get<1>(*bestSolution) << endl;

	return 0;
}