#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

#include "plant_propagation.h"

using namespace std;

// http://www.sfu.ca/~ssurjano/schwef.html
double schwefel(const Member& member) {
	double d = member.size();

	double x = member[0];
	double y = member[1];
	return 418.9829 * d - ((x * sin(sqrt(abs(x)))) + (y * sin(sqrt(abs(y)))));
};

// from paper
double minMaxFitness(double v, double min, double max)
{
	return (max - v) / (max - min);
}

int main()
{
	string funName = "schwefel";

	vector<CoordBound> coordinateBounds{
		CoordBound(-500.0, 500.0),
		CoordBound(-500.0, 500.0)
	};

	// Generate a population P = {p_i, i = 1, ...., m}
	const size_t m = 50;

	// iteration g and max number of iterations
	const size_t gMax = 100;

	// max numbers of runners per plant
	const size_t nMax = 5;

	const size_t nRuns = 10;

	vector<MemberWithValue> results;
	results.reserve(nRuns);

	auto knownOptimumMember = Member{ 420.9687, 420.9687 };
	double knownOptimum = schwefel(knownOptimumMember);

	cout << "runPlantPropagation(" << funName << ") with parameters: " << endl;
	cout << "\t" << "m:\t\t " << m << endl;
	cout << "\t" << "g_max:\t\t " << gMax << endl;
	cout << "\t" << "n_max:\t\t " << nMax << endl;
	cout << "\t" << "dims:\t\t " << coordinateBounds.size() << endl;
	cout << "\t" << "bounds:\t\t " << printBounds(coordinateBounds) << endl;
	cout << "\t" << "knownOptimum:\t " << printMember(knownOptimumMember) << endl;
	cout << "\t" << "globalMinimum:\t " << knownOptimum << endl;

	for (size_t i = 0; i < nRuns; ++i) {

		auto finalPopulation = runPlantPropagation1(
			coordinateBounds,
			m,
			gMax,
			nMax,
			minMaxFitness,
			schwefel
		);

		auto bestSolution = finalPopulation[0];

		double objectiveValue = schwefel(bestSolution);

		cout << "run: " << i << ", objectiveValue=" << objectiveValue << endl;

		results.push_back(MemberWithValue(bestSolution, objectiveValue));
	}

	auto bestSolution = min_element(results.begin(), results.end(), compareMemberWithValueLower);

	cout << "--------" << endl;
	cout << "bestSolution: " << printMember(get<0>(*bestSolution)) << endl;
	cout << "bestObjective: " << get<1>(*bestSolution) << endl;

	return 0;
}