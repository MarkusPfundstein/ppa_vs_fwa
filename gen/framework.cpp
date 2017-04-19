#include <iostream>
#include <algorithm>
#include <tuple>
#include "framework.h"
#include "ppa.h"
#include "fwa.h"

map<string, function<Population(Parameters*)>> ALGO_MAP = {
	{ "ppa", runPPA },
	{ "fwa", runFWA }
};

map <string, function<double(const Member&)>> OBJ_MAP = {
	{ "rosenbrock", rosenbrock },
	{ "griewank", griewank },
	{ "schwefel", schwefel },
	{ "easom", easom },
	{ "ackleys_path", ackleys_path }
};

Population run(Parameters *ps)
{
	return ALGO_MAP.at(ps->algorithm)(ps);
}

int runExperiments(size_t nRuns, Parameters *ps)
{
	if (ps->objectiveFunctionName.size() > 0) {
		ps->objectiveFunction = OBJ_MAP.at(ps->objectiveFunctionName);
	}

	cout << "\t" << "algorithm:\t\t " << ps->algorithm << endl;
	cout << "\t" << "objectiveFunction:\t " << ps->objectiveFunctionName << endl;
	cout << "\t" << "initialSize:\t\t " << ps->initialSize << endl;
	cout << "\t" << "maxGenerations:\t\t " << ps->maxGenerations << endl;
	cout << "\t" << "dimensions:\t\t " << ps->coordinateBounds.size() << endl;
	cout << "\t" << "bounds[0]:\t\t " << printBound(ps->coordinateBounds.front()) << endl;
	if (ps->knownOptimum != nullptr) {
		cout << "\t" << "knownOptimum:\t\t " << printMember(*ps->knownOptimum) << endl;
		cout << "\t" << "globalMinimum:\t\t " << ps->objectiveFunction(*ps->knownOptimum) << endl;
	}

	vector<MemberWithValue> results;
	results.reserve(nRuns);

	for (size_t i = 0; i < nRuns; ++i) {

		auto finalPopulation = run(ps);

		auto bestSolution = finalPopulation[0];

		double objectiveValue = ps->objectiveFunction(bestSolution);

		cout << "run: " << (i + 1) << ", objectiveValue=" << objectiveValue << endl;

		results.push_back(MemberWithValue(bestSolution, objectiveValue));
	}

	auto bestSolution = min_element(results.begin(), results.end(), compareMemberWithValueLower);

	cout << "--------" << endl;
	cout << "bestSolution: " << printMember(get<0>(*bestSolution)) << endl;
	cout << "bestObjective: " << get<1>(*bestSolution) << endl;

	return 0;
}
