#include <iostream>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <fstream>
#include "framework.h"
#include "ppa.h"
#include "fwa.h"

map<string, function<Population(Parameters*, ValueCollector&)>> ALGO_MAP = {
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

Population run(Parameters *ps, ValueCollector &vc, double *timeTakenMs)
{
	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
	
	Population p = ALGO_MAP.at(ps->algorithm)(ps, vc);

	chrono::steady_clock::time_point end = chrono::steady_clock::now();

	*timeTakenMs = chrono::duration_cast<chrono::microseconds>(end - begin).count() / 1000.0;

	return p;
}

string writeCollectorData(const ValueCollector &collector, double timeTakenMs, const Parameters& ps)
{
	stringstream ss;

	const size_t n = collector.bestMembersInGeneration.size() - 1;
	ss << "{" << endl;
	ss << "\t\"bestMembersInGeneration\": {" << endl;
	ss << "\t\t\"coordinates\": [" << endl;
	size_t i = 0;
	for (auto d : collector.bestMembersInGeneration) {
		ss << "\t\t\t[" << printMember(d) << "]";
		if (i++ < n) {
			ss << ",";
		}
		ss << endl;
	}
	ss << "\t\t], " << endl;
	ss << "\t\t\"values\": [" << endl;
	i = 0;
	for (auto d : collector.bestMembersInGeneration) {
		ss << "\t\t\t" << ps.objectiveFunction(d);
		if (i++ < n) {
			ss << ",";
		}
		ss << endl;
	}
	ss << "\t\t]" << endl;
	ss << "\t}," << endl;
	ss << "\t\"timeTaken\": " << timeTakenMs << endl;
	ss << "}" << endl;

	return ss.str();
}

void writeJsonStart(ostream &out)
{
	out << "[" << endl; 
}

void writeJsonEnd(ostream &out)
{
	out << "]" << endl;
}

void writeJsonOutput(const string& s, ostream &out, size_t run, size_t runMax)
{
	out << s;
	if (run < (runMax - 1)) {
		out << ", ";
	}
	out << endl;
}


int runExperiments(size_t nRuns, Parameters *ps, string writeValuesPath)
{
	if (ps->objectiveFunctionName.size() > 0) {
		ps->objectiveFunction = OBJ_MAP.at(ps->objectiveFunctionName);
	}

	cout << "algorithm:\t\t " << ps->algorithm << endl;
	cout << "objectiveFunction:\t " << ps->objectiveFunctionName << endl;
	cout << "initialSize:\t\t " << ps->initialSize << endl;
	cout << "maxGenerations:\t\t " << ps->maxGenerations << endl;
	cout << "dimensions:\t\t " << ps->coordinateBounds.size() << endl;
	cout << "bounds[0]:\t\t (" << printBound(ps->coordinateBounds.front()) << ")" << endl;
	cout << "initBounds[0]:\t\t (" << printBound(ps->initBounds.front()) << ")" << endl;
	if (ps->knownOptimum != nullptr) {
		cout << "knownOptimum:\t\t (" << printMember(*ps->knownOptimum) << ")" << endl;
		cout << "globalMinimum:\t\t " << ps->objectiveFunction(*ps->knownOptimum) << endl;
	}
	cout << ps->printParameters() << endl;

	ofstream out;
	if (writeValuesPath != "") {
		out.open(writeValuesPath, ostream::out);
		if (!out.is_open()) {
			cerr << "error opening: " << writeValuesPath << endl;
			return -1;
		}
		cout << "write intermediate values to: " << writeValuesPath << endl;
		writeJsonStart(out);
	}

	cout << "start experiment, runs: " << nRuns << endl;

	vector<MemberWithValue> results;
	results.reserve(nRuns);

	for (size_t i = 0; i < nRuns; ++i) {

		ValueCollector collector;
		collector.bestMembersInGeneration.reserve(ps->maxGenerations);
		
		double timeTakenMs;

		auto finalPopulation = run(ps, collector, &timeTakenMs);

		if (out.is_open()) {
			writeJsonOutput(writeCollectorData(collector, timeTakenMs, *ps), out, i, nRuns);
		}

		auto bestSolution = finalPopulation[0];

		double objectiveValue = ps->objectiveFunction(bestSolution);

		cout << "run: " << (i + 1) << ", objectiveValue=" << objectiveValue << ", ms: " << timeTakenMs <<  endl;

		results.push_back(MemberWithValue(bestSolution, objectiveValue));
	}

	if (out.is_open()) {
		writeJsonEnd(out);
		out.close();
	}

	auto bestSolution = min_element(results.begin(), results.end(), compareMemberWithValueLower);

	cout << "--------" << endl;
	cout << "bestSolution: (" << printMember(get<0>(*bestSolution)) << ")" << endl;
	cout << "bestObjective: " << get<1>(*bestSolution) << endl;

	return 0;
}
