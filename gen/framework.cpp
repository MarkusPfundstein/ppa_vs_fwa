#include <iostream>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include "framework.h"
#include "ppa.h"
#include "fwa.h"
#include "ppa-sbppa.h"

map<string, function<Population(Parameters*, ValueCollector&)>> ALGO_MAP = {
	{ "ppa", runPPA },
	{ "fwa", runFWA },
	{ "ppa-sbppa", runPPA_sbPPA },
	{ "ppalevy", runPPALevy }
};

map <string, function<double(const Member&)>> OBJ_MAP = {
	{ "rosenbrock", rosenbrock },
	{ "griewank", griewank },
	{ "schwefel", schwefel1_2 },
	{ "schwefel7", schwefel7 },
	{ "easom", easom },
	{ "ackley", ackleys_path },
	{ "michalewicz12", michalewicz12 },
	{ "sphere", sphere },
	{ "rastrigrin", rastrigrin },
	{ "ellipse", ellipse },
	{ "cigar", cigar },
	{ "tablet", tablet },
	{ "schwefelFWA", schwefelFWA }
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
	return "{}";

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

void writeJsonStart(ostream &out, Parameters *ps)
{
	out << "{" << endl;
	out << "\t\"params\":{" << endl;

	out << "\t\t\"algorithm\":\"" << ps->algorithm << "\"," << endl;
	out << "\t\t\"objectiveFunction\":\"" << ps->objectiveFunctionName << "\"," << endl;
	out << "\t\t\"initialSize\": " << ps->initialSize << "," << endl;
	out << "\t\t\"maxGenerations\":" << ps->maxGenerations << "," << endl;
	out << "\t\t\"maxFevals\":" << ps->maxFunctionEvaluations << "," << endl;
	out << "\t\t\"shift\":" << ps->shiftValue << "," << endl;
	out << "\t\t\"dimensions\":" << ps->coordinateBounds.size() << "," << endl;
	out << "\t\t\"bounds[0]\": [" << printBound(ps->coordinateBounds.front()) << "]" << "," << endl;
	out << "\t\t\"initBounds[0]\": [" << printBound(ps->initBounds.front()) << "]";
	if (ps->knownOptimum != nullptr) {
		out << "," << endl;
		out << "\t\t\"knownOptimum\": [" << printMember(*ps->knownOptimum) << "]" << "," << endl;
		out << "\t\t\"globalMinimum\": " << ps->objectiveFunction(*ps->knownOptimum) << endl;
	}
	else {
		out << endl;
	}
	out << "\t}," << endl;		// end params
	out << "\t\"runs\" : [" << endl; 
}

void writeJsonEnd(ostream &out, double overallBestValue, Member bestMember, double mean, double sd)
{
	out << "\t]," << endl; // end "runs"
	out << "\t\"stats\": {" << endl;
	out << "\t\t\"mean\": " << mean << "," <<  endl;
	out << "\t\t\"sd\": " << sd << "," << endl;
	out << "\t\t\"best\": " << overallBestValue << "," << endl;
	out << "\t\t\"bestMember\": [";
	out << printMember(bestMember);
	out << "]" << endl;
	out << "\t}" << endl;  // end "stats
	out << "}" << endl;
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
	srand((unsigned int)time(NULL));


	if (ps->objectiveFunctionName.size() > 0) {
		ps->objectiveFunction = OBJ_MAP.at(ps->objectiveFunctionName);
	}

	cout << "algorithm:\t\t " << ps->algorithm << endl;
	cout << "objectiveFunction:\t " << ps->objectiveFunctionName << endl;
	cout << "initialSize:\t\t " << ps->initialSize << endl;
	cout << "maxGenerations:\t\t " << ps->maxGenerations << endl;
	cout << "maxFevals:\t\t " << ps->maxFunctionEvaluations << endl;
	cout << "shiftValue:\t\t " << ps->shiftValue << endl;
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
		writeJsonStart(out, ps);
	}

	setShift(ps->shiftValue);
	cout << "getShift: " << getShift() << endl;

	cout << "start experiment, runs: " << nRuns << endl;

	vector<MemberWithValue> results;
	results.reserve(nRuns);

	double errorSum = 0.0;
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
		printf("#%d:\t%.9e\n", (i + 1), objectiveValue);
		cout << "fevals:      " << collector.numberFunctionEvaluations << "/" << ps->maxFunctionEvaluations << endl;
		cout << "generations: " << collector.numberGenerations << "/" << ps->maxGenerations << endl;

		errorSum += objectiveValue;

		results.push_back(MemberWithValue(bestSolution, objectiveValue));
	}

	double average = errorSum / nRuns;
	double quadratic = 0.0;
	for (size_t i = 0; i < nRuns; i++)
	{
		quadratic += (average - get<1>(results[i])) * (average - get<1>(results[i]));
	}
	double deviation = sqrt(quadratic / nRuns);

	auto bestSolution = min_element(results.begin(), results.end(), compareMemberWithValueLower);

	cout << "--------" << endl;
	cout << "bestSolution: (" << printMember(get<0>(*bestSolution)) << ")" << endl;
	cout << "bestObjective: " << get<1>(*bestSolution) << endl;
	printf("Error is %.9e +/- %.9e\n", average, deviation);


	if (out.is_open()) {
		writeJsonEnd(out, get<1>(*bestSolution), get<0>(*bestSolution), average, deviation);
		out.close();
	}


	//cout << "enter something to exit... " << endl;
	//char temps[256];
	//cin >> temp;

	return 0;
}
