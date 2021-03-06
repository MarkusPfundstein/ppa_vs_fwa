#include "optionparser.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <memory>
#include "framework.h"

using namespace std;

template<typename Out>
void split(const string &s, char delim, Out result) {
	stringstream ss;
	ss.str(s);
	string item;
	while (std::getline(ss, item, delim)) {
		*(result++) = item;
	}
}


vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, back_inserter(elems));
	return elems;
}

string getStringArg(option::Option* options, int enumName, string defaultValue = "")
{
	auto opt = options[enumName];
	if (!opt) {
		return defaultValue;
	}
	return string(static_cast<option::Option>(opt).arg);
}

int getNumericArg(option::Option* options, int enumName, int defaultValue = -1)
{
	string arg = getStringArg(options, enumName);
	if (arg == "") {
		return defaultValue;
	}
	return stoi(arg);
}

//--algorithm = ppa --function = rosenbrock --initial - size = 30 --dimensions = 10
enum  optionIndex { 
	UNKNOWN, 
	HELP, 
	ALGORITHM, 
	FUNCTION, 
	INITIAL_SIZE, 
	DIMENSIONS,
	NUMBER_RUNS,
	MAX_GENERATIONS,
	MAX_FEVALS,
	UNI_MIN_BOUND,
	UNI_MAX_BOUND,
	INIT_MIN_BOUND,
	INIT_MAX_BOUND,
	KNOWN_OPTIMUM,
	PPA_NMAX,
	FWA_AMAX,
	FWA_MAX_SPARKS,
	WRITE_VALUES,

	SHIFT
};

const option::Descriptor usage[] = {
	{ UNKNOWN, 0,"", "",        Arg::Unknown, "USAGE: solver-sop [options]\n\n"
	"Options:" },
	{ HELP,    0,"", "help",    Arg::None,    "  \t--help  \tPrint usage and exit." },
	{ WRITE_VALUES,    0, "",  "write-values",     Arg::Optional,     "  --write-values \tfile in which we write intermediate values" },

	/* GENERAL */
	{ ALGORITHM,       0, "a", "algorithm",        Arg::NonEmpty, "  -a <arg>, \t--algorithm=<arg>  \tAlgorithm to use." },
	{ FUNCTION,        0, "f", "function",         Arg::NonEmpty, "  -f <arg>, \t--function=<arg>   \tSOP function to use." },
	{ INITIAL_SIZE,    0, "n", "initial-size",     Arg::Numeric,  "  -n <arg>, \t--initial-size=<arg>  \tinitial size of population." },
	{ DIMENSIONS,      0, "a", "dimensions",       Arg::Numeric,  "  -d <arg>, \t--dimensions=<arg>  \thow many dimensions. function must support it" },
	{ NUMBER_RUNS,     0, "t", "number-runs",      Arg::Numeric,  "  -t <arg>, \t--number-runs=<arg>  \thow often do we want to run the experiment?" },
	{ MAX_GENERATIONS, 0, "g", "max-generations",  Arg::Numeric,  "  -g <arg>, \t--max-generations=<arg>  \thow often do we want to reproduce?" },
	{ MAX_FEVALS,      0,   "", "max-fevals",      Arg::Numeric,  "  --max-fevals" },
	{ UNI_MIN_BOUND,   0, "u", "min-bound",        Arg::Numeric,  "  -u <arg>, \t--min-bound=<arg>  \tound for uniform coordinate system" },
	{ UNI_MAX_BOUND,   0, "v", "max-bound",        Arg::Numeric,  "  -v <arg>, \t--max-bound=<arg>  \tbound for uniform coordinate system" },
	{ KNOWN_OPTIMUM,   0, "o", "known-optimum",    Arg::Optional, "  -o x1,x2,...\t--known-optimum=x1,x2,..." },
	{ INIT_MIN_BOUND,  0, "",  "init-min-bound",   Arg::Numeric,  " --init-min-bound=<arg>" },
	{ INIT_MAX_BOUND,  0, "",  "init-max-bound",   Arg::Numeric,  " --init-max-bound=<arg>" },

	/* PPA SPECIFIC */
	{ PPA_NMAX,        0, "",  "ppa-nmax",         Arg::Numeric,  " --ppa-nmax=<arg>\tmax runners per plant" },
	
	/* FWA SPECIFIC */
	{ FWA_AMAX,        0, "",  "fwa-amax",         Arg::Numeric,  " --fwa-amax=<arg>\tmax amplitude per explosion" },
	{ FWA_MAX_SPARKS,  0, "",  "fwa-max-sparks",   Arg::Numeric,  " --fwa-max-sparks=<arg>\tmax sparks per firework" },

	{ SHIFT,  0, "",   "origin-shift",   Arg::Numeric,  " --origin-shift=<arg>\twanna shift the origin? put a number here" },

	{ UNKNOWN, 0,"", "",        Arg::None, "no idea what you mean" },
	{ 0, 0, 0, 0, 0, 0 } 
};

bool getParameters(option::Option *options, Parameters **p)
{
	string algorithm = getStringArg(options, ALGORITHM);
	string objectiveFunction = getStringArg(options, FUNCTION);
	int dimensions = getNumericArg(options, DIMENSIONS, -1);
	int initialSize = getNumericArg(options, INITIAL_SIZE, 30);
	int maxGenerations = getNumericArg(options, MAX_GENERATIONS, 30);
	int maxFunctionEvals = getNumericArg(options, MAX_FEVALS, 10000 * dimensions);
	int minBound = getNumericArg(options, UNI_MIN_BOUND, -1);
	int maxBound = getNumericArg(options, UNI_MAX_BOUND, -1);
	int initMinBound = getNumericArg(options, INIT_MIN_BOUND, minBound);
	int initMaxBound = getNumericArg(options, INIT_MAX_BOUND, maxBound);
	int ppaNmax = getNumericArg(options, PPA_NMAX, 5);
	int fwaAmax = getNumericArg(options, FWA_AMAX, 40);
	int fwaMaxSparks = getNumericArg(options, FWA_MAX_SPARKS, 50);
	int shift = getNumericArg(options, SHIFT, 0);	

	string knownOptimumString = getStringArg(options, KNOWN_OPTIMUM);

	if (dimensions == -1) {
		cerr << "dimension must be set" << endl;
		return false;
	}
	if (minBound == -1 || maxBound == -1) {
		cerr << "uni-min-bound and uni-max-bound must both be set" << endl;
		return false;
	}

	if (algorithm == "ppa") {
		auto *ppa = new PPA();
		ppa->nMax = ppaNmax;

		*p = static_cast<Parameters*>(ppa);
	}
	else if (algorithm == "ppalevy") {
		auto *params = new PPA();
		params->nMax = ppaNmax;

		*p = static_cast<Parameters*>(params);
	}
	else if (algorithm == "ppalevy-birds") {
		auto *params = new PPA();
		params->nMax = ppaNmax;

		*p = static_cast<Parameters*>(params);
	}
	else if (algorithm == "ppa-sbppa") {
		auto *params = new PPA_sbPPA();

		*p = static_cast<Parameters*>(params);
	}
	else if (algorithm == "fwa") {
		auto *fwa = new FWA();
		fwa->Amax = fwaAmax;
		fwa->maxSparks = fwaMaxSparks;
		*p = static_cast<Parameters*>(fwa);
	}
	else {
		//throw new exception("unknown algorithm");
		cerr << "[cli.cpp] unknown algorithm: " << algorithm << endl;
		return false;
	}
	// from here on we must always jump to errorAndDeleteP

	if (knownOptimumString.size() > 0) {
		// brat , what a shitty code
		
		auto tokens = split(knownOptimumString, ',');
		if (tokens.size() < static_cast<size_t>(dimensions)) {
			cerr << "knownOptimum dimensionality must be equal to dimensions" << endl;
			goto errorAndDeleteP;
		}
		if (tokens.size() > static_cast<size_t>(dimensions)) {
			cout << "WARNING: knownOptimum.dim > dimensionality" << endl;
		}

		auto *member = new Member();
		for (auto &t : tokens) {
			double f = stod(t);
			member->push_back(f);
		}
		(*p)->knownOptimum = member;
	}

	// standard assignments
	(*p)->algorithm = algorithm;
	(*p)->objectiveFunctionName = objectiveFunction;
	(*p)->initialSize = initialSize;
	(*p)->coordinateBounds = createUniformCoordinateBounds(dimensions, minBound, maxBound);
	(*p)->initBounds = createUniformCoordinateBounds(dimensions, initMinBound, initMaxBound);
	(*p)->maxGenerations = maxGenerations;
	(*p)->maxFunctionEvaluations = maxFunctionEvals;
	(*p)->shiftValue = shift;
	
	return true;

errorAndDeleteP: 
	delete (*p);
	return false;
}

int main(int argc, char* argv[])
{
	cout.precision(9);

	argc -= (argc > 0); argv += (argc > 0); // skip program name argv[0] if present
	option::Stats stats(usage, argc, argv);
	option::Option* options = new option::Option[stats.options_max];
	option::Option* buffer = new option::Option[stats.buffer_max];

	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
		return 1;

	if (options[HELP] || argc == 0)
	{
		return 0;
	}

	int numberRuns = getNumericArg(options, NUMBER_RUNS, 10);
	string writeValuesPath = getStringArg(options, WRITE_VALUES, "");

	Parameters *p;
	if (getParameters(options, &p)) {
		runExperiments(numberRuns, p, writeValuesPath);
		// yeah yeah I know. should be in deconstructor etc.. just hacking atm
		if (p->knownOptimum) {
			delete p->knownOptimum;
		}
		delete p;
	}

	delete[] options;
	delete[] buffer;

	return 0;
}
