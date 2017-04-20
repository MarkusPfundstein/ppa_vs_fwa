#pragma once

#include <string>
#include <map>
#include <functional>
#include <string>
#include <sstream>
#include "population.h"
#include "functions.h"

struct ValueCollector {
	vector<Member> bestMembersInGeneration;
};

struct Parameters {
	// set bounds (min, max) for coordinate system. determines dimension of problem.
	vector<CoordBound> coordinateBounds;

	// name of algorithm
	string algorithm;

	// initial size of population. 
	size_t initialSize;

	// how many generations do we want to go on before we stop?
	size_t maxGenerations;

	// name of the objective function
	string objectiveFunctionName;
	
	// the objective function we want to evaluate
	function<double(const Member&)> objectiveFunction;

	// optional: Member of population that is known to be optimum. 
	// if set we can control to stop when we are in certain area of optimum.
	Member *knownOptimum = nullptr;

	virtual string printParameters() = 0;

	Parameters(string alg) : algorithm(alg) {}
};

/* algorithm specific */
struct PPA : public Parameters {
	// maximum number of runners per plant
	size_t nMax = 5;

	// function to normalize objective function value into range (0, 1). Defaults to minMaxFitness
	function<double(double, double, double)> fitnessFunction = minMaxFitness;

	string printParameters()
	{
		stringstream ss;
		ss << "nMax:\t\t\t " << nMax;
		return ss.str();
	}

	PPA() : Parameters("ppa") {}
};

struct FWA : public Parameters {
	// maximum amplitude per explosion (area in which we generate new sparks)
	double Amax = 40.0;
	
	// maximum number of sparks that an explosion can generate
	size_t maxSparks = 50;
	// parameters to control upper and lower bound for spark generation
	double a = 0.04;
	double b = 0.8;

	// maximum number of gaussian mutations that we want to generate per iteration
	size_t gaussianMutations = 5;

	// distance function to calculate distance between Members. Defaults to euclideanDistance
	function<double(const Member&, const Member&)> distanceFunction = euclideanDistance;

	string printParameters()
	{
		stringstream ss;
		ss << "Amax:\t\t\t " << Amax << endl;
		ss << "maxSparks:\t\t " << maxSparks << endl;
		ss << "a:\t\t\t " << a << endl;
		ss << "b:\t\t\t " << b << endl;
		ss << "gaussianM:\t\t " << gaussianMutations << endl;
		return ss.str();
	}

	FWA() : Parameters("fwa") {}
};

int runExperiments(size_t n, Parameters *ps, string writeValuesPath = "");

template <typename T> 
const T& castParameters(Parameters *ps)
{
	return *static_cast<T*>(ps);
}

