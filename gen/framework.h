#pragma once

#include <string>
#include <map>
#include <functional>
#include "population.h"
#include "functions.h"

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
};

/* algorithm specific */
struct PPA : public Parameters {
	// maximum number of runners per plant
	size_t nMax;

	// function to normalize objective function value into range (0, 1). Defaults to minMaxFitness
	function<double(double, double, double)> fitnessFunction;

	PPA() : Parameters() {
		algorithm = "ppa";
		fitnessFunction = minMaxFitness;
	}
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
	function<double(const Member&, const Member&)> distanceFunction;

	FWA() : Parameters() {
		algorithm = "fwa";
		distanceFunction = euclideanDistance;
	}
};

Population run(Parameters *params);

int runExperiments(size_t n, Parameters *ps);

template <typename T> 
const T& castParameters(Parameters *ps)
{
	return *static_cast<T*>(ps);
}

