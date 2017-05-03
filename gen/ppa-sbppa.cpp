#include <algorithm>
#include <unordered_set>
#include <iostream>

#include "ppa-sbppa.h"
#include "population.h"
#include "rng.h"



Population runPPA_sbPPA(Parameters *ps, ValueCollector &vc)
{
	throw new runtime_error("ppa-sbppa sucks");
#if 0
	const PPA_sbPPA &params = castParameters<PPA_sbPPA>(ps);
	


	const size_t NP = params.initialSize; // population size
	const size_t n = params.coordinateBounds.size();	// dimensions
	int evals = 0;	  // objective function evaluations so far

	Population Pbest;
	Pbest.reserve(NP);

	// do NP trial runs -> create NP random populations, evaluate fitness for each and put best one into Pbest
	for (size_t r = 0; r < NP; ++r) {
		auto P = createRandomPopulation(params.initialSize, params.initBounds);
		auto objectiveValues = evalObjectiveFunctionForPopulation(P, params.objectiveFunction/*, params.fitnessFunction, normalizeTan*/, &evals);
		auto m = *max_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);
		Pbest.push_back(get<0>(m));
	}

	size_t g;
	for (g = 0; g < params.maxGenerations && (size_t)evals < params.maxFunctionEvaluations; ++g) {		

		auto objectiveValues = evalObjectiveFunctionForPopulation(Pbest, params.objectiveFunction, &evals);
		auto m = *max_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);
		std::cout << "best: " << get<1>(m) << std::endl;

		for (size_t i = 0; i < NP; ++i) {
			// if Poiss(lambda)_i > 0.05)
			if (P(params.Lambda, params.A) > params.poissonThreshold) {
				for (size_t j = 0; j < n; ++j) {
					// do we update coordinate at j?
					if (fRand(0, 1) < params.PR) {
						
						double xj = Pbest[i][j];
						auto boundj = params.coordinateBounds[j];

						double aj = get<0>(boundj);
						double bj = get<1>(boundj);

						const double phi = fRand(aj, bj); // random variable within search space
						double step = phi * L(xj - phi, n);
						double newX = xj + step;
						if (newX < aj) {
							newX = aj;
						}
						if (newX > bj) {
							newX = bj;
						}
						Pbest[i][j] = newX;
						
					}
				}
			}
			else {
				auto F = minmax_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

				const double Fmin = get<1>(*F.first);
				const double Fmax = get<1>(*F.second);

				vector<MemberWithValue> fitnesses;
				fitnesses.reserve(objectiveValues.size());

				double v = get<1>(objectiveValues[i]);
				double fx = Fmin == Fmax ? 0.5 : minMaxFitness(v, Fmin, Fmax);
				double ni = normalizeTan(fx);

				for (size_t j = 0; j < n; ++j) {
					double drj = 2.0 * (1.0 - ni) * (fRand(0, 1) - 0.5);
					// each drj will be in (-1, 1)

					double xj = Pbest[i][j];
					auto boundj = params.coordinateBounds[j];

					double aj = get<0>(boundj);
					double bj = get<1>(boundj);

					// update solution
					// TO-DO: Validate mapping back to solution space. I think this here is a bit crude
					double newX = xj + (bj - aj) * drj;
					if (newX < aj) {
						newX = aj;
					}
					if (newX > bj) {
						newX = bj;
					}
					Pbest[i][j] = newX;
				}
			}
		}
		
	}

	vc.numberFunctionEvaluations = evals;
	vc.numberGenerations = g;

	return Pbest;
#endif
}

