#include <algorithm>
#include <unordered_set>
#include <iostream>

#include "ppa-sbppa.h"
#include "population.h"
#include "rng.h"

const double PI = 3.14159265358979323846;

double P(double Lambda, const unsigned NP)
{
       // need to adjust by 1 to not consider 0 agents
	int k = rand() % NP + 1;
	return poisson(Lambda, k);
}

double sigmaU(double beta)
{
//	sigma = (gamma(1 + beta)*sin(pi*beta / 2) / (gamma((1 + beta) / 2)*beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);


	//double numerator = tgamma(1.0 + beta) * sin(PI * beta / 2.0);
	//double denominator = tgamma((1.0 + beta) / 2.0) * beta * pow(2, (beta - 1) / 2.0);
	double numerator = tgamma(1.0 + beta) * sin(PI * beta / 2.0);
	double denominator = tgamma((1.0 + beta) / 2.0) * beta * pow(2, (beta - 1) / 2.0);
	
	return pow(numerator / denominator, 1.0 / beta);
}

double Levy(double beta = 1.5)
{
	double sigma = sigmaU(beta);
	
	// matlab
	double u = sampleGaussian(0, 1) * sigma;
	double v = sampleGaussian(0, 1);

	// paper 
	//double u = sampleGaussian(0, sigma * sigma);
	//double v = sampleGaussian(0, 1);

	double step = u / pow(abs(v), 1.0 / beta);
	return step;
}

Population runPPA_sbPPA(Parameters *ps, ValueCollector &vc)
{
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

	vector<double> pDist;
	pDist.reserve(NP);

	for (size_t i = 0; i < NP; ++i) {
		// Matlab: 	//0.11	0.19	0.220	0.183	0.122	0.0679659986039712	0.0323	0.013	0.0049	0.00166
		// 3.3 with NP = 10 gives me same results as matlab code
		auto d = poisson(3.3, i + 1);       
		pDist.push_back(d);
	}

	for (g = 0; g < params.maxGenerations && (size_t)evals < params.maxFunctionEvaluations; ++g) {		
		auto objectiveValues = evalObjectiveFunctionForPopulation(Pbest, params.objectiveFunction, &evals);

		/*
		if (g % 1000 == 0) {

			auto bestSolution = *min_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);
			std::cout << "best: " << get<1>(bestSolution) << std::endl;

		}
		*/
		for (size_t i = 0; i < NP; ++i) {
			// if Poiss(lambda)_i > 0.05)
			//     p2 = poisspdf(p1,npop/3);  <- poisspdf(k, lamda)
			auto &X = Pbest[i];
			Member Xnew(X);
			// EXPLORATION
			if (pDist[i] > params.poissonThreshold) {
				for (size_t j = 0; j < n; ++j) {
					// do we update coordinate at j?
					if (fRand(0, 1) < params.PR) {
						
						double xj = X[j];
						auto boundj = params.coordinateBounds[j];

						double aj = get<0>(boundj);
						double bj = get<1>(boundj);

		
						//0.01*step*(pop(i,s)-((rand(1)*b(s))+a(s)));
						double step = Levy();
						double stepSize = 0.01 * step * (xj - ((fRand(0, 1) * bj + aj)));
						double newX = xj + stepSize;
						if (newX < aj) {
							newX = aj;
						}
						if (newX > bj) {
							newX = bj;
						}
						Xnew[j] = newX;
					}
				}
			}

			// EXPLOITATION
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

					double xj = X[j];
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
					Xnew[j] = newX;
				}
			}

			double valCur = params.objectiveFunction(X);
			double valNew = params.objectiveFunction(Xnew);
			if (valNew < valCur) {
				Pbest[i] = Xnew;
			}
		}
		
	}

	vc.numberFunctionEvaluations = evals;
	vc.numberGenerations = g;

	return Pbest;
}

