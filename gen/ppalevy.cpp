#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include "rng.h"
#include "ppa.h"

using namespace std;

const double PI = 3.14159265358979323846;

double Poisson(double Lambda, const unsigned NP)
{
	// need to adjust by 1 to not consider 0 agents
	int k = rand() % NP + 1;
	return poisson(Lambda, k);
}

double sigmaU(double beta)
{
	if (!(1.0 <= beta && beta <= 2.0)) {
		throw new runtime_error("beta out of range [1, 2]");
	}

	double numerator = tgamma(1.0 + beta) * sin(PI * beta * 0.5);
	double denominator = tgamma(0.5 * (1.0 + beta)) * beta * pow(2, 0.5 * (beta - 1));

	return (numerator / denominator) / beta;
}

double Levy(double val, double scale)
{
	double beta = 1.5;

	double sigma = sigmaU(beta);

	double u = sampleGaussian(0, sigma * sigma);
	double v = sampleGaussian(0, 1);

	double s = u / pow(abs(v), 1.0 / beta);
	return s * scale;
}

// O(2n)
static vector<MemberWithValue> calculateFitnessForPopulation2(
	const Population& population,
	function<double(const Member&)> objF,
	function<double(double, double, double)> f,
	function<double(double)> normF,
	int *evals)
{
	// O(n)
	auto objectiveValues = evalObjectiveFunctionForPopulation(population, objF, evals);
	//for (auto f : objectiveValues) cout << get<1>(f) << endl;

	auto minMaxMembers = minmax_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

	const double min = get<1>(*minMaxMembers.first);
	const double max = get<1>(*minMaxMembers.second);

	vector<MemberWithValue> fitnesses;
	fitnesses.reserve(objectiveValues.size());

	// O(n)
	for (const auto &m : objectiveValues) {
		double fval = min == max ? 0.5 : f(get<1>(m), min, max);
		auto member = get<0>(m);
		fitnesses.push_back(MemberWithValue(member, normF(fval)));
	}

	return fitnesses;
}

vector<Member> computeDisplacedSeeds2(double k, const MemberWithValue &m, const vector<CoordBound> &bounds)
{
	vector<Member> newRunners;

	Member member = get<0>(m);
	double ni = get<1>(m);

	// n_r = number of runners to generate for solution i in current population
	size_t nr = static_cast<size_t>(ceil(k * ni * fRand(0, 1)));

	for (size_t r = 0; r < nr; ++r) {
		Member newMember;
		newMember.reserve(member.size());

		for (size_t j = 0; j < member.size(); ++j) {
			double xj = member[j];

			if (fRand(0, 1) < 0.8) {

				auto boundj = bounds[j];

				double aj = get<0>(boundj);
				double bj = get<1>(boundj);

				// update solution
				// TO-DO: Validate mapping back to solution space. I think this here is a bit crude
				double phi = fRand(aj, bj);
				double newxj = xj + Levy(xj - phi, 40);		// 40 = amplitude
				if (newxj < aj) {
					newxj = aj;
				}
				if (newxj > bj) {
					newxj = bj;
				}
				newMember.push_back(newxj);
			}
			else {
				newMember.push_back(xj);
			}


		}
		newRunners.push_back(newMember);
	}

	return newRunners;
}

vector<Member> computeNewRunners2(double nMax, const MemberWithValue &m, const vector<CoordBound> &bounds)
{
	vector<Member> newRunners;

	Member member = get<0>(m);
	double ni = get<1>(m);

	// n_r = number of runners to generate for solution i in current population
	size_t nr = static_cast<size_t>(ceil(nMax * ni * fRand(0, 1)));

	for (size_t r = 0; r < nr; ++r) {
		Member newMember;
		newMember.reserve(member.size());

		for (size_t j = 0; j < member.size(); ++j) {
			// distance is inversely proportional to fitness
			double drj = 2.0 * (1.0 - ni) * (fRand(0, 1) - 0.5);
			// each drj will be in (-1, 1)

			double xj = member[j];
			auto boundj = bounds[j];

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

			newMember.push_back(newX);
		}
		newRunners.push_back(newMember);
	}

	return newRunners;
}

Population runPPALevy(Parameters *ps, ValueCollector &vc)
{
	const PPA &params = castParameters<PPA>(ps);

	auto P = createRandomPopulation(params.initialSize, params.initBounds);

	int evals = 0;
	size_t g = 0;
	for (; g < params.maxGenerations && (size_t)evals < params.maxFunctionEvaluations; ++g) {
		// compute N_i = f(p_i) forall p in P
		auto N = calculateFitnessForPopulation2(
			P,
			params.objectiveFunction,
			params.fitnessFunction,
			normalizeTan,
			&evals
		);

		// sort P in descending order of N
		// Note: Sorting is only possible in ascending direction. Hence we have to reverse afterwards
		stable_sort(
			N.begin(),
			N.end(),
			compareMemberWithValueLower
		);
		reverse(N.begin(), N.end());

		auto best = get<0>(N.front());

		vc.bestMembersInGeneration.push_back(best);


		// create new Population phi
		Population phi;
		phi.push_back(best);

		// if we are *not* in the last iteration, we create new runners
		if (g < params.maxGenerations - 1) {
			for (size_t i = 0; i < params.initialSize; ++i) {
				// r_i <- set of runners where both the size of the set and the distance for each
				//        runner (individually) is proportional to the fitness N_i
				auto R = computeNewRunners2(params.nMax, N[i], params.coordinateBounds);

				// phi = phi UNION r
				// Note: we don't do union here. (see erase below)
				for (auto r : R) {
					phi.push_back(r);
				}

				
				// compute birds arriving and doing a levy
				unsigned k = iRand(1, 10);
				//if (Poisson(1.1, k) > 0.05) {
				auto R2 = computeDisplacedSeeds2(k, N[i], params.coordinateBounds);
				for (auto r : R2) {
					phi.push_back(r);
				}
				
			}

			// P <- phi
			P = phi;
		}
		// if we are in the last iteration, copy members in (sorted) N into P
		else {
			P.erase(P.begin(), P.end());
			for (auto n : N) {
				P.push_back(get<0>(n));
			}
		}
	}

	vc.numberFunctionEvaluations = evals;
	vc.numberGenerations = g;

	return P;
}
