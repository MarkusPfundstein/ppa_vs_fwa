#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include "common.h"
#include "plant_propagation.h"

using namespace std;

double normalizeTan(double f) {
	return 0.5 * (tanh(4.0 * f - 2.0) + 1.0);
}


vector<MemberWithValue> calculateFitnessForPopulation(
	const Population& population,
	function<double(const Member&)> objF,
	function<double(double, double, double)> f,
	function<double(double)> normF)
{
	auto objectiveValues = evalObjectiveFunctionForPopulation(population, objF);
	//for (auto f : objectiveValues) cout << get<1>(f) << endl;

	auto minMaxMembers = minmax_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

	const double min = get<1>(*minMaxMembers.first);
	const double max = get<1>(*minMaxMembers.second);

	vector<MemberWithValue> fitnesses;
	fitnesses.reserve(objectiveValues.size());

	for (const auto &m : objectiveValues) {
		double fval = min == max ? 0.5 : f(get<1>(m), min, max);
		auto member = get<0>(m);
		fitnesses.push_back(MemberWithValue(member, normF(fval)));
	}

	return fitnesses;
}

vector<double> normalizeFitness(const vector<double> &fitnesses, function<double(double)> f)
{
	vector<double> normalizedFitnesses;
	normalizedFitnesses.reserve(fitnesses.size());

	transform(
		fitnesses.begin(),
		fitnesses.end(),
		back_inserter(normalizedFitnesses),
		f
	);

	return normalizedFitnesses;
}

vector<Member> computeNewRunners(double nMax, const MemberWithValue &m, const vector<CoordBound> &bounds)
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

Population runPlantPropagation1(
	const vector<CoordBound> coordinateBounds,
	const size_t m,
	const size_t gMax,
	const size_t nMax,
	function<double(double, double, double)> fitnessFunction,
	function<double(const Member&)> objectiveFunction
)
{
	auto P = createRandomPopulation(m, coordinateBounds);

	for (size_t g = 0; g < gMax; ++g) {
		// compute N_i = f(p_i) forall p in P
		auto N = calculateFitnessForPopulation(
			P,
			objectiveFunction,
			fitnessFunction,
			normalizeTan
		);

		// sort P in descending order of N
		// Note: Sorting is only possible in ascending direction. Hence we have to reverse afterwards
		stable_sort(
			N.begin(),
			N.end(),
			compareMemberWithValueLower
		);
		reverse(N.begin(), N.end());

		// create new Population phi
		Population phi;

		// if we are *not* in the last iteration, we create new runners
		if (g < gMax - 1) {
			for (size_t i = 0; i < m; ++i) {
				// r_i <- set of runners where both the size of the set and the distance for each
				//        runner (individually) is proportional to the fitness N_i
				auto R = computeNewRunners(nMax, N[i], coordinateBounds);

				// phi = phi UNION r
				// Note: we don't do union here. (see erase below)
				for (auto r : R) {
					phi.push_back(r);
				}
				// Note: erase not really necessary but can be used to integrate the union operation in the paper
				//phi.erase(unique(phi.begin(), phi.end(), compareMember), phi.end());
			}

			// P <- phi
			P.erase(P.begin(), P.end());
			copy(phi.begin(), phi.end(), back_inserter(P));
		}
		// if we are in the last iteration, copy members in (sorted) N into P
		else {
			P.erase(P.begin(), P.end());
			for (auto n : N) {
				P.push_back(get<0>(n));
			}
		}
	}

	return P;
}
