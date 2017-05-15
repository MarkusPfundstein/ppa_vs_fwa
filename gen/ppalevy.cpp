#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include "rng.h"
#include "ppa.h"

using namespace std;

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

vector<Member> computeNewRunners2(double nMax, const MemberWithValue &m, const vector<CoordBound> &bounds)
{
	vector<Member> newRunners;

	Member member = get<0>(m);
	double ni = get<1>(m);

	// n_r = number of runners to generate for solution i in current population
	size_t nr = static_cast<size_t>(ceil(nMax * ni * fRand(0, 1)));
	if (nr < 1) {
		nr = 1;
	}

	for (size_t r = 0; r < nr; ++r) {
		Member newMember;
		newMember.reserve(member.size());

		int fireworksDirections = rand() % bounds.size() + 1;


		// try to not pertubate *every* dimension
		vector<int> flags;
		for (size_t i = 0; i < bounds.size(); ++i) {
			flags.push_back(0);
		}
		int fireworksSelectedIndex = 0;
		for (int k = 0; k < fireworksDirections; k++) {
			do
			{
				fireworksSelectedIndex = rand() % bounds.size();
			} while (flags[fireworksSelectedIndex]);

			flags[fireworksSelectedIndex] = 1;
		}

		for (size_t j = 0; j < member.size(); ++j) {
			if (flags[j] == 1) {
				// distance is inversely proportional to fitness
				double drj = 2.0 * (1.0 - ni) * (fRand(0, 1) - 0.5);
				// each drj will be in (-1, 1)

				double xj = member[j];
				auto boundj = bounds[j];

				double aj = get<0>(boundj);
				double bj = get<1>(boundj);

				// update solution
				// TO-DO: Validate mapping back to solution space. I think this here is a bit crude
				double newxj = xj + (bj - aj) * drj;
				if (newxj < aj) {
					newxj = aj;
				}
				if (newxj > bj) {
					newxj = bj;
				}

				newMember.push_back(newxj);
			}
			else {
				newMember.push_back(member[j]);
			}
		}
		newRunners.push_back(newMember);
	}

	return newRunners;
}

Population runPPALevy(Parameters *ps, ValueCollector &vc)
{

	std::cout << "new levy with k-select" << std::endl;
	const PPA &params = castParameters<PPA>(ps);

	auto P = createRandomPopulation(params.initialSize, params.initBounds);

	int evals = 0;


	const unsigned A = 10;  // number birds


	size_t g = 0;
	for (; g < params.maxGenerations && (size_t)evals < params.maxFunctionEvaluations; ++g) {
		// sort P in descending order of N
		// Note: Sorting is only possible in ascending direction. Hence we have to reverse afterwards

		auto N = calculateFitnessForPopulation2(P, params.objectiveFunction, params.fitnessFunction, normalizeTan, &evals);

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
			int Ac = A;
			

			vector<MemberWithValue> survivors = {};
			/* TAKE k BEST & n-k random , k < n*/
			const int k = 3;
			for (size_t i = 0; i < k; ++i) {
				survivors.push_back(MemberWithValue(get<0>(N[i]), get<1>(N[i])));
			}
			vector<int> indicesSeen = {};
			for (size_t i = k; i < params.initialSize; ++i) {
				// random index between 3 and N.size
				int randIndex;
				do {
					randIndex = iRand(k, N.size() - 1);
					// randIndex not yet found
					if (find(indicesSeen.begin(), indicesSeen.end(), randIndex) == indicesSeen.end()) {
						indicesSeen.push_back(randIndex);
						break;
					}
				} while (1);

				auto &n = N[randIndex];
				survivors.push_back(MemberWithValue(get<0>(n), get<1>(n)));
			}


			for (size_t i = 0; i < survivors.size(); ++i) {
				// r_i <- set of runners where both the size of the set and the distance for each
				//        runner (individually) is proportional to the fitness N_i
				auto R = computeNewRunners2(params.nMax, survivors[i], params.coordinateBounds);

				// phi = phi UNION r
				// Note: we don't do union here. (see erase below)
				for (auto r : R) {
					phi.push_back(r);
				}
			}


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
