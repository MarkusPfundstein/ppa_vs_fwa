#include <algorithm>
#include <unordered_set>
#include <iostream>

#include "fwa.h"
#include "population.h"
#include "rng.h"

// fwa consists of:
// - explosion operator
//   how do we evaluate explosion of sparks. 
//   three parameters:
//   - explosion strength & explosion amplitude
//     fireworks with better fitness produce more sparks with smaller amplitude (like short runners)
//     fireworks with less fitness produce less sparks but larger amplitude and are used to explore the feasible space (like long-runners)
//   - displacement operation
//     ensure diversification
// - gaussian mutation operator
//   generates new, mutated sparks between best and selected fireworks
// - mapping rule
//   how do we get sparks outside of the search space back into it?
// - selection strategy
//   how do we select best sparks for further generations?


// Algorithm:
// 1. randomly generat fireworks in the feasible space
// 2. calculate fitness value of each firework according to the fitness function. The number of sparks is calculated based on immune concentration
//    theory and the fireworks with better fitness produce more sparks
// 3. considering fireworks phenomena in the real world and the landscape of the functions, the fireworks generate sparks within a certain
//    amplitude. The explosion amplitude is determined by the fitness value of a firework. The explosion amplitude is smaller for higher fitness
//    and vice versa. Each spark represents a solution in the feasible space. To keep the diversity of the population, mutation operation is 
//    needed and Gaussian mutation is one of them.
// 4. Calculate the best fitness value. If the terminal condition is met, stop the algorithm. Otherwise, continue with the iteration process. The
//    best sparks and the selected sparks form a new population.

#define EPS 1e-38

void boundCheck(double *x, const CoordBound &b)
{
	const double LBOUND = get<0>(b);
	const double UBOUND = get<1>(b);
	if (*x < LBOUND || *x > UBOUND)
	{
		double absoluteValue = abs(*x);
		while (absoluteValue > 0)
			absoluteValue -= (UBOUND - LBOUND);
		*x = absoluteValue + UBOUND;
	}
}


Member generateGaussianSpark(const Member& xi, const vector<CoordBound> &bounds)
{
	int fireworksDirections = rand() % bounds.size() + 1;

	Member xj(xi);

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

	double g = sampleGaussian_1_1();

	for (size_t k = 0; k < bounds.size(); ++k) {
		if (flags[k]) {
			double xjk = xj[k];

			xjk *= g;
			boundCheck(&xjk, bounds[k]);

			xj[k] = xjk;
		}
	}

	return xj;
}


Population generateGaussianSparks(const Population& oldPop, int mMuts, const vector<CoordBound> &bounds)
{
	Population newSparks;
	newSparks.reserve(mMuts);

	// in C code rand() % N
	auto dimIndices = randomIndices(oldPop.size());
	// calculate gaussian sparks
	for (int i = 0; i < mMuts; ++i) {
		const auto& selectedSpark = oldPop[dimIndices[i]];

		newSparks.push_back(generateGaussianSpark(selectedSpark, bounds));
	}

	return newSparks;
}

Member generateNewSpark(const Member& xi, const double A, const vector<CoordBound> &bounds)
{
	int fireworksDirections = rand() % bounds.size() + 1;

	Member xj(xi);

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

	double h = A * fRand(-1, 1);

	for (size_t k = 0; k < bounds.size(); ++k) {
		if (flags[k]) {
			double xjk = xj[k];

			xjk += h;
			boundCheck(&xjk, bounds[k]);

			xj[k] = xjk;
		}
	}	

	return xj;
}

Population generateNewSparksFromExplosion(const Member& parentSpark, const size_t s, const double A, const vector<CoordBound> &bounds)
{
	Population newSparks;

	for (size_t i = 0; i < s; ++i) {
		newSparks.push_back(generateNewSpark(parentSpark, A, bounds));
	}

	return newSparks;
}

const Population calculateNewSparks(
	const vector<MemberWithValue> &objectiveValues, 
	const size_t m, 
	const double a, 
	const double b,
	const double Amax,
	const vector<CoordBound> &bounds
)
{
	Population newSparks;

	const auto minMax = minmax_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

	const double ymin = get<1>(*minMax.first);
	const double ymax = get<1>(*minMax.second);
	
	const double am = a * m;
	const double bm = b * m;

	// reduce objective values to get sum of all
	double genSparkSum = 0.0;
	double amplitudeSum = 0.0;

	for (const auto &it : objectiveValues) {
		const double fval = get<1>(it);
		genSparkSum += (ymax - fval);
		amplitudeSum += (fval - ymin);
	}
	vector<int> numFireworks;
	vector<double> amplitudes;

	for (const auto &it : objectiveValues) {
		// Eq. 2, calculate relative fitness
		const double fval = get<1>(it);
		//cout << printMember(member) << endl;
		const double s = m * ((ymax - fval + EPS) / (genSparkSum + EPS));

		// Eq. 3, ensure upper and lower bounds of sparks
		int sr;
		// minimum sparks
		if (s < am) {
			sr = (int)round(am);
		}
		// maximum sparks
		else if (s > bm) {
			sr = (int)round(bm);
		}
		else {
			sr = (int)round(s);
		}
		numFireworks.push_back(sr);
		//sr = m;
		// Eq. 4 amplitude
		double A = Amax * ((fval - ymin + EPS) / (amplitudeSum + EPS));
		if (isnan(A)) {
			cout << "fval: " << fval << ", yim: " << ymin << ", amplitudeSum: " << amplitudeSum << endl;
			cout << "A ERROR: " << A << endl;
			exit(1);
		}

		amplitudes.push_back(A);
	}

	for (size_t i = 0; i < objectiveValues.size(); ++i) {
		const auto &member = get<0>(objectiveValues[i]);
		const int sr = numFireworks[i];
		const double A = amplitudes[i];
		auto generatedSparks = generateNewSparksFromExplosion(member, static_cast<size_t>(sr), A, bounds);
		for (const auto& s : generatedSparks) {
			newSparks.push_back(s);
		}
	}

	return newSparks;
}

double computeRxi(const Member &xi, const vector<MemberWithValue>& K, function<double(const Member&, const Member&)> distanceFunc)
{
	double currentSum = 0.0;

	// calculate sum of distance for current spark
	for (const auto &xj : K) {
		currentSum += distanceFunc(xi, get<0>(xj));
	}
	return currentSum;
}

Population selectSparksForNextGeneration(vector<MemberWithValue>& K, const size_t n, function<double(const Member&, const Member&)> distanceFunc)
{
	vector<MemberWithValue> probs;
	probs.reserve(static_cast<size_t>(pow(K.size(), 2)));
	// for every spark, calculate distance to every other spark
	double R_total = 0.0;
	for (const auto &xi : K) {
		R_total += computeRxi(get<0>(xi), K, distanceFunc);
	}

	for (const auto &xi : K) {
		const double R_xi = computeRxi(get<0>(xi), K, distanceFunc);

		const double P_xi = R_xi / R_total;

		probs.push_back(MemberWithValue(get<0>(xi), P_xi));
	}

	Population newPopulation;
	newPopulation.reserve(n);
	newPopulation.push_back(get<0>(K[0]));

	for (size_t i = 1; i < n; i++)
	{
		double possibilityRand = (double)rand() / RAND_MAX;
		double currentpossibility = 0.0;
		size_t possibilityIndex;
		for (size_t j = 0; j < n; j++)
		{
			currentpossibility += get<1>(probs[j]);
			if (currentpossibility >= possibilityRand)
			{
				possibilityIndex = j;
				break;
			}
			//Be careful for the situation when 0.9999 >= 1.0000
			if (j == n - 1) possibilityIndex = j;
		}

		newPopulation.push_back(get<0>(K[possibilityIndex]));
	}

	return newPopulation;
}

Population runFWA(Parameters *ps, ValueCollector &vc)
{
	const FWA& fwa = castParameters<FWA>(ps);

	auto P = createRandomPopulation(fwa.initialSize, fwa.initBounds);

	int evals = 0;
	size_t g = 0;
	for (; g < fwa.maxGenerations && (size_t)evals < fwa.maxFunctionEvaluations; ++g) {
		// set of n fireworks at n locations
		auto objectiveValues = evalObjectiveFunctionForPopulation(P, fwa.objectiveFunction, &evals);
		if (g == fwa.maxGenerations - 1) {
			// lower is better
			stable_sort(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

			P = Population();
			for (const auto &l : objectiveValues) {
				P.push_back(move(get<0>(l)));
			}

			vc.bestMembersInGeneration.push_back(P[0]);
		} else {
			// for each firework x_i do:
			//     calculate number of sparks s_i that firework x_i yields (Eq. 3)
			//     obtain locations of s_i sparks of firework x_i (Alg. 1)
			auto newSparks = calculateNewSparks(objectiveValues, fwa.maxSparks, fwa.a, fwa.b, fwa.Amax, fwa.coordinateBounds);

			// for k = 1:mMuts do:
			//		randomly select a firework x_i
			//      generate a specific gaussian spark for the firework (Alg. 2)
			auto gaussianSparks = generateGaussianSparks(P, fwa.gaussianMutations, fwa.coordinateBounds);

			
			auto objectiveValuesNewSparks = evalObjectiveFunctionForPopulation(newSparks, fwa.objectiveFunction, &evals);
			auto objectiveValuesGaussianSparks = evalObjectiveFunctionForPopulation(gaussianSparks, fwa.objectiveFunction, &evals);

			vector<MemberWithValue> allObjectiveValues(objectiveValues);// (objectiveValues);
			copy(objectiveValuesNewSparks.begin(), objectiveValuesNewSparks.end(), back_inserter(allObjectiveValues));
			copy(objectiveValuesGaussianSparks.begin(), objectiveValuesGaussianSparks.end(), back_inserter(allObjectiveValues));

			// select the best location and keep it for next explosion generation
			auto bestMemberIt = min_element(allObjectiveValues.begin(), allObjectiveValues.end(), compareMemberWithValueLower);
			auto bestLocation = get<0>(*bestMemberIt);
			auto bestValue = get<1>(*bestMemberIt);
			allObjectiveValues.erase(bestMemberIt);
			allObjectiveValues.insert(allObjectiveValues.begin(), MemberWithValue(bestLocation, bestValue));

			vc.bestMembersInGeneration.push_back(bestLocation);

			// initial size
			P = selectSparksForNextGeneration(allObjectiveValues, fwa.initialSize, fwa.distanceFunction);
		}
	}

	vc.numberFunctionEvaluations = evals;
	vc.numberGenerations = g;

	return P;
}