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

Member generateGaussianSpark(const Member& xi, const vector<CoordBound> &bounds)
{
	// alg. 1
	// initialize the location of the spark (xj = xi)
	Member xj(xi);

	auto dimIndices = randomIndices(xj.size());

	// calculate the displacement g = Gaussian(1, 1);
	const double g = sampleGaussian_1_1();

	// for each dimension x^j_k \in { pre-selected z dimensions of xj } do:
	const size_t z = static_cast<size_t>(round(xj.size() * fRand(0, 1)));
	for (size_t i = 0; i < z; ++i) {
		const size_t k = dimIndices[i];
		double &xjk = xj[k];

		xjk *= g;

		// out of bounds? get bound for dimension and check
		const auto& bound = bounds[k];
		const double xmink = get<0>(bound);
		const double xmaxk = get<1>(bound);

		if (xjk < xmink || xjk > xmaxk) {
			// bounds should actually be ints I guess
			xjk = xmink + static_cast<int>(abs(xjk)) % static_cast<int>((xmaxk - xmink));
		}
	}

	return xj;
}

Population generateGaussianSparks(const Population& oldPop, int mMuts, const vector<CoordBound> &bounds)
{
	Population newSparks;
	newSparks.reserve(mMuts);

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
	// alg. 1
	// initialize the location of the spark (xj = xi)
	Member xj(xi);

	auto dimIndices = randomIndices(xj.size());

	// calculate the displacement h = A_i * rand(0, 1)
	const double h = A * fRand(-1, 1);

	// for each dimension x^j_k \in { pre-selected z dimensions of xj } do:
	const size_t z = static_cast<size_t>(round(xj.size() * fRand(0, 1)));
	for (size_t i = 0; i < z; ++i) {
		const size_t k = dimIndices[i];
		double &xjk = xj[k];

		xjk += h;

		// out of bounds? get bound for dimension and check
		const auto& bound = bounds[k];
		const double xmink = get<0>(bound);
		const double xmaxk = get<1>(bound);

		if (xjk < xmink || xjk > xmaxk) {
			// bounds should actually be ints I guess
			xjk = xmink + static_cast<int>(abs(xjk)) % static_cast<int>((xmaxk - xmink));
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
	
	const double eps = numeric_limits<double>::epsilon();

	const double am = a * m;
	const double bm = b * m;

	// reduce objective values to get sum of all
	double total_yMaxMinusF = 0.0;
	double total_fMinus_yMin = 0.0;

	for (const auto &it : objectiveValues) {
		const double fval = get<1>(it);
		total_yMaxMinusF += (ymax - fval);
		total_fMinus_yMin += (fval - ymin);
	}

	for (const auto &it : objectiveValues) {
		// Eq. 2, calculate relative fitness
		const double fval = get<1>(it);
		const auto& member = get<0>(it);
		const double s = m * (((ymax - fval) + eps) / (total_yMaxMinusF + eps));

		// Eq. 3, ensure upper and lower bounds of sparks
		double sr; 
		// minimum sparks
		if (s < am) {
			sr = round(am);
		}
		// maximum sparks
		else if (s > bm) {
			sr = round(bm);
		}
		else {
			sr = round(s);
		}

		// Eq. 4 amplitude
		const double A = Amax * (((fval - ymin) + eps) / (total_fMinus_yMin + eps));

		auto generatedSparks = generateNewSparksFromExplosion(member, static_cast<size_t>(sr), A, bounds);
		for (const auto& s : generatedSparks) {
			newSparks.push_back(move(s));
		}
	}

	return newSparks;
}

double computeRxi(const Member &xi, const Population& K, const Member& bestLocation, function<double(const Member&, const Member&)> distanceFunc)
{
	double currentSum = 0.0;

	// calculate sum of distance for current spark
	for (const auto &xj : K) {
		if (compareMemberEquals(xj, bestLocation)) {
			continue;
		}
		currentSum += distanceFunc(xi, xj);
	}
	return currentSum;
}

Population selectSparksForNextGeneration(const Population& K, const Member& bestLocation, size_t n, function<double(const Member&, const Member&)> distanceFunc)
{
	vector<MemberWithValue> vecs;
	vecs.reserve(static_cast<size_t>(pow(K.size(), 2)));
	// for every spark, calculate distance to every other spark
	double R_total = 0.0;
	for (const auto &xi : K) {
		if (compareMemberEquals(xi, bestLocation)) {
			continue;
		}

		R_total += computeRxi(xi, K, bestLocation, distanceFunc);
	}

	for (const auto &xi : K) {
		if (compareMemberEquals(xi, bestLocation)) {
			continue;
		}

		const double R_xi = computeRxi(xi, K, bestLocation, distanceFunc);

		const double P_xi = R_xi / R_total;

		vecs.push_back(MemberWithValue(Member(xi), P_xi));
	}

	// sort
	stable_sort(
		vecs.begin(),
		vecs.end(),
		compareMemberWithValueLower
	);
	
	// vecs.erase(unique(vecs.begin(), vecs.end(), compareMemberWithValueSameValue), vecs.end());
	
	// pick (n - 1) values
	Population newPopulation;
	newPopulation.reserve(n - 1);
	while (n-- > 1) {
		
		const double r = fRand(0.0, 1.0);
		//cout << "r: " << r << endl;

		double s = 0.0;
		for (size_t i = 0; i < vecs.size(); ++i) {
			if (i == vecs.size() - 1) {
				newPopulation.push_back(get<0>(vecs[i]));
			}
			else {
				const auto& xi = vecs[i];
				const auto& xi2 = vecs[i + 1];

				s += get<1>(xi);
				if (s <= r && r < s + get<1>(xi2)) {
					const auto& n = get<0>(xi2);
					newPopulation.push_back(n);
					break;
				}
			}
		}
	}
	
	return newPopulation;
}

Population runFWA(Parameters *ps, ValueCollector &vc)
{
	const FWA& fwa = castParameters<FWA>(ps);

	auto P = createRandomPopulation(fwa.initialSize, fwa.coordinateBounds);

	for (size_t g = 0; g < fwa.maxGenerations; ++g) {
		// set of n fireworks at n locations
		auto objectiveValues = evalObjectiveFunctionForPopulation(P, fwa.objectiveFunction);
		if (g == fwa.maxGenerations - 1) {
			// lower is better
			stable_sort(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

			P = Population();
			for (const auto &l : objectiveValues) {
				P.push_back(move(get<0>(l)));
			}

			vc.bestMembersInGeneration.push_back(P[0]);

			return P;
		}
		else {
			// for each firework x_i do:
			//     calculate number of sparks s_i that firework x_i yields (Eq. 3)
			//     obtain locations of s_i sparks of firework x_i (Alg. 1)
			auto newSparks = calculateNewSparks(objectiveValues, fwa.maxSparks, fwa.a, fwa.b, fwa.Amax, fwa.coordinateBounds);

			// for k = 1:mMuts do:
			//		randomly select a firework x_i
			//      generate a specific gaussian spark for the firework (Alg. 2)
			auto gaussianSparks = generateGaussianSparks(P, fwa.gaussianMutations, fwa.coordinateBounds);

			// select the best location and keep it for next explosion generation
			auto bestLocation = get<0>(*min_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower));

			vc.bestMembersInGeneration.push_back(bestLocation);

			// randomly select (n - 1) locations from the two types of sparks and the current fireworks
			Population allSparks(newSparks);
			copy(gaussianSparks.begin(), gaussianSparks.end(), back_inserter(allSparks));
			copy(P.begin(), P.end(), back_inserter(allSparks));

			// remove duplicates (urks) and bestLocation
			Population setSparks;
			setSparks.reserve(allSparks.size());
			for (const auto &s1 : allSparks) {
				if (compareMemberEquals(s1, bestLocation)) continue;

				bool isInIt = false;
				for (const auto &s2 : setSparks) {
					if (compareMemberEquals(s1, s2)) {
						isInIt = true;
						break;
					}
				}
				if (!isInIt) {
					setSparks.push_back(s1);
				}
			}

			P = selectSparksForNextGeneration(setSparks, bestLocation, fwa.initialSize, fwa.distanceFunction);
			P.push_back(bestLocation);
		}
	}

	return P;
}