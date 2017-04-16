#include <algorithm>

#include <iostream>

#include "fwa.h"
#include "population.h"

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



const int calculateNewSparks(
	const vector<MemberWithValue> &objectiveValues, 
	const size_t m, 
	const double a, 
	const double b,
	const double Amax
)
{
	const auto minMax = minmax_element(objectiveValues.begin(), objectiveValues.end(), compareMemberWithValueLower);

	const double ymin = get<1>(*minMax.first);
	const double ymax = get<1>(*minMax.second);
	
	const double eps = numeric_limits<double>::epsilon();

	// reduce objective values to get sum of all
	double total_yMaxMinusF = 0.0;
	double total_fMinus_yMin = 0.0;
	for (const auto &it : objectiveValues) {
		const double current = get<1>(it);
		total_yMaxMinusF += (ymax - current);
		total_fMinus_yMin += (current - ymin);
	}

	const double am = a * m;
	const double bm = b * m;
	cout << "am: " << am << endl;
	cout << "bm: " << bm << endl;
	for (const auto &it : objectiveValues) {
		// Eq. 2, calculate relative fitness
		const double current = get<1>(it);
		const double s = m * (((ymax - current) + eps) / (total_yMaxMinusF + eps));

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
		const double A = Amax * (((current - ymin) + eps) / (total_fMinus_yMin + eps));

		cout << "--------------" << endl;
		cout << "member: " << printMember(get<0>(it)) << endl;
		cout << "\tfval: " << current << endl;
		cout << "\t\ts: " << s << endl;
		cout << "\t\tsr: " << sr << endl;
		cout << "\t\tA: " << A << endl;
	}

	return 0;
}

Population fwa(
	const vector<CoordBound> coordinateBounds,
	const size_t n,
	const size_t gMax,
	function<double(double, double, double)> fitnessFunction,
	function<double(const Member&)> objectiveFunction
)
{
	auto P = createRandomPopulation(n, coordinateBounds);

	// temporary set one member to opt
	P[0] = Member({ 1, 1 });

	// n = 5, m = 50, a = 0.04, b = 0.8, Amax = 50, mMuts = 5
	// n = initial number of fireworks
	// m = total number of sparks         (like nMax in PPA?)
	// Amax = maximum explosion amplitude
	// mMuts = number of mutation sparks
	// a, b = some constants
	
	const double Amax = 40.0;
	const size_t m = 50;       
	const size_t mMuts = 5;
	const double a = 0.04;
	const double b = 0.8;

	for (size_t g = 0; g < gMax; ++g) {

		cout << printPopulation(P);

		auto objectiveValues = evalObjectiveFunctionForPopulation(P, objectiveFunction);

		calculateNewSparks(objectiveValues, m, a, b, Amax);

	}

	return P;
}