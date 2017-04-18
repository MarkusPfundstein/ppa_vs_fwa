#include <chrono>
#include <functional>
#include <random>
#include <algorithm>

#include "rng.h"

using namespace std;

// global seed
// TO-DO: Make possible to fix
const auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
// our random value generator
mt19937 mt = mt19937(static_cast<unsigned int>(seed));

auto randomDouble = bind(uniform_real_distribution<double>(0, RAND_MAX), mt);

auto randomInt = bind(uniform_int_distribution<int>(0, RAND_MAX), mt);

double sampleGaussian_1_1()
{
	static auto gaussianDist = normal_distribution<>(1.0, 1.0);


	return gaussianDist(mt);
}

double fRand(double min, double max)
{
	auto f = randomDouble() / RAND_MAX;
	return min + f * (max - min);
}

int iRand(int min, int max)
{
	auto f = randomInt();
	return f % (max - min + 1) + min;
}

// select randomly z dimensions -> generate array of indices.. shuffle them
vector<int> randomIndices(const size_t size)
{
	vector<int> dimIndices;
	dimIndices.reserve(size);
	for (size_t i = 0; i < size; ++i) {
		dimIndices.push_back(i);
	}
	random_shuffle(dimIndices.begin(), dimIndices.end(), [](int i) { return randomInt() % i; });

	return dimIndices;
}

