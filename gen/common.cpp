#include <chrono>
#include <functional>
#include <random>

#include "common.h"

using namespace std;

// global seed
// TO-DO: Make possible to fix
const auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
// our random value generator
auto randomDouble = bind(uniform_real_distribution<double>(0, RAND_MAX), mt19937(static_cast<unsigned int>(seed)));

double fRand(double fMin, double fMax)
{
	auto f = randomDouble() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

