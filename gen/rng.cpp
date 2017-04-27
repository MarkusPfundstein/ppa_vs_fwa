#include <chrono>
#include <functional>
#include <random>
#include <algorithm>

#include "rng.h"

using namespace std;
const double PI = 3.141592653589793238462643383279502884;

#define NO_RAND 0
// global seed
// TO-DO: Make possible to fix
#if NO_RAND
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
	double f = randomDouble() / RAND_MAX;
	return min + f * (max - min);
}

int iRand(int min, int max)
{
	int f = randomInt();
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
#else
double sampleGaussian_1_1()
{
	double mu = 1.0;
	double sigma = 1.0;
	double x1 = 0.0;
	double x2, y;

	//Warning without outlet
	while (x1 == 0.0)
	{
		x1 = (double)rand() / RAND_MAX;
	}

	x2 = (double)rand() / RAND_MAX;

	//double log (double); 以e为底的对数
	y = sqrt(-2 * log(x1)) * cos(2 * PI * x2);
	double temp = mu + y * sigma;

	//_finite returns 0 if the argument is infinite or a NAN
	if (!_finite(temp))
	{
		//printf("X1 means to zero in function randGaussian.\n");
		//system("pause");
		exit(0);
	}
	return mu + y * sigma;

}

double fRand(double min, double max)
{
	return min + (double)rand() / RAND_MAX * (max - min);
}

int iRand(int min, int max)
{
	return min + (double)rand() / RAND_MAX * (max - min);
}

// select randomly z dimensions -> generate array of indices.. shuffle them
vector<int> randomIndices(const size_t size)
{
	vector<int> dimIndices;
	dimIndices.reserve(size);
	for (size_t i = 0; i < size; ++i) {
		dimIndices.push_back(i);
	}
	random_shuffle(dimIndices.begin(), dimIndices.end(), [](int i) { return rand() % i; });

	return dimIndices;
}

#endif
