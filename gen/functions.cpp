#include "functions.h"

// https://en.wikipedia.org/wiki/Rosenbrock_function
double rosenbrock2d(const Member& member)
{
	double a = 1;
	double b = 100.0;

	double x = member[0];
	double y = member[1];
	return pow((a - x), 2) + b * pow((y - pow(x, 2)), 2);
};

double griewank(const Member& member) {
	// 1/4000 * SUM_i(x_i^2) - PROD_i(cos(x/sqrt(i)) + 1

	double sum = 0.0;
	for (double xi : member) {
		sum += ( pow(xi, 2) / 4000.0 );
	}
	double prod = 1.0;

	int i = 1;
	for (double xi : member) {
		prod *= cos(xi / sqrt(i));
		++i;
	}

	return sum - prod + 1.0;
}

// http://www.sfu.ca/~ssurjano/schwef.html
double schwefel2d(const Member& member) {
	double d = member.size();

	double x = member[0];
	double y = member[1];
	return 418.9829 * d - ((x * sin(sqrt(abs(x)))) + (y * sin(sqrt(abs(y)))));
};

double minMaxFitness(double v, double min, double max)
{
	return (max - v) / (max - min);
}

double euclideanDistance(const Member& m1, const Member& m2)
{
	double sum = 0.0;
	for (size_t i = 0; i < m1.size(); ++i) {
		sum += pow(m1[i] - m2[i], 2);
	}

	return sqrt(sum);
};