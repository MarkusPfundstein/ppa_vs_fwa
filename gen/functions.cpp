#include "functions.h"

const double PI = 3.141592653589793238462643383279502884;

double SHIFT = 0;

void setShift(double newShift)
{
	SHIFT = newShift;
}

double getShift() {
	return SHIFT;
}

template <typename RT, typename VT>
RT foldLeft(const vector<VT> &xs, VT s, function<RT(RT, VT)> f)
{
	RT sx = s;
	for (const auto& x : xs) {
		sx = f(sx, x);
	}
	return sx;
}

template <typename T>
T sum(const vector<T> &xs, function<T(T)> f)
{
	return foldLeft<T, T>(xs, 0, [&](auto acc, auto x) { return acc + f(x); });
}

template <typename T>
T prod(const vector<T> &xs, function<T(T)> f)
{
	return foldLeft<T, T>(xs, 1, [&](auto acc, auto x) { return acc * f(x); });
}

/* FUNCTIONS FOR FINAL TESTS */

double rosenbrock(const Member& member) {
	double sum = 0.0;

	// shift value -> new origin (o, o, o)
	// unshifted o = 1
	const double o = SHIFT;

	const size_t n = member.size();
	for (size_t i = 0; i < (n - 1); ++i) {
		const double x_1 = member[i + 1] - o + 1;
		const double x = member[i] - o + 1;

		const double inner1 = x_1 - pow(x, 2);
		const double part1 = 100.0 * pow(inner1, 2);
		const double part2 = pow(x - 1.0, 2);

		sum += (part1 + part2);
	}

	return sum;
}

double griewank(const Member& member) {
	// 1/4000 * SUM_i(x_i^2) - PROD_i(cos(x/sqrt(i)) + 1

	const double S = sum<double>(member, [](auto x) { return (pow(x + SHIFT, 2) / 4000.0); });

	int i = 1;
	double prod = 1.0;
	for (double xi : member) {
		prod *= cos((xi + SHIFT) / sqrt(i));
		++i;
	}

	return S - prod + 1.0;
}

double ackleys_path(const Member& member)
{
	const double a = 20;
	const double b = 0.2;
	const double c = 2 * PI;
	const double n = member.size();

	const double o = SHIFT;

	const double s1 = sum<double>(member, [o](auto x) { return pow(x + o, 2); });
	const double s2 = sum<double>(member, [o, c](auto x) { return cos(c * (x + o)); });

	const double e1 = exp(-b * sqrt(s1 / n));
	const double e2 = exp(s2 / n);
	
	return -a * e1 - e2 + a + exp(1);
}

double sphere(const Member &X)
{
	return sum<double>(X, [](auto x) { return (x + SHIFT) * (x + SHIFT); });
}

double rastrigrin(const Member &X)
{
	return sum<double>(X, [](auto x) {
		return (x + SHIFT) * (x + SHIFT) - 10.0 * cos(2 * PI * (x + SHIFT)) + 10;
	});
}

double ellipse(const Member &X)
{
	const double D = X.size();
	
	double s = 0.0;

	for (int i = 1; i <= D; ++i) {
		const double x = X[i - 1] + SHIFT;
		const double exponent = 4.0 * (i - 1.0) / (D - 1.0);

		s += pow(10, exponent) * x * x;
	}

	return s;
}

double cigar(const Member& X)
{
	double x0 = X[0] + SHIFT;

	double s = x0 * x0;

	for (int i = 1; i < X.size(); ++i) {
		s += pow(10, 4) * (X[i] + SHIFT) * (X[i] + SHIFT);
	}

	return s;
}

double tablet(const Member& X)
{
	double x0 = X[0] + SHIFT;

	double s = pow(10, 4) * x0 * x0;

	for (int i = 1; i < X.size(); ++i) {
		s += (X[i] + SHIFT) * (X[i] + SHIFT);
	}

	return s;
}

double schwefelFWA(const Member& X)
{
	double x0 = X[0] + SHIFT;

	return sum<double>(X, [x0](auto xx) {
		const double x = xx + SHIFT;

		return pow(x0 - x * x, 2) + pow(x - 1.0, 2);
	});
}

/* END OF FUNCTIONS FINAL TEST */

double schwefel1_2(const Member& x)
{
	double f = 0;
	double a;
	for (size_t i = 0; i < x.size(); i++)
	{
		a = 0.0;
		for (size_t j = 0; j <= i; j++) {
			a += x[j] + SHIFT;
		}
		f += a * a;
	}
	return f;
}

double schwefel7(const Member& member)
{
	const double c = 418.98288727243369;

	const double s = sum<double>(member, [](auto x) { return (x + SHIFT) * sin(sqrt(abs(x + SHIFT)));  });

	return c * member.size() - s;
}

double easom(const Member& member)
{
	const double x1 = member[0] + SHIFT;
	const double x2 = member[1] + SHIFT;

	return -cos(x1) * cos(x2) * exp(-(pow(x1 - PI, 2) + pow(x2 - PI, 2)));
}

double michalewicz12(const Member& y)
{
	const double m = 10.0;
	double s = 0.0;
	for (int i = 0; i < y.size(); ++i) {
		const double x = y[i] + SHIFT;
		s += sin(x) * pow(sin(i * x * x / PI), 2.0 * m);
	}

	return -s;
}



/* others */

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


double normalizeTan(double f) {
	return 0.5 * (tanh(4.0 * f - 2.0) + 1.0);
	//return tanh(f - 0.5) + 0.5;
}