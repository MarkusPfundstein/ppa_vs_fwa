#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <random>

using namespace std;

double fRand(double min, double max);
int iRand(int min, int max);
vector<int> randomIndices(const size_t size);

double sampleGaussian(double mu, double sigma);
double sampleGaussian_1_1();

double poisson(double Lambda, unsigned k, double t = 1.0);

#endif