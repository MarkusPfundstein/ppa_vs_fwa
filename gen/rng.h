#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <random>

using namespace std;

double fRand(double min, double max);
int iRand(int min, int max);
vector<int> randomIndices(const size_t size);

double sampleGaussian_1_1();

#endif