#pragma once

#include "population.h"

/* objective function */
double rosenbrock2d(const Member& member);
double schwefel2d(const Member& member);

/* fitness normalizers */
double minMaxFitness(double v, double min, double max);

/* distance functions */
double euclideanDistance(const Member&, const Member&);