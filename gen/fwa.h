#ifndef FIREWORKS_H
#define FIREWORKS_H

#include "population.h"

Population fwa(
	const vector<CoordBound> coordinateBounds,
	const size_t m,
	const size_t gMax,
	function<double(double, double, double)> fitnessFunction,
	function<double(const Member&)> objectiveFunction
);

#endif
