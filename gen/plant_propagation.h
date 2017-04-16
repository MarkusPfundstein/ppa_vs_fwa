#ifndef PLANT_PROPAGATION_H
#define PLANT_PROPAGATION_H

#include "population.h"

Population runPlantPropagation1(
	const vector<CoordBound> coordinateBounds,
	const size_t m,
	const size_t gMax,
	const size_t nMax,
	function<double(double, double, double)> fitnessFunction,
	function<double(const Member&)> objectiveFunction
);

#endif