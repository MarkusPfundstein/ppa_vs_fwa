#ifndef PLANT_PROPAGATION_H
#define PLANT_PROPAGATION_H

#include "population.h"
#include <functional>

Population runPlantPropagation1(
	const std::vector<CoordBound> coordinateBounds,
	const size_t m,
	const size_t gMax,
	const size_t nMax,
	std::function<double(double, double, double)> fitnessFunction,
	std::function<double(const Member&)> objectiveFunction
);

#endif