#pragma once

#include "population.h"


// http://www.geatbx.com/docu/fcnindex-01.html#P150_6749

/* objective function */
double rosenbrock(const Member& member);
double griewank(const Member& member);
double schwefel1_2(const Member& member);
double schwefel7(const Member& member);
double easom(const Member& member);
double ackleys_path(const Member& member);
double michalewicz12(const Member& x);

/* fitness normalizers */
double minMaxFitness(double v, double min, double max);

/* distance functions */
double euclideanDistance(const Member&, const Member&);


double normalizeTan(double f);

