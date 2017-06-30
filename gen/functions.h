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

double sphere(const Member &X);
double rastrigrin(const Member &X);
double ellipse(const Member &X);
double cigar(const Member& X);
double tablet(const Member& X);
double schwefelFWA(const Member& X);

/* fitness normalizers */
double minMaxFitness(double v, double min, double max);

/* distance functions */
double euclideanDistance(const Member&, const Member&);


double normalizeTan(double f);



void setShift(double);
double getShift();
