#ifndef POPULATION_H
#define POPULATION_H


#include <tuple>
#include <vector>
#include <string>
#include <functional>

using namespace std;

typedef double COORDBOUND_TYPE;

// Bounds for each dimension (in paper called a, b)
typedef tuple<COORDBOUND_TYPE, COORDBOUND_TYPE> CoordBound;

// Member of population -> Coordinate
typedef vector<double> Member;
typedef vector<Member> Population;

// tuple that combines a Coordinate with a double -> e.g. fitness value
typedef tuple<Member, double> MemberWithValue;

bool compareMemberWithValueLower(const MemberWithValue& v, const MemberWithValue& v2);
bool compareMemberWithValueSameValue(const MemberWithValue& v, const MemberWithValue& v2);
bool compareMemberEquals(const Member & m1, const Member& m2);

vector<MemberWithValue> evalObjectiveFunctionForPopulation(const Population &population, function<double(const Member&)> f, int *evals);

Population createRandomPopulation(size_t np, const vector<CoordBound> &bounds);

vector<CoordBound> createUniformCoordinateBounds(size_t n, COORDBOUND_TYPE min, COORDBOUND_TYPE max);

string printBound(const CoordBound& bound);
string printBounds(const vector<CoordBound> &bounds);
string printMember(const Member& m);
string printPopulation(const Population& population);

#endif
