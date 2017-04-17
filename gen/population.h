#ifndef POPULATION_H
#define POPULATION_H


#include <tuple>
#include <vector>
#include <string>
#include <functional>

using namespace std;

// Bounds for each dimension (in paper called a, b)
typedef std::tuple<double, double> CoordBound;

// Member of population -> Coordinate
typedef std::vector<double> Member;
typedef std::vector<Member> Population;

// tuple that combines a Coordinate with a double -> e.g. fitness value
typedef std::tuple<Member, double> MemberWithValue;

bool compareMemberWithValueLower(const MemberWithValue& v, const MemberWithValue& v2);
bool compareMemberWithValueSameValue(const MemberWithValue& v, const MemberWithValue& v2);
bool compareMemberEquals(const Member & m1, const Member& m2);

vector<MemberWithValue> evalObjectiveFunctionForPopulation(const Population &population, function<double(const Member&)> f);

Population createRandomPopulation(size_t np, const std::vector<CoordBound> &bounds);

std::string printBounds(const std::vector<CoordBound> &bounds);
std::string printMember(const Member& m);
std::string printPopulation(const Population& population);

#endif
