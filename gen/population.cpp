#include <sstream>
#include "population.h"
#include "common.h"

using namespace std;


bool compareMemberWithValueLower(const MemberWithValue& v, const MemberWithValue& v2) {
	return get<1>(v) < get<1>(v2);
}

bool compareMember(const Member & m1, const Member& m2)
{
	// http://stackoverflow.com/a/17341
	double eps = std::numeric_limits<double>::epsilon();
	for (size_t i = 0; i < m1.size(); ++i) {
		if (::abs(m1[i] - m2[i]) >= eps) {
			return false;
		}
	}
	return true;
}

Member createRandomInstance(const vector<CoordBound> &bounds)
{
	Member data;
	data.reserve(bounds.size());

	for (auto bound : bounds) {

		const double min = get<0>(bound);
		const double max = get<1>(bound);

		data.push_back(fRand(min, max));
	}

	return data;
}

Population createRandomPopulation(size_t np, const vector<CoordBound> &bounds) {
	Population pop;

	for (size_t i = 0; i < np; ++i) {
		pop.push_back(createRandomInstance(bounds));
	}

	return pop;
}

std::string printBounds(const std::vector<CoordBound> &bounds)
{
	stringstream ss;
	for (size_t i = 0; i < bounds.size(); ++i) {
		const auto &b = bounds[i];
		ss << "(" << get<0>(b) << ", " << get<1>(b) << ")";
		if (i < bounds.size() - 1) {
			ss << ", ";
		}
	}
	return ss.str();
}

std::string printMember(const Member& m)
{
	stringstream ss;
	ss << "(";
	for (size_t i = 0; i < m.size(); ++i) {
		ss << m[i];
		if (i < m.size() - 1) {
			ss << ",";
		}
	}
	ss << ")";
	return ss.str();
}

std::string printMemberWithValue(const MemberWithValue& m)
{
	return "";
}