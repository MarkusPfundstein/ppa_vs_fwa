#ifndef COMMON_H
#define COMMON_H

double fRand(double min, double max);
int iRand(int, int);

template<class BidiIter>
BidiIter fisherYatesShuffle(BidiIter begin, BidiIter end, size_t num_random) {
	size_t left = distance(begin, end);
	while (num_random--) {
		BidiIter r = begin;
		advance(r, iRand(0, RAND_MAX) % left);
		swap(*begin, *r);
		++begin;
		--left;
	}
	return begin;
}

#endif