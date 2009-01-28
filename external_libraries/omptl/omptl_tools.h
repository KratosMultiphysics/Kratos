// Copyright (C) 2006 Fokko Beekhof
// Email contact: Fokko.Beekhof@cui.unige.ch

// The OMPTL library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef OMPTL_TOOLS_H
#define OMPTL_TOOLS_H

#include <utility>
#include <vector>
#include <cassert>
#include <algorithm>

namespace omptl
{

// Log of the number of operations that is expected to run faster in a single
// thread.
const unsigned int C = 8;

template<typename Iterator>
bool _linear_serial_is_faster(Iterator first, Iterator last,
			     const unsigned int P)
{
	const unsigned int N = ::std::distance(first, last);
	// Approximation: (log2(N) - 1) <= l < =  log2(N)
	unsigned int l = 1 << (8*sizeof(unsigned int) - 1);
	while(l & N)
		l /= 2;

	return (N < P) && (l < C);
}

template<typename Iterator>
bool _logn_serial_is_faster(Iterator first, Iterator last,
			    const unsigned int P)
{
	const unsigned int N = ::std::distance(first, last);
	// Approximation: (log2(N) - 1) <= l < =  log2(N)
	unsigned int l = 1 << (8*sizeof(unsigned int) - 1);
	while(l & N)
		l /= 2;

	return (N < P) && (l < (1 << C));
}

template<typename Iterator>
bool _nlogn_serial_is_faster(Iterator first, Iterator last,
			    const unsigned int P)
{
	const unsigned int N = ::std::distance(first, last);

	// Approximation: (log2(N) - 1) <= l < =  log2(N)
	unsigned int l = 1 << (8*sizeof(unsigned int) - 1);
	while(l & N)
		l /= 2;

	return (l < P) && (l*N < (1 << C));
}

template<typename Iterator1, typename Iterator2>
void _copy_partitions(
		const ::std::pair<Iterator1, Iterator1> *source_partitions,
		Iterator2 first, Iterator2 *dest_partitions,
		const unsigned int P)
{
	for (unsigned int i = 0; i < P; ++i)
	{
		dest_partitions[i] = first;

		// The last "advance is very important, it may create space
		// if it is an InsertIterator or something like that.
		::std::advance(first, ::std::distance(
						source_partitions[i].first,
						source_partitions[i].second) );
	}
}

// Divide a given range into P partitions
template<typename Iterator>
void _partition_range(Iterator first, Iterator last,
		::std::pair<Iterator, Iterator> *partitions,
		const unsigned int P)
{
	typedef ::std::pair<Iterator, Iterator> Partition;

	const unsigned int N = ::std::distance(first, last);
	const unsigned int range = N / P + ((N%P)? 1 : 0);

	// All but last partition have same range
	Iterator currentLast = first;
	::std::advance(currentLast, range);
	for (unsigned int i = 0; i < P - 1; ++i)
	{
		partitions[i] = Partition(first, currentLast);
		::std::advance(first, range);
		::std::advance(currentLast, range);
	}

	// Last range may be shorter
	partitions[P - 1] = Partition(first, last);
}

// Given a range, re-arrange the items such that all elements smaller than
// the pivot precede all other values. Returns an Iterator to the first
// element not smaller than the pivot.
template<typename Iterator>
Iterator _stable_pivot_range(Iterator first, Iterator last,
			     const typename Iterator::value_type pivot)
{
	Iterator pivotIt = last;
	while (first < last)
	{
		if (*first < pivot)
			++first;
		else
		{
			Iterator high = first;
			while ( (++high < last) && !(*high < pivot) )
				/* nop */;
			if (high < last)
				::std::iter_swap(first, last);
			first = pivotIt = ++high;
		}
	}

	return pivotIt;
}

template<typename Iterator>
Iterator _pivot_range(Iterator first, Iterator last,
	const typename Iterator::value_type pivot)
{
	while (first < last)
	{
		if (*first < pivot)
			++first;
		else
		{
			while ( (first < --last) && !(*last < pivot) )
				/* nop */;
			::std::iter_swap(first, last);
		}
	}

	return last;
}

template<typename Iterator>
void _partition_range_by_pivots(Iterator first, Iterator last,
	const ::std::vector<typename Iterator::value_type> &pivots,
	::std::pair<Iterator, Iterator> *partitions,
	const unsigned int P)
{
	typedef ::std::pair<Iterator, Iterator> Partition;
	Iterator ptable[P];
	typename Iterator::value_type pvts[pivots.size()];

	::std::vector<Iterator> borders;

	bool used[pivots.size()];
	::std::fill(&used[0], used + pivots.size(), false);

	borders.push_back(first);
	borders.push_back(last);
	partitions[0].first	= first;
	partitions[0].second	= last;
	unsigned int p = 1;
	for (/* nop */; (1 << p) <= (int)P; ++p)
	{
		const int PROC = (1 << p);
		const int PIVOTS = (1 << (p-1));
		OMPTL_ASSERT(PIVOTS <= (int)pivots.size());

		int t;
		#pragma omp parallel for default(shared) private(t)
		for (t = 0; t < PIVOTS; ++t)
		{
			const int index = int(P / PROC) +
						2 * t * int(P / PROC) - 1;
			OMPTL_ASSERT(index < (int)pivots.size());
			OMPTL_ASSERT(!used[index]);
			used[index] = true;
			pvts[t] = pivots[index];
/*::std::cout << "pvts T: " << t << " --> " << index <<
	" " << pvts[t] << ::std::endl;*/
		}

		#pragma omp parallel for default(shared) private(t)
		for (t = 0; t < PIVOTS; ++t)
			ptable[t] = _pivot_range(partitions[t].first,
						 partitions[t].second,
						 pvts[t]);

		for (t = 0; t < PIVOTS; ++t)
		{
// ::std::cout << "ADD: " << ::std::distance(first, ptable[t]) << ::std::endl;
			borders.push_back(ptable[t]);
		}

		::std::sort(borders.begin(), borders.end());

		#pragma omp parallel for default(shared) private(t)
		for (t = 0; t < (int)borders.size() - 1; ++t)
		{
			partitions[t].first	= borders[t];
			partitions[t].second	= borders[t + 1];
		}

/*::std::cout << "PASS: " << p << ::std::endl;
		for (t = 0; t < (1 << p); ++t)
::std::cout << t << ": " << ::std::distance(first, partitions[t].first)
		<< " - " << ::std::distance(first, partitions[t].second)
		<< ::std::endl;*/
	}

	for (unsigned int i = 0; i < pivots.size(); ++i)
		if(!used[i])
			pvts[i] = pivots[i];

	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < (int)pivots.size(); ++t)
		if (!used[t])
			ptable[t] = _pivot_range(partitions[t].first,
						partitions[t].second,
						pvts[t]);


	for (unsigned int i = 0; i < pivots.size(); ++i)
	{
		if (!used[i])
		{
// ::std::cout << "LAST ADD: " << ::std::distance(first, ptable[i]) << ::std::endl;
			borders.push_back(ptable[i]);
		}
	}

	::std::sort(borders.begin(), borders.end());

	OMPTL_ASSERT(borders.size() - 1 == P);
 	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < (int)P; ++t)
	{
		partitions[t].first	= borders[t];
		partitions[t].second	= borders[t + 1];
	}

// ::std::cout << "LAST: " << p << ::std::endl;
// 	for (t = 0; t < P; ++t)
// ::std::cout << t << ": " << ::std::distance(first, partitions[t].first)
// 	<< " - " << ::std::distance(first, partitions[t].second)
// 	<< ::std::endl;

}

template<typename Iterator>
void _partition_range_stable_by_pivots(Iterator first, Iterator last,
	const ::std::vector<typename Iterator::value_type> &pivots,
	::std::pair<Iterator, Iterator> *partitions,
	const unsigned int P)
{
	typedef ::std::pair<Iterator, Iterator> Partition;

	Iterator start = first;
	for (unsigned int i = 0; i < P - 1; ++i)
	{
		Iterator low = start;

		while (low < last)
		{
			// Find a value not lower than the pivot.
			while( (*low < pivots[i]) && (low < last) )
				::std::advance(low, 1);

			// Entire range scanned ?
			if (low == last) break;

			// Find a value lower than the pivot, starting from
			// low, working our way up.
			Iterator high = low;
			::std::advance(high, 1);
			while( !(*high < pivots[i]) && (high < last) )
				::std::advance(high, 1);

			// Entire range scanned ?
			if (high == last) break;

			// Swap values
			OMPTL_ASSERT( !(*low < pivots[i]) && (*high < pivots[i]) );
			::std::iter_swap(low, high);
		}

		partitions[i] = Partition(start, low);
		start = low;
	}
	partitions[P - 1] = Partition(start, last);
}

/*
 * The sample ratio is used to sample more data. This way, the pivots can be
 * chosen more wisely, which is our only guarantee we can generate partitions
 * of equal size.
 */
template<typename RandomAccessIterator>
void _find_pivots(RandomAccessIterator first, RandomAccessIterator last,
	::std::vector<typename RandomAccessIterator::value_type> &pivots,
	const unsigned int P,
	unsigned int SAMPLE_RATIO = 10)
{
	OMPTL_ASSERT(SAMPLE_RATIO > 0);
	const unsigned int N = ::std::distance(first, last);
	OMPTL_ASSERT(N > P);

	// Adjust the constant. Erm.
	while (SAMPLE_RATIO * (P + 1) > N)
		SAMPLE_RATIO /= 2;

	pivots.clear();
	pivots.reserve(P - 1);

	::std::vector<typename RandomAccessIterator::value_type> samples;
	const unsigned int NSAMPLES = SAMPLE_RATIO * P + SAMPLE_RATIO;
	samples.reserve(NSAMPLES);

	for (unsigned int i = 0; i < NSAMPLES; ++i)
	{
		const unsigned int offset = i * (N-1) / (NSAMPLES - 1);
		OMPTL_ASSERT(offset < N);
		samples.push_back(*(first + offset));
// std::cout << "offset: " << offset << " sample: " << samples[i] << std::endl;
	}
	OMPTL_ASSERT(samples.size() == NSAMPLES);

	// Sort samples to create relative ordering in data
	::std::sort(samples.begin(), samples.end());

	// Take pivots from sampled data
	for (unsigned int i = 1; i < P; ++i)
	{
		pivots.push_back(samples[i * samples.size() / P]);
/*std::cout << "pivot: " << i << " idx: " << (i * samples.size() / P)
	<< " " << pivots[i-1] << std::endl;*/
	}
	OMPTL_ASSERT(pivots.size() == P - 1);
}

}  // namespace omptl

#endif /* OMPTL_TOOLS_H */
