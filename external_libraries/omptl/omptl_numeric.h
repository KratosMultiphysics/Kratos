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

#ifndef OMPTL_NUMERIC_H
#define OMPTL_NUMERIC_H

#include <utility>
#include <numeric>

#include <omptl_algorithm>
#include <omptl_tools.h>

namespace omptl
{

template <class InputIterator, class T>
T accumulate(InputIterator first, InputIterator last, T init, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::accumulate(first, last, init);

	::std::pair<InputIterator, InputIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	T results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::accumulate(partitions[t].first,
					 partitions[t].second, T(0));

	return ::std::accumulate(results, results + P, init);
}

template <class InputIterator, class T, class BinaryFunction>
T accumulate(InputIterator first, InputIterator last, T init,
             BinaryFunction binary_op, unsigned int P)
{
	return ::std::accumulate(first, last, init, binary_op);
}

template <class InputIterator, class OutputIterator, class BinaryFunction>
OutputIterator adjacent_difference(InputIterator first, InputIterator last,
                                   OutputIterator result,
				   BinaryFunction binary_op, unsigned int P)
{
	return ::std::adjacent_difference(first, last, result, binary_op);
}

// TODO
template <class InputIterator1, class InputIterator2, class T,
          class BinaryFunction1, class BinaryFunction2>
T inner_product(InputIterator1 first1, InputIterator1 last1,
                InputIterator2 first2, T init,
		BinaryFunction1 binary_op1, BinaryFunction2 binary_op2,
		unsigned int P)
{
	return ::std::inner_product(first1, last1, first2, init,
				    binary_op1, binary_op2);
}


template <class InputIterator, class OutputIterator, class BinaryOperation>
OutputIterator partial_sum(InputIterator first, InputIterator last,
                           OutputIterator result, BinaryOperation binary_op,
			   unsigned int P)
{
	return ::std::partial_sum(first, last, result, binary_op);
}

} // namespace omptl

#endif /* OMPTL_NUMERIC_H */
