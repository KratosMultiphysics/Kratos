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


#ifndef OMPTL_ALGORITHM_H
#define OMPTL_ALGORITHM_H

#ifndef OMPTL_ALGORITHM
#warning <omptl_algorithm> not included. Do not include <omptl_algorithm.h>!
#include <omptl_algorithm>
#endif /* OMPTL_ALGORITHM */

#include <omptl_tools.h>
#include <omptl_numeric>

#include <iterator>

namespace omptl
{

/*
 * Not (yet) paralellized due to data dependance.
 */
template <class ForwardIterator>
ForwardIterator adjacent_find(ForwardIterator first, ForwardIterator last,
			      unsigned int P)
{
	return ::std::adjacent_find(first, last);
}

/*
 * Not (yet) paralellized due to data dependance.
 */
template <class ForwardIterator, class BinaryPredicate>
ForwardIterator adjacent_find(ForwardIterator first, ForwardIterator last,
                              BinaryPredicate binary_pred, unsigned int P)
{
	return ::std::adjacent_find(first, last, binary_pred);
}

template <class ForwardIterator, class T, class StrictWeakOrdering>
bool binary_search(ForwardIterator first, ForwardIterator last, const T& value,
                   StrictWeakOrdering comp, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::binary_search(first, last, value, comp);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	bool results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::binary_search(partitions[t].first,
						  partitions[t].second,
						  value, comp);

	return ::std::count(static_cast<bool *>(results), results + P, false);
}

template <class IteratorIn, class IteratorOut,
	  class IteratorInTag, class IteratorOutTag>
IteratorOut _copy(IteratorIn first, IteratorIn last, IteratorOut result,
		unsigned int P, IteratorInTag, IteratorOutTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::copy(first, last, result);

	::std::pair<IteratorIn, IteratorIn> source_partitions[P];
	::omptl::_partition_range(first, last, source_partitions, P);

	IteratorOut dest_partitions[P];
	::omptl::_copy_partitions(source_partitions, result, dest_partitions,P);

	IteratorOut results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::copy(source_partitions[t].first,
					     source_partitions[t].second,
					     dest_partitions[t]);

	return results[P - 1];
}

template <class InputIterator, class OutputIterator, class InputIteratorTag>
OutputIterator _copy(InputIterator first, InputIterator last,
		    OutputIterator result, unsigned int P,
		    InputIteratorTag, ::std::output_iterator_tag)
{
	return ::std::copy(first, last, result);
}

template <class InputIterator, class OutputIterator, class OutputIteratorTag>
OutputIterator _copy(InputIterator first, InputIterator last,
		    OutputIterator result, unsigned int P,
		    ::std::input_iterator_tag, OutputIteratorTag)
{
	return ::std::copy(first, last, result);
}

template <class InputIterator, class OutputIterator>
OutputIterator copy(InputIterator first, InputIterator last,
		    OutputIterator result, unsigned int P)
{
	return ::omptl::_copy(first, last, result, P,
	typename ::std::iterator_traits<InputIterator>::iterator_category(),
	typename ::std::iterator_traits<OutputIterator>::iterator_category());
}

template <class BidirectionalIterator1, class BidirectionalIterator2>
BidirectionalIterator2 copy_backward(BidirectionalIterator1 first,
                                     BidirectionalIterator1 last,
                                     BidirectionalIterator2 result,
				     unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::copy_backward(first, last, result);

	::std::pair<BidirectionalIterator1, BidirectionalIterator1>
							source_partitions[P];
	::omptl::_partition_range(first, last, source_partitions, P);

	BidirectionalIterator2 dest_partitions[P];
	::omptl::_copy_partitions(source_partitions, result,
				  dest_partitions, P);

	BidirectionalIterator2 results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::copy_backward(source_partitions[t].first,
						  source_partitions[t].second,
						  dest_partitions[t]);

	return results[P - 1];
}

template <class Iterator, class EqualityComparable, class IteratorTag>
typename ::std::iterator_traits<Iterator>::difference_type
_count(Iterator first, Iterator last, const EqualityComparable& value,
	unsigned int P, ::std::input_iterator_tag)
{
	return ::std::count(first, last, value);
}

template <class Iterator, class EqualityComparable, class IteratorTag>
typename ::std::iterator_traits<Iterator>::difference_type
_count(Iterator first, Iterator last,
      const EqualityComparable& value, unsigned int P, IteratorTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::count(first, last, value);

	::std::pair<Iterator, Iterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	typedef typename ::std::iterator_traits<Iterator>::difference_type Tdif;
	Tdif results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::count(partitions[t].first,
					 partitions[t].second, value);

	return ::std::accumulate(static_cast<Tdif *>(results),
					 results + P, 0);
}

template <class InputIterator, class EqualityComparable>
typename ::std::iterator_traits<InputIterator>::difference_type
count(InputIterator first, InputIterator last,
      const EqualityComparable& value, unsigned int P)
{
	return ::omptl::_count(first, last, value, P,
	typename ::std::iterator_traits<InputIterator>::iterator_category());
}

template <class InputIterator, class EqualityComparable, class Size>
void count(InputIterator first, InputIterator last,
           const EqualityComparable& value, Size& n, unsigned int P)
{
	n = ::omptl::count(first, last, value, P);
}


template <class InputIterator, class Predicate>
typename ::std::iterator_traits<InputIterator>::difference_type
_count_if(InputIterator first, InputIterator last, Predicate pred,
	 unsigned int P, ::std::input_iterator_tag)
{
	return ::std::count_if(first, last, pred);
}

template <class Iterator, class Predicate, class IteratorTag>
typename ::std::iterator_traits<Iterator>::difference_type
_count_if(Iterator first, Iterator last, Predicate pred,
	 unsigned int P, IteratorTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::count_if(first, last, pred);

	::std::pair<Iterator, Iterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	typedef typename ::std::iterator_traits<Iterator>::difference_type Tdif;
	Tdif results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::count_if(partitions[t].first,
					 partitions[t].second, pred);

	return ::omptl::accumulate<Tdif>(&results[0], results + P, Tdif(0));
}

template <class InputIterator, class Predicate>
typename ::std::iterator_traits<InputIterator>::difference_type
count_if(InputIterator first, InputIterator last,
	 Predicate pred, unsigned int P)
{
	return ::omptl::_count_if(first, last, pred,
	typename ::std::iterator_traits<InputIterator>::iterator_category());
}

template <class InputIterator, class Predicate, class Size>
void count_if(InputIterator first, InputIterator last,
              Predicate pred, Size& n, unsigned int P)
{
	n = ::omptl::count_if(first, last, pred, P);
}

template <class Iterator1, class Iterator2, class BinaryPredicate,
	  class Iterator1Tag, class Iterator2Tag>
bool _equal(Iterator1 first1, Iterator1 last1,
	    Iterator2 first2, BinaryPredicate binary_pred,
	    unsigned int P, Iterator1Tag, Iterator2Tag)
{
	if (_linear_serial_is_faster(first1, last1, P))
		return ::std::equal(first1, last1, first2, binary_pred);

	::std::pair<Iterator1, Iterator1> source_partitions[P];
	::omptl::_partition_range(first1, last1, source_partitions, P);

	Iterator2 dest_partitions[P];
	::omptl::_copy_partitions(source_partitions, first2,
				  dest_partitions, P);

	bool results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::equal(source_partitions[t].first,
					     source_partitions[t].second,
					     dest_partitions[t], binary_pred);

	return ::std::count(static_cast<bool *>(results), results + P, false);
}

template <class InputIterator1, class InputIterator2,
          class BinaryPredicate, class InputIterator1Tag>
bool _equal(InputIterator1 first1, InputIterator1 last1,
           InputIterator2 first2, BinaryPredicate binary_pred, unsigned int P,
	  InputIterator1Tag, ::std::input_iterator_tag)
{
	return ::std::equal(first1, last1, first2, binary_pred);
}

template <class InputIterator1, class InputIterator2,
          class BinaryPredicate, class InputIterator2Tag>
bool _equal(InputIterator1 first1, InputIterator1 last1,
            InputIterator2 first2, BinaryPredicate binary_pred, unsigned int P,
	    ::std::input_iterator_tag, InputIterator2Tag)
{
	return ::std::equal(first1, last1, first2, binary_pred);
}

template <class InputIterator1, class InputIterator2,
          class BinaryPredicate>
bool equal(InputIterator1 first1, InputIterator1 last1,
           InputIterator2 first2, BinaryPredicate binary_pred, unsigned int P)
{
	return ::omptl::_equal(first1, last1, first2, binary_pred,
	typename ::std::iterator_traits<InputIterator1>::iterator_category(),
	typename ::std::iterator_traits<InputIterator2>::iterator_category());
}

//TODO
template <class ForwardIterator, class T, class StrictWeakOrdering>
::std::pair<ForwardIterator, ForwardIterator>
equal_range(ForwardIterator first, ForwardIterator last, const T& value,
            StrictWeakOrdering comp, unsigned int P)
{
	return ::std::equal_range(first, last, value, comp);
}

template <class ForwardIterator, class T>
void fill(ForwardIterator first, ForwardIterator last,
	  const T& value, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		::std::fill(first, last, value);

	::std::pair<ForwardIterator, ForwardIterator> source_partitions[P];
	::omptl::_partition_range(first, last, source_partitions, P);

	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::fill(source_partitions[t].first,
			    source_partitions[t].second, value);
}

template <class Iterator, class Size, class T, class IteratorTag>
Iterator _fill_n(Iterator first, Size n, const T& value,
			unsigned int P, IteratorTag)
{
	const Size range = (n / P) + ( (n % P) ? 1 : 0 );
	Size ranges[P];
	::std::fill(&ranges[0], P - 1, range);
	ranges[P - 1] = n - (P - 1) * range;

	Iterator source_partitions[P];
	source_partitions[0] = first;
	for (unsigned int i = 1; i < P; ++i)
	{
		source_partitions[i] = source_partitions[i - 1];
		::std::advance(source_partitions[i], range);
	}

	Iterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::fill_n(source_partitions[t], ranges[t],
					   value);

	return results[P - 1];
}

template <class OutputIterator, class Size, class T>
OutputIterator _fill_n(OutputIterator first, Size n, const T& value,
			unsigned int P, ::std::output_iterator_tag)
{
	return ::std::fill_n(first, n, value);
}

template <class OutputIterator, class Size, class T>
OutputIterator fill_n(OutputIterator first, Size n,
		      const T& value, unsigned int P)
{
	return ::omptl::_fill_n(first, n, value, P,
	typename ::std::iterator_traits<OutputIterator>::iterator_category());
}

template<class InputIterator, class EqualityComparable>
InputIterator _find(InputIterator first, InputIterator last,
                   const EqualityComparable& value, unsigned int P,
		   ::std::input_iterator_tag)
{
	return std::find(first, last, value);
}

template<class Iterator, class EqualityComparable, class IteratorTag>
Iterator _find(Iterator first, Iterator last, const EqualityComparable& value,
		unsigned int P, IteratorTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::find(first, last, value);

	::std::pair<Iterator, Iterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	Iterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::find(partitions[t].first,
					 partitions[t].second, value);

	Iterator *result = ::std::find_if(results, results + P,
			::std::bind2nd(::std::not_equal_to<Iterator>(), last) );

	if (result != results + P)
		return *result;

	return last;
}

template<class InputIterator, class EqualityComparable>
InputIterator find(InputIterator first, InputIterator last,
                   const EqualityComparable& value, unsigned int P)
{
	return ::omptl::_find(first, last, value, P,
		::std::iterator_traits<InputIterator>::iterator_category());

}

template<class InputIterator, class Predicate>
InputIterator _find_if( InputIterator first, InputIterator last,
			Predicate pred, unsigned int P,
			::std::input_iterator_tag)
{
	return ::std::find_if(first, last, pred);
}

template<class Iterator, class Predicate, class IteratorTag>
Iterator _find_if( Iterator first, Iterator last, Predicate pred,
			unsigned int P, IteratorTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::find_if(first, last, pred);

	::std::pair<Iterator, Iterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	Iterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::find_if(partitions[t].first,
					 partitions[t].second, pred);

	Iterator *result = ::std::find_if(results, results + P,
		::std::bind2nd(::std::not_equal_to<Iterator>(), last) );

	if (result != results + P)
		return *result;

	return last;
}

template<class InputIterator, class Predicate>
InputIterator find_if(InputIterator first, InputIterator last,
                      Predicate pred, unsigned int P)
{
	return ::omptl::_find_if(first, last, pred, P,
	typename ::std::iterator_traits<InputIterator>::iterator_category());
}

// TODO
template <class ForwardIterator1, class ForwardIterator2,
          class BinaryPredicate>
ForwardIterator1 find_end(ForwardIterator1 first1, ForwardIterator1 last1,
			  ForwardIterator2 first2, ForwardIterator2 last2,
			  BinaryPredicate comp, unsigned int P)
{
	return ::std::find_end(first1, last1, first2, last2, comp);
}

template <class InputIterator, class ForwardIterator, class BinaryPredicate>
InputIterator _find_first_of(InputIterator first1, InputIterator last1,
			ForwardIterator first2, ForwardIterator last2,
			BinaryPredicate comp, unsigned int P,
			::std::input_iterator_tag)
{
	return ::std::find_first_of(first1, last1, first2, last2, comp);
}

// find_first_of suffers from a loss of efficiency, and potentially a loss of
// performance when executed in parallel!
template <class Iterator, class ForwardIterator, class BinaryPredicate,
	  class IteratorTag>
Iterator _find_first_of(Iterator first1, Iterator last1,
			ForwardIterator first2, ForwardIterator last2,
			BinaryPredicate comp, unsigned int P, IteratorTag)
{
	if (_linear_serial_is_faster(first1, last2, P))
		return ::std::find_first_of(first1, last1,
					    first2, last2, comp);

	::std::pair<Iterator, Iterator> partitions[P];
	::omptl::_partition_range(first1, last2, partitions, P);

	Iterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::find_first_of(partitions[t].first, partitions[t].second,
					first2, last2, comp);

	Iterator *result = ::std::find_if(results, results + P,
		::std::bind2nd(::std::not_equal_to<Iterator>(), last1) );

	if (result != results + P)
		return *result;

	return last1;
}

template <class InputIterator, class ForwardIterator, class BinaryPredicate>
InputIterator find_first_of(InputIterator first1, InputIterator last1,
			ForwardIterator first2, ForwardIterator last2,
			BinaryPredicate comp, unsigned int P)
{
	return ::omptl::_find_first_of(first1, last1, first2, last2, comp,
	typename ::std::iterator_traits<InputIterator>::iterator_category() );
}

template <class InputIterator, class UnaryFunction>
UnaryFunction _for_each(InputIterator first, InputIterator last,
		   UnaryFunction f, unsigned int P, ::std::input_iterator_tag)
{
	return ::std::for_each(first, last, f);
}

template <class Iterator, class UnaryFunction, class IteratorTag>
UnaryFunction _for_each(Iterator first, Iterator last,
		       UnaryFunction f, unsigned int P, IteratorTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::for_each(first, last, f);

	::std::pair<Iterator, Iterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::for_each(partitions[t].first, partitions[t].second, f);

	return f;
}

template <class InputIterator, class UnaryFunction>
UnaryFunction for_each(InputIterator first, InputIterator last,
		       UnaryFunction f, unsigned int P)
{
	return ::omptl::_for_each(first, last, f, P,
	typename ::std::iterator_traits<InputIterator>::iterator_category() );
}

template <class ForwardIterator, class Generator>
void generate(ForwardIterator first, ForwardIterator last, Generator gen)
{
	::std::generate(first, last, gen);
}

template <class ForwardIterator, class Generator>
void par_generate(ForwardIterator first, ForwardIterator last,
		  Generator gen, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
	{
		::std::generate(first, last, gen);
		return;
	}

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::generate(partitions[t].first, partitions[t].second, gen);
}

template <class RandomAccessIterator, class StrictWeakOrdering>
void push_heap(RandomAccessIterator first, RandomAccessIterator last,
               StrictWeakOrdering comp, unsigned int P)
{
	return ::std::push_heap(first, last, comp);
}

template <class RandomAccessIterator, class StrictWeakOrdering>
inline void pop_heap(RandomAccessIterator first, RandomAccessIterator last,
                     StrictWeakOrdering comp, unsigned int P)
{
	return ::std::pop_heap(first, last, comp);
}

template <class RandomAccessIterator, class StrictWeakOrdering>
void make_heap(RandomAccessIterator first, RandomAccessIterator last,
               StrictWeakOrdering comp, unsigned int P)
{
	return ::std::make_heap(first, last, comp);
}

template <class RandomAccessIterator, class StrictWeakOrdering>
void sort_heap(RandomAccessIterator first, RandomAccessIterator last,
               StrictWeakOrdering comp, unsigned int P)
{
	return ::std::sort_heap(first, last, comp);
}

template <class InputIterator1, class InputIterator2, class StrictWeakOrdering,
	  class Iterator1Tag>
bool _includes(InputIterator1 first1, InputIterator1 last1,
              InputIterator2 first2, InputIterator2 last2,
              StrictWeakOrdering comp, unsigned int P,
		Iterator1Tag, ::std::input_iterator_tag)
{
	return ::std::includes(first1, last1, first2, last2, comp);
}

template <class InputIterator1, class InputIterator2, class StrictWeakOrdering,
	  class Iterator2Tag>
bool _includes(InputIterator1 first1, InputIterator1 last1,
              InputIterator2 first2, InputIterator2 last2,
              StrictWeakOrdering comp, unsigned int P,
		::std::input_iterator_tag, Iterator2Tag)
{
	return ::std::includes(first1, last1, first2, last2, comp);
}

template <class Iterator1, class Iterator2, class StrictWeakOrdering,
	  class Iterator1Tag, class Iterator2Tag>
bool _includes(Iterator1 first1, Iterator1 last1,
              Iterator2 first2, Iterator2 last2,
              StrictWeakOrdering comp, unsigned int P,
		Iterator1Tag, Iterator2Tag)
{
	if (_linear_serial_is_faster(first2, last2, P))
		return ::std::includes(first1, last1, first2, last2, comp);

	::std::pair<Iterator2, Iterator2> partitions[P];
	::omptl::_partition_range(first2, last2, partitions, P);

	bool results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::includes(first1, last1,
						partitions[t].first,
						partitions[t].second, comp);

	return (::std::count(&results[0], results + P, true) == P);
}

template <class InputIterator1, class InputIterator2, class StrictWeakOrdering>
bool includes(InputIterator1 first1, InputIterator1 last1,
              InputIterator2 first2, InputIterator2 last2,
              StrictWeakOrdering comp, unsigned int P)
{
	return ::omptl::_includes(first1, last1, first2, last2, comp,
	typename ::std::iterator_traits<InputIterator1>::iterator_category(),
	typename ::std::iterator_traits<InputIterator2>::iterator_category());
}

template <class InputIterator1, class InputIterator2, class BinaryPredicate>
bool lexicographical_compare(InputIterator1 first1, InputIterator1 last1,
                             InputIterator2 first2, InputIterator2 last2,
                             BinaryPredicate comp, unsigned int P)
{
	return ::std::lexicographical_compare(first1, last1,
					      first2, last2, comp);
}

template <class ForwardIterator, class T, class StrictWeakOrdering>
ForwardIterator lower_bound(ForwardIterator first, ForwardIterator last,
                            const T& value, StrictWeakOrdering comp,
			    unsigned int P)
{
	if (_logn_serial_is_faster(first, last, P))
		return ::std::lower_bound(first, last, value, comp);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	ForwardIterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::lower_bound(partitions[t].first,
						partitions[t].second,
						value, comp);

	ForwardIterator *result = ::std::find_if(results, results + P,
		::std::bind2nd(::std::not_equal_to<ForwardIterator>(), last) );

	if (result != results + P)
		return *result;

	return last;
}

// Not parallelized, dependencies between data.
template <class InputIterator1, class InputIterator2, class OutputIterator,
          class StrictWeakOrdering>
OutputIterator merge(InputIterator1 first1, InputIterator1 last1,
                     InputIterator2 first2, InputIterator2 last2,
                     OutputIterator result,
		     StrictWeakOrdering comp, unsigned int P)
{
	return ::std::merge(first1, last1, first2, last2, result, comp);
}

template <class ForwardIterator, class BinaryPredicate>
ForwardIterator min_element(ForwardIterator first, ForwardIterator last,
                            BinaryPredicate comp, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::min_element(first, last, comp);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	ForwardIterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::min_element(partitions[t].first,
					 	partitions[t].second, comp);

	// There should be a better way. Is there a way to compare
	// dereferenced iterators ?
	ForwardIterator result = results[0];
	for (unsigned int i = 1; i < P; ++i)
		if (comp(results[i], result)) // i.e.: (results[i] < result)
			result = results[i];

	return result;
}

template <class ForwardIterator, class BinaryPredicate>
ForwardIterator max_element(ForwardIterator first, ForwardIterator last,
                            BinaryPredicate comp, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::max_element(first, last, comp);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	ForwardIterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::max_element(partitions[t].first,
					 	partitions[t].second, comp);

	// There should be a better way. Is there a way to compare
	// dereferenced iterators ?
	ForwardIterator result = results[0];
	for (unsigned int i = 1; i < P; ++i)
		if (comp(result, results[i])) // i.e.: (result < results[i])
			result = results[i];

	return result;
}

template <class InputIterator1, class InputIterator2, class BinaryPredicate,
	  class Iterator1Tag>
::std::pair<InputIterator1, InputIterator2>
_mismatch(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
         BinaryPredicate binary_pred, unsigned int P, Iterator1Tag,
	 ::std::input_iterator_tag)
{
	return ::std::mismatch(first1, last1, first2, binary_pred);
}

template <class InputIterator1, class InputIterator2, class BinaryPredicate,
	  class Iterator2Tag>
::std::pair<InputIterator1, InputIterator2>
_mismatch(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
         BinaryPredicate binary_pred, unsigned int P, ::std::input_iterator_tag,
	 Iterator2Tag)
{
	return ::std::mismatch(first1, last1, first2, binary_pred);
}

template <class Iterator1, class Iterator2, class BinaryPredicate,
	  class Iterator1Tag, class Iterator2Tag>
::std::pair<Iterator1, Iterator2>
_mismatch(Iterator1 first1, Iterator1 last1, Iterator2 first2,
         BinaryPredicate binary_pred, unsigned int P,
	 Iterator1Tag, Iterator2Tag)
{
	if (_linear_serial_is_faster(first1, last1, P))
		return ::std::mismatch(first1, last1, first2, binary_pred);

	::std::pair<Iterator1, Iterator1> source_partitions[P];
	::omptl::_partition_range(first1, last1, source_partitions, P);

	Iterator2 dest_partitions[P];
	::omptl::_copy_partitions(source_partitions, first2,
				  dest_partitions, P);

	::std::pair<Iterator1, Iterator2> results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::mismatch(source_partitions[t].first,
					     source_partitions[t].second,
					     dest_partitions[t], binary_pred);

	// This could have been done more elegantly with select1st
	for (unsigned int i = 0; i < P - 1; ++i)
		if (results[i].first != source_partitions[i].second)
			return results[i];

	return results[P - 1];
}

template <class InputIterator1, class InputIterator2, class BinaryPredicate>
::std::pair<InputIterator1, InputIterator2>
mismatch(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
         BinaryPredicate binary_pred, unsigned int P)
{
	return ::omptl::_mismatch(first1, last1, first2, binary_pred, P,
	typename ::std::iterator_traits<InputIterator1>::iterator_category(),
	typename ::std::iterator_traits<InputIterator2>::iterator_category());
}

// TODO How can this be parallelized ?
template <class RandomAccessIterator, class StrictWeakOrdering>
void nth_element(RandomAccessIterator first, RandomAccessIterator nth,
                 RandomAccessIterator last,
		 StrictWeakOrdering comp, unsigned int P)
{
	::std::nth_element(first, nth, last, comp);
}

template <class RandomAccessIterator, class StrictWeakOrdering>
void partial_sort(RandomAccessIterator first,
                  RandomAccessIterator middle,
                  RandomAccessIterator last,
                  StrictWeakOrdering comp, unsigned int P)
{
	::omptl::_pivot_range(first, last, *middle);
	::omptl::sort(first, middle, P);
}

// Not parallelized due to dependencies.
template <class InputIterator, class RandomAccessIterator>
RandomAccessIterator
partial_sort_copy(InputIterator first, InputIterator last,
                  RandomAccessIterator result_first,
                  RandomAccessIterator result_last, unsigned int P)
{
	return ::std::partial_sort_copy(first, last,
					result_first, result_last);
}

// Not parallelized due to dependencies.
template <class InputIterator, class RandomAccessIterator,
          class StrictWeakOrdering>
RandomAccessIterator
partial_sort_copy(InputIterator first, InputIterator last,
                  RandomAccessIterator result_first,
                  RandomAccessIterator result_last, StrictWeakOrdering comp,
		  unsigned int P)
{
	return ::std::partial_sort_copy(first, last,
					result_first, result_last, comp);
}

// Not (yet) parallelized, not straightforward due to possible dependencies
// between subtasks.
template <class ForwardIterator, class Predicate>
ForwardIterator partition(ForwardIterator first, ForwardIterator last,
			  Predicate pred, unsigned int P)
{
	return ::std::partition(first, last, pred);
}

// Not (yet) parallelized, not straightforward due to possible dependencies
// between subtasks.
template <class ForwardIterator, class Predicate>
ForwardIterator stable_partition(ForwardIterator first, ForwardIterator last,
                                 Predicate pred, unsigned int P)
{
	return ::std::stable_partition(first, last, pred);
}

template <class BidirectionalIterator, class StrictWeakOrdering>
bool next_permutation(BidirectionalIterator first, BidirectionalIterator last,
                      StrictWeakOrdering comp, unsigned int P)
{
	return ::std::next_permutation(first, last, comp);
}

template <class BidirectionalIterator, class StrictWeakOrdering>
bool prev_permutation(BidirectionalIterator first, BidirectionalIterator last,
                      StrictWeakOrdering comp, unsigned int P)
{
	return ::std::prev_permutation(first, last, comp);
}

template <class RandomAccessIterator, class RandomNumberGenerator>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
                    unsigned int P)
{
	::std::random_shuffle(first, last);
}

template <class RandomAccessIterator, class RandomNumberGenerator>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last)
{
	::std::random_shuffle(first, last);
}


template <class RandomAccessIterator, class RandomNumberGenerator>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
                    RandomNumberGenerator& rgen, unsigned int P)
{
	::std::random_shuffle(first, last, rgen);
}

template <class RandomAccessIterator, class RandomNumberGenerator>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
                    RandomNumberGenerator& rgen)
{
	::std::random_shuffle(first, last, rgen);
}

// Not (yet) parallelized, not straightforward due to possible dependencies
// between subtasks.
template <class ForwardIterator, class T>
ForwardIterator remove(ForwardIterator first, ForwardIterator last,
                       const T& value, unsigned int P)
{
	return ::std::remove(first, last, value);
}

// Not (yet) parallelized, not straightforward due to possible dependencies
// between subtasks.
template <class ForwardIterator, class Predicate>
ForwardIterator remove_if(ForwardIterator first, ForwardIterator last,
                          Predicate pred, unsigned int P)
{
	return ::std::remove_if(first, last, pred);
}

// Not parallelized due to possible complications with OutputIterators.
// No par_remove_copy exists due to possible dependencies between subtasks.
template <class InputIterator, class OutputIterator, class T>
OutputIterator remove_copy(InputIterator first, InputIterator last,
                           OutputIterator result, const T& value,
			   unsigned int P)
{
	return ::std::remove_copy(first, last, result, value);
}

// Not parallelized due to possible complications with OutputIterators.
// No par_remove_copy_if exists due to possible dependencies between subtasks.
template <class InputIterator, class OutputIterator, class Predicate>
OutputIterator remove_copy_if(InputIterator first, InputIterator last,
                              OutputIterator result, Predicate pred,
			      unsigned int P)
{
	return ::std::remove_copy(first, last, result, pred);
}

template <class ForwardIterator, class T>
void replace(ForwardIterator first, ForwardIterator last, const T& old_value,
             const T& new_value, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		::std::replace(first, last, old_value, new_value);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::replace(partitions[t].first, partitions[t].second,
				old_value, new_value);
}

// TODO
template <class InputIterator, class OutputIterator, class T>
OutputIterator replace_copy(InputIterator first, InputIterator last,
                            OutputIterator result, const T& old_value,
                            const T& new_value, unsigned int P)
{
	::std::replace_copy(first, last, result, old_value, new_value);
}

// TODO
template <class InputIterator, class OutputIterator, class Predicate, class T>
OutputIterator replace_copy_if(InputIterator first, InputIterator last,
                               OutputIterator result, Predicate pred,
                               const T& new_value, unsigned int P)
{
	::std::replace_copy_if(first, last, result, pred, new_value);
}

template <class ForwardIterator, class Predicate, class T>
void replace_if(ForwardIterator first, ForwardIterator last, Predicate pred,
                const T& new_value, unsigned int P)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::replace_if(first, last, pred, new_value);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::replace_if(partitions[t].first, partitions[t].second,
				pred, new_value);
}

// TODO
template <class BidirectionalIterator>
void reverse(BidirectionalIterator first, BidirectionalIterator last,
	     unsigned int P = omp_get_max_threads())
{
	::std::reverse(first, last);
}

// TODO
template <class BidirectionalIterator, class OutputIterator>
OutputIterator reverse_copy(BidirectionalIterator first,
			    BidirectionalIterator last,
			    OutputIterator result,
			    unsigned int P = omp_get_max_threads())
{
	return ::std::reverse_copy(first, last, result);
}

// TODO
template <class ForwardIterator>
ForwardIterator rotate( ForwardIterator first, ForwardIterator middle,
			ForwardIterator last,
			unsigned int P = omp_get_max_threads())
{
	return ::std::rotate(first, middle, last);
}

// TODO
template <class ForwardIterator, class OutputIterator>
OutputIterator rotate_copy(ForwardIterator first, ForwardIterator middle,
                           ForwardIterator last, OutputIterator result,
			   unsigned int P = omp_get_max_threads())
{
	return ::std::rotate(first, middle, last, result);
}

template <class ForwardIterator1, class ForwardIterator2,
	  class BinaryPredicate>
ForwardIterator1 search(ForwardIterator1 first1, ForwardIterator1 last1,
                        ForwardIterator2 first2, ForwardIterator2 last2,
                        BinaryPredicate binary_pred, unsigned int P)
{
	if (_linear_serial_is_faster(first1, last1, P))
		return ::std::search(first1, last1, first2, last2,
					 binary_pred);

	::std::pair<ForwardIterator1, ForwardIterator1> source_partitions[P];
	::omptl::_partition_range(first1, last1, source_partitions, P);

	ForwardIterator1 results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::search(source_partitions[t].first,
					   source_partitions[t].second,
					   first2, last2, binary_pred);

	ForwardIterator1 *result = ::std::find_if(results, results + P,
		::std::bind2nd(::std::not_equal_to<ForwardIterator1>(),
				last1));

	if (result != results + P)
		return *result;

	return last1;
}

// TODO
template <class ForwardIterator, class Integer,
          class T, class BinaryPredicate>
ForwardIterator search_n(ForwardIterator first, ForwardIterator last,
                         Integer count, const T& value,
                         BinaryPredicate binary_pred, unsigned int P)
{
	::std::search_n(first, last, count, value, binary_pred);
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class StrictWeakOrdering>
OutputIterator set_difference(InputIterator1 first1, InputIterator1 last1,
				InputIterator2 first2, InputIterator2 last2,
				OutputIterator result, StrictWeakOrdering comp,
				unsigned int P)
{
	return ::std::set_difference(first1, last1,
				     first2, last2, result, comp);
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class StrictWeakOrdering>
OutputIterator set_intersection(InputIterator1 first1, InputIterator1 last1,
				InputIterator2 first2, InputIterator2 last2,
				OutputIterator result, StrictWeakOrdering comp,
			 	unsigned int P)
{
	return ::std::set_intersection( first1, last1,
					first2, last2, result, comp);
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class StrictWeakOrdering>
OutputIterator set_symmetric_difference(InputIterator1 first1,
					InputIterator1 last1,
					InputIterator2 first2,
					InputIterator2 last2,
					OutputIterator result,
					StrictWeakOrdering comp,
			 		unsigned int P)
{
	return ::std::set_symmetric_difference( first1, last1,
						first2, last2, result, comp);
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class StrictWeakOrdering>
OutputIterator set_union(InputIterator1 first1, InputIterator1 last1,
			 InputIterator2 first2, InputIterator2 last2,
			 OutputIterator result, StrictWeakOrdering comp,
			 unsigned int P)
{
	return ::std::set_union(first1, last1, first2, last2, result, comp);
}

template<typename RandomAccessIterator>
void sort(RandomAccessIterator first, RandomAccessIterator last, unsigned int P)
{
	if (::omptl::_nlogn_serial_is_faster(first, last, P))
	{
		::std::sort(first, last);
		return;
	}

	// Generate pivots
	::std::vector<typename RandomAccessIterator::value_type> pivots;
	_find_pivots(first, last, pivots, P);

	// Sort sufficiently to respect pivot order
	::std::pair<RandomAccessIterator, RandomAccessIterator> partitions[P];
	::omptl::_partition_range_by_pivots(first, last, pivots,
					    partitions, P);

	// Sort
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::sort(partitions[t].first, partitions[t].second);
}

template<typename RandomAccessIterator>
void stable_sort(RandomAccessIterator first, RandomAccessIterator last,
		 unsigned int P)
{
	if (::omptl::_nlogn_serial_is_faster(first, last, P))
	{
		::std::stable_sort(first, last);
		return;
	}

	// Generate pivots
	::std::vector<typename RandomAccessIterator::value_type> pivots;
	_find_pivots(first, last, pivots, P);

	// Sort sufficiently to respect pivot order
	::std::pair<RandomAccessIterator, RandomAccessIterator> partitions[P];
	::omptl::_partition_range_stable_by_pivots(first, last, pivots,
						   partitions, P);

	// Sort
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		::std::stable_sort(partitions[t].first, partitions[t].second);
}

template <class ForwardIterator1, class ForwardIterator2>
ForwardIterator2 swap_ranges(ForwardIterator1 first1, ForwardIterator1 last1,
                             ForwardIterator2 first2, unsigned int P)
{
	if (_linear_serial_is_faster(first1, last1, P))
		return ::std::swap_ranges(first1, last1, first2);

	::std::pair<ForwardIterator1, ForwardIterator1> source_partitions[P];
	::omptl::_partition_range(first1, last1, source_partitions, P);

	ForwardIterator2 dest_partitions[P];
	::omptl::_copy_partitions(source_partitions, first2,
				  dest_partitions, P);

	ForwardIterator2 results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::swap_ranges(source_partitions[t].first,
					     source_partitions[t].second,
					     dest_partitions[t]);

	return results[P - 1];
}

template <class InputIterator, class OutputIterator, class UnaryFunction,
	 class IteratorInTag>
OutputIterator _transform(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryFunction op,
			 unsigned int P, IteratorInTag,
			 ::std::output_iterator_tag)
{
	return ::std::transform(first, last, result, op);
}

template <class InputIterator, class OutputIterator, class UnaryFunction,
	 class IteratorOutTag>
OutputIterator _transform(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryFunction op,
			 unsigned int P, ::std::input_iterator_tag,
			 IteratorOutTag)
{
	return ::std::transform(first, last, result, op);
}

template <class IteratorIn, class IteratorOut, class UnaryFunction,
	 class IteratorInTag, class IteratorOutTag>
IteratorOut _transform(IteratorIn first, IteratorIn last,
			IteratorOut result, UnaryFunction op,
			unsigned int P, IteratorInTag, IteratorOutTag)
{
	if (_linear_serial_is_faster(first, last, P))
		return ::std::transform(first, last, result, op);

	::std::pair<IteratorIn, IteratorIn> source_partitions[P];
	::omptl::_partition_range(first, last, source_partitions, P);

	IteratorOut dest_partitions[P];
	::omptl::_copy_partitions(source_partitions, result,
				  dest_partitions, P);

	IteratorOut results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::transform(source_partitions[t].first,
					     source_partitions[t].second,
					     dest_partitions[t], op);

	return results[P - 1];
}

template <class InputIterator, class OutputIterator, class UnaryFunction>
OutputIterator transform(InputIterator first, InputIterator last,
                         OutputIterator result, UnaryFunction op,
			 unsigned int P)
{
	return ::omptl::_transform(first, last, result, op, P,
	typename ::std::iterator_traits<InputIterator>::iterator_category(),
	typename ::std::iterator_traits<OutputIterator>::iterator_category());
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class BinaryFunction, class Iterator2Tag, class IteratorOutTag>
OutputIterator _transform(InputIterator1 first1, InputIterator1 last1,
			InputIterator2 first2, OutputIterator result,
			BinaryFunction binary_op, unsigned int P,
			::std::input_iterator_tag, Iterator2Tag, IteratorOutTag)
{
	return ::std::transform(first1, last1, first2, result, binary_op);
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class BinaryFunction, class Iterator1Tag, class IteratorOutTag>
OutputIterator _transform(InputIterator1 first1, InputIterator1 last1,
			InputIterator2 first2, OutputIterator result,
			BinaryFunction binary_op, unsigned int P,
			Iterator1Tag, ::std::input_iterator_tag, IteratorOutTag)
{
	return ::std::transform(first1, last1, first2, result, binary_op);
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class BinaryFunction, class Iterator1Tag, class Iterator2Tag>
OutputIterator _transform(InputIterator1 first1, InputIterator1 last1,
			InputIterator2 first2, OutputIterator result,
			BinaryFunction binary_op, unsigned int P,
			Iterator1Tag, Iterator2Tag, ::std::output_iterator_tag)
{
	return ::std::transform(first1, last1, first2, result, binary_op);
}

template <class Iterator1, class Iterator2, class IteratorOut,
          class BinaryFunction, class Iterator1Tag, class Iterator2Tag,
	  class IteratorOutTag>
IteratorOut _transform(Iterator1 first1, Iterator1 last1,
			Iterator2 first2, IteratorOut result,
			BinaryFunction binary_op, unsigned int P,
			Iterator1Tag, Iterator2Tag, IteratorOutTag)
{
	if (_linear_serial_is_faster(first1, last1, P))
		return ::std::transform(first1, last1, first2,
					 result, binary_op);

	::std::pair<Iterator1, Iterator1> source_partitions1[P];
	::omptl::_partition_range(first1, last1, source_partitions1, P);

	Iterator2 source_partitions2[P];
	::omptl::_copy_partitions(source_partitions1, first2,
				 source_partitions2 , P);

	IteratorOut dest_partitions[P];
	::omptl::_copy_partitions(source_partitions1, result,
				  dest_partitions, P);

	IteratorOut results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::transform(source_partitions1[t].first,
					      source_partitions1[t].second,
					      source_partitions2[t],
					      dest_partitions[t], binary_op);

	return results[P - 1];
}

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class BinaryFunction>
OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                         InputIterator2 first2, OutputIterator result,
                         BinaryFunction binary_op, unsigned int P)
{
	return ::omptl::_transform(first1, last1, first2, result, binary_op, P,
	typename ::std::iterator_traits<InputIterator1>::iterator_category(),
	typename ::std::iterator_traits<InputIterator2>::iterator_category(),
	typename ::std::iterator_traits<OutputIterator>::iterator_category());
}

template <class ForwardIterator, class BinaryPredicate>
ForwardIterator unique(ForwardIterator first, ForwardIterator last,
                       BinaryPredicate binary_pred =
		       ::std::equal_to<typename ForwardIterator::value_type>(),
		       unsigned int P = omp_get_max_threads())
{
	return ::std::unique_copy(first, last, binary_pred);
}

template <class InputIterator, class OutputIterator, class BinaryPredicate>
OutputIterator unique_copy(InputIterator first, InputIterator last,
                           OutputIterator result,
                           BinaryPredicate binary_pred =
		       ::std::equal_to<typename InputIterator::value_type>(),
		       unsigned int P = omp_get_max_threads())
{
	return ::std::unique_copy(first, last, result, binary_pred);
}

template <class ForwardIterator, class T, class StrictWeakOrdering>
ForwardIterator upper_bound(ForwardIterator first, ForwardIterator last,
                            const T& value, StrictWeakOrdering comp,
			    unsigned int P)
{
	if (_logn_serial_is_faster(first, last, P))
		return ::std::upper_bound(first, last, value, comp);

	::std::pair<ForwardIterator, ForwardIterator> partitions[P];
	::omptl::_partition_range(first, last, partitions, P);

	ForwardIterator results[P];
	int t;
	#pragma omp parallel for default(shared) private(t)
	for (t = 0; t < int(P); ++t)
		results[t] = ::std::upper_bound(partitions[t].first,
					 	partitions[t].second,
						value, comp);

	// There has to be a better way...
	for (unsigned int i = P - 1; i > 0; --i)
		if (results[i] != last)
			return results[i];

	return results[0];
}

} /* namespace omptl */

#endif /* OMPTL_ALGORITHM_H */
