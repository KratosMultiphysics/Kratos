#include <vector>

#include <omptl>
#include <omptl_numeric>
#include <omptl_algorithm>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <omp.h>

const unsigned int N = 100 * (1 << 20);

int main (int argc, char * const argv[])
{
	// Number of threads is derived from environment
	// variable "OMP_NUM_THREADS"
	std::cout << N << std::endl;
	std::cout << "Threads: " << omp_get_max_threads() << std::endl;

	std::vector<unsigned int> v1(N);

	omptl::generate(v1.begin(), v1.end(), std::rand);
	omptl::sort(v1.begin(), v1.end());
	//omptl::random_shuffle(v1.begin(), v1.end() );

	std::vector<unsigned int> v2(N);
	omptl::copy(v1.begin(), v1.end(), v2.begin());
	omptl::for_each(v2.begin(), v2.end(), std::sqrt<unsigned int>);
	std::cout << "Nr 3's: " << omptl::count(v2.begin(), v2.end(), 3)
		<< std::endl;
	std::cout << "Sum: "
		<< omptl::accumulate(v2.begin(), v2.end(), 0) << std::endl;

	std::cout << *v1.begin() << std::endl;

	return 0;
}

