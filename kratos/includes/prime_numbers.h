//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//

#if !defined(KRATOS_PRIME_NUMBERS_H_INCLUDED )
#define  KRATOS_PRIME_NUMBERS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <array>

// Project includes
#include "includes/exception.h"

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Gives a prime number before or after given number.
  /** Has an array of precalculated value up to 1e8 and calculates the rest.
  */
  class PrimeNumbers
    {
		///@name Constants
		///@{

		static constexpr std::size_t mNumberOfPreCalculatedPrimes = 100000;

		static constexpr std::size_t mLargestPreCalculatedPrime = 1299689;

		///@}
    public:
		///@name Type Definitions
		///@{

		///@}
		///@name Life Cycle
      ///@{

      /// Default constructor.
		PrimeNumbers() {}

	  /// Copy constructor.
	  PrimeNumbers(PrimeNumbers const& rOther) = delete;
	  
	  /// Destructor.
      virtual ~PrimeNumbers(){}


      ///@}
      ///@name Operators
      ///@{
	  
	  /// Assignment operator.
	  PrimeNumbers& operator=(PrimeNumbers const& rOther) = delete;

	  std::size_t operator[](std::size_t Index) {
		  return GetPreCalculatedPrime(Index);
	  }

	  ///@}
      ///@name Operations
      ///@{


      ///@}
      ///@name Access
      ///@{

	  static std::size_t GetNumberOfPreCalculatedPrimes() {
		  return mNumberOfPreCalculatedPrimes;
	  }
	  
	  static std::size_t LargestPreCalculatedPrime() {
		  return mLargestPreCalculatedPrime;
	  }

	  static std::size_t GetPreCalculatedPrime(std::size_t Index) {
		  KRATOS_DEBUG_ERROR_IF(Index >= mNumberOfPreCalculatedPrimes) << "Index " << Index 
			  << " is larger than Number of precalculated primes " << mNumberOfPreCalculatedPrimes << std::endl;
		  return mPrecalculatedPrimes[Index];
	  }

      ///@}
      ///@name Inquiry
      ///@{

	  static bool IsPrime(std::size_t TheNumber);


      ///@}
      ///@name Input and output
      ///@{


      ///@}

    private:
      ///@name Static Member Variables
      ///@{

		static const std::array<std::size_t, mNumberOfPreCalculatedPrimes> mPrecalculatedPrimes;

      ///@}
      ///@name Private Operations
      ///@{

		static std::size_t FindGreaterEqualPrecalculatedPrimeNumber(std::size_t TheNumber);

      ///@}
 
    }; // Class PrimeNumbers

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PRIME_NUMBERS_H_INCLUDED  defined
