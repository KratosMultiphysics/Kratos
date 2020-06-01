//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Quirin Aumann
//

#if !defined(KRATOS_COMPLEX_SORT_UTILITY_H_INCLUDED)
#define KRATOS_COMPLEX_SORT_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
namespace ComplexSortUtility
{
    /**
     * @brief Sorts and pairs complex conjugate numbers
     * @param rV A vector of complex numbers to be paired
     * @param tol Tolerance. Defaults to 100*EPS
     * @details This utility pairs complex conjugate values and sorts them with increasing real part.
     *      In each pair, the negative complex part comes first; real numbers are sorted to the end of the vector.
     *      An error is thrown, if the complex numbers cannot be paired in complex conjugates.
     * @author Quirin Aumann
     */
    template <typename ComplexVectorType>
    void PairComplexConjugates(ComplexVectorType& rV, const double tol = std::numeric_limits<double>::epsilon()*100)
    {
        using complex = std::complex<double>;

        KRATOS_TRY
        KRATOS_ERROR_IF( (tol >= 1.) || (tol < 0.) ) << "Invalid tolerance provided\n";

        //find real entries and put them at the end of the vector by sorting imaginary part descending
        std::sort(rV.begin(), rV.end(), [](complex a, complex b) {
            return std::abs(std::imag(a)) > std::abs(std::imag(b));
        });

        //count real entries
        const size_t n_real = std::count_if(rV.begin(), rV.end(), [tol](complex x) {
            return std::abs(std::imag(x)) < tol;
        });

        if( n_real > 0 )
        {
            //sort real entries
            std::sort(rV.end()-n_real, rV.end(), [](complex a, complex b) {
                return std::real(a) < std::real(b);
            });
        }

        const size_t n_cplx = rV.size() - n_real;
        KRATOS_ERROR_IF( (n_cplx % 2) != 0 ) << "Complex numbers cannot be paired (odd number of complex numbers provided)\n";

        //sort remaing complex entries by real part
        std::sort(rV.begin(), rV.end()-n_real, [](complex a, complex b) {
            return std::real(a) < std::real(b);
        });

        //check if real parts occur in pairs
        for( size_t i=0; i<n_cplx; i+=2)
        {
            KRATOS_ERROR_IF( std::abs(std::real(rV(i)) - std::real(rV(i+1))) > tol )
                << "Complex numbers " << rV(i) << "," << rV(i+1) << " cannot be paired (real parts not in pairs): difference "
                << std::abs(std::real(rV(i)) - std::real(rV(i+1))) << " at tolerance " << tol << "\n";
        }

        //now build slices with identical real parts and check inside them for complex conjugate pairs,
        size_t i = 0;
        while( i<n_cplx )
        {
            //count number of identical real parts
            size_t j = 1;
            while( (std::abs(std::real(rV(i)) - std::real(rV(i+j))) < tol) )
            {
                j++;
                //break if end of array is reached
                if( i+j == rV.size() )
                    break;
            }

            //if two identical real parts are found, check for the imaginary part
            if( j == 2 )
            {
                //compare for conjugate imaginary parts
                if( std::abs(std::imag(rV(i)) + std::imag(rV(i+1))) < tol )
                {
                    //make sure negative imaginary part comes first
                    if( std::imag(rV(i)) > 0 )
                    {
                        rV(i) = std::conj(rV(i));
                        rV(i+1) = std::conj(rV(i+1));
                    }
                }
                else
                {
                    KRATOS_ERROR << "Complex numbers " << rV(i) << "," << rV(i+1) << " cannot be paired (imaginary parts not conjugate): difference "
                        << std::abs(std::real(rV(i)) - std::real(rV(i+1))) << " at tolerance " << tol << "\n";
                }
            }
            //if more identical real parts are found, loop over the imaginary parts
            else if( (j % 2) == 0 )
            {
                //sort by absolute imaginary part
                std::sort(rV.begin()+i, rV.begin()+i+j, [](complex a, complex b) {
                    return std::abs(std::imag(a)) < std::abs(std::imag(b));
                });

                size_t ii = i;
                while( ii<j )
                {
                    //count number of identical absolute imaginary parts
                    size_t jj = 1;
                    while( std::abs(std::abs(std::imag(rV(ii))) - std::abs(std::imag(rV(ii+jj)))) < tol )
                    {
                        jj++;
                        //break if end of array is reached
                        if( ii+jj == rV.size() )
                            break;
                    }

                    KRATOS_ERROR_IF( (jj % 2) != 0 )
                        << "Complex numbers cannot be paired (imaginary parts not in pairs)\n";

                    for( size_t k = 0; k<jj+1; k+=2 )
                    {
                        //break if end of array is reached
                        if( ii+k == rV.size() )
                            break;
                        //compare for conjugate imaginary parts
                        if( std::abs(std::imag(rV(ii+k)) + std::imag(rV(ii+k+1))) < tol )
                        {
                            //make sure negative imaginary part comes first
                            if( std::imag(rV(ii+k)) > 0 )
                            {
                                rV(ii+k) = std::conj(rV(ii+k));
                                rV(ii+k+1) = std::conj(rV(ii+k+1));
                            }
                        }
                        else
                        {
                            KRATOS_ERROR << "Complex numbers " << rV(ii+k) << "," << rV(ii+k+1) << " cannot be paired (imaginary parts not conjugate): difference "
                                << std::abs(std::imag(rV(ii+k)) + std::imag(rV(ii+k+1))) << " at tolerance " << tol << "\n";
                        }

                    }
                    ii += jj;
                }

            }
            else
            {
                KRATOS_ERROR << "Complex numbers cannot be paired (imaginary parts not in pairs)\n";
            }

            i += j;
        }

        KRATOS_CATCH("")
    }

} // namespace ComplexSortUtility

} // namespace Kratos

#endif