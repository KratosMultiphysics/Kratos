/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-11-10 14:23:33 $
//   Revision:            $Revision: 1.7 $
//
//


#if !defined(KRATOS_UBLAS_SPACE_H_INCLUDED )
#define  KRATOS_UBLAS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>
#include <numeric>

#ifdef _OPENMP
#include "omp.h"
//#include "omptl"
//#include "omptl_algorithm"
//#include "omptl_numeric"
#endif


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"


namespace Kratos
{
#ifdef _OPENMP
    // The function object multiplies an element by a Factor

    template <class Type>
    class MultValueNoAdd
    {
    private:
        Type Factor; // The value to multiply by
    public:
        // Constructor initializes the value to multiply by

        MultValueNoAdd(const Type& _Val) : Factor(_Val)
        {
        }

        // The function call for the element to be multiplied

        inline Type operator () (const Type& elem) const
        {
            return elem * Factor;
        }
    };

    template <class Type>
    class MultAndAddValue : public std::binary_function<Type, Type, Type>
    {
    private:
        Type Factor; // The value to multiply by
    public:
        // Constructor initializes the value to multiply by

        MultAndAddValue(const Type& _Val) : Factor(_Val)
        {
        }

        // The function call for the element to be multiplied

        inline Type operator () (const Type& elem1, const Type& elem2) const
        {
            return elem1 * Factor + elem2;
        }
    };
#endif

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.

    /** Detail class definition.
     */
    template<class TDataType, class TMatrixType, class TVectorType>
    class UblasSpace
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of UblasSpace
        KRATOS_CLASS_POINTER_DEFINITION(UblasSpace);

        typedef TDataType DataType;

        typedef TMatrixType MatrixType;

        typedef TVectorType VectorType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef typename boost::shared_ptr< TMatrixType > MatrixPointerType;
        typedef typename boost::shared_ptr< TVectorType > VectorPointerType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        UblasSpace()
        {
        }

        /// Destructor.

        virtual ~UblasSpace()
        {
        }


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        static MatrixPointerType CreateEmptyMatrixPointer()
        {
            return MatrixPointerType(new TMatrixType());
        }

        static VectorPointerType CreateEmptyVectorPointer()
        {
            return VectorPointerType(new TVectorType());
        }

        /// return size of vector rV

        static IndexType Size(VectorType const& rV)
        {
            return rV.size();
        }

        /// return number of rows of rM

        static IndexType Size1(MatrixType const& rM)
        {
            return rM.size1();
        }

        /// return number of columns of rM

        static IndexType Size2(MatrixType const& rM)
        {
            return rM.size2();
        }

        /// rXi = rMij

        static void GetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
        {
            rX = column(rM, j);
        }


        ///////////////////////////////// TODO: Take a close look to this method!!!!!!!!!!!!!!!!!!!!!!!!!
        /// rMij = rXi

        static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
        {
            rX = row(rM, j);
        }

        /// rY = rX

        static void Copy(MatrixType const& rX, MatrixType& rY)
        {
            rY.assign(rX);
        }

        /// rY = rX

        static void Copy(VectorType const& rX, VectorType& rY)
        {
#ifndef _OPENMP
            rY.assign(rX);
#else
            //omptl::copy(rX.begin(), rX.end(), rY.begin());
            const int size = rX.size();
            if (rY.size() != size)
                rY.resize(size, false);
            
#pragma omp parallel for 
            for (int i = 0; i < size; i++)
                rY[i] = rX[i];
#endif
        }

        /// rX * rY

        static TDataType Dot(VectorType const& rX, VectorType const& rY)
        {
#ifndef _OPENMP
            return inner_prod(rX, rY);
#else
            vector<unsigned int> partition;
            int number_of_threads = omp_get_max_threads();
            CreatePartition(number_of_threads, rX.size(), partition);

            vector< TDataType > partial_results(number_of_threads);

            int i;
#pragma omp parallel for private(i)
            for (i = 0; i < number_of_threads; i++)
            {
                partial_results[i] = std::inner_product(rX.data().begin() + partition[i],
                        rX.data().begin() + partition[i + 1],
                        rY.data().begin() + partition[i],
                        TDataType());
            }

            double total = TDataType();
            for (int i = 0; i < number_of_threads; i++)
                total += partial_results[i];
            return total;
#endif
        }


        /// ||rX||2

        static double TwoNorm(VectorType const& rX)
        {
            return sqrt(Dot(rX, rX));
        }

        static void Mult(Matrix& rA, VectorType& rX, VectorType& rY)
        {
            axpy_prod(rA, rX, rY, true);
        }

        static void Mult(compressed_matrix<TDataType>& rA, VectorType& rX, VectorType& rY)
        {
#ifndef _OPENMP
            axpy_prod(rA, rX, rY, true);
#else
            ParallelProductNoAdd(rA, rX, rY);
#endif
        }

        static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
        {
            axpy_prod(rX, rA, rY, true);
        } // rY = rAT * rX

        static inline SizeType GraphDegree(IndexType i, TMatrixType& A)
        {
            typename MatrixType::iterator1 a_iterator = A.begin1();
            std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            return ( std::distance(a_iterator.begin(), a_iterator.end()));
#else
            return ( std::distance(begin(a_iterator, boost::numeric::ublas::iterator1_tag()),
                    end(a_iterator, boost::numeric::ublas::iterator1_tag())));
#endif
        }

        static inline void GraphNeighbors(IndexType i, TMatrixType& A, std::vector<IndexType>& neighbors)
        {
            neighbors.clear();
            typename MatrixType::iterator1 a_iterator = A.begin1();
            std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (typename MatrixType::iterator2 row_iterator = a_iterator.begin();
                    row_iterator != a_iterator.end(); ++row_iterator)
            {
#else
            for (typename MatrixType::iterator2 row_iterator = begin(a_iterator,
                    boost::numeric::ublas::iterator1_tag());
                    row_iterator != end(a_iterator,
                    boost::numeric::ublas::iterator1_tag()); ++row_iterator)
            {
#endif
                neighbors.push_back(row_iterator.index2());
            }
        }


        //********************************************************************
        //checks if a multiplication is needed and tries to do otherwise

        static void InplaceMult(VectorType& rX, const double A)
        {

            if (A == 1.00)
            {
            } else if (A == -1.00)
            {
#ifndef _OPENMP
                typename VectorType::iterator x_iterator = rX.begin();
                typename VectorType::iterator end_iterator = rX.end();
                while (x_iterator != end_iterator)
                {
                    *x_iterator = -*x_iterator;
                    x_iterator++;

                }
#else
                //omptl::transform(rX.begin(), rX.end(), rX.begin(), std::negate<double>());
                const int size = rX.size();
                
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] = -rX[i];

#endif
            } else
            {
#ifndef _OPENMP
                rX *= A;
#else
                //                 rX *= A;
                //omptl::transform(rX.begin(), rX.end(), rX.begin(), MultValueNoAdd<double>(A));
                const int size = rX.size();

#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] *= A;
#endif
            }
        }

        //********************************************************************
        //checks if a multiplication is needed and tries to do otherwise
        //ATTENTION it is assumed no aliasing between rX and rY
        // X = A*y;

        static void Assign(VectorType& rX, const double A, const VectorType& rY)
        {
#ifndef _OPENMP
            if (A == 1.00)
                noalias(rX) = rY;
            else if (A == -1.00)
                noalias(rX) = -rY;
            else
                noalias(rX) = A*rY;
#else
            const int size = rY.size();
            if (rX.size() != size)
                rX.resize(size, false);

            if (A == 1.00)
            {
                //omptl::copy(rY.begin(), rY.end(), rX.begin());
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] = rY[i];
            } else if (A == -1.00)
            {
                //omptl::transform(rY.begin(), rY.end(), rX.begin(), std::negate<double>());
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] = -rY[i];
            } else
            {
                //                noalias(rX) = A*rY;
                //omptl::transform(rY.begin(), rY.end(), rX.begin(), MultValueNoAdd<double>(A)); //todo
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] = A * rY[i];
            }
#endif
        }

        //********************************************************************
        //checks if a multiplication is needed and tries to do otherwise
        //ATTENTION it is assumed no aliasing between rX and rY
        // X += A*y;

        static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
        {
#ifndef _OPENMP
            if (A == 1.00)
                noalias(rX) += rY;
            else if (A == -1.00)
                noalias(rX) -= rY;
            else
                noalias(rX) += A*rY;
#else
            const int size = rY.size();
            if (rX.size() != size)
                rX.resize(size, false);

            if (A == 1.00)
            {
                //omptl::transform(rY.data().begin(), rY.data().end(), rX.data().begin(), rX.data().begin(), std::plus<double>());
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] += rY[i];
            } else if (A == -1.00)
            {
                //omptl::transform(rY.data().begin(), rY.data().end(), rX.data().begin(), rX.data().begin(), std::minus<double>());
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] -= rY[i];
            } else
            {
                //omptl::transform(rY.data().begin(), rY.data().end(), rX.data().begin(), rX.data().begin(), MultAndAddValue<double>(A));
#pragma omp parallel for 
                for (int i = 0; i < size; i++)
                    rX[i] += A * rY[i];
            }


#endif

        }

        //********************************************************************

        static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
        {
            Assign(rZ, A, rX); //rZ = A*rX
            UnaliasedAdd(rZ, B, rY); //rZ += B*rY
        }

        static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
        {
            InplaceMult(rY, B);
            UnaliasedAdd(rY, A, rX);
        }


        /// rA[i] * rX

        static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
        {
            return inner_prod(row(rA, i), rX);
        }


        /// rX = A

        static void Set(VectorType& rX, TDataType A)
        {
            std::fill(rX.begin(), rX.end(), A);
        }

        static void Resize(MatrixType& rA, SizeType m, SizeType n)
        {
            rA.resize(m, n, false);
        }

        static void Resize(VectorType& rX, SizeType n)
        {
            rX.resize(n, false);
        }

        static void Clear(MatrixPointerType& pA)
        {
            pA->clear();
            pA->resize(0, 0, false);
        }

        static void Clear(VectorPointerType& pX)
        {
            pX->clear();
            pX->resize(0, false);
        }

        /*	static void Clear(MatrixType& rA)
                {rA.clear();}
	
        static void Clear(VectorType& rX) {rX.clear();}*/

        template<class TOtherMatrixType>
        inline static void ClearData(TOtherMatrixType& rA)
        {
            rA.clear();
        }

        inline static void ClearData(compressed_matrix<TDataType>& rA)
        {
            rA.clear();
            //    	rA.value_data() = unbounded_array<TDataType>();
            //if(rA.non_zeros() != 0) rA.value_data() = unbounded_array<TDataType>();
        }

        inline static void ClearData(VectorType& rX)
        {
            rX = VectorType();
        }

        template<class TOtherMatrixType>
        inline static void ResizeData(TOtherMatrixType& rA, SizeType m)
        {
            rA.resize(m, false);
            //            std::fill(rA.begin(), rA.end(), TDataType());
#ifndef _OPENMP
            std::fill(rA.begin(), rA.end(), TDataType());
#else
            //omptl::fill(rA.begin(), rA.end(), TDataType());
            ParallelFill(rA.begin(), rA.end(), TDataType());

#endif
        }

        inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m)
        {
            rA.value_data().resize(m);
#ifndef _OPENMP
            std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#else
            //omptl::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
            ParallelFill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#endif
        }

        inline static void ResizeData(VectorType& rX, SizeType m)
        {
            rX.resize(m, false);
#ifndef _OPENMP
            std::fill(rX.begin(), rX.end(), TDataType());
#else
            //omptl::fill(rX.begin(), rX.end(), TDataType());
            ParallelFill(rX.begin(), rX.end(), TDataType());
#endif
        }

        template<class TOtherMatrixType>
        inline static void SetToZero(TOtherMatrixType& rA)
        {
#ifndef _OPENMP
            std::fill(rA.begin(), rA.end(), TDataType());
#else
            //omptl::fill(rA.begin(), rA.end(), TDataType());
            ParallelFill(rA.begin(), rA.end(), TDataType());
#endif
        }

        inline static void SetToZero(compressed_matrix<TDataType>& rA)
        {
#ifndef _OPENMP
            std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#else
            //omptl::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
            ParallelFill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#endif
        }

        inline static void SetToZero(VectorType& rX)
        {
#ifndef _OPENMP
            std::fill(rX.begin(), rX.end(), TDataType());
#else
            //omptl::fill(rX.begin(), rX.end(), TDataType());
            ParallelFill(rX.begin(), rX.end(), TDataType());
#endif
        }

        template<class TOtherMatrixType, class TEquationIdVectorType>
        inline static void AssembleLHS(
        MatrixType& A,
        TOtherMatrixType& LHS_Contribution,
        TEquationIdVectorType& EquationId
        )
        {
            unsigned int system_size = A.size1();
            unsigned int local_size = LHS_Contribution.size1();

            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                unsigned int i_global = EquationId[i_local];
                if (i_global < system_size)
                {
                    for (unsigned int j_local = 0; j_local < local_size; j_local++)
                    {
                        unsigned int j_global = EquationId[j_local];
                        if (j_global < system_size)
                            A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                    }
                }
            }
        }


        //        static void GatherLocalValues(Vector& global_indices, )
        //		{
        //		axpy_prod(rA, rX, rY, true);
        //		}




        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.

        virtual std::string Info() const
        {
            return "UBlasSpace";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "UBlasSpace";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
        {
        }

        //***********************************************************************

        inline static bool IsDistributed()
        {
            return false;
        }

        //***********************************************************************

        inline static double GetValue(const VectorType& x, std::size_t I)
        {
            return x[I];
        }
        //***********************************************************************

        static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, double* pValues)
        {
            KRATOS_TRY

            for (std::size_t i = 0; i < IndexArray.size(); i++)
                pValues[i] = x[IndexArray[i]];

            KRATOS_CATCH("")
        }
        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{


        ///@}
        ///@name Private Operators
        ///@{

#ifdef _OPENMP
        //y += A*x in parallel

        static void ParallelProductNoAdd(const MatrixType& A, const VectorType& in, VectorType& out)
        {
            //create partition
            vector<unsigned int> partition;
            unsigned int number_of_threads = omp_get_max_threads();
            unsigned int number_of_initialized_rows = A.filled1() - 1;
            CreatePartition(number_of_threads, number_of_initialized_rows, partition);
            //parallel loop
#pragma omp parallel
            {
                int thread_id = omp_get_thread_num();
                int number_of_rows = partition[thread_id + 1] - partition[thread_id];
                typename compressed_matrix<TDataType>::index_array_type::const_iterator row_iter_begin = A.index1_data().begin() + partition[thread_id];
                typename compressed_matrix<TDataType>::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
                typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;
                //                  typename VectorType::iterator output_vec_begin = out.begin()+partition[thread_id];


                partial_product_no_add(number_of_rows,
                        row_iter_begin,
                        index_2_begin,
                        value_begin,
                        in,
                        partition[thread_id],
                        out
                        );
            }
        }

        static void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
        {
            partitions.resize(number_of_threads + 1);
            int partition_size = number_of_rows / number_of_threads;
            partitions[0] = 0;
            partitions[number_of_threads] = number_of_rows;
            for (unsigned int i = 1; i < number_of_threads; i++)
                partitions[i] = partitions[i - 1] + partition_size;
        }

        template <class TIterartorType>
        static void ParallelFill(TIterartorType Begin, TIterartorType End, TDataType const& Value)
        {
#ifndef _OPENMP
            std::fill(Begin, End, Value);
#else
            int size = End-Begin;
            //#pragma omp parallel for 
            for (int i = 0; i < size; i++)
            {
                *(Begin+i) = Value;
            }
#endif

        }

        /**
         * calculates partial product resetting to Zero the output before
         */
        static void partial_product_no_add(
                int number_of_rows,
                typename compressed_matrix<TDataType>::index_array_type::const_iterator row_begin,
                typename compressed_matrix<TDataType>::index_array_type::const_iterator index2_begin,
                typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin,
                const VectorType& input_vec,
                unsigned int output_begin_index,
                VectorType& output_vec
                //                 typename VectorType::iterator output_vec_begin
                )
        {
            int row_size;
            int kkk = output_begin_index;
            typename MatrixType::index_array_type::const_iterator row_it = row_begin;
            for (int k = 0; k < number_of_rows; k++)
            {
                row_size = *(row_it + 1)-*row_it;
                row_it++;
                TDataType t = TDataType();

                for (int i = 0; i < row_size; i++)
                    t += *value_begin++ * (input_vec[*index2_begin++]);

                output_vec[kkk++] = t;
                //                 *output_vec_begin++ = t;

            }
        }
#endif


        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        UblasSpace & operator=(UblasSpace const& rOther);

        /// Copy constructor.
        UblasSpace(UblasSpace const& rOther);


        ///@}

    }; // Class UblasSpace



    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    //   inline std::istream& operator >> (std::istream& rIStream,
    // 				    UblasSpace& rThis);

    //   /// output stream function
    //   inline std::ostream& operator << (std::ostream& rOStream,
    // 				    const UblasSpace& rThis)
    //     {
    //       rThis.PrintInfo(rOStream);
    //       rOStream << std::endl;
    //       rThis.PrintData(rOStream);

    //       return rOStream;
    //     }
    ///@}


} // namespace Kratos.

#endif // KRATOS_UBLAS_SPACE_H_INCLUDED  defined 


