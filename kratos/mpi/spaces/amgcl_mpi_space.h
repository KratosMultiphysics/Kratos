//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef KRATOS_AMGCL_MPI_SPACE_H_INCLUDED
#define KRATOS_AMGCL_MPI_SPACE_H_INCLUDED

// System includes

// External includes
#include "amgcl/mpi/distributed_matrix.hpp"

// Project includes
#include "includes/define.h"
#include "mpi/utilities/amgcl_mpi_dof_updater.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<class TMatrixType, class TVectorType>
class AmgclMPISpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AmgclMPISpace
    KRATOS_CLASS_POINTER_DEFINITION(AmgclMPISpace);

    typedef double DataType;

    typedef TMatrixType MatrixType;

    typedef TVectorType VectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef typename Kratos::shared_ptr< TMatrixType > MatrixPointerType;
    typedef typename Kratos::shared_ptr< TVectorType > VectorPointerType;

    typedef AmgclMPIDofUpdater< AmgclMPISpace<TMatrixType,TVectorType> > DofUpdaterType;
    typedef typename DofUpdaterType::UniquePointer DofUpdaterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AmgclMPISpace(){}

    /// Destructor.
    virtual ~AmgclMPISpace(){}

    ///@}
    ///@name Operations
    ///@{

    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static VectorPointerType CreateEmptyVectorPointer()
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// return size of vector rV
    static IndexType Size(VectorType const& rV)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// return number of rows of rM
    static IndexType Size1(MatrixType const& rM)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// return number of columns of rM
    static IndexType Size2(MatrixType const& rM)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// rY = rX
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// rY = rX
    static void Copy(VectorType const& rX, VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// rX * rY
    static double Dot(VectorType& rX, VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /// ||rX||2
    static double TwoNorm(VectorType const& rX)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static void Mult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    // rY = rAT * rX
    static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static void InplaceMult(VectorType& rX, const double A)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;
    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;
    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    // rZ = (A * rX) + (B * rY)
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    // rY = (A * rX) + (B * rY)
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static void Clear(MatrixPointerType& pA)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static void Clear(VectorPointerType& pX)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    template<class TOtherMatrixType>
    static void SetToZero(TOtherMatrixType& rA)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static void SetToZero(VectorType& rX)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static constexpr bool IsDistributed()
    {
        return true;
    }

    //***********************************************************************

    static double GetValue(const VectorType& x, std::size_t I)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }
    //***********************************************************************

    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, double* pValues)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(const char* pFileName, /*const*/ TOtherMatrixType& rM, const bool Symmetric)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    template< class VectorType >
    static bool WriteMatrixMarketVector(const char* pFileName, const VectorType& rV)
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    static DofUpdaterPointerType CreateDofUpdater()
    {
        DofUpdaterType tmp;
        return tmp.Create();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
std::stringstream buffer;
    buffer << "AmgclMPISpace" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "AmgclMPISpace";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

private:
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AmgclMPISpace& operator=(AmgclMPISpace const& rOther) = delete;

    /// Copy constructor.
    AmgclMPISpace(AmgclMPISpace const& rOther) = delete;

    ///@}

}; // Class AmgclMPISpace

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_AMGCL_MPI_SPACE_H_INCLUDED  defined


