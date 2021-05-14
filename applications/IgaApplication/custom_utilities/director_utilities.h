//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_DIRECTOR_UTILITIES_H_INCLUDED )
#define  KRATOS_DIRECTOR_UTILITIES_H_INCLUDED

// System includes
#include "includes/define.h"

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    class KRATOS_API(IGA_APPLICATION) DirectorUtilities
    {
    public:
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

        typedef typename SparseSpaceType::MatrixType SparseMatrixType;
        typedef typename LocalSpaceType::VectorType DenseVectorType;

        /// Constructor
        DirectorUtilities(
            ModelPart & rModelPart,
            Parameters JsonParameters);

        void ComputeDirectors();

    private:
        ModelPart& mrModelPart;

        const Parameters mParameters;

    }; // Class DirectorUtilities

}  // namespace Kratos.

#endif // KRATOS_DIRECTOR_UTILITIES_H_INCLUDED  defined