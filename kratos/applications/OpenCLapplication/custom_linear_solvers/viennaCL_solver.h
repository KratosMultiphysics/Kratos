/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
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


#if !defined(KRATOS_GPU_VIENNACL_SOLVER_H_INCLUDED )
#define  KRATOS_GPU_VIENNACL_SOLVER_H_INCLUDED

//#define  VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK

// System includes 

//#include <ctime>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"

// Must be set if you want to use ViennaCL algorithms on ublas objects
#define VIENNACL_HAVE_UBLAS 1

// uncomment to enable experimental double precision support with ATI Stream SDK:
//#define VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK

//
// ViennaCL includes
//
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
//#include "viennacl/linalg/bicgstab_tuned.hpp"
#include "viennacl/linalg/gmres.hpp"


namespace Kratos
{

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
    enum OpenCLPrecision
    {
        Single, Double
    };

    enum OpenCLSolverType
    {
        CG, BiCGStab, GMRES
    };

    enum OpenCLPreconditionerType
    {
        NoPreconditioner, ILU
    };

    class ViennaCLSolver : public IterativeSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>, Preconditioner<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector> >, Reorderer<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector> > >
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of ViennaCLSolver
        typedef boost::shared_ptr<ViennaCLSolver> Pointer;

        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

        typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

        typedef Reorderer<SparseSpaceType, DenseSpaceType > ReordererType;

        typedef IterativeSolver<SparseSpaceType, DenseSpaceType, Preconditioner<SparseSpaceType, DenseSpaceType>, ReordererType > BaseType;

        typedef SparseSpaceType::MatrixType SparseMatrixType;

        typedef SparseSpaceType::VectorType VectorType;

        typedef DenseSpaceType::MatrixType DenseMatrixType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        ViennaCLSolver()
        {
        }

        ViennaCLSolver(double NewTolerance,
                unsigned int NewMaxIterationsNumber,
                OpenCLPrecision precision,
                OpenCLSolverType solver_type,
                OpenCLPreconditionerType preconditioner_type

                ) : BaseType(NewTolerance, NewMaxIterationsNumber)
        {
            havePreconditioner = false;
            maxIter = NewMaxIterationsNumber;
            tol = NewTolerance;

            mprecision = precision;
            msolver_type = solver_type;
            mpreconditioner_type = preconditioner_type;

            mentries_per_row = 10;
            mdrop_tolerance = 1e-3;
        }

        /// Copy constructor.

        ViennaCLSolver(const ViennaCLSolver& Other) : BaseType(Other)
        {
        }

        /// Destructor.

        virtual ~ViennaCLSolver()
        {
        }


        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.

        ViennaCLSolver & operator=(const ViennaCLSolver& Other)
        {
            BaseType::operator=(Other);
            return *this;
        }


        ///@}
        ///@name Operations
        ///@{

        void SetILUEntriesPerRow(unsigned int entries)
        {
            mentries_per_row = entries;
        }

        void SetILUDropTolerance(double tol)
        {
            mdrop_tolerance = tol;
        }

        /** Normal solve method.
            Solves the linear system Ax=b and puts the result on SystemVector& rX.
            rX is also th initial guess for iterative methods.
            @param rA. System matrix
            @param rX. Solution vector. it's also the initial
            guess for iterative linear solvers.
            @param rB. Right hand side vector.
         */
        bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
        {
            if (IsNotConsistent(rA, rX, rB))
                return false;

            if (mprecision == Single)
            {
                //create ViennaCL data structure
                typedef float scalar_type;
                viennacl::compressed_matrix<scalar_type,8u> gpu_A;
                viennacl::vector<scalar_type> gpu_B(rX.size());
                viennacl::vector<scalar_type> gpu_X(rX.size());

                copy(rB.begin(), rB.end(), gpu_B.begin());
                copy(rA, gpu_A);

                //solve the linear system of equations using ViennaCL's OpenCL implementation
                if (mpreconditioner_type == NoPreconditioner)
                {
                    if (msolver_type == CG)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == BiCGStab)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == GMRES)
                    {
                        viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                }
                else if (mpreconditioner_type == ILU)
                {
                    viennacl::linalg::ilut_precond< viennacl::compressed_matrix<scalar_type,8u> > vcl_ilut(gpu_A, viennacl::linalg::ilut_tag(mentries_per_row, mdrop_tolerance));

                    if (msolver_type == CG)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver, vcl_ilut);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == BiCGStab)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver, vcl_ilut);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == GMRES)
                    {
                        viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver, vcl_ilut);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }

                }

                //copy back to CPU
                copy(gpu_X.begin(), gpu_X.end(), rX.begin());
            } else if (mprecision == Double)
            {
                //create ViennaCL data structure
                typedef double scalar_type;
                viennacl::compressed_matrix<scalar_type,8u> gpu_A;
                viennacl::vector<scalar_type> gpu_B(rX.size());
                viennacl::vector<scalar_type> gpu_X(rX.size());

                copy(rB.begin(), rB.end(), gpu_B.begin());
                copy(rA, gpu_A);
                //
                //solve the linear system of equations using ViennaCL's OpenCL implementation
                if (mpreconditioner_type == NoPreconditioner)
                {
                    if (msolver_type == CG)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == BiCGStab)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
//                        gpu_X = solve_tuned(gpu_A, gpu_B, custom_solver);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver);
                       BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == GMRES)
                    {
                        viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                } else if (mpreconditioner_type == ILU)
                {
                    viennacl::linalg::ilut_precond< viennacl::compressed_matrix<scalar_type,8u> > vcl_ilut(gpu_A, viennacl::linalg::ilut_tag(mentries_per_row, mdrop_tolerance));

                    if (msolver_type == CG)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver, vcl_ilut);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == BiCGStab)
                    {
                        viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver, vcl_ilut);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }
                    if (msolver_type == GMRES)
                    {
                        viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                        gpu_X = solve(gpu_A, gpu_B, custom_solver, vcl_ilut);
                        BaseType::mIterationsNumber = custom_solver.iters();
                        BaseType::mResidualNorm = custom_solver.error();
                    }

                }
                //
                // 	    //copy back to CPU
                copy(gpu_X.begin(), gpu_X.end(), rX.begin());
            }



            bool is_solved = true;
            return is_solved;
        }

        /** Multi solve method for solving a set of linear systems with same coefficient matrix.
            Solves the linear system Ax=b and puts the result on SystemVector& rX.
            rX is also th initial guess for iterative methods.
            @param rA. System matrix
            @param rX. Solution vector. it's also the initial
            guess for iterative linear solvers.
            @param rB. Right hand side vector.
         */
        bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
        {
            bool is_solved = false;

            return is_solved;
        }

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Return information about this object.

        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "GPU ViennaCL Biconjugate gradient stabilized linear solver ";
            return buffer.str();
        }

        /// Print information about this object.

        void PrintInfo(std::ostream& OStream) const
        {
            OStream << "GPU ViennaCL Biconjugate gradient stabilized linear solver with ";
        }

        /// Print object's data.

        void PrintData(std::ostream& OStream) const
        {
            BaseType::PrintData(OStream);
        }


        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{
        //This var control preconditioner
        bool havePreconditioner;
        size_t maxIter;
        double tol;

        unsigned int mentries_per_row;
        double mdrop_tolerance;

        OpenCLPrecision mprecision;
        OpenCLSolverType msolver_type;
        OpenCLPreconditionerType mpreconditioner_type;


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


        ///@}

    }; // Class ViennaCLSolver

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& IStream,
            ViennaCLSolver& rThis)
    {
        return IStream;
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& OStream,
            const ViennaCLSolver& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);

        return OStream;
    }
    ///@}


} // namespace Kratos.

#endif // KRATOS_GPU_VIENNACL_SOLVER_H_INCLUDED  defined 


