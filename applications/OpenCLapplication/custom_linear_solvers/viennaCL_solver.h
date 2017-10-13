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
#define VIENNACL_WITH_OPENCL 1

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
#include "viennacl/linalg/gmres.hpp"

#include <amgcl/amgcl.hpp>
#include "amgcl/level_viennacl.hpp"
#include "amgcl/operations_viennacl.hpp"
#include <amgcl/aggr_plain.hpp>
#include <amgcl/interp_aggr.hpp>


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
    NoPreconditioner, ILU, AMG_DAMPED_JACOBI, AMG_SPAI0
};

class ViennaCLSolver : public IterativeSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>, Preconditioner<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector> >, Reorderer<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector> > >
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ViennaCLSolver
    KRATOS_CLASS_POINTER_DEFINITION( ViennaCLSolver );

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
	mndof = 1;
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
            viennacl::compressed_matrix<scalar_type,1u> gpu_A;
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
                viennacl::linalg::ilut_precond< viennacl::compressed_matrix<scalar_type,1u> > vcl_ilut(gpu_A, viennacl::linalg::ilut_tag(mentries_per_row, mdrop_tolerance));

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
 else if (mpreconditioner_type == AMG_DAMPED_JACOBI)
            {
		std::vector<scalar_type> data_float;
		data_float.resize(rA.value_data().size());
		copy(rA.value_data().begin(),rA.value_data().end(),data_float.begin());
		//for (unsigned int i=0; i!=rA.value_data().size(); i++)
		//	data_float[i]=(float)(rA.value_data()[i]);	

		amgcl::sparse::matrix_map<scalar_type, std::size_t> Aamgcl(
				    rA.size1(), 
				    rA.size2(), 
				    reinterpret_cast<const std::size_t*>(rA.index1_data().begin()),
				    reinterpret_cast<const std::size_t*>(rA.index2_data().begin()),
				    //rA.value_data().begin()
				    data_float.data()
				    );

		viennacl::hyb_matrix<scalar_type> gpu_A;
		viennacl::copy(amgcl::sparse::viennacl_map(Aamgcl), gpu_A);


                typedef amgcl::solver<
				scalar_type, int,
				amgcl::interp::aggregation<amgcl::aggr::plain>,
				amgcl::level::viennacl<amgcl::GPU_MATRIX_HYB, amgcl::relax::damped_jacobi>
				> AMG;
				
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		//std::cout<<"detected dofs per node: "<< mndof <<'\n'; //for debugging purposes
		//clock_t start = std::clock(); //for debugging purposes
		AMG amg((Aamgcl), prm);
		//double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; //for debugging purposes
		//std::cout<<"preconditioner time: "<< duration <<'\n'; //for debugging purposes

		amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg);

				
                if (msolver_type == CG)
                {
                    viennacl::linalg::cg_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }



                if (msolver_type == BiCGStab)
                {
                    viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }
                if (msolver_type == GMRES)
                {
                    viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }

            }
	    else if (mpreconditioner_type == AMG_SPAI0)
            {

		std::vector<scalar_type> data_float;
		data_float.resize(rA.value_data().size());
		copy(rA.value_data().begin(),rA.value_data().end(),data_float.begin());
		
		amgcl::sparse::matrix_map<scalar_type, std::size_t> Aamgcl(
				    rA.size1(), 
				    rA.size2(), 
				    reinterpret_cast<const std::size_t*>(rA.index1_data().begin()),
				    reinterpret_cast<const std::size_t*>(rA.index2_data().begin()),
				    //rA.value_data().begin()
				    data_float.data()
				    );

		viennacl::hyb_matrix<scalar_type> gpu_A;
		viennacl::copy(amgcl::sparse::viennacl_map(Aamgcl), gpu_A);

                typedef amgcl::solver<
				scalar_type, int,
				amgcl::interp::aggregation<amgcl::aggr::plain>,
				amgcl::level::viennacl<amgcl::GPU_MATRIX_HYB, amgcl::relax::spai0>
				> AMG;
				
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		//clock_t start = std::clock(); //for debugging purposes
		AMG amg((Aamgcl), prm);
		//double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; //for debugging purposes
		//std::cout<<"preconditioner time: "<< duration <<'\n'; //for debugging purposes

		amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg);

				
                if (msolver_type == CG)
                {
                    viennacl::linalg::cg_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }



                if (msolver_type == BiCGStab)
                {
                    viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }
                if (msolver_type == GMRES)
                {
                    viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }

            }

            //copy back to CPU
            copy(gpu_X.begin(), gpu_X.end(), rX.begin());
        }
        else if (mprecision == Double)
        {
            //create ViennaCL data structure
            typedef double scalar_type;
            //typedef float scalar_type;
            
            viennacl::vector<scalar_type> gpu_B(rX.size());
            viennacl::vector<scalar_type> gpu_X(rX.size());

            copy(rB.begin(), rB.end(), gpu_B.begin());
	    
            //
            //solve the linear system of equations using ViennaCL's OpenCL implementation
            if (mpreconditioner_type == NoPreconditioner)
            {
		viennacl::compressed_matrix<scalar_type> gpu_A;
            	copy(rA, gpu_A);

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
            }
            else if (mpreconditioner_type == ILU)
            {
		viennacl::compressed_matrix<scalar_type> gpu_A;
            	copy(rA, gpu_A);

                viennacl::linalg::ilut_precond< viennacl::compressed_matrix<scalar_type,1u> > vcl_ilut(gpu_A, viennacl::linalg::ilut_tag(mentries_per_row, mdrop_tolerance));

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
            else if (mpreconditioner_type == AMG_DAMPED_JACOBI)
            {	
		amgcl::sparse::matrix_map<scalar_type, std::size_t> Aamgcl(
				    rA.size1(), 
				    rA.size2(), 
				    reinterpret_cast<const std::size_t*>(rA.index1_data().begin()),
				    reinterpret_cast<const std::size_t*>(rA.index2_data().begin()),
				    rA.value_data().begin()
				    );

		viennacl::hyb_matrix<scalar_type> gpu_A;
		viennacl::copy(amgcl::sparse::viennacl_map(Aamgcl), gpu_A);


                typedef amgcl::solver<
				scalar_type, int,
				amgcl::interp::aggregation<amgcl::aggr::plain>,
				amgcl::level::viennacl<amgcl::GPU_MATRIX_HYB, amgcl::relax::damped_jacobi>
				> AMG;
				
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		//std::cout<<"detected dofs per node: "<< mndof <<'\n'; //for debugging purposes
		//clock_t start = std::clock(); //for debugging purposes
		AMG amg((Aamgcl), prm);
		//double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; //for debugging purposes
		//std::cout<<"preconditioner time: "<< duration <<'\n'; //for debugging purposes

		amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg);

				
                if (msolver_type == CG)
                {
                    viennacl::linalg::cg_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }



                if (msolver_type == BiCGStab)
                {
                    viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }
                if (msolver_type == GMRES)
                {
                    viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }

            }
	    else if (mpreconditioner_type == AMG_SPAI0)
            {
		amgcl::sparse::matrix_map<scalar_type, std::size_t> Aamgcl(
				    rA.size1(), 
				    rA.size2(), 
				    reinterpret_cast<const std::size_t*>(rA.index1_data().begin()),
				    reinterpret_cast<const std::size_t*>(rA.index2_data().begin()),
				    rA.value_data().begin()
				    );

		viennacl::hyb_matrix<scalar_type> gpu_A;
		viennacl::copy(amgcl::sparse::viennacl_map(Aamgcl), gpu_A);

                typedef amgcl::solver<
				scalar_type, int,
				amgcl::interp::aggregation<amgcl::aggr::plain>,
				amgcl::level::viennacl<amgcl::GPU_MATRIX_HYB, amgcl::relax::spai0>
				> AMG;
				
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		//clock_t start = std::clock(); //for debugging purposes
		AMG amg((Aamgcl), prm);
		//double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; //for debugging purposes
		//std::cout<<"preconditioner time: "<< duration <<'\n'; //for debugging purposes

		amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg);

				
                if (msolver_type == CG)
                {
                    viennacl::linalg::cg_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }



                if (msolver_type == BiCGStab)
                {
                    viennacl::linalg::bicgstab_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
                    BaseType::mIterationsNumber = custom_solver.iters();
                    BaseType::mResidualNorm = custom_solver.error();
                }
                if (msolver_type == GMRES)
                {
                    viennacl::linalg::gmres_tag custom_solver(tol, maxIter);
                    gpu_X = viennacl::linalg::solve(gpu_A,
 						gpu_B,
 						custom_solver,
 						amgcl::make_viennacl_precond< viennacl::vector<scalar_type> >(amg));
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

	/** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
        int old_ndof = -1;
		unsigned int old_node_id = rdof_set.begin()->Id();
		int ndof=0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
		{
			
			if(it->EquationId() < rA.size1() )
			{
				unsigned int id = it->Id();
				if(id != old_node_id)
				{
					old_node_id = id;
					if(old_ndof == -1) old_ndof = ndof;
					else if(old_ndof != ndof) //if it is different than the block size is 1
					{
						old_ndof = -1;
						break;
					}
					
					ndof=1;
				}
				else
				{
					ndof++;
				}
			}
		}
		
		if(old_ndof == -1) 
			mndof = 1;
		else
			mndof = ndof;
			
		KRATOS_WATCH(mndof);


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

	int mndof;
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


