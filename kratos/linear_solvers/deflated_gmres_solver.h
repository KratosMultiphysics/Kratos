//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_DEFLATED_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_DEFLATED_GMRES_SOLVER_H_INCLUDED
// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
// External includes
// Project includes
#include "includes/define.h"
#include "reorderer.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "linear_solvers/iterative_solver.h"
#include <boost/numeric/ublas/vector.hpp>
#include "utilities/openmp_utils.h"

//#define NO_PRECOND 

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
/** This solver is designed for the solution of mixed U-P problems.
 * It uses a block structure diving the matrix in UU PP UP PU blocks
 * and uses "standard" linear solvers for the different blocks as well as a GMRES for the outer part
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DeflatedGMRESSolver :
    public IterativeSolver<TSparseSpaceType, TDenseSpaceType,TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of DeflatedGMRESSolver
    KRATOS_CLASS_POINTER_DEFINITION (DeflatedGMRESSolver);
    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;    
    typedef typename TSparseSpaceType::VectorType VectorType;
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
    typedef typename TDenseSpaceType::VectorType DenseVectorType;
    typedef std::size_t  SizeType;
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    DeflatedGMRESSolver (typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer pred_solver,
                         double NewMaxTolerance,
                         unsigned int NewMaxIterationsNumber,
                         unsigned int m, unsigned int max_reduced_size
                        ) : BaseType (NewMaxTolerance, NewMaxIterationsNumber)
    {
        //saving the linear solvers to be used in the solution process
        //mpsolver_UU_block = psolver_UU_block;
        //mpsolver_PP_block = psolver_PP_block;
	//this is the solver used at the prediction step before entering the GMRES loop... can be direct or iterative
	mPred_solver = pred_solver;
        mBlocksAreAllocated = false;
        mis_initialized = false;
        mm = m;
	mmax_reduced_size=max_reduced_size;
	KRATOS_WATCH("Quasi-deflated solver created")
	std::cout<<"Krylov space size is"<< mm<<std::endl;
	std::cout<<"Maximum deflated matrix size is"<< mmax_reduced_size<<std::endl;	
	myfile.open("iterations.txt");
    }
    /// Copy constructor.
    DeflatedGMRESSolver (const DeflatedGMRESSolver& Other)
    {
        KRATOS_THROW_ERROR (std::logic_error,"copy constructor not correctly implemented","");
    }
    /// Destructor.
    ~DeflatedGMRESSolver() override {}
    ///@}
    ///@name Operators
    ///@{
    /// Assignment operator.
    DeflatedGMRESSolver& operator= (const DeflatedGMRESSolver& Other)
    {
        return *this;
    }
    ///@}
    ///@name Operations
    ///@{
    /** This function is designed to be called as few times as possible. It creates the data structures
     * that only depend on the connectivity of the matrix (and not on its coefficients)
     * so that the memory can be allocated once and expensive operations can be done only when strictly
     * needed
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void Initialize (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
	if (mBlocksAreAllocated == true)
	{
	    
	    mis_initialized = true;
	}
	else
	{
	  std::cout << "linear solver intialization is deferred to the moment at which blocks are available" << std::endl;
	}
    }
    /** This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void InitializeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {     
        //copy to local matrices
        if (mBlocksAreAllocated == false)
        {
            FillBlockMatrices (true, rA, mK, mG, mD, mS);
            mBlocksAreAllocated = true;
        }
        else
        {
            FillBlockMatrices (false, rA, mK, mG, mD, mS);
            mBlocksAreAllocated = true;
        }
        
        if(mis_initialized == false) this->Initialize(rA,rX,rB);
        
    }

    /** This function actually performs the solution work, eventually taking advantage of what was done before in the
     * Initialize and InitializeSolutionStep functions.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void PerformSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        unsigned int m = mm;
        unsigned int max_iter = BaseType::GetMaxIterationsNumber();
        double tol = BaseType::GetTolerance();
        gmres_solve (rA,rX,rB,m,max_iter,tol);
    } 
    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void FinalizeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        
    }
    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        mK.clear();
        mG.clear();
        mD.clear();
        mS.clear();
        mBlocksAreAllocated = false;        
	mPred_solver->Clear();
	mu.clear();
	mp.clear();
	mru.clear();
	mrp.clear();
        mis_initialized = false;
    }

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rVectorx is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mis_initialized == false)
            this->Initialize (rA,rX,rB);

        this->InitializeSolutionStep (rA,rX,rB);

        this->PerformSolutionStep (rA,rX,rB);

        this->FinalizeSolutionStep (rA,rX,rB);

        return false;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rVectorx is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     @param rB. Right hand side vector.
    */
    bool Solve (SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
    }

    /** Eigenvalue and eigenvector solve method for derived eigensolvers */
     void Solve (SparseMatrixType& K,
                         SparseMatrixType& M,
                         DenseVectorType& Eigenvalues,
                         DenseMatrixType& Eigenvectors) override
    {}

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
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
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        //count pressure dofs
        unsigned int n_pressure_dofs = 0;
        unsigned int tot_active_dofs = 0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
		{
		  if (it->EquationId() < rA.size1())
		  {
		      tot_active_dofs += 1;
		      if (it->GetVariable().Key() == PRESSURE)
			  n_pressure_dofs += 1;
		  }
		}

        if (tot_active_dofs != rA.size1() )
            KRATOS_THROW_ERROR (std::logic_error,"total system size does not coincide with the free dof map","");

        //resize arrays as needed
        mpressure_indices.resize (n_pressure_dofs,false);

        unsigned int other_dof_size = tot_active_dofs - n_pressure_dofs;
        mother_indices.resize (other_dof_size,false);
        mglobal_to_local_indexing.resize (tot_active_dofs,false);
        mis_pressure_block.resize (tot_active_dofs,false);
        //construct aux_lists as needed
        //"other_counter[i]" i will contain the position in the global system of the i-th NON-pressure node
        //"pressure_counter[i]" will contain the in the global system of the i-th NON-pressure node
        //
        //mglobal_to_local_indexing[i] will contain the position in the local blocks of the
        unsigned int pressure_counter = 0;
        unsigned int other_counter = 0;
        unsigned int global_pos = 0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
        {
            if (it->EquationId() < rA.size1())
            {
                if (it->GetVariable().Key() == PRESSURE)
                {
                    mpressure_indices[pressure_counter] = global_pos;
                    mglobal_to_local_indexing[global_pos] = pressure_counter;
                    mis_pressure_block[global_pos] = true;
                    pressure_counter++;
                }
                else
                {
                    mother_indices[other_counter] = global_pos;
                    mglobal_to_local_indexing[global_pos] = other_counter;
                    mis_pressure_block[global_pos] = false;
                    other_counter++;
                }
                global_pos++;
            }
        }
    }
/*
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
        //count pressure dofs
        unsigned int n_pressure_dofs = 0;
        unsigned int tot_active_dofs = 0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
//             if (it->IsFixed() != true)
            {
                tot_active_dofs += 1;
                if (it->GetVariable().Key() == PRESSURE)
                    n_pressure_dofs += 1;
            }
	//KRATOS_WATCH(rA.size1())
	//KRATOS_WATCH(tot_active_dofs)
//         if (tot_active_dofs != rA.size1() )
//             KRATOS_THROW_ERROR (std::logic_error,"total system size does not coincide with the free dof map","");

        //resize arrays as needed
        mpressure_indices.resize (n_pressure_dofs,false);

        unsigned int other_dof_size = tot_active_dofs - n_pressure_dofs;
        mother_indices.resize (other_dof_size,false);
        mglobal_to_local_indexing.resize (tot_active_dofs,false);
        mis_pressure_block.resize (tot_active_dofs,false);
        //construct aux_lists as needed
        //"other_counter[i]" i will contain the position in the global system of the i-th NON-pressure node
        //"pressure_counter[i]" will contain the in the global system of the i-th NON-pressure node
        //
        //mglobal_to_local_indexing[i] will contain the position in the local blocks of the
        unsigned int pressure_counter = 0;
        unsigned int other_counter = 0;
        unsigned int global_pos = 0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
        {
//             if (it->IsFixed() != true)
//             {
                if (it->GetVariable().Key() == PRESSURE)
                {
                    mpressure_indices[pressure_counter] = global_pos;
                    mglobal_to_local_indexing[global_pos] = pressure_counter;
                    mis_pressure_block[global_pos] = true;
                    pressure_counter++;
                }
                else
                {
                    mother_indices[other_counter] = global_pos;
                    mglobal_to_local_indexing[global_pos] = other_counter;
                    mis_pressure_block[global_pos] = false;
                    other_counter++;
                }
                global_pos++;
//             }
        }
    }
*/
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
    std::string Info() const override
    {
        return "Linear solver";
    }
    /// Print information about this object.
    void PrintInfo (std::ostream& rOStream) const override
    {
        rOStream << "Linear solver";
    }
    /// Print object's data.
    void PrintData (std::ostream& rOStream) const override
    {
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
    ///this function generates the subblocks of matrix A
    ///as A = ( K G ) u
    ///       ( D S ) p
    /// subblocks are allocated or nor depending on the value of "need_allocation"
    void FillBlockMatrices (bool need_allocation, SparseMatrixType& rA, SparseMatrixType& K, SparseMatrixType& G, SparseMatrixType& D, SparseMatrixType& S )
    {
        KRATOS_TRY
	KRATOS_WATCH("FILLING BLOCK MATRICES")
        //get access to A data
        const std::size_t* index1 = rA.index1_data().begin();
        const std::size_t* index2 = rA.index2_data().begin();
        const double*	   values = rA.value_data().begin();

        SparseMatrixType L(mpressure_indices.size(),mpressure_indices.size() );

        if (need_allocation == true)
        {
            K.clear();
            G.clear();
            D.clear();
            S.clear();
            L.clear();

            //do allocation
            K.resize (mother_indices.size()   ,mother_indices.size() );
            G.resize (mother_indices.size()   ,mpressure_indices.size() );
            D.resize (mpressure_indices.size(),mother_indices.size() );
            S.resize (mpressure_indices.size(),mpressure_indices.size() );
	    
	    mrp.resize(mpressure_indices.size() );
	    mru.resize(mother_indices.size() );
	    mp.resize(mpressure_indices.size());
	    mu.resize(mother_indices.size());


            //KRATOS_WATCH (mglobal_to_local_indexing);
            //allocate the blocks by push_back
            for (unsigned int i=0; i<rA.size1(); i++)
            {
                unsigned int row_begin = index1[i];
                unsigned int row_end   = index1[i+1];
                unsigned int local_row_id = mglobal_to_local_indexing[i];

                if ( mis_pressure_block[i] == false) //either K or G
                {
                    for (unsigned int j=row_begin; j<row_end; j++)
                    {
                        unsigned int col_index = index2[j];
                        double value = values[j];
                        unsigned int local_col_id = mglobal_to_local_indexing[col_index];
                        if (mis_pressure_block[col_index] == false) //K block
                            K.push_back ( local_row_id, local_col_id, value);
                        else //G block
                            G.push_back ( local_row_id, local_col_id, value);
                    }
                }
                else //either D or S
                {
                    for (unsigned int j=row_begin; j<row_end; j++)
                    {
                        unsigned int col_index = index2[j];
                        double value = values[j];
                        unsigned int local_col_id = mglobal_to_local_indexing[col_index];
                        if (mis_pressure_block[col_index] == false) //D block
                            D.push_back ( local_row_id, local_col_id, value);
                        else //S block
                            L.push_back ( local_row_id, local_col_id, value);
                    }
                }
            }
	
            S = L;           
            VectorType diagK (mother_indices.size() );
            ComputeDiagonalByLumping (K,diagK);

           

        }
        else //allocation is not needed so only do copying
        {
            for (unsigned int i=0; i<rA.size1(); i++)
            {
                unsigned int row_begin = index1[i];
                unsigned int row_end   = index1[i+1];
                unsigned int local_row_id =  mglobal_to_local_indexing[i];
                if ( mis_pressure_block[i] == false ) //either K or G
                {
                    for (unsigned int j=row_begin; j<row_end; j++)
                    {
                        unsigned int col_index = index2[j];
                        double value = values[j];
                        unsigned int local_col_id = mglobal_to_local_indexing[col_index];
                        if (mis_pressure_block[col_index] == false) //K block
                            K( local_row_id, local_col_id) = value;
                        else //G block
                            G( local_row_id, local_col_id) = value;
                    }
                }
                else //either D or S
                {
                    for (unsigned int j=row_begin; j<row_end; j++)
                    {
                        unsigned int col_index = index2[j];
                        double value = values[j];
                        unsigned int local_col_id = mglobal_to_local_indexing[col_index];
                        if (mis_pressure_block[col_index] == false) //D block
                            D( local_row_id, local_col_id) = value;
                        else //S block
                            L( local_row_id, local_col_id) = value;
                    }
                }
            }

	S = L;

            VectorType diagK (mother_indices.size() );
            ComputeDiagonalByLumping (K,diagK);

         

        }





        KRATOS_CATCH ("")
    }
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
    /// A counted pointer to the reorderer object.
    //typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpsolver_UU_block;
    //typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpsolver_PP_block;
    typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mPred_solver;

    unsigned int mm;
    unsigned int mmax_reduced_size;
    bool mBlocksAreAllocated;
    bool mis_initialized;
    DenseVector<unsigned int> mpressure_indices;
    DenseVector<unsigned int> mother_indices;
    DenseVector<int> mglobal_to_local_indexing;
    DenseVector<int> mis_pressure_block;
    SparseMatrixType mK;
    SparseMatrixType mG;
    SparseMatrixType mD;
    SparseMatrixType mS;
    
    VectorType mrp;
	VectorType mru;
	VectorType mp;
	VectorType mu;
	
    std::ofstream myfile;
    ///@}
    ///@name Private Operators
    ///@{
    inline void GeneratePlaneRotation (const double &dx, const double &dy, double &cs, double &sn)
    {
        if (dy == 0.0)
        {
            cs = 1.0;
            sn = 0.0;
        }
        else if (dx == 0.0)
        {
            cs = 0.0;
            sn = 1.0;
        }
        else
        {
            const double rnorm = 1.0/sqrt (dx*dx + dy*dy);
            cs = fabs (dx) * rnorm;
            sn = cs * dy / dx;
        }
    }

    inline void ApplyPlaneRotation (double &dx, double &dy, const double &cs, const double &sn)
    {
        double temp  =  cs * dx + sn * dy;
        dy = cs * dy - sn * dx;
        dx = temp;
    }

    void Update (VectorType& y, VectorType& x, int k, Matrix& h, VectorType& s, std::vector< VectorType >& V)
    {
        for (unsigned int i=0; i<s.size(); i++)
            y[i] = s[i];
        /*		for(unsigned int i=s.size(); i<y.size(); i++)
        			y[i] = 0.0;*/
        // Backsolve:
        for (int i = k; i >= 0; --i)
        {
            y (i) /= h (i,i);
            for (int j = i - 1; j >= 0; --j)
                y (j) -= h (j,i) * y (i);
        }
        //create new search dir
        for (int j = 0; j <= k; ++j)
            TSparseSpaceType::UnaliasedAdd (x, y[j], V[j]);   // x +=  y(j)* V[j];
    }

    int gmres_solve ( SparseMatrixType& A,
                      VectorType& x,
                      const VectorType& b,
                      unsigned int& m,
                      unsigned int& max_iter,
                      double& tol)
    {
        const unsigned int dim = A.size1();
        if (m == 0)
            KRATOS_THROW_ERROR (std::logic_error,"the dimension of the GMRES krylov space can not be set to zero. Please change the value of m","")
            if (m > max_iter)
                m = max_iter;
	//KRATOS_WATCH("Krylov space size")
	//KRATOS_WATCH(m)
        VectorType s (m+1), sn (m+1), w (dim), r (dim), y (m+1);
	/////////THINGS NECESSARY FOR DEFLATION///////////////////////////////////////////////////////////////////////////
	SparseMatrixType S_deflated;											//
	//THIS ww is the matrix W in vector form									//
	std::vector<int> W;												//	
															//
	DeflationUtils::ConstructW(mmax_reduced_size, mS, W, S_deflated);						//
	DeflationUtils::FillDeflatedMatrix(mS, W, S_deflated);								//	
	int red_dim=S_deflated.size1();											//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	VectorType  cs (m+1);
        Matrix  H (m+1, m+1);
        int restart = 0;
	int p_dim=mS.size1();
	VectorType output(p_dim);//WT*lambda	
	VectorType temp (dim,0.0);

        double normb = TSparseSpaceType::TwoNorm (b);
	/*KRATOS_WATCH(normb);*/
        if (normb < 1e-16) //ARBITRARY SMALL NUMBER!
        {
            normb = 1e-16;
        }
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//get the residual
        //r = b - Ax
        TSparseSpaceType::Mult (A,x,r);
        TSparseSpaceType::ScaleAndAdd (1.00, b, -1.00, r); //r = b - r	
	
	//CHECKING IF THE MATRIX IS INVERTIBLE!!!! If it is not (i.e. if S_deflated*Identity=0, we add a number to diagonal)
	//CheckDeflatedMatrix(S_deflated);	

#ifndef NO_PRECOND	
	KRATOS_WATCH("SOLVING DEFLATED PRESSURE")
	SolveDeflatedPressure( output, r, S_deflated, W);		
	KRATOS_WATCH("SOLVED DEFLATED PRESSURE")
	//update x: by modifying its part corresponding to pressure	
	WritePPart (temp, output);		
	TSparseSpaceType::ScaleAndAdd(1.00, temp, 1.00, x);	
	
	TSparseSpaceType::Mult (A,x,r);
        TSparseSpaceType::ScaleAndAdd (1.00, b, -1.00, r); //r = b - r
	//KRATOS_WATCH(r)
#endif
	
	
        const double rel_tol = tol*normb;
        double beta = TSparseSpaceType::TwoNorm (r);
        if (beta <= rel_tol)   //finalize!
        {
            tol = beta / normb;
            max_iter = 0;
            return 0;
        }
        unsigned int j;
        int err = 0;
        std::vector< VectorType > V (m+1);
        for (j = 0; j <= m; ++j)
            V[j].resize (dim,false);
        j = 1;
        while (j <= max_iter)
        {
            TSparseSpaceType::Assign (V[0], 1.0/beta, r); //V[0] = r /(T)beta;
            TSparseSpaceType::SetToZero (s);
            s[0] = beta;
            for (unsigned int i = 0; (i < m) && (j <= max_iter); ++i, ++j)
            {
                TSparseSpaceType::Mult (A,V[i],w); //w = A*V[i];	

                for (unsigned int k = 0; k <= i; k++)
                {
                    H (k, i) = TSparseSpaceType::Dot (V[k], w);
                    w -= H (k, i) * V[k];
                }
#ifndef NO_PRECOND
		Modify_w(  w,  W, dim, p_dim, red_dim);
#endif

                const double normw = TSparseSpaceType::TwoNorm (w);
                H (i+1, i) = normw;                
                // This breakdown is a good one ...
                if (normw == 0)
                    TSparseSpaceType::Copy (V[i+1], w); //V[i+1] = w;
                else
                    TSparseSpaceType::Assign (V[i+1], 1.0/normw, w); //V[i+1] = w / normw;
                for (unsigned int k = 0; k < i; k++)
                    ApplyPlaneRotation (H (k,i), H (k+1,i), cs (k), sn (k) );
                GeneratePlaneRotation (H (i,i), H (i+1,i), cs (i), sn (i) );
                ApplyPlaneRotation (H (i,i), H (i+1,i), cs (i), sn (i) );
                ApplyPlaneRotation (s (i), s (i+1), cs (i), sn (i) );
                beta = fabs (s (i+1) );
		std::cout << "iter = " <<  j << "  estimated res ratio = " << beta << std::endl;
                //KRATOS_WATCH (beta);
                if (beta <= rel_tol)
                {
                    this->Update (y, x, i, H, s, V);
		   //WRITE THE NUMBER OF ITERATIONS INTO A FILE
		   myfile <<j<<"\n";
		   

                    return 0;
                }
		//IF WE SURPASS THE MAX ITERATION NUMBER WE WILL ALSO PRINT IT TO FILE
		else if (j>=max_iter)
			myfile <<j<<"\n";
            }
	    
            this->Update (y,x, m - 1, H, s, V);	 
	    
            //r = b - Ax
            TSparseSpaceType::Mult (A,x,r);
	    TSparseSpaceType::ScaleAndAdd (1.00, b, -1.00, r); //r = b - r	    
            beta = TSparseSpaceType::TwoNorm (r);
            
	    std::cout << "number of iterations at convergence = " << j << std::endl;
            if (beta < rel_tol)
            {
                return 0;
            }
            ++restart;
        }
        err = 1;
        return err;
    }


    void CheckDeflatedMatrix(SparseMatrixType& S_deflated)
    {
    std::size_t reduced_size = S_deflated.size1();

    VectorType identity(reduced_size,1.0);
    VectorType res(reduced_size,0.0);
	TSparseSpaceType::Mult (S_deflated,identity,res);

KRATOS_WATCH(res)
KRATOS_WATCH(norm_2(res))	
    }

    //FUNCTION THAT SOLVES THE SYSTEM WTLW*lambda=WT*r or WTLW*d_lambda=WT*w.. in the first case output=W*lambda, in the second ouput=W*d_lambda
    //void SolveDeflatedPressure( VectorType& output, VectorType& r, SparseMatrixType& S_deflated, std::vector<int>& W, LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType>& Factorization)
    void SolveDeflatedPressure( VectorType& output, VectorType& r, SparseMatrixType& S_deflated, std::vector<int>& W)
    {
	///////////////////////////////////////////////////////////
	// put here deflation i.e. solve for WTLWp=WTr fixing w
	//ww is a deflation matrix W written in a vector format	
	//extracted the part of the residual corresponding to the pressure - r_p
	VectorType rp;	
	//get the lower part of the residual vector, corresponding to pressure dofs
	GetPPart (r, rp);
	std::size_t reduced_size = S_deflated.size1();
	//std::size_t full_size = mS.size1();
	VectorType WT_rp(reduced_size), lambda(reduced_size);
	
	//w_T_r is the residual multiplied by the WT
	DeflationUtils::ApplyWtranspose(W, rp, WT_rp);	
	mPred_solver->Solve(S_deflated, lambda, WT_rp);	
	KRATOS_WATCH(norm_2(lambda) );

	DeflationUtils::ApplyW(W, lambda, output);       
	///////////////////////////////////////////////////////////////////////////////////////////////

    }
    //makes w orthogonal to W
    void Modify_w( VectorType& w, std::vector<int>& W, std::size_t full_glob_size, std::size_t full_size, std::size_t reduced_size)
    {
    //std::size_t full_glob_size = A.size1();
    //std::size_t reduced_size = S_deflated.size1();
    //std::size_t full_size = mS.size1();
	

    VectorType wp (full_size); 
    VectorType WT_wp(reduced_size); 
    VectorType W_WT_wp (full_size); 
    VectorType temp (full_glob_size);     

    GetPPart(w,wp);
    VectorType ModulusSquared(reduced_size);  
    VectorType identity(full_size,1.0);
    DeflationUtils::ApplyWtranspose(W, identity, ModulusSquared);

    DeflationUtils::ApplyWtranspose(W, wp, WT_wp);
    //scale down Wt_wp with the squared modulus
	
    for(unsigned int i=0; i<reduced_size; i++)
	{
	WT_wp[i] /= ModulusSquared[i];
	
	}

    DeflationUtils::ApplyW(W, WT_wp, W_WT_wp);      
        
    wp-=W_WT_wp;
    WritePPart(w,wp);
    }

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to u-dofs
    void GetUPart (const VectorType& rtot, VectorType& ru)
    {
        if (ru.size() != mother_indices.size() )
            ru.resize (mother_indices.size(), false);
        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(ru.size()); i++)
            ru[i] = rtot[mother_indices[i]];
    }

    //this function extracts from a vector which has the size of the
    //overall r, the part that corresponds to p-dofs
    void GetPPart (const VectorType& rtot, VectorType& rp)
    {
        if (rp.size() != mpressure_indices.size() )
            rp.resize (mpressure_indices.size(), false);
        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(rp.size()); i++)
            rp[i] = rtot[mpressure_indices[i]];
    }

    void WriteUPart (VectorType& rtot, const VectorType& ru)
    {
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(ru.size()); i++)
            rtot[mother_indices[i]] = ru[i];
    }

    void WritePPart (VectorType& rtot, const VectorType& rp)
    {
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(rp.size()); i++)
            rtot[mpressure_indices[i]] = rp[i];
    }

    void ComputeDiagonalByLumping (SparseMatrixType& A,VectorType& diagA)
    {
        if (diagA.size() != A.size1() )
            diagA.resize (A.size1() );
        //get access to A data
        const std::size_t* index1 = A.index1_data().begin();
//        const std::size_t* index2 = A.index2_data().begin();
        const double*	   values = A.value_data().begin();

        #pragma omp parallel for
        for (int i=0; i< static_cast<int>(A.size1()); i++)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i+1];
            double temp = 0.0;
            for (unsigned int j=row_begin; j<row_end; j++)
                temp += values[j]*values[j];

            diagA[i] = sqrt(temp);
        }
    }

    double CheckMatrix (SparseMatrixType& A)
    {
        //get access to A data
        const std::size_t* index1 = A.index1_data().begin();
        const std::size_t* index2 = A.index2_data().begin();
        const double*	   values = A.value_data().begin();
        double norm = 0.0;
        for (unsigned int i=0; i<A.size1(); i++)
        {
            unsigned int row_begin = index1[i];
            unsigned int row_end   = index1[i+1];
            if (row_end - row_begin == 0)
                std::cout << "line " << i << " has no elements" << std::endl;
            //KRATOS_THROW_ERROR(std::logic_error, "line found with no entries on line ",i)
            for (unsigned int j=row_begin; j<row_end; j++)
            {
                if (index2[j]>A.size2() )
                    KRATOS_THROW_ERROR (std::logic_error, "array above size of A","")
                    norm += values[j]*values[j];
            }
        }
        return sqrt (norm);
    }

    

    /// Helper function for Sytem matrix functions
    void SortCols (
        std::vector<unsigned int>& ColList,
        std::size_t& NumCols)
    {
        bool swap = true;
        unsigned int d = NumCols;
        int temp;
        while ( swap || d > 1 )
        {
            swap = false;
            d = (d+1) /2;
            for ( unsigned int i=0; i< (NumCols - d); i++)
                if ( ColList[i+d] < ColList[i] )
                {
                    temp = ColList[i+d];
                    ColList[i+d] = ColList[i];
                    ColList[i] = temp;
                    swap = true;
                }
        }
    }

   
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
}; // Class DeflatedGMRESSolver
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  DeflatedGMRESSolver<TSparseSpaceType, TDenseSpaceType,TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}
/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DeflatedGMRESSolver<TSparseSpaceType, TDenseSpaceType,TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo (rOStream);
    rOStream << std::endl;
    rThis.PrintData (rOStream);
    return rOStream;
}
///@}
}  // namespace Kratos.
#endif // KRATOS_DEFLATED_GMRES_SOLVER_H_INCLUDED  defined
