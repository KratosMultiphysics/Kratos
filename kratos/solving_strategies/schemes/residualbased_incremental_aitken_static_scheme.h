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


#if !defined( KRATOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H_INCLUDED )
#define  KRATOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

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

/// A scheme for the solution of a problem using Aitken iterations.
/** Aitken's method intends to improve convergence by introducing a relaxation factor in the update of the solution after each iteration.
  */
template< class TSparseSpace,class TDenseSpace >
class ResidualBasedIncrementalAitkenStaticScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedIncrementalAitkenStaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedIncrementalAitkenStaticScheme);

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    //typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    //typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** @param DefaultOmega Default relaxation factor to use in the first iteration, where Aitken's factor cannot be computed. Use a value between 0 and 1.
      */
    ResidualBasedIncrementalAitkenStaticScheme(double DefaultOmega):
        mDefaultOmega(DefaultOmega),
        mOldOmega(DefaultOmega)
    {}

    /// Destructor.
    ~ResidualBasedIncrementalAitkenStaticScheme() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Initialize the iteration counter at the begining of each solution step
    /**
      * @param r_model_part The problem's ModelPart
      * @param A System matrix
      * @param Dx Solution vector (containing the increment of the unknowns obtained in the present iteration)
      * @param b Right hand side vector
      */
    void InitializeSolutionStep(ModelPart &r_model_part,
                                        TSystemMatrixType &A,
                                        TSystemVectorType &Dx,
                                        TSystemVectorType &b) override
    {
        BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);
        if (TSparseSpace::Size(mPreviousDx) != TSparseSpace::Size(Dx)) {
            TSparseSpace::Resize(mPreviousDx, TSparseSpace::Size(Dx));
        }
        TSparseSpace::SetToZero(mPreviousDx);
        mIterationCounter = 0;
    }

    /// Increase the iteration counter at the begining of each iteration
    /**
      * @param r_model_part The problem's ModelPart
      * @param A System matrix
      * @param Dx Solution vector (containing the increment of the unknowns obtained in the present iteration)
      * @param b Right hand side vector
      */
    void InitializeNonLinIteration(ModelPart &r_model_part,
                                           TSystemMatrixType &A,
                                           TSystemVectorType &Dx,
                                           TSystemVectorType &b) override
    {
        BaseType::InitializeNonLinIteration(r_model_part,A,Dx,b);
    }


    /// Update the degrees of freedom of the problem using Aitken's accelerator
    /**
      * @param r_model_part The problem's ModelPart
      * @param A System matrix
      * @param Dx Solution vector (containing the increment of the unknowns obtained in the present iteration)
      * @param b Right hand side vector
      */
    void Update(ModelPart &r_model_part,
                        DofsArrayType &rDofSet,
                        TSystemMatrixType &A,
                        TSystemVectorType &Dx,
                        TSystemVectorType &b) override
    {
        mIterationCounter++;

        // Compute relaxation factor
        double Omega;

        if (mIterationCounter > 1)
        {
            double Num = 0.0;
            double Den = 0.0;

            for (unsigned int i = 0; i < Dx.size(); i++)
            {
                double Diff = Dx[i] - mPreviousDx[i];
                Num += mPreviousDx[i] * Diff;
                Den += Diff * Diff;
            }

            Omega = - mOldOmega * Num / Den;
        }
        else
        {
            // Initialize the process using min(DefaultOmega,Omega_from_last_step)
            if (mOldOmega < mDefaultOmega)
                Omega = mOldOmega;
            else
                Omega = mDefaultOmega;
        }

        //KRATOS_WATCH(Omega);

        // Update using relaxation factor
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() += Omega * Dx[i_dof->EquationId()];
            }
        }

        // Store results for next iteration
        boost::numeric::ublas::noalias(mPreviousDx) = Dx;
        mOldOmega = Omega;
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

    /// Keeps track of the iteration number to ensure that Aiken's accelerator is not used until there are enough previous steps to use it.
    unsigned int mIterationCounter;

    /// Initial value for the Aitken's Omega coefficient.
    /** This is the value used to initialize Aitken's method for the first time. For additional solutions after the first step, the minimum between
      * this value and the final value of omega in the previous step (stored in mOldOmega) is used.
      */
    const double mDefaultOmega;

    /// Value of Aitken's Omega during the previous iteration.
    double mOldOmega;

    /// Vector of solution variations obtained in the previous iteration.
    TSystemVectorType mPreviousDx;

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

    /// Assignment operator.
    ResidualBasedIncrementalAitkenStaticScheme& operator=(ResidualBasedIncrementalAitkenStaticScheme const& rOther) {}

    /// Copy constructor.
    ResidualBasedIncrementalAitkenStaticScheme(ResidualBasedIncrementalAitkenStaticScheme const& rOther) {}


    ///@}

}; // Class ResidualBasedIncrementalAitkenStaticScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H_INCLUDED  defined
