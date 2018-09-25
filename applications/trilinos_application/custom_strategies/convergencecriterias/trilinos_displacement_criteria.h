//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_TRILINOS_DISPLACEMENT_CRITERIA )
#define  KRATOS_TRILINOS_DISPLACEMENT_CRITERIA

//  System includes

//  External includes

//  Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "Epetra_MpiComm.h"

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

// Short class definition
// Detail class definition

// URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

// URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

// URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

// URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


// URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

// URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

// URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

// URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


template<class TSparseSpace,
         class TDenseSpace
         >
class TrilinosDisplacementCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( TrilinosDisplacementCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}
    ///@name Life Cycle

    ///@{

    // * Constructor.

    TrilinosDisplacementCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRatioTolerance = NewRatioTolerance;
        mAlwaysConvergedNorm = AlwaysConvergedNorm;

        //mActualizeRHSIsNeeded = false;
    }

    // * Destructor.

    virtual ~TrilinosDisplacementCriteria() {}


    ///@}
    ///@name Operators

    ///@{

    // Criterias that need to be called after getting the solution
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            TDataType mFinalCorrectionNorm = TSparseSpace::TwoNorm(Dx);

            TDataType ratio = 0.00;

            mReferenceDispNorm = CalculateReferenceNorm(rDofSet, rModelPart);

            ratio = mFinalCorrectionNorm/mReferenceDispNorm;

            double aaa = SparseSpaceType::Size(Dx);

            double AbsoluteNorm = (mFinalCorrectionNorm/sqrt(aaa));

            if(rModelPart.GetCommunicator().MyPID() == 0) //print performed only by the first processor
                std::cout << "DISPLACEMENT CRITERIA :: obtained tol = " << ratio << ";  expected ratio = " << mRatioTolerance << "absolute tol = " << AbsoluteNorm << std::endl;

            if ( ratio <= mRatioTolerance || AbsoluteNorm<mAlwaysConvergedNorm ) // || (mFinalCorrectionNorm/x.size())<=1e-7)
            {
                if(rModelPart.GetCommunicator().MyPID() == 0) //print performed only by the first processor
                    std::cout <<"convergence is achieved" <<std::endl;

                return true;
            }
            else
            {
                return false;
            }
        }
        else //in this case all the displacements are imposed!
        {
            return true;
        }
    }

    void Initialize(
        ModelPart& r_model_part
    ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )  override {}



    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
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
    TDataType mRatioTolerance;
    TDataType mAlwaysConvergedNorm;

    // Epetra_MpiComm& mrComm;


    TDataType mReferenceDispNorm;
    ///@}
    ///@name Private Operators
    ///@{

    TDataType CalculateReferenceNorm(DofsArrayType& rDofSet, ModelPart& rModelPart)
    {
        TDataType local_ReferenceDispNorm = TDataType();
        TDataType temp;

        const double rank = rModelPart.GetCommunicator().MyPID(); // double because I want to compare with PARTITION_INDEX

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree() && (i_dof->GetSolutionStepValue(PARTITION_INDEX) == rank))
            {
                temp = i_dof->GetSolutionStepValue();
                local_ReferenceDispNorm += temp*temp;
            }
        }

        //perform the sum between all of the nodes
        TDataType ReferenceDispNorm = local_ReferenceDispNorm;
        rModelPart.GetCommunicator().SumAll(ReferenceDispNorm);


        ReferenceDispNorm = std::sqrt(ReferenceDispNorm);
        return ReferenceDispNorm;
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

}; //  Class ClassName

///@}

///@name Type Definitions
///@{


///@}

}  //  namespace Kratos.

#endif //  KRATOS_TRILINOS_DISPLACEMENT_CRITERIA  defined
