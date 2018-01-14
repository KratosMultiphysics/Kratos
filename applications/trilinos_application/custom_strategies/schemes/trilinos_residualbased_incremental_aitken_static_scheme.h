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

#ifndef KRATOS_TRILINOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H
#define KRATOS_TRILINOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

// Application includes
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"


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
class TrilinosResidualBasedIncrementalAitkenStaticScheme : public TrilinosResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosResidualBasedIncrementalAitkenStaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualBasedIncrementalAitkenStaticScheme);

    typedef TrilinosResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef Kratos::shared_ptr < TSystemVectorType > SystemVectorPointerType;

    //typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    //typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** @param DefaultOmega Default relaxation factor to use in the first iteration, where Aitken's factor cannot be computed. Use a value between 0 and 1.
      */
    TrilinosResidualBasedIncrementalAitkenStaticScheme(double DefaultOmega):
        TrilinosResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>(),
        mDefaultOmega(DefaultOmega),
        mOldOmega(DefaultOmega)
    {}

    /// Destructor.
    virtual ~TrilinosResidualBasedIncrementalAitkenStaticScheme() {}


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
    virtual void InitializeSolutionStep(ModelPart &r_model_part,
                                        TSystemMatrixType &A,
                                        TSystemVectorType &Dx,
                                        TSystemVectorType &b)
    {
        BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);
        mIterationCounter = 0;
    }

    /// Increase the iteration counter at the begining of each iteration
    /**
      * @param r_model_part The problem's ModelPart
      * @param A System matrix
      * @param Dx Solution vector (containing the increment of the unknowns obtained in the present iteration)
      * @param b Right hand side vector
      */
    virtual void InitializeNonLinIteration(ModelPart &r_model_part,
                                           TSystemMatrixType &A,
                                           TSystemVectorType &Dx,
                                           TSystemVectorType &b)
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
    virtual void Update(ModelPart &r_model_part,
                        DofsArrayType &rDofSet,
                        TSystemMatrixType &A,
                        TSystemVectorType &Dx,
                        TSystemVectorType &b)
    {
        KRATOS_TRY;

        mIterationCounter++;
        double Omega;

        if (mIterationCounter > 1)
        {
            // Compute relaxation factor
            Omega = this->CalculateOmega(Dx);
        }
        else
        {
            // Initialize the process using min(DefaultOmega,Omega_from_last_step)
            if (mOldOmega < mDefaultOmega)
                Omega = mOldOmega;
            else
                Omega = mDefaultOmega;
        }


//        KRATOS_WATCH(Omega);

        this->UpdateWithRelaxation(rDofSet,Dx,Omega);

        // Store results for next iteration
        SystemVectorPointerType Tmp = SystemVectorPointerType( new TSystemVectorType(Dx));
        mpPreviousDx.swap(Tmp);
        mOldOmega = Omega;

        KRATOS_CATCH("");
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

    /// Calculate the value of the Aitken relaxation factor for the current iteration
    double CalculateOmega(TSystemVectorType& Dx)
    {
        KRATOS_TRY;

        double Num = 0.0;
        double Den = 0.0;

        int MyLength = Dx.MyLength();

        // Safety check, to be sure that we didn't move nodes along processors
        if (MyLength != mpPreviousDx->MyLength())
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Unexpected error in Trilinos Aitken iterations: the Dx vector has a different size than in previous iteration.","");
        }

//        int NumVectors = Dx.NumVectors();

        // These point to data held by the Epetra_Vector, do not delete[] here.
        double* pDxValues;
        int DxLDA;
        double* pOldDxValues;
        int OldDxLDA;

        Dx.ExtractView(&pDxValues,&DxLDA);
        mpPreviousDx->ExtractView(&pOldDxValues,&OldDxLDA);

        for (int i = 0; i < MyLength; i++)
        {
            double Diff = pDxValues[i] - pOldDxValues[i];
            Num += pOldDxValues[i] * Diff;
            Den += Diff * Diff;
        }

        Dx.Comm().Barrier();

        // Sum values across processors
        double SendBuff[2] = {Num,Den};
        double RecvBuff[2];

        Dx.Comm().SumAll(&SendBuff[0],&RecvBuff[0],2);

        Num = RecvBuff[0];
        Den = RecvBuff[1];

        double Omega = - mOldOmega * Num / Den;
        return Omega;

        KRATOS_CATCH("");
    }

    void UpdateWithRelaxation(DofsArrayType& rDofSet,
                              TSystemVectorType& Dx,
                              const double Omega)
    {
        KRATOS_TRY;

        if (!this->DofImporterIsInitialized())
            this->InitializeDofImporter(rDofSet,Dx);

        Kratos::shared_ptr<Epetra_Import> pImporter = this->pGetImporter();

        int system_size = TSparseSpace::Size(Dx);

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp( pImporter->TargetMap() );

        //importing in the new temp vector the values
        int ierr = temp.Import(Dx,*pImporter,Insert);
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

        double* temp_values; //DO NOT make delete of this one!!
        temp.ExtractView( &temp_values );

        Dx.Comm().Barrier();

        //performing the update
        typename DofsArrayType::iterator dof_begin = rDofSet.begin();
        for(unsigned int iii=0; iii<rDofSet.size(); iii++)
        {
            int global_id = (dof_begin+iii)->EquationId();
            if(global_id < system_size)
            {
                double aaa = temp[pImporter->TargetMap().LID(global_id)];
                (dof_begin+iii)->GetSolutionStepValue() += Omega * aaa;
            }
        }

        KRATOS_CATCH("");
    }

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
    SystemVectorPointerType mpPreviousDx;

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
    TrilinosResidualBasedIncrementalAitkenStaticScheme& operator=(TrilinosResidualBasedIncrementalAitkenStaticScheme const& rOther) {}

    /// Copy constructor.
    TrilinosResidualBasedIncrementalAitkenStaticScheme(TrilinosResidualBasedIncrementalAitkenStaticScheme const& rOther) {}


    ///@}

}; // Class TrilinosResidualBasedIncrementalAitkenStaticScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_TRILINOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H
