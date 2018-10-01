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


#ifndef KRATOS_TRILINOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H
#define KRATOS_TRILINOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

// Application includes
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
class TrilinosResidualBasedIncrementalAitkenStaticScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosResidualBasedIncrementalAitkenStaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualBasedIncrementalAitkenStaticScheme);

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

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
        ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>(),
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
    void InitializeSolutionStep(ModelPart &r_model_part,
                                        TSystemMatrixType &A,
                                        TSystemVectorType &Dx,
                                        TSystemVectorType &b) override
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

    void Clear() override
    {
        mpDofImporter.reset();
        mImporterIsInitialized = false;
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
        auto dof_begin = rDofSet.begin();
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


    virtual void InitializeDofImporter(DofsArrayType& rDofSet,
                                       TSystemVectorType& Dx)
    {
        int system_size = TSparseSpace::Size(Dx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        //filling the array with the global ids
        int counter = 0;
        for(auto i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            int id = i_dof->EquationId();
            if( id < system_size )
            {
                index_array[counter] = id;
                counter += 1;
            }
        }

        std::sort(index_array.begin(),index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(),index_array.end());
        index_array.resize(NewEnd-index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        Dx.Comm().SumAll(&tot_update_dofs,&check_size,1);
        KRATOS_ERROR_IF( (check_size < system_size) &&  (Dx.Comm().MyPID() == 0) )
            << "Dof count is not correct. There are less dofs then expected." << std::endl
            << "Expected number of active dofs = " << system_size << " dofs found = " << check_size << std::endl;

        //defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,Dx.Comm() );

        //defining the importer class
        Kratos::shared_ptr<Epetra_Import> pDofImporter = Kratos::make_shared<Epetra_Import>(dof_update_map,Dx.Map());
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    /// Get pointer Epetra_Import instance that can be used to import values from Dx to the owner of each Dof.
    /**
     * @note Important: always check that the Importer is initialized before calling using
     * DofImporterIsInitialized or initialize it with InitializeDofImporter.
     * @return Importer
     */
    Kratos::shared_ptr<Epetra_Import> pGetImporter()
    {
        return mpDofImporter;
    }

    ///@}
    ///@name Protected Inquiry
    ///@{

    bool DofImporterIsInitialized() const
    {
        return mImporterIsInitialized;
    }

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

    bool mImporterIsInitialized = false;

    Kratos::shared_ptr<Epetra_Import> mpDofImporter = nullptr;

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
