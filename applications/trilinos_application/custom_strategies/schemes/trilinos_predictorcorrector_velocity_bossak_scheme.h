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


#if !defined(KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME )
#define  KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME


/* System includes */


/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "../../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

namespace Kratos
{




/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

  This class provides the implementation of the basic tasks that are needed by the solution strategy.
  It is intended to be the place for tailoring the solution strategies to problem specific tasks.

        Detail class definition.

          \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

                \URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

                  \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

                        \URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


                                \URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

                                  \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

                                        \URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

                                          \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class TrilinosPredictorCorrectorVelocityBossakScheme
    : public ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace, TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(TrilinosPredictorCorrectorVelocityBossakScheme );

    typedef ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */

    TrilinosPredictorCorrectorVelocityBossakScheme(double NewAlphaBossak, double MoveMeshStrategy)
        : ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy)
    {

        std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme" << std::endl;
    }

    /** Destructor.
     */
    virtual ~TrilinosPredictorCorrectorVelocityBossakScheme()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    /**
            Performing the update of the solution.
     */

    //***************************************************************************

    virtual void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dv,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY

        BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

        this->AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    void BasicUpdateOperations(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY;

        if ( !this->DofImporterIsInitialized() )
            this->InitializeDofImporter(rDofSet,Dx);

        int system_size = TSparseSpace::Size(Dx);

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp( mpDofImporter->TargetMap() );

        //importing in the new temp vector the values
        int ierr = temp.Import(Dx,*mpDofImporter,Insert) ;
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

        double* temp_values;
        temp.ExtractView( &temp_values );

        b.Comm().Barrier();

        //performing the update
        typename DofsArrayType::iterator dof_begin = rDofSet.begin();
        for(unsigned int iii=0; iii<rDofSet.size(); iii++)
        {
            int global_id = (dof_begin+iii)->EquationId();
            if(global_id < system_size)
            {
                double aaa = temp[mpDofImporter->TargetMap().LID(global_id)];
                (dof_begin+iii)->GetSolutionStepValue() += aaa;
            }
        }

        //removing unnecessary memory
//			delete [] temp_values; //DO NOT DELETE THIS!!

        KRATOS_CATCH("")
    }

    bool DofImporterIsInitialized()
    {
        return mImporterIsInitialized;
    }

    virtual void Clear()
    {
        mpDofImporter.reset();
        mImporterIsInitialized = false;
    }


    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    virtual void InitializeDofImporter(DofsArrayType& rDofSet,
                                       TSystemVectorType& Dx)
    {
        int system_size = TSparseSpace::Size(Dx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        //filling the array with the global ids
        int counter = 0;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
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
        if ( (check_size < system_size) &&  (Dx.Comm().MyPID() == 0) )
        {
            std::stringstream Msg;
            Msg << "Dof count is not correct. There are less dofs then expected." << std::endl;
            Msg << "Expected number of active dofs = " << system_size << " dofs found = " << check_size << std::endl;
            KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"")
        }

        //defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,Dx.Comm() );

        //defining the importer class
        Kratos::shared_ptr<Epetra_Import> pDofImporter = Kratos::make_shared<Epetra_Import>(dof_update_map,Dx.Map());
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    bool mImporterIsInitialized;
    
    Kratos::shared_ptr<Epetra_Import> mpDofImporter;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME  defined */



