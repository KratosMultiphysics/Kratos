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

/* *********************************************************
*
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-11-11 14:03:41 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_TRILINOS_STATIC_SCHEME )
#define  KRATOS_TRILINOS_STATIC_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"
#include "Epetra_Import.h"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

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
class TrilinosResidualBasedIncrementalUpdateStaticScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( TrilinosResidualBasedIncrementalUpdateStaticScheme);

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

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
    TrilinosResidualBasedIncrementalUpdateStaticScheme():
        Scheme<TSparseSpace,TDenseSpace>(),
        mImporterIsInitialized(false)
    {}

    /** Destructor.
    */
    virtual ~TrilinosResidualBasedIncrementalUpdateStaticScheme() {}


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
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY;

        if (!DofImporterIsInitialized())
            this->InitializeDofImporter(rDofSet,Dx);

        int system_size = TSparseSpace::Size1(A);

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp( mpDofImporter->TargetMap() );

        //importing in the new temp vector the values
        int ierr = temp.Import(Dx,*mpDofImporter,Insert);
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

        double* temp_values; //DO NOT make delete of this one!!
        temp.ExtractView( &temp_values );

        Dx.Comm().Barrier();


// ModelPart::NodesContainerType::iterator node_it = r_model_part.Nodes().find(2756);
// std::cout << A.Comm().MyPID() << " node 2756 " << node_it->FastGetSolutionStepValue(PARTITION_INDEX) << " disp_x id " << node_it->pGetDof(DISPLACEMENT_X)->EquationId() << std::endl;

// std::cout << "rank=" << A.Comm().MyPID() << "dof with id 117 "<< rDofSet.find(117) << std::endl;

        //performing the update
        typename DofsArrayType::iterator dof_begin = rDofSet.begin();
        for(unsigned int iii=0; iii<rDofSet.size(); iii++)
        {
            int global_id = (dof_begin+iii)->EquationId();
            if(global_id < system_size)
            {
                if((dof_begin+iii)->IsFree())
                {
                    double aaa = temp[mpDofImporter->TargetMap().LID(global_id)];
                    /*		if(global_id == 117) std::cout << "rank = " << b.Comm().MyPID() << " global id" << global_id << "local num" << iii << " map num " << dof_update_map.LID(global_id) << " value= " << aaa << std::endl;*/
                    (dof_begin+iii)->GetSolutionStepValue() += aaa;
                }
            }
        }

//                        delete [] temp_values;  //deleting this is WRONG! do not do it!!

        KRATOS_CATCH("")
    }



    //void Predict(
    //	const String& ElementGroupName,
    //	DofsArrayType& rDofSet,
    //	TSystemMatrixType& A,
    //	TSystemVectorType& Dx,
    //	TSystemVectorType& b,
    //	ProcessInfo& CurrentProcessInfo
    //	)
    //{
    //	double CurrentTime = CurrentProcessInfo.GetCurrentTime();
    //	double DeltaTime = CurrentProcessInfo.GetDeltaTime();
    //	double OldTime = CurrentTime - DeltaTime;
    //
    //	int i;
    //	typename DofsArrayType::iterator it2;
    //
    //	//predicting variables
    //	for (it2=rDofSet.begin();it2 != rDofSet.end(); ++it2)
    //	{
    //		// N.B. fixed values are not predicted!!
    //		if ( !(*it2)->IsFixed()  )
    //		{
    //			const Dof& X = *(*it2);
    //
    //			mpModel->Value(*it2) =  mpModel->Value(X.GetVariable(), X, OldTime);
    //		}
    //	}
    //
    //}
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce
    the selected time integration scheme. It "asks" the matrix needed to the element and
    performs the operations needed to introduce the seected time integration scheme.

    this function calculates at the same time the contribution to the LHS and to the RHS
    of the system
    */
    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    )
    {
        KRATOS_TRY
        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        (rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    //***************************************************************************
    //***************************************************************************
    virtual void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************
    virtual void Calculate_LHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);
        (rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    /** functions totally analogous to the precedent but applied to
    the "condition" objects
    */
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        (rCurrentCondition)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        (rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);
        KRATOS_CATCH("")
    }

    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        (rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);
        KRATOS_CATCH("")
    }


    /*@} */
    /**@name Operations */
    /*@{ */


    virtual void Clear()
    {
        mpDofImporter.reset();
        mImporterIsInitialized = false;
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    bool DofImporterIsInitialized()
    {
        return mImporterIsInitialized;
    }

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
        boost::shared_ptr<Epetra_Import> pDofImporter( new Epetra_Import(dof_update_map,Dx.Map()) );
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }


    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /// Get pointer Epetra_Import instance that can be used to import values from Dx to the owner of each Dof.
    /**
     * @note Important: always check that the Importer is initialized before calling using
     * DofImporterIsInitialized or initialize it with InitializeDofImporter.
     * @return Importer
     */
    boost::shared_ptr<Epetra_Import> pGetImporter()
    {
        return mpDofImporter;
    }

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

    boost::shared_ptr<Epetra_Import> mpDofImporter;

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

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_STANDARD_STATIC_SCHEME  defined */
