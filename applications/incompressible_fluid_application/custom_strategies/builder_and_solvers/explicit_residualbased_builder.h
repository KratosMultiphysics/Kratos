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
*   Last Modified by:    $Author: kkazem $
*   Date:                $Date: 2008-11-19 16:12:53 $
*   Revision:            $Revision: 1.10 $
*
* ***********************************************************/


#if !defined(KRATOS_EXPLICIT_RESIDUAL_BASED_BUILDER )
#define    KRATOS_EXPLICIT_RESIDUAL_BASED_BUILDER


/* System includes */
#include <set>
// #include <omp.h>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#else
#include <ctime>
#endif


/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

#include "includes/model_part.h"
#include "containers/array_1d.h"
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

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

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
         class TDenseSpace , //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ExplicitResidualBasedBuilder
    : public ResidualBasedEliminationBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ExplicitResidualBasedBuilder );


    typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ExplicitResidualBasedBuilder(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedEliminationBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
    {

        std::cout << "using the ExplicitResidualBasedBuilder builder and solver " << std::endl;

    }


    /** Destructor.
    */
    virtual ~ExplicitResidualBasedBuilder() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    //**************************************************************************
    //**************************************************************************

    //**************************************************************************
    //**************************************************************************
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        ModelPart::ElementsContainerType::iterator elem_bg = r_model_part.ElementsBegin();
        int n_elems = r_model_part.Elements().size();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        #pragma omp parallel for firstprivate(n_elems, elem_bg)
        for( int ii=0; ii<n_elems; ++ii)
        {
            //calculate min_dt
            ModelPart::ElementsContainerType::iterator it = elem_bg + ii;

            Element::GeometryType& geom = it->GetGeometry();
            double air_water = it->GetValue(IS_WATER_ELEMENT);

            unsigned int nodes_num = geom.size();
            unsigned int dim = it->GetGeometry().WorkingSpaceDimension();

            //calculate elemental Right Hand Side Contribution
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
            Element::EquationIdVectorType EquationId;
            pScheme->Calculate_RHS_Contribution(*(it.base()),RHS_Contribution,EquationId,CurrentProcessInfo);

            //add RHS_Elemental to its nodes

            unsigned int type_ind = dim+1;

            unsigned int rhs_size = RHS_Contribution.size();
            unsigned int air_water_size = type_ind*nodes_num;
            if(rhs_size != air_water_size)
                type_ind = dim;


            for (unsigned int i = 0; i <geom.size(); i++)
            {
                unsigned int index = i*type_ind;

                geom[i].SetLock();
                array_1d<double,3>& node_rhs_vel = geom[i].FastGetSolutionStepValue(RHS);
                double& node_rhs_water_p = geom[i].FastGetSolutionStepValue(RHS_WATER);
                double& node_rhs_air_p = geom[i].FastGetSolutionStepValue(RHS_AIR);


                //add velocity rhs
                for(unsigned int kk=0; kk<dim; kk++)
                    node_rhs_vel[kk] += RHS_Contribution[index+kk];

                //add pressure rhs
                if( nodes_num == (dim+1) )
                {

                    if( air_water== 1.0)
                        node_rhs_water_p += RHS_Contribution[index+dim];

                    else if( air_water== 0.0)
                        node_rhs_air_p += RHS_Contribution[index+dim];


// 				    else
// 					      KRATOS_WATCH("55555555555555555555 neither air nor water!!! 5555555555555555555");
                }
                geom[i].UnSetLock();
            }

            // loop for the rest of shell nodes
            if(nodes_num == dim)
            {
                WeakPointerVector< Node < 3 > >& neighb = it->GetValue(NEIGHBOUR_NODES);
                unsigned int ngh_num=0;

                for (unsigned int ind = 0; ind < 3; ind++)
                {

                    if (neighb[ind].Id() != geom[ind].Id())
                    {
                        unsigned int ngh_index = (3 + ngh_num)*3 ;

                        neighb[ind].SetLock();
                        array_1d<double,3>& ngh_rhs_vel = neighb[ind].FastGetSolutionStepValue(RHS);//deine temlate dim

                        for(unsigned int kk=0; kk<dim; kk++)
                            ngh_rhs_vel[kk] += RHS_Contribution[ngh_index+kk];

                        neighb[ind].UnSetLock();
                        ngh_num++;

                    }

                }
// KRATOS_WATCH("INSIDE EXPLICIT RESIDUALBASED BUILDER FILLING NEIGHBOR RHS");
            }

        }
        KRATOS_WATCH("INSIDE EXPLICIT RESIDUALBASED BUILDER FILLING AFTER ADDING ELEMENT");

        // assemble all conditions
// 			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
// 			{
// 				//calculate elemental contribution
// 				pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
//
// 			    Condition::GeometryType& geom = (*it)->GetGeometry();
// 			    unsigned int nodes_num = geom.size();
// 			    unsigned int dim = (*it)->GetGeometry().WorkingSpaceDimension();
//
// 				for (unsigned int i = 0; i <geom.size(); i++)
// 				    {
// 					unsigned int index = i*dim;
//
// 				      array_1d<double,3>& node_rhs_vel = geom[i].FastGetSolutionStepValue(RHS);//deine temlate dim
//
//
// 				      //add velocity rhs
// 				      for(unsigned int kk=0; kk<dim; kk++)
// 					    node_rhs_vel[kk] += RHS_Contribution[index+kk];
//
// 				    }
//
// 				//assemble the elemental contribution
// 				//AssembleRHS(b,RHS_Contribution,EquationId);
// 			}


// 			KRATOS_WATCH("44444444444444444444444");

// 			RHS_Contribution.resize(0,false);

        // assemble all conditions
// 			for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
// 			{
// 				//calculate elemental contribution
// 				pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);
//
// 				Condition::GeometryType& geom = (*it)->GetGeometry();
// 				unsigned int nodes_num = geom.size();
// 				unsigned int dim = (*it)->GetGeometry().WorkingSpaceDimension();
//
// 				    for (unsigned int i = 0; i <geom.size(); i++)
// 					{
// 					    unsigned int index = i*dim;
//
// 					  array_1d<double,3>& node_rhs_vel = geom[i].FastGetSolutionStepValue(RHS);//deine temlate dim
//
//
// 					  //add velocity rhs
// 					  for(unsigned int kk=0; kk<dim; kk++)
// 						node_rhs_vel[kk] += RHS_Contribution[index+kk];
//
// 					}
//
// 				//assemble the elemental contribution
// 				//AssembleRHS(b,RHS_Contribution,EquationId);
// 			}



// 			  #ifdef _OPENMP
// 			  double stop_prod = omp_get_wtime();
// 			  std::cout << "Time for calculating Calculate_Elements_RHS_and_Add  = " << stop_prod - start_prod << std::endl;
// 			  #endif

        //conditions are calculated serial
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
        Element::EquationIdVectorType EquationId;

        for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);

            if(RHS_Contribution.size() != 0)
            {
                Condition::GeometryType& geom = (*it)->GetGeometry();
                //unsigned int nodes_num = geom.size();
                unsigned int dim = (*it)->GetGeometry().WorkingSpaceDimension();

                for (unsigned int i = 0; i <geom.size(); i++)
                {
                    unsigned int index = i*dim;

                    array_1d<double,3>& node_rhs_vel = geom[i].FastGetSolutionStepValue(RHS);//deine temlate dim


                    //add velocity rhs
                    for(unsigned int kk=0; kk<dim; kk++)
                        node_rhs_vel[kk] += RHS_Contribution[index+kk];

                }
            }
            /*KRATOS_WATCH(RHS_Contribution);*/
            //assemble the elemental contribution
            //AssembleRHS(b,RHS_Contribution,EquationId);
        }
        KRATOS_WATCH("END OF EXPLICIT RESIDUALBASED BUILDER ");

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
    }


    //**************************************************************************
    //**************************************************************************
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        //refresh RHS to have the correct reactions

    }

    //**************************************************************************
    //**************************************************************************
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {}

    //**************************************************************************
    //**************************************************************************
    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {}

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
    */
    void Clear()
    {
// 			this->mDofSet = DofsArrayType();

// 			if(this->mpReactionsVector != NULL)
// 				TSparseSpace::Clear( (this->mpReactionsVector) );
// 			this->mReactionsVector = TSystemVectorType();

        if (this->GetEchoLevel()>0)
        {

            KRATOS_WATCH("ExplicitResidualBasedBuilder Clear Function called");
        }
    }
    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme,
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY
        KRATOS_WATCH("Explicit ResizeAndInitializeVectors");

        if(pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0,0) );
            pA.swap(pNewA);
        }
        if(pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0) );
            pDx.swap(pNewDx);
        }
        if(pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0) );
            pb.swap(pNewb);
        }
        if(BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0) );
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& A  = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b  = *pb;


        A.resize(1,1,false);
        Dx.resize(1,false);
        b.resize(1,false);

        //resizing the system vectors and matrix


        KRATOS_CATCH("")

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
    //******************************************************************************************
    //******************************************************************************************
    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(int i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }

    /*@} */
    /**@name Member Variables */
    /*@{ */

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

}; /* Class ResidualBasedEliminationBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

