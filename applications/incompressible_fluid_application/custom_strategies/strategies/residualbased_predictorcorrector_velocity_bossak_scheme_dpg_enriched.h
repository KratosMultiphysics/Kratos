/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
 *   Last Modified by:    $Author: Kazem $
 *   Date:                $Date: 2008-07-25 14:48:17 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_DPG_ENRICHED )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_DPG_ENRICHED


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/coordinate_transformation_utilities.h"
/* NOTE: the include below uses an ugly relative path because this file is also used within the trilinos application.
 * The correct way to do this would be to either formalize the dependence (so that the trilinos application explicitly
 * includes this folder in its include path) or invert it (the incompressible fluid depending on the trilinos application),
 * which is probably the better option for Kratos. Until either of those happen, we need this patch. JC
 */
#include "../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"


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
class ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched : public ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched);

    typedef ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace> BaseType;

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
    ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched
           (double NewAlphaBossak, double MoveMeshStrategy,unsigned int DomainSize)
        : ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace, TDenseSpace>(NewAlphaBossak,MoveMeshStrategy),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0)
    {

        std::cout << "using the ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched" << std::endl;
    }

    /** Destructor.
     */
    virtual ~ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched()
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

        virtual void Update(ModelPart& r_model_part,
                            DofsArrayType& rDofSet,
                            TSystemMatrixType& A,
                            TSystemVectorType& Dv,
                            TSystemVectorType& b) override
        {
            KRATOS_TRY;

            mRotationTool.RotateVelocities(r_model_part);

            BaseType::BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            mRotationTool.RecoverVelocities(r_model_part);

            BaseType::AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            KRATOS_CATCH("")
        }
    //*************************************************************************************
    //*************************************************************************************

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);

        //loop over elements to copy teh enriched pressure to AUX
	ModelPart::ElementsContainerType::iterator elem_bg = r_model_part.ElementsBegin();
	int n_elems = r_model_part.Elements().size();
	for ( int jj=0; jj<n_elems; ++jj)
	      {
		const double enriche_pr = (elem_bg + jj)->GetValue(PRESSUREAUX);
	        (elem_bg + jj)->SetValue(AUX_INDEX,enriche_pr);
	      }

    }
    //***************************************************************************

    /** this function is designed to be called in the builder and solver
    to introduce
    the selected time integration scheme. It "asks" the matrix needed to
    the element and
    performs the operations needed to introduce the seected time
    integration scheme.

      this function calculates at the same time the contribution to the
    LHS and to the RHS
      of the system
     */
    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY
        int k = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
        //basic operations for the element considered
        (rCurrentElement)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

//std::cout << rCurrentElement->Id() << " RHS = " << RHS_Contribution << std::endl;
        (rCurrentElement)->CalculateMassMatrix(BaseType::mMass[k], CurrentProcessInfo);
        (rCurrentElement)->CalculateLocalVelocityContribution(BaseType::mDamp[k], RHS_Contribution, CurrentProcessInfo);

        (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (statics is already included)

        this->AddDynamicsToLHS(LHS_Contribution, BaseType::mDamp[k], BaseType::mMass[k], CurrentProcessInfo);
        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, BaseType::mDamp[k], BaseType::mMass[k], CurrentProcessInfo);

	//In case of enriched element save the last row of LHS and the last column of RHS and condense the matrices
	const unsigned int Dim = (rCurrentElement)->GetGeometry().WorkingSpaceDimension();
	const unsigned int NumberOfNodes = (rCurrentElement)->GetGeometry().size();
	const unsigned int enriched_size = (Dim + 1) * NumberOfNodes + 1;

        if( LHS_Contribution.size1() == enriched_size )
	{
	  //save the last row
	  Vector& enriched_vec = (rCurrentElement)->GetValue(GAPS);
	  enriched_vec.resize(enriched_size + 1, false );

	  for(unsigned int ii = 0; ii < enriched_size; ++ii)
	    enriched_vec[ii] = LHS_Contribution(enriched_size-1,ii);
	  //add RHS contribution
	  enriched_vec[enriched_size] = RHS_Contribution[enriched_size-1];

	  //Condense
	  const unsigned int local_size = (Dim + 1) * NumberOfNodes;
	 // LocalSystemMatrixType aux_LHS;
      boost::numeric::ublas::bounded_matrix <double, 16,16> aux_LHS;
	  //LocalSystemVectorType aux_RHS;
	  array_1d<double,16>  aux_RHS;
	  //aux_LHS.resize(local_size,local_size,false);
	  //aux_RHS.resize(local_size,false);

	  double inv_D = 0.0;
	  if( LHS_Contribution(enriched_size-1,enriched_size-1) != 0.0)
	    inv_D = 1.0 / LHS_Contribution(enriched_size-1,enriched_size-1);
	  else
	  {
          KRATOS_ERROR << "!!!!ENRICHED ELEMENT MUST HAVE NON-ZERO DIAGONAL MEMBER!!!!" << std::endl
                       << "Id: " << rCurrentElement->Id() << " Nodes:" << std::endl
                       << (rCurrentElement)->GetGeometry()[0] << std::endl
                       << (rCurrentElement)->GetGeometry()[1] << std::endl
                       << (rCurrentElement)->GetGeometry()[2] << std::endl
                       << (rCurrentElement)->GetGeometry()[3] << std::endl
                       << "LHS contribution: " << LHS_Contribution << std::endl;
	  }

	  double inv_D_f_2 = RHS_Contribution[enriched_size-1] * inv_D;

	  for(unsigned int ii = 0; ii < local_size; ++ii)
	   {
	     //Condensed RHS
	     int last_index = enriched_size-1;
	     aux_RHS[ii] = RHS_Contribution[ii] - LHS_Contribution(ii,last_index) * inv_D_f_2;

	     for(unsigned int jj = 0; jj < local_size; ++jj)
		aux_LHS(ii,jj) = LHS_Contribution(ii,jj) - LHS_Contribution(ii,last_index) * inv_D * LHS_Contribution(last_index,jj);

	   }

	  //Replace
	  LHS_Contribution.resize(local_size,local_size,false);
	  RHS_Contribution.resize(local_size,false);

	  LHS_Contribution = aux_LHS;
	  RHS_Contribution = aux_RHS;
	}

	// If there is a slip condition, apply it on a rotated system of coordinates
	mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

	mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

        KRATOS_CATCH("")

    }
    //*************************************************************************************
    //*************************************************************************************
    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        int k = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        (rCurrentElement)->CalculateMassMatrix(BaseType::mMass[k], CurrentProcessInfo);

        (rCurrentElement)->CalculateLocalVelocityContribution(BaseType::mDamp[k], RHS_Contribution, CurrentProcessInfo);

        (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, BaseType::mDamp[k], BaseType::mMass[k], CurrentProcessInfo);

	//In case of enriched element save the last row of LHS and the last column of RHS and condense the matrices
	const unsigned int Dim = (rCurrentElement)->GetGeometry().WorkingSpaceDimension();
	const unsigned int NumberOfNodes = (rCurrentElement)->GetGeometry().size();
	const unsigned int enriched_size = (Dim + 1) * NumberOfNodes + 1;

        if( RHS_Contribution.size() == enriched_size )
	{
	  const unsigned int local_size = enriched_size-1;
	  LocalSystemVectorType aux_RHS;
	  aux_RHS.resize(local_size,false);

	  for(unsigned int ii = 0; ii < local_size; ++ii)
	   {
	     //copy RHS
	     aux_RHS[ii] = RHS_Contribution[ii];
	   }
	  //Replace
	  RHS_Contribution.resize(local_size,false);

	  RHS_Contribution = aux_RHS;

	 }
	  // If there is a slip condition, apply it on a rotated system of coordinates
	  mRotationTool.Rotate(RHS_Contribution,rCurrentElement->GetGeometry());
	  mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());

    }
    //*************************************************************************************
    //*************************************************************************************
        /** functions totally analogous to the precedent but applied to
        the "condition" objects
         */
        virtual void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                            LocalSystemMatrixType& LHS_Contribution,
                                                            LocalSystemVectorType& RHS_Contribution,
                                                            Element::EquationIdVectorType& EquationId,
                                                            ProcessInfo& CurrentProcessInfo) override
        {
            KRATOS_TRY
            int k = OpenMPUtils::ThisThread();

            //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
            (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->CalculateMassMatrix(this->mMass[k], CurrentProcessInfo);
            //(rCurrentCondition)->CalculateDampingMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(this->mDamp[k], RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);


            this->AddDynamicsToLHS(LHS_Contribution, this->mDamp[k], this->mMass[k], CurrentProcessInfo);

            this->AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, this->mDamp[k], this->mMass[k], CurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }

        virtual void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                                          LocalSystemVectorType& RHS_Contribution,
                                                          Element::EquationIdVectorType& EquationId,
                                                          ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            int k = OpenMPUtils::ThisThread();

            //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
            //Initializing the non linear iteration for the current condition
            (rCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentCondition)->CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
            (rCurrentCondition)->CalculateMassMatrix(this->mMass[k],rCurrentProcessInfo);
            //(rCurrentCondition)->CalculateDampingMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(this->mDamp[k], RHS_Contribution,rCurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId,rCurrentProcessInfo);

            //adding the dynamic contributions (static is already included)
            this->AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, this->mDamp[k], this->mMass[k],rCurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("");
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
        CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool;

    /*@} */
    /**@name Member Variables */
    /*@{ */
    /*		Matrix mMass;
                    Matrix mDamp;

                    Vector mvel;
                    Vector macc;
                    Vector maccold;

                    DofsVectorType mElementalDofList;
     */

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

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */

