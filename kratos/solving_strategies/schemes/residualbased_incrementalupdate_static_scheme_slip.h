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


#ifndef KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SCHEME_SLIP_H
#define KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SCHEME_SLIP_H

/* System includes */


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
//#include "includes/variables.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "utilities/coordinate_transformation_utilities.h"

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

/// Scheme for the solution of problems involving a slip condition.
/** This Scheme is a reimplementation of ResidualBasedIncrementalUpdateStaticScheme that can be used to
  * apply slip conditions along the edges of the model. The problem is solved using rotated coordinate axes on
  * the nodes where this condition will be applied, with the first coordinate of the rotated system aligned with the
  * normal to the boundary on each of these nodes.
  */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ResidualBasedIncrementalUpdateStaticSchemeSlip : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticSchemeSlip);

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> RotationToolType;
    typedef typename CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double>::Pointer RotationToolPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The pseudo static scheme (parameters)
     * @param ThisParameters Configuration parameters
     */
    explicit ResidualBasedIncrementalUpdateStaticSchemeSlip(Parameters ThisParameters)
        : BaseType(),
          mpRotationTool(Kratos::make_shared<RotationToolType>(ThisParameters["domain_size"].GetInt(),ThisParameters["block_size"].GetInt(),IS_STRUCTURE,0.0))
    {}

    /// Constructor.
    /** @param DomainSize Number of spatial dimensions (2 or 3).
      * @param BlockSize Number of matrix and vector rows associated to each node. Only the first DomainSize rows will be rotated.
      */
    explicit ResidualBasedIncrementalUpdateStaticSchemeSlip(unsigned int DomainSize,
                                                   unsigned int BlockSize):
        BaseType(),
        mpRotationTool(Kratos::make_shared<RotationToolType>(DomainSize,BlockSize,IS_STRUCTURE,0.0))
    {}

    /// Constructor providing a custom rotation tool.

    /** @param pRotationTool a pointer to the helper class that manages DOF rotation.
      */
    explicit ResidualBasedIncrementalUpdateStaticSchemeSlip(RotationToolPointerType pRotationTool):
        BaseType(),
        mpRotationTool(pRotationTool)
    {}

    /// Destructor.
    ~ResidualBasedIncrementalUpdateStaticSchemeSlip() override{}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Update the degrees of freedom after a solution iteration.
    void Update(ModelPart& r_model_part,
                        DofsArrayType& rDofSet,
                        TSystemMatrixType& A,
                        TSystemVectorType& Dx,
                        TSystemVectorType& b) override
    {
        KRATOS_TRY;

        mpRotationTool->RotateVelocities(r_model_part);

        BaseType::Update(r_model_part,rDofSet,A,Dx,b);

        mpRotationTool->RecoverVelocities(r_model_part);

        KRATOS_CATCH("");
    }


    /// Obtain an element's local contribution to the system and apply slip conditions if needed.
    void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                              LocalSystemMatrixType& LHS_Contribution,
                                              LocalSystemVectorType& RHS_Contribution,
                                              Element::EquationIdVectorType& EquationId,
                                              ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        BaseType::CalculateSystemContributions(rCurrentElement,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

        mpRotationTool->Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
        mpRotationTool->ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

        KRATOS_CATCH("");
    }

    /// Obtain an element's local contribution to the RHS and apply slip conditions if needed.
    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                            LocalSystemVectorType& RHS_Contribution,
                                            Element::EquationIdVectorType& EquationId,
                                            ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        BaseType::Calculate_RHS_Contribution(rCurrentElement,RHS_Contribution,EquationId,CurrentProcessInfo);

        mpRotationTool->Rotate(RHS_Contribution,rCurrentElement->GetGeometry());
        mpRotationTool->ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());

        KRATOS_CATCH("");
    }

    /// Obtain an element's local contribution to the system matrix and apply slip conditions if needed.
    void Calculate_LHS_Contribution(Element::Pointer rCurrentElement,
                                            LocalSystemMatrixType& LHS_Contribution,
                                            Element::EquationIdVectorType& EquationId,
                                            ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        BaseType::Calculate_LHS_Contribution(rCurrentElement,LHS_Contribution,EquationId,CurrentProcessInfo);

        LocalSystemVectorType Temp = ZeroVector(LHS_Contribution.size1());
        mpRotationTool->Rotate(LHS_Contribution,Temp,rCurrentElement->GetGeometry());
        mpRotationTool->ApplySlipCondition(LHS_Contribution,Temp,rCurrentElement->GetGeometry());

        KRATOS_CATCH("");
    }


    /// Obtain a condition's local contribution to the system and apply slip conditions if needed.
    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                        LocalSystemMatrixType& LHS_Contribution,
                                                        LocalSystemVectorType& RHS_Contribution,
                                                        Element::EquationIdVectorType& EquationId,
                                                        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        BaseType::Condition_CalculateSystemContributions(rCurrentCondition,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

        mpRotationTool->Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
        mpRotationTool->ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

        KRATOS_CATCH("");
    }

    /// Obtain a condition's local contribution to the RHS and apply slip conditions if needed.
    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                                      LocalSystemVectorType& RHS_Contribution,
                                                      Element::EquationIdVectorType& EquationId,
                                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        BaseType::Condition_Calculate_RHS_Contribution(rCurrentCondition,RHS_Contribution,EquationId,CurrentProcessInfo);

        mpRotationTool->Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
        mpRotationTool->ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());

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

    /// Rotation tool instance
    RotationToolPointerType mpRotationTool;

    ///@}
    ///@name Serialization
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

}; // class

///@}
///@name Type Definitions
///@{

///@}

///@} // group

}  // namespace Kratos

#endif // KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SCHEME_SLIP_H
