//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined(KRATOS_GENERAL_U_PL_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_GENERAL_U_PL_DIFF_ORDER_CONDITION_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) GeneralUPlDiffOrderCondition : public Condition
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( GeneralUPlDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    GeneralUPlDiffOrderCondition();

    // Constructor 1
    GeneralUPlDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor 2
    GeneralUPlDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~GeneralUPlDiffOrderCondition();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rConditionDofList,const ProcessInfo& rCurrentProcessInfo ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo ) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ConditionVariables
    {
        //Variables at all integration points
        Matrix NuContainer;
        Matrix NpContainer;
        GeometryType::JacobiansType JContainer;

        //Variables at each integration point
        Vector Nu; //Contains the displacement shape functions at every node
        Vector Np; //Contains the pressure shape functions at every node
        double IntegrationCoefficient;

        //Imposed condition at all nodes
        Vector ConditionVector;
    };

    // Member Variables

    IntegrationMethod mThisIntegrationMethod;

    Geometry< Node >::Pointer mpPressureGeometry;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag);

    void InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateKinematics(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual void CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables);

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables);

    virtual void CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }

}; // class GeneralUPlDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_GENERAL_U_PL_DIFF_ORDER_CONDITION_H_INCLUDED defined
