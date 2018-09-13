//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_GENERAL_U_PW_DIFF_ORDER_CONDITION_H_INCLUDED )
#define  KRATOS_GENERAL_U_PW_DIFF_ORDER_CONDITION_H_INCLUDED

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

class KRATOS_API(POROMECHANICS_APPLICATION) GeneralUPwDiffOrderCondition : public Condition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( GeneralUPwDiffOrderCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    GeneralUPwDiffOrderCondition();
    
    // Constructor 1
    GeneralUPwDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    GeneralUPwDiffOrderCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~GeneralUPwDiffOrderCondition();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
    void Initialize() override;
 
    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;
    
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo ) override;
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;
    
    void EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo ) override;

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
    
    Geometry< Node<3> >::Pointer mpPressureGeometry;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
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
    
}; // class GeneralUPwDiffOrderCondition.

} // namespace Kratos.

#endif // KRATOS_GENERAL_U_PW_DIFF_ORDER_CONDITION_H_INCLUDED defined 
