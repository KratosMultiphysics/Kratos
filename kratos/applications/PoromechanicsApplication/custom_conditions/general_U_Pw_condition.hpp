//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_GENERAL_U_PW_CONDITION_H_INCLUDED )
#define  KRATOS_GENERAL_U_PW_CONDITION_H_INCLUDED

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

class GeneralUPwCondition : public Condition
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( GeneralUPwCondition );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    GeneralUPwCondition();
    
    // Constructor 1
    GeneralUPwCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    
    // Constructor 2
    GeneralUPwCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~GeneralUPwCondition();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;
 
    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );

    void EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   

    struct ConditionVariables
    {
        //Variables at all integration points
        Matrix NContainer;
        GeometryType::JacobiansType JContainer;

        //Variables at each integration point
        Vector Np; //Contains the shape functions at every node
        double IntegrationCoefficient;
        
        //Imposed condition at all nodes
        Vector ConditionVector;
        
        //FIC condition variables
        double BiotModulusInverse;
        double NewmarkCoefficient;
        double ElementLength;
    };
    
    // Member Variables
    
    IntegrationMethod mThisIntegrationMethod;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
                                bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag);

    void InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateKinematics(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual void CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber);

    virtual void CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight);

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables);

    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables);
    
    virtual void CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }
    
}; // class GeneralUPwCondition.

} // namespace Kratos.

#endif // KRATOS_GENERAL_U_PW_CONDITION_H_INCLUDED defined 
