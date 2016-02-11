//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_SMALL_STRAIN_U_PW_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_U_PW_FIC_ELEMENT_H_INCLUDED

/* Project includes */
#include "includes/serializer.h"
#include "custom_elements/small_strain_U_Pw_element.hpp"

namespace Kratos
{

class SmallStrainUPwFICElement : public SmallStrainUPwElement
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUPwFICElement );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallStrainUPwFICElement();
    
    // Constructor 1
    SmallStrainUPwFICElement(IndexType NewId, GeometryType::Pointer pGeometry);
    
    // Constructor 2
    SmallStrainUPwFICElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    
    // Destructor
    virtual ~SmallStrainUPwFICElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
    void Initialize();
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);
    
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    GeometryType::ShapeFunctionsGradientsType mDN_DXContainer; //Contains the matrices with the shape functions derivatives for every integration point
    Vector mdetJContainer;
    Matrix mExtrapolationMatrix;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateExtrapolationMatrix(Matrix& rExtrapolationMatrix, const unsigned int dimension);
    
    
    void InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);
    
    void ExtrapolateShapeFunctionsDerivatives(ElementalVariables& rVariables);
    
    void CalculateKinematics(ElementalVariables& rVariables, unsigned int PointNumber);
    
    void CalculateStrainDerivativeTerm(ElementalVariables& rVariables);
    
    
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables);
    
    void CalculateAndAddStrainDerivativeMatrix(MatrixType& rLeftHandSideMatrix,ElementalVariables& rVariables);
    
    void CalculateAndAddStressDerivativeMatrix(MatrixType& rLeftHandSideMatrix,ElementalVariables& rVariables);
    
    void CalculateStressDerivativeTerm(Matrix& rStressDerivativeTerm, const ElementalVariables& rVariables);
    
    void CalculateAndAddPressureDerivativeMatrix(MatrixType& rLeftHandSideMatrix,ElementalVariables& rVariables);
    
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementalVariables& rVariables);
    
    void CalculateAndAddStrainDerivativeFlow(VectorType& rRightHandSideVector,ElementalVariables& rVariables);
    
    void CalculateAndAddStressDerivativeFlow(VectorType& rRightHandSideVector,ElementalVariables& rVariables);
    
    void CalculateAndAddPressureDerivativeFlow(VectorType& rRightHandSideVector,ElementalVariables& rVariables);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallStrainUPwElement )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallStrainUPwElement )
    }
    
    
}; // Class SmallStrainUPwFICElement

} // namespace Kratos

#endif // KRATOS_SMALL_STRAIN_U_PW_FIC_ELEMENT_H_INCLUDED  defined 
