//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_U_PW_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

// Application includes
#include "custom_utilities/element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwElement : public Element
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwElement );
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwElement(IndexType NewId = 0) : Element( NewId ) {}

    /// Constructor using an array of nodes
    UPwElement(IndexType NewId, const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}
    
    /// Constructor using Geometry
    UPwElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element( NewId, pGeometry ) {}
    
    /// Constructor using Properties
    UPwElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element( NewId, pGeometry, pProperties ) 
    {
        mThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    }

    /// Copy Constructor
    UPwElement(UPwElement const& rOther) : Element(rOther), mThisIntegrationMethod(rOther.mThisIntegrationMethod), mConstitutiveLawVector(rOther.mConstitutiveLawVector) {}

    /// Destructor
    virtual ~UPwElement() {}

    /// Assignment operator.
    UPwElement & operator=(UPwElement const& rOther)
    {
        Element::operator=(rOther);

        mThisIntegrationMethod = rOther.mThisIntegrationMethod;

        mConstitutiveLawVector.clear();
        mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size());

        for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
        {
            mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
        }

        return *this;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;
    
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;
    
    int Check(const ProcessInfo& rCurrentProcessInfo);
    
    void Initialize();
    
    void GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );
    
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo );
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo );
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
    
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void GetValuesVector(Vector& rValues, int Step = 0);
    
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0);
    
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo);
    
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo );
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    /// Member Variables
    
    GeometryData::IntegrationMethod mThisIntegrationMethod;
    
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    virtual void CaculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo );
    
    virtual void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    virtual void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );
    
    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight);
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }


}; // Class UPwElement

} // namespace Kratos

#endif // KRATOS_U_PW_ELEMENT_H_INCLUDED  defined 
