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


#if !defined(KRATOS_U_PL_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PL_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

// Application includes
#include "custom_utilities/poro_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlElement : public Element
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlElement );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPlElement(IndexType NewId = 0) : Element( NewId ) {}

    /// Constructor using an array of nodes
    UPlElement(IndexType NewId, const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPlElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element( NewId, pGeometry ) {}

    /// Constructor using Properties
    UPlElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element( NewId, pGeometry, pProperties )
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    /// Destructor
    virtual ~UPlElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo ) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable, const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable, const std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<double >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void AddExplicitContribution(const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double, 3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


protected:

    /// Member Variables

    GeometryData::IntegrationMethod mThisIntegrationMethod;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    Matrix mIntrinsicPermeability;

    std::vector<double> mImposedZStrainVector; /// The vector containing the imposed z strains (for 2.5D element: 2D geom with 3D CL)

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight);


    virtual void CalculateFluxResidual (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateMixBodyForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateNegInternalForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateExplicitContributions (VectorType& rFluxResidual, VectorType& rBodyForce, VectorType& rNegInternalForces, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateLumpedMassMatrix( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateDampingMatrixWithLumpedMass( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateInertialForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateDampingForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    /// Assignment operator.
    UPlElement & operator=(UPlElement const& rOther);

    /// Copy constructor.
    UPlElement(UPlElement const& rOther);


}; // Class UPlElement

} // namespace Kratos

#endif // KRATOS_U_PL_ELEMENT_H_INCLUDED  defined
