// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_U_PW_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_U_PW_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "custom_retention/retention_law.h"
#include "custom_retention/retention_law_factory.h"

// Application includes
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwBaseElement : public Element
{

public:

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwBaseElement );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwBaseElement(IndexType NewId = 0) : Element( NewId ) {}

    /// Constructor using an array of nodes
    UPwBaseElement(IndexType NewId, const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwBaseElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element( NewId, pGeometry ) {}

    /// Constructor using Properties
    UPwBaseElement(IndexType NewId,
                   GeometryType::Pointer pGeometry,
                   PropertiesType::Pointer pProperties) : Element( NewId, pGeometry, pProperties ) {}

    /// Destructor
    virtual ~UPwBaseElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, 
                    const ProcessInfo& rCurrentProcessInfo) const override;

    void ResetConstitutiveLaw() override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                                const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo ) const override;

    void CalculateMassMatrix( MatrixType& rMassMatrix,
                              const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateDampingMatrix( MatrixType& rDampingMatrix,
                                 const ProcessInfo& rCurrentProcessInfo ) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                      const std::vector<double>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      const std::vector<Matrix>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "U-Pw Base class Element #" << this->Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "U-Pw Base class Element #" << this->Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    std::vector<RetentionLaw::Pointer> mRetentionLawVector;

    std::vector<Vector> mStressVector;
    std::vector<Vector> mStateVariablesFinalized;
    bool mIsInitialised = false;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                                   const ProcessInfo& CurrentProcessInfo );

    virtual void CalculateAll( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& CurrentProcessInfo,
                               const bool CalculateStiffnessMatrixFlag,
                               const bool CalculateResidualVectorFlag);

    virtual double CalculateIntegrationCoefficient( const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                                    const IndexType& PointNumber,
                                                    const double& detJ);

    /**
     * @brief This methods gives us a matrix with the increment of displacement
     * @param DeltaDisplacement The matrix containing the increment of displacements
     * @return DeltaDisplacement: The matrix containing the increment of displacements
     */
    Matrix& CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const;

    double CalculateDerivativesOnInitialConfiguration(const GeometryType& Geometry,
                                                      Matrix& DNu_DX0,
                                                      const IndexType& GPoint,
                                                      IntegrationMethod ThisIntegrationMethod) const;

    virtual unsigned int GetNumberOfDOF() const;
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
    UPwBaseElement & operator=(UPwBaseElement const& rOther);

    /// Copy constructor.
    UPwBaseElement(UPwBaseElement const& rOther);


}; // Class UPwBaseElement

} // namespace Kratos

#endif // KRATOS_GEO_U_PW_ELEMENT_H_INCLUDED  defined
