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


#if !defined(KRATOS_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwSmallStrainElement : public UPwElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwSmallStrainElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPwElement<TDim,TNumNodes>::mConstitutiveLawVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwSmallStrainElement(IndexType NewId = 0) : UPwElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPwSmallStrainElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPwElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPwElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPwElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwSmallStrainElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        ///Properties variables
        double DynamicViscosityInverse;
        double FluidDensity;
        double Density;
        double BiotCoefficient;
        double BiotModulusInverse;
        BoundedMatrix<double,TDim, TDim> PermeabilityMatrix;

        ///ProcessInfo variables
        double VelocityCoefficient;
        double DtPressureCoefficient;

        ///Nodal variables
        array_1d<double,TNumNodes> PressureVector;
        array_1d<double,TNumNodes> DtPressureVector;
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        array_1d<double,TNumNodes*TDim> VelocityVector;
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;

        ///General elemental variables
        Vector VoigtVector;

        ///Variables computed at each GP
        Matrix B;
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu;
        array_1d<double,TDim> BodyAcceleration;
        double IntegrationCoefficient;
        ///Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix F;
        double detF;
        ///Auxiliary Variables
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes> UPMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes*TDim> PUMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrix;
        Matrix UVoigtMatrix;
        BoundedMatrix<double,TNumNodes,TDim> PDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
        array_1d<double,TNumNodes> PVector;
    };

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint);

    void ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize);


    void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo ) override;


    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT);


    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);


    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    UPwSmallStrainElement & operator=(UPwSmallStrainElement const& rOther);

    /// Copy constructor.
    UPwSmallStrainElement(UPwSmallStrainElement const& rOther);

}; // Class UPwSmallStrainElement

} // namespace Kratos

#endif // KRATOS_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
