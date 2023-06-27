//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti
//                   Lorena Casallas
//                   Ignasi de Pouplana
//


#if !defined(KRATOS_U_PW_PG_SMALL_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_PG_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_Pg_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{UPwPgElement

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwPgSmallStrainElement : public UPwPgElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwPgSmallStrainElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwPgElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPwPgElement<TDim,TNumNodes>::mConstitutiveLawVector;
    using UPwPgElement<TDim,TNumNodes>::mIntrinsicPermeability;
    using UPwPgElement<TDim,TNumNodes>::mImposedZStrainVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwPgSmallStrainElement(IndexType NewId = 0) : UPwPgElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPwPgSmallStrainElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPwPgElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwPgSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPwPgElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwPgSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPwPgElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwPgSmallStrainElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        ///Properties variables
        double WaterDynamicViscosity;
        double GasDynamicViscosity;
        double FluidDensity;
        double GasDensity;
        double Density;
        double BiotCoefficient;
        double BiotModulusInverse;
        double SolidCompressibilityCoeff;
        double FluidCompressibilityCoeff;
        double GasCompressibilityCoeff;
        double Porosity
        double HenryCoefficient;
        bool AddGasDiffusion = false;

        //Element state variables
        double Sw;
        double dSwdpc;
        double ipCapilarPressure;
        
        // Water retention curve parameters
        int WaterSaturationLaw;
        double PoreSizeFactor;
        double ResidualWaterSaturation;
        double GasEntryPressure;

        ///ProcessInfo variables
        double VelocityCoefficient;
        double DtPressureCoefficient;

        ///Nodal variables
        array_1d<double,TNumNodes> PressureVector;
        array_1d<double,TNumNodes> GasPressureVector;
        array_1d<double,TNumNodes> DtPressureVector;
        array_1d<double,TNumNodes> DtGasPressureVector;
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
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes> UPwMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes> UPgMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes*TDim> PwUMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes*TDim> PgUMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> NpNpT;
        BoundedMatrix<double,TNumNodes,TNumNodes> PermMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PwPwMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PgPgMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PwPgMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PgPwMatrix;
        Matrix UVoigtMatrix;
        BoundedMatrix<double,TNumNodes,TDim> PDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
        array_1d<double,TNumNodes> PVector;
        array_1d<double,TNumNodes> CapilarPressureVector;
    };

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SaveGPGradPressure(Matrix& rGradPressureContainer, const array_1d<double,TDim>& GradPressure, const unsigned int& GPoint);

    void SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint);

    void ExtrapolateGPValues(const Matrix& GradPressureContainer, const Matrix& StressContainer, const unsigned int& VoigtSize);

    void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT);

    void CalculateKinematics(Matrix& rGradNpT,
                                Matrix& rB,
                                Vector& rStrainVector,
                                const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                const array_1d<double,TNumNodes*TDim>& DisplacementVector,
                                const unsigned int& GPoint);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void GetCouplingCompressibilityCoefficients(double Cwu, double Cgu, ElementVariables& rVariables);

    void CalculateWaterSaturationDegree(ElementVariables& rVariables);

    void GetCompressibilityCoefficients(double Cww, double Cwg, double Cgw, double Cgg, const ElementVariables& Variables);

    double EffectiveSaturation(double Sw, double Swr);

    double WaterRelativePermeability(double Se, const ElementVariables& Variables);

    double GasRelativePermeability(double Se, const ElementVariables& Variables);

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);


    void CalculateFluxResidual (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMixBodyForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateNegInternalForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateExplicitContributions (VectorType& rFluxResidual, VectorType& rBodyForce, VectorType& rNegInternalForces, const ProcessInfo& rCurrentProcessInfo) override;

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
    UPwPgSmallStrainElement & operator=(UPwPgSmallStrainElement const& rOther);

    /// Copy constructor.
    UPwPgSmallStrainElement(UPwPgSmallStrainElement const& rOther);

}; // Class UPwPgSmallStrainElement

} // namespace Kratos

#endif // KRATOS_U_PW_PG_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
