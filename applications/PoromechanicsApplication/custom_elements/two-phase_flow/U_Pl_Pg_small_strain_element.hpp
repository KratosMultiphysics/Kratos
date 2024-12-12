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
//                   Ignasi de Pouplana
//                   Xavier Tort
//                   Lorena Casallas
//


#if !defined(KRATOS_U_PL_PG_SMALL_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PL_PG_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/two-phase_flow/U_Pl_Pg_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlPgSmallStrainElement : public UPlPgElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlPgSmallStrainElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlPgElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPlPgElement<TDim,TNumNodes>::mConstitutiveLawVector;
    using UPlPgElement<TDim,TNumNodes>::mSaturationLawVector;
    using UPlPgElement<TDim,TNumNodes>::mIntrinsicPermeability;
    using UPlPgElement<TDim,TNumNodes>::mImposedZStrainVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPlPgSmallStrainElement(IndexType NewId = 0) : UPlPgElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPlPgSmallStrainElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPlPgElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPlPgSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPlPgElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPlPgSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPlPgElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPlPgSmallStrainElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

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
        double LiquidDynamicViscosity;
        double GasDynamicViscosity;
        double LiquidDensity;
        double GasDensity;
        double SolidDensity;
        double Density;
        double BiotCoefficient;
        double BiotModulusInverse;
        double SolidCompressibilityCoeff;
        double LiquidCompressibilityCoeff;
        double GasCompressibilityCoeff;
        double Porosity;
        double GasHenrySolubilityCoeff;
        double GasMolarWeight;
        double Temperature;
        double TortuosityCoeff;
        double GasDiffusionCoeff;
        double LiquidDensityRefPressure;

        ///ProcessInfo variables
        double VelocityCoefficient;
        double DtLiquidPressureCoefficient;
        double DtGasPressureCoefficient;

        ///Nodal variables
        array_1d<double,TNumNodes> LiquidPressureVector;
        array_1d<double,TNumNodes> GasPressureVector;
        array_1d<double,TNumNodes> DtLiquidPressureVector;
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
        //Saturation Law parameters
        double Sl;
        double dSldPc;
        double krl;
        double krg;
        ///Auxiliary Variables
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes> UPMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes> UPMatrixAux;
        BoundedMatrix<double,TNumNodes,TNumNodes*TDim> PUMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrixAux;
        Matrix UVoigtMatrix;
        BoundedMatrix<double,TNumNodes,TDim> PDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
        array_1d<double,TNumNodes> PVector;
        array_1d<double,TNumNodes> PVectorAux;
    };

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLumpedMassMatrix( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateAndAddMassMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);


    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddHydromechanicalCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddHydromechanicalCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);


    void SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint);

    void ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize);

    void InitializeElementVariables(ElementVariables& rVariables, ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    SaturationLaw::Parameters& rSaturationParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& rCurrentProcessInfo);

    void CalculateKinematics(Matrix& rGradNpT,
                                Matrix& rB,
                                Vector& rStrainVector,
                                const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                const array_1d<double,TNumNodes*TDim>& DisplacementVector,
                                const unsigned int& GPoint);

    void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT);

    void GetHydromechanicalCouplingCoefficients(double& Clu, double& Cgu, ElementVariables& rVariables);

    void GetCompressibilityCoefficients(double& Cll, double& Clg, double& Cgl, double& Cgg, const ElementVariables& Variables);

    void UpdateDensity(ElementVariables& Variables);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        typedef UPlPgElement<TDim,TNumNodes> BaseElement;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseElement )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlPgElement<TDim,TNumNodes> BaseElement;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseElement )
    }

    /// Assignment operator.
    UPlPgSmallStrainElement & operator=(UPlPgSmallStrainElement const& rOther);

    /// Copy constructor.
    UPlPgSmallStrainElement(UPlPgSmallStrainElement const& rOther);

}; // Class UPlPgSmallStrainElement

} // namespace Kratos

#endif // KRATOS_U_PL_PG_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
