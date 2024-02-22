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


#if !defined(KRATOS_U_PL_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PL_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/one-phase_flow/U_Pl_element.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlSmallStrainInterfaceElement : public UPlElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlSmallStrainInterfaceElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPlElement<TDim,TNumNodes>::mConstitutiveLawVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPlSmallStrainInterfaceElement(IndexType NewId = 0) : UPlElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPlSmallStrainInterfaceElement(IndexType NewId, const NodesArrayType& ThisNodes) : UPlElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPlSmallStrainInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPlElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPlSmallStrainInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPlElement<TDim,TNumNodes>( NewId, pGeometry, pProperties )
    {
        /// Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
    }

    /// Destructor
    ~UPlSmallStrainInterfaceElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct SFGradAuxVariables
    {
        array_1d<double,TDim> GlobalCoordinatesGradients;
        array_1d<double,TDim> LocalCoordinatesGradients;

        BoundedMatrix<double,TNumNodes,TDim-1> ShapeFunctionsNaturalGradientsMatrix;
        BoundedMatrix<double,TDim-1,TDim-1> LocalCoordinatesGradientsMatrix;
        BoundedMatrix<double,TDim-1,TDim-1> LocalCoordinatesGradientsInvMatrix;
        BoundedMatrix<double,TNumNodes,TDim-1> ShapeFunctionsGradientsMatrix;
    };

    struct InterfaceElementVariables
    {
        ///Properties variables
        double DynamicViscosityInverse;
        double LiquidDensity;
        double Density;
        double BiotCoefficient;
        double BiotModulusInverse;

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
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        array_1d<double,TDim> VoigtVector;

        ///Variables computed at each GP

        ///Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix F;
        double detF;
        ///Auxiliary Variables
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu;
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix;
        array_1d<double,TDim> BodyAcceleration;
        double IntegrationCoefficient;
        double JointWidth;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes*TDim> UMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TNumNodes> UPMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes*TDim> PUMatrix;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrix;
        BoundedMatrix<double,TDim,TDim> DimMatrix;
        BoundedMatrix<double,TNumNodes*TDim,TDim> UDimMatrix;
        BoundedMatrix<double,TNumNodes,TDim> PDimMatrix;
        array_1d<double,TNumNodes*TDim> UVector;
        array_1d<double,TNumNodes> PVector;
    };

    /// Member Variables

    std::vector<double> mInitialGap;
    std::vector<bool> mIsOpen;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateInitialGap(const GeometryType& Geom);


    void ExtrapolateGPValues (const std::vector<double>& JointWidthContainer);


    void CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo ) override;


    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void InitializeElementVariables(InterfaceElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& rCurrentProcessInfo);

    void CalculateRotationMatrix(BoundedMatrix<double,TDim,TDim>& rRotationMatrix, const GeometryType& Geom);

    void CalculateJointWidth(double& rJointWidth,const double& NormalRelDisp,const double& InitialJointWidth,const unsigned int& GPoint);

    void CheckAndCalculateJointWidth(double& rJointWidth,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    double& rNormalRelDisp,const double& InitialJointWidth,const unsigned int& GPoint);

    template< class TMatrixType >
    void CalculateShapeFunctionsGradients(TMatrixType& rGradNpT, SFGradAuxVariables& rAuxVariables,const Matrix& Jacobian,
                                            const BoundedMatrix<double,TDim,TDim>& RotationMatrix,
                                            const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint);


    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);


    void CalculateAndAddRHS(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);


    void CalculateFluxResidual (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMixBodyForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateNegInternalForce (VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateExplicitContributions (VectorType& rFluxResidual, VectorType& rBodyForce, VectorType& rNegInternalForces, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLumpedMassMatrix( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo ) override;


    void InterpolateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues );

    template< class TValueType >
    void InterpolateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        typedef UPlElement<TDim,TNumNodes> BaseElement;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseElement )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlElement<TDim,TNumNodes> BaseElement;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseElement )
    }

    /// Assignment operator.
    UPlSmallStrainInterfaceElement & operator=(UPlSmallStrainInterfaceElement const& rOther);

    /// Copy constructor.
    UPlSmallStrainInterfaceElement(UPlSmallStrainInterfaceElement const& rOther);

}; // Class UPlSmallStrainInterfaceElement

} // namespace Kratos

#endif // KRATOS_U_PL_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED  defined
