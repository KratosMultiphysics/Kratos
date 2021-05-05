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

#if !defined(KRATOS_GEO_U_PW_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_U_PW_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainInterfaceElement : public UPwBaseElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwSmallStrainInterfaceElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwBaseElement<TDim,TNumNodes>::mConstitutiveLawVector;
    using UPwBaseElement<TDim,TNumNodes>::mRetentionLawVector;
    using UPwBaseElement<TDim,TNumNodes>::mStressVector;
    using UPwBaseElement<TDim,TNumNodes>::mStateVariablesFinalized;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwSmallStrainInterfaceElement(IndexType NewId = 0) : UPwBaseElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPwSmallStrainInterfaceElement(IndexType NewId,
                                   const NodesArrayType& ThisNodes) : UPwBaseElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwSmallStrainInterfaceElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry) : UPwBaseElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwSmallStrainInterfaceElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
                                   : UPwBaseElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwSmallStrainInterfaceElement() override {}

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                     std::vector<double>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                      std::vector<array_1d<double,3>>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;


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
        bool IgnoreUndrained;
        double DynamicViscosityInverse;
        double FluidDensity;
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

        ///Retention Law parameters
        double FluidPressure;
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;
    };

    /// Member Variables

    std::vector<double> mInitialGap;
    std::vector<bool> mIsOpen;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void ModifyInactiveElementStress(const double &JointWidth, Vector &StressVector);

    void CalculateOnLobattoIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                            std::vector<array_1d<double,3>>& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnLobattoIntegrationPoints(const Variable<Matrix>& rVariable,
                                             std::vector<Matrix>& rOutput,
                                             const ProcessInfo& rCurrentProcessInfo);

    void CalculateInitialGap(const GeometryType& Geom);

    void ExtrapolateGPValues (const std::vector<double>& JointWidthContainer);

    void CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                           const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& CurrentProcessInfo,
                       const bool CalculateStiffnessMatrixFlag,
                       const bool CalculateResidualVectorFlag) override;

    void InitializeElementVariables(InterfaceElementVariables& rVariables,
                                    ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType& Geom,
                                    const PropertiesType& Prop,
                                    const ProcessInfo& CurrentProcessInfo);

    void CalculateRotationMatrix(BoundedMatrix<double,TDim,TDim>& rRotationMatrix, const GeometryType& Geom);

    void CalculateJointWidth(double& rJointWidth,
                             const double& NormalRelDisp,
                             const double& MinimumJointWidth,
                             const unsigned int& GPoint);

    void CheckAndCalculateJointWidth(double& rJointWidth,
                                    ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    double& rNormalRelDisp,
                                    const double& MinimumJointWidth,
                                    const unsigned int& GPoint);

    template< class TMatrixType >
    void CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,
                                          SFGradAuxVariables& rAuxVariables,
                                          const Matrix& Jacobian,
                                          const BoundedMatrix<double,TDim,TDim>& RotationMatrix,
                                          const Matrix& DN_De,
                                          const Matrix& Ncontainer,
                                          const double& JointWidth,
                                          const unsigned int& GPoint);


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

    void InterpolateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues );

    void UpdateElementalVariableStressVector(InterfaceElementVariables &rVariables, const unsigned int &PointNumber);

    void UpdateElementalVariableStressVector(Vector &StressVector, const unsigned int &PointNumber);

    void UpdateStressVector(const InterfaceElementVariables &rVariables, const unsigned int &PointNumber);

    void UpdateStressVector(const Vector &StressVector, const unsigned int &PointNumber);

    template< class TValueType >
    void InterpolateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues );

    void SetRetentionParameters(const InterfaceElementVariables &rVariables,
                                RetentionLaw::Parameters &rRetentionParameters);
    double CalculateFluidPressure( const InterfaceElementVariables &rVariables, const unsigned int &PointNumber );
    
    double CalculateBulkModulus(const Matrix &ConstitutiveMatrix);

    void InitializeBiotCoefficients( InterfaceElementVariables &rVariables,
                                     const double &BulkModulus );

    void CalculateRetentionResponse( InterfaceElementVariables& rVariables,
                                     RetentionLaw::Parameters& rRetentionParameters,
                                     const unsigned int &GPoint );

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
    UPwSmallStrainInterfaceElement & operator=(UPwSmallStrainInterfaceElement const& rOther);

    /// Copy constructor.
    UPwSmallStrainInterfaceElement(UPwSmallStrainInterfaceElement const& rOther);

}; // Class UPwSmallStrainInterfaceElement

} // namespace Kratos

#endif // KRATOS_GEO_U_PW_SMALL_STRAIN_INTERFACE_ELEMENT_H_INCLUDED  defined
