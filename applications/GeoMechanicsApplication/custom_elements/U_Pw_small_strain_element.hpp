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

#if !defined(KRATOS_GEO_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/comparison_utilities.hpp"
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainElement : public UPwBaseElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwSmallStrainElement );

    typedef std::size_t IndexType;
    typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    /// The definition of the sizetype
    typedef std::size_t SizeType;
    using UPwBaseElement<TDim,TNumNodes>::mConstitutiveLawVector;
    using UPwBaseElement<TDim,TNumNodes>::mRetentionLawVector;
    using UPwBaseElement<TDim,TNumNodes>::mStressVector;
    using UPwBaseElement<TDim,TNumNodes>::mStateVariablesFinalized;
    using UPwBaseElement<TDim,TNumNodes>::mIsInitialised;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwSmallStrainElement(IndexType NewId = 0) : UPwBaseElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPwSmallStrainElement(IndexType NewId,
                          const NodesArrayType& ThisNodes) : UPwBaseElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwSmallStrainElement(IndexType NewId,
                          GeometryType::Pointer pGeometry) : UPwBaseElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwSmallStrainElement(IndexType NewId,
                          GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties) : UPwBaseElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwSmallStrainElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create( IndexType NewId,
                             NodesArrayType const& ThisNodes,
                             PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties ) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                      std::vector<array_1d<double,3>>& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "U-Pw small strain Element #" << this->Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "U-Pw small strain Element #" << this->Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        ///Properties variables
        bool IgnoreUndrained;
        bool ConsiderGeometricStiffness;
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
        BoundedMatrix<double, TDim, TNumNodes*TDim> Nu;
        array_1d<double, TDim> BodyAcceleration;
        double IntegrationCoefficient;

        ///Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix F;
        double detF;
        Vector detJContainer;
        Matrix NContainer;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

        ///Retention Law parameters
        double FluidPressure;
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;

        // needed for updated Lagrangian:
        double detJ0;

        //Auxiliary Variables
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

    void SaveGPStress( Matrix &rStressContainer,
                       const Vector &StressVector,
                       const unsigned int &GPoint );

    void ExtrapolateGPValues( const Matrix &StressContainer,
                              const unsigned int &VoigtSize );


    void CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                           const ProcessInfo& CurrentProcessInfo ) override;


    void CalculateAll( MatrixType &rLeftHandSideMatrix,
                       VectorType &rRightHandSideVector,
                       const ProcessInfo& CurrentProcessInfo,
                       const bool CalculateStiffnessMatrixFlag,
                       const bool CalculateResidualVectorFlag ) override;

    void InitializeElementVariables( ElementVariables &rVariables,
                                     const ProcessInfo &CurrentProcessInfo );

    void SetConstitutiveParameters(ElementVariables &rVariables,
                                   ConstitutiveLaw::Parameters &rConstitutiveParameters);

    void SetRetentionParameters(const ElementVariables &rVariables,
                                RetentionLaw::Parameters &rRetentionParameters);

    virtual void CalculateKinematics( ElementVariables &rVariables, const unsigned int &PointNumber );

    void InitializeBiotCoefficients( ElementVariables &rVariables,
                                     const double &BulkModulus );

    void CalculateBMatrix( Matrix &rB,
                           const Matrix &GradNpT );


    virtual void CalculateAndAddLHS(MatrixType &rLeftHandSideMatrix, ElementVariables &rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType &rLeftHandSideMatrix, ElementVariables &rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType &rLeftHandSideMatrix, ElementVariables &rVariables);

    void CalculateAndAddCompressibilityMatrix(MatrixType &rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddPermeabilityMatrix(MatrixType &rLeftHandSideMatrix, ElementVariables &rVariables);

    virtual void CalculateAndAddRHS(VectorType &rRightHandSideVector, ElementVariables &rVariables);

    void CalculateAndAddStiffnessForce(VectorType &rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType &rRightHandSideVector, ElementVariables &rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables &rVariables);

    void CalculateAndAddCompressibilityFlow(VectorType &rRightHandSideVector, ElementVariables &rVariables);

    void CalculateAndAddPermeabilityFlow(VectorType &rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType &rRightHandSideVector, ElementVariables &rVariables);

    void UpdateElementalVariableStressVector(ElementVariables &rVariables, const unsigned int &PointNumber);
    void UpdateElementalVariableStressVector(Vector &StressVector, const unsigned int &PointNumber);
    void UpdateStressVector(const ElementVariables &rVariables, const unsigned int &PointNumber);
    void UpdateStressVector(const Vector &StressVector, const unsigned int &PointNumber);

    double CalculateBulkModulus(const Matrix &ConstitutiveMatrix);

    virtual void CalculateCauchyAlmansiStrain( ElementVariables &rVariables );
    virtual void CalculateCauchyGreenStrain( ElementVariables &rVariables );
    virtual void CalculateCauchyStrain( ElementVariables &rVariables );
    virtual void CalculateStrain( ElementVariables &rVariables );

    void InitializeNodalVariables( ElementVariables &rVariables );
    void InitializeProperties( ElementVariables &rVariables );
    double CalculateFluidPressure( const ElementVariables &rVariables, const unsigned int &PointNumber );

    void CalculateRetentionResponse( ElementVariables &rVariables,
                                     RetentionLaw::Parameters &rRetentionParameters,
                                     const unsigned int &GPoint );

    void CalculateExtrapolationMatrix(BoundedMatrix<double,TNumNodes,TNumNodes> &rExtrapolationMatrix);

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

    // Assignment operator.
    UPwSmallStrainElement & operator=(UPwSmallStrainElement const& rOther);

    // Copy constructor.
    UPwSmallStrainElement(UPwSmallStrainElement const& rOther);

}; // Class UPwSmallStrainElement

} // namespace Kratos

#endif // KRATOS_GEO_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
