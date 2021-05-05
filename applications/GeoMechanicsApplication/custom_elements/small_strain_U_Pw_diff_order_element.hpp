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

#if !defined(KRATOS_GEO_SMALL_STRAIN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_SMALL_STRAIN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED

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
#include "custom_utilities/comparison_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUPwDiffOrderElement : public Element
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SmallStrainUPwDiffOrderElement );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    SmallStrainUPwDiffOrderElement();

    // Constructor 1
    SmallStrainUPwDiffOrderElement(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor 2
    SmallStrainUPwDiffOrderElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    virtual ~SmallStrainUPwDiffOrderElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList, 
                    const ProcessInfo& rCurrentProcessInfo) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void CalculateDampingMatrix( MatrixType& rDampingMatrix,
                                 const ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                      const std::vector<double>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      const std::vector<Matrix>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
                                     std::vector<int>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double> &rVariable,
                                      std::vector<double> &rOutput,
                                      const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector> &rVariable,
                                      std::vector<Vector> &rOutput,
                                      const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>> &rVariable,
                                      std::vector<array_1d<double,3>> &rOutput,
                                      const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix> &rVariable,
                                      std::vector< Matrix> &rOutput,
                                      const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo ) override;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "U-Pw small strain different order Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "U-Pw small strain different order Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        //Variables at all integration points
        Matrix NuContainer;
        Matrix NpContainer;
        GeometryType::ShapeFunctionsGradientsType DNu_DXContainer;
        GeometryType::ShapeFunctionsGradientsType DNp_DXContainer;
        Vector detJuContainer;

        //Variables at each integration point
        Vector Nu; //Contains the displacement shape functions at every node
        Vector Np; //Contains the pressure shape functions at every node
        Matrix DNu_DX; //Contains the global derivatives of the displacement shape functions
        Matrix DNp_DX; //Contains the global derivatives of the pressure shape functions
        Matrix B;
        double IntegrationCoefficient;
        Vector StrainVector;
        Matrix ConstitutiveMatrix;
        Vector StressVector; //It is the "Effective Stress Vector": sigma'_ij = sigma_ij + alpha*pw*delta_ij

        //Variables needed for consistency with the general constitutive law
        double detF;
        Matrix F;

        // needed for updated Lagrangian:
        double detJ0;  // displacement
        double detJp0; // pore pressure

        //Nodal variables
        Vector BodyAcceleration;
        Vector DisplacementVector;
        Vector VelocityVector;
        Vector PressureVector;
        Vector PressureDtVector;

        ///Retention Law parameters
        double FluidPressure;
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;

        //Properties and processinfo variables
        bool IgnoreUndrained;
        bool ConsiderGeometricStiffness;
        double BiotCoefficient;
        double BiotModulusInverse;
        double DynamicViscosityInverse;
        Matrix IntrinsicPermeability;
        double NewmarkCoefficient1;
        double NewmarkCoefficient2;
    };

    // Member Variables

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    std::vector<RetentionLaw::Pointer> mRetentionLawVector;

    GeometryType::Pointer  mpPressureGeometry;
    std::vector<Vector> mStressVector;
    std::vector<Vector> mStressVectorFinalized;
    std::vector<Vector> mStateVariablesFinalized;
    bool mIsInitialised = false;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    virtual void CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                                   const ProcessInfo& CurrentProcessInfo );

    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void InitializeNodalVariables(ElementVariables& rVariables);

    void InitializeProperties(ElementVariables& rVariables);

    void InitializeBiotCoefficients( ElementVariables& rVariables, const double &BulkModulus );

    virtual void CalculateKinematics(ElementVariables& rVariables, const unsigned int &PointNumber);

    virtual void CalculateKinematicsOnInitialConfiguration(ElementVariables& rVariables, unsigned int PointNumber);

    double CalculateDerivativesOnInitialConfiguration(const GeometryType& Geometry,
                                                      Matrix& DN_DX,
                                                      const IndexType& PointNumber,
                                                      IntegrationMethod ThisIntegrationMethod) const;

    void UpdateElementalVariableStressVector(ElementVariables& rVariables, unsigned int PointNumber);

    void UpdateStressVector(const ElementVariables& rVariables, unsigned int PointNumber);

    void SetConstitutiveParameters(ElementVariables& rVariables,
                                   ConstitutiveLaw::Parameters& rConstitutiveParameters);

    virtual void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, double detJ, double weight);

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

    double CalculateBulkModulus(const Matrix &ConstitutiveMatrix);

    void CalculateBMatrix( Matrix& rB,
                           const Matrix& DNu_DX );

    void AssignPressureToIntermediateNodes();

    void AssembleUBlockMatrix(Matrix &rLeftHandSideMatrix,
                              const Matrix &StiffnessMatrix) const;

    virtual void CalculateCauchyAlmansiStrain( ElementVariables& rVariables );
    virtual void CalculateCauchyGreenStrain( ElementVariables& rVariables );
    virtual void CalculateCauchyStrain( ElementVariables& rVariables );
    virtual void CalculateStrain( ElementVariables& rVariables );

    double CalculateFluidPressure( const ElementVariables &rVariables, const unsigned int &PointNumber );

    void SetRetentionParameters(const ElementVariables& rVariables,
                                RetentionLaw::Parameters& rRetentionParameters);

    void CalculateRetentionResponse( ElementVariables &rVariables,
                                     RetentionLaw::Parameters &rRetentionParameters,
                                     const unsigned int &GPoint );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }


    // Private Operations

    template < class TValueType >
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType> &Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }

}; // Class SmallStrainUPwDiffOrderElement

} // namespace Kratos

#endif // KRATOS_GEO_SMALL_STRAIN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED  defined
