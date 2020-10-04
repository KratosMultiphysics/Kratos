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

// Application includes
#include "custom_utilities/comparison_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUPwDiffOrderElement : public Element
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUPwDiffOrderElement );

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

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void GetValueOnIntegrationPoints(const Variable<int>& rVariable,
                                     std::vector<int>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                     std::vector<double>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                     std::vector<Vector>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                     std::vector<array_1d<double,3>>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                     std::vector<Matrix>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                     std::vector<ConstitutiveLaw::Pointer>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
        double detJ0;
        Matrix J0;
        Matrix InvJ0;

        //Nodal variables
        Vector BodyAcceleration;
        Vector DisplacementVector;
        Vector VelocityVector;
        Vector PressureVector;
        Vector PressureDtVector;

        //Properties and processinfo variables
        bool IgnoreUndrained;
        double BiotCoefficient;
        double BiotModulusInverse;
        double DynamicViscosity;
        Matrix IntrinsicPermeability;
        double NewmarkCoefficient1;
        double NewmarkCoefficient2;
    };

    // Member Variables

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    //Geometry< Node<3> >::Pointer mpPressureGeometry;
    GeometryType::Pointer  mpPressureGeometry;
    std::vector<Vector> mStressVector;
    std::vector<Vector> mStateVariablesFinalized;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void InitializeNodalVariables(ElementVariables& rVariables);

    void InitializeProperties(ElementVariables& rVariables);

    void InitializeBiotCoefficients( ElementVariables& rVariables, const double &BulkModulus );

    void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

    void UpdateElementalVariableStressVector(ElementVariables& rVariables, unsigned int PointNumber);

    void UpdateStressVector(const ElementVariables& rVariables, unsigned int PointNumber);

    void SetElementalVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters);

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, double detJ, double weight);

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

    // convention factor for liquid pressure
    const double signFactor = 1.0;

}; // Class SmallStrainUPwDiffOrderElement

} // namespace Kratos

#endif // KRATOS_GEO_SMALL_STRAIN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED  defined
