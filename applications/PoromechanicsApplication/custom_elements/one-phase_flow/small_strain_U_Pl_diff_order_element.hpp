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


#if !defined(KRATOS_SMALL_STRAIN_U_PL_DIFF_ORDER_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_U_PL_DIFF_ORDER_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

// Application includes
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/poro_element_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SmallStrainUPlDiffOrderElement : public Element
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SmallStrainUPlDiffOrderElement );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallStrainUPlDiffOrderElement();

    // Constructor 1
    SmallStrainUPlDiffOrderElement(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor 2
    SmallStrainUPlDiffOrderElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    virtual ~SmallStrainUPlDiffOrderElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable, const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable, const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable, const std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementalVariables
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
        Matrix GradNpT; //Contains the global derivatives of the pressure shape functions
        Matrix B;
        double IntegrationCoefficient;
        Vector StrainVector;
        Matrix ConstitutiveMatrix;
        Vector StressVector; //It is the "Effective Stress Vector": sigma'_ij = sigma_ij + alpha*pl*delta_ij

        //Variables needed for consistency with the general constitutive law
        double detF;
        Matrix F;

        //Nodal variables
        Vector BodyAcceleration;
        Vector DisplacementVector;
        Vector VelocityVector;
        Vector PressureVector;
        Vector PressureDtVector;

        //Properties and processinfo variables
        double BiotCoefficient;
        double BiotModulusInverse;
        double DynamicViscosity;
        double NewmarkCoefficient1;
        double NewmarkCoefficient2;
    };

    // Member Variables

    GeometryData::IntegrationMethod mThisIntegrationMethod;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    Geometry< Node >::Pointer mpPressureGeometry;

    Matrix mIntrinsicPermeability;
    
    std::vector<double> mImposedZStrainVector; /// The vector containing the imposed z strains (for 2.5D element: 2D geom with 3D CL)

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag);

    void InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void InitializeNodalVariables (ElementalVariables& rVariables);

    void InitializeProperties (ElementalVariables& rVariables);

    void CalculateKinematics(ElementalVariables& rVariables, unsigned int PointNumber);

    void SetElementalVariables(ElementalVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters);

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, double detJ, double weight);


    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables);

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables);

    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables);


    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables);

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
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType>& Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }

}; // Class SmallStrainUPlDiffOrderElement

} // namespace Kratos

#endif // KRATOS_SMALL_STRAIN_U_PL_DIFF_ORDER_ELEMENT_H_INCLUDED  defined
