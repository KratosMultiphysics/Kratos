// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//

// System includes

// External includes


// Project includes
#include "custom_elements/generic_total_lagrangian_femdem_element.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "fem_to_dem_application_variables.h"


namespace Kratos
{
namespace
{
class LargeDisplacementKinematics
{
public:
    LargeDisplacementKinematics(Element::GeometryType const& rGeom) : mrGeom(rGeom) { mThisIntegrationMethod = rGeom.GetDefaultIntegrationMethod();}
    LargeDisplacementKinematics(Element::GeometryType const& rGeom, Element::GeometryType::IntegrationMethod ThisIntegrationMethod) : mrGeom(rGeom), mThisIntegrationMethod(ThisIntegrationMethod) {}

    void DeformationGradient(std::size_t IntegrationPoint, Matrix& rF)
    {
        if (IntegrationPoint != mCurrentIntegrationPoint)
            Recalculate(IntegrationPoint);
        GeometryUtils::DeformationGradient(mJ, mInvJ0, rF);
    }

    void ShapeFunctionsGradients(std::size_t IntegrationPoint, Matrix& rDN_DX0)
    {
        if (IntegrationPoint != mCurrentIntegrationPoint)
            Recalculate(IntegrationPoint);
        Matrix const& rDN_De = mrGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod)[IntegrationPoint];
        GeometryUtils::ShapeFunctionsGradients(rDN_De, mInvJ0, rDN_DX0);
    }

    double DetJ0(std::size_t IntegrationPoint)
    {
        if (IntegrationPoint != mCurrentIntegrationPoint)
            Recalculate(IntegrationPoint);
        return mDetJ0;
    }

private:
    void Recalculate(std::size_t IntegrationPoint)
    {
        // Calculate mInvJ0 and mDetJ0.
        GeometryUtils::JacobianOnInitialConfiguration(
            mrGeom, mrGeom.IntegrationPoints(mThisIntegrationMethod)[IntegrationPoint], mJ);
        MathUtils<double>::InvertMatrix(mJ, mInvJ0, mDetJ0);
        // Calculate mJ.
        mrGeom.Jacobian(mJ, IntegrationPoint, mThisIntegrationMethod);
        // Update current integration point.
        mCurrentIntegrationPoint = IntegrationPoint;
    }

    Element::GeometryType const& mrGeom; /// Currently geometry
    Element::GeometryType::IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods
    Matrix mJ;
    Matrix mInvJ0;
    double mDetJ0 = 0.0;
    std::size_t mCurrentIntegrationPoint = -1;
};

void CalculateGreenLagrangeStrainSensitivity(Matrix const& rF,
                                             Matrix const& rF_Deriv,
                                             Matrix& rE_Deriv)
{
    rE_Deriv = 0.5 * (prod(trans(rF_Deriv), rF) + prod(trans(rF), rF_Deriv));
}
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::GenericTotalLagrangianFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseSolidElement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
    if (mThresholds.size() != NumberOfEdges)
        mThresholds.resize(NumberOfEdges);
    noalias(mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (mDamages.size() != NumberOfEdges)
        mDamages.resize(NumberOfEdges);
    noalias(mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::GenericTotalLagrangianFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
    if (mThresholds.size() != NumberOfEdges)
        mThresholds.resize(NumberOfEdges);
    noalias(mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (mDamages.size() != NumberOfEdges)
        mDamages.resize(NumberOfEdges);
    noalias(mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianFemDemElement>(NewId, GetGeometry().Create( ThisNodes ), pProperties);
}

//************************************************************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianFemDemElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::~GenericTotalLagrangianFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Pointer p_new_elem = Kratos::make_intrusive<GenericTotalLagrangianFemDemElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::InitializeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
    if (this->GetValue(RECOMPUTE_NEIGHBOURS)) {
        this->ComputeEdgeNeighbours(rCurrentProcessInfo);
        this->SetValue(RECOMPUTE_NEIGHBOURS, false);
    }

    this->InitializeInternalVariablesAfterMapping();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters cl_values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options=cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_values.SetStrainVector(this_constitutive_variables.StrainVector);
    cl_values.SetStressVector(this_constitutive_variables.StressVector);
    cl_values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // Reading integration points
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this->GetStressMeasure());
    }
    this->SetValue(STRESS_VECTOR, this_constitutive_variables.StressVector);
    this->SetValue(STRAIN_VECTOR, this_constitutive_variables.StrainVector);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes   = this->GetGeometry().size();
    const SizeType dimension         = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size           = GetStrainSize();
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if (CalculateStiffnessMatrixFlag == true) { // Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) { // Calculation of the matrix is required
        if (rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize(mat_size, false);
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    // Reading integration points
    const auto& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options = cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    cl_values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;
    const double characteristic_length = this->CalculateCharacteristicLength(this);

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this->GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = this->GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && this->GetProperties().Has(THICKNESS))
            int_to_reference_weight *= this->GetProperties()[THICKNESS];

        bool is_damaging = false;
        Vector damages_edges = ZeroVector(NumberOfEdges);

        if (yield_surface != "Elastic") {
            // Loop over edges of the element...
            Vector average_stress_edge(VoigtSize);
            Vector average_strain_edge(VoigtSize);

            for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
                noalias(average_stress_edge) = this_constitutive_variables.StressVector;
                noalias(average_strain_edge) = this_constitutive_variables.StrainVector;
                this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
                this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);
 
                damages_edges[edge] = mDamages[edge];
                double threshold = mThresholds[edge];
                
                this->IntegrateStressDamageMechanics(threshold, damages_edges[edge], average_strain_edge, 
                                                     average_stress_edge, edge, characteristic_length, cl_values, 
                                                     is_damaging);

            } // Loop over edges
        }

        // Calculate the elemental Damage...
        const double damage_element = this->CalculateElementalDamage(damages_edges);

        Vector r_strain_vector;
        this->CalculateGreenLagrangeStrainVector(r_strain_vector, this_kinematic_variables.F);
        const Vector& r_stress_vector = this_constitutive_variables.StressVector;
        const Vector& r_integrated_stress_vector = (1.0 - damage_element)*r_stress_vector;

        if (CalculateStiffnessMatrixFlag == true) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            if (is_damaging == true && norm_2(r_strain_vector) > tolerance) {
                Matrix tangent_tensor;
                if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 0) {
                    tangent_tensor = this_constitutive_variables.D;
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 1) {
                    tangent_tensor = (1.0 - damage_element) * this_constitutive_variables.D;
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 2) {
                    this->CalculateTangentTensor(tangent_tensor, r_strain_vector, r_integrated_stress_vector, this_kinematic_variables.F, this_constitutive_variables.D, cl_values);
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 3) {
                    this->CalculateTangentTensorSecondOrder(tangent_tensor, r_strain_vector, r_integrated_stress_vector, this_kinematic_variables.F, this_constitutive_variables.D, cl_values);
                }
                this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, tangent_tensor, int_to_reference_weight);                
            } else {
                this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, (1.0 - damage_element)*this_constitutive_variables.D, int_to_reference_weight);
            }
            /* Geometric stiffness matrix */
            this->CalculateAndAddKg(rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight);
        }

        if (CalculateResidualVectorFlag == true) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, r_integrated_stress_vector, int_to_reference_weight);
        }
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo 
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes   = this->GetGeometry().size();
    const SizeType dimension         = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size           = GetStrainSize();
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options = cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // If strain has to be computed inside of the constitutive law with PK2
    cl_values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;
    const double characteristic_length = this->CalculateCharacteristicLength(this);

    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this->GetStressMeasure());

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(cl_values, GetStressMeasure());

        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->FinalizeSolutionStep( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = this->GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && this->GetProperties().Has(THICKNESS))
            int_to_reference_weight *= this->GetProperties()[THICKNESS];

        bool is_damaging = false;
        if (yield_surface != "Elastic") {
            // Loop over edges of the element...
            Vector average_stress_edge(VoigtSize);
            Vector average_strain_edge(VoigtSize);

            for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
                noalias(average_stress_edge) = this_constitutive_variables.StressVector;
                noalias(average_strain_edge) = this_constitutive_variables.StrainVector;
                this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
                this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

                this->IntegrateStressDamageMechanics(mThresholds[edge], mDamages[edge], average_strain_edge, 
                                                     average_stress_edge, edge, characteristic_length, cl_values, 
                                                     is_damaging);
            } // Loop over edges
        }

        // Calculate the elemental Damage...
        mDamage = this->CalculateElementalDamage(mDamages);

        if (mDamage >= 0.98) {
            this->Set(ACTIVE, false);
            mDamage = 0.98;
            // We set a "flag" to generate the DEM 
            rCurrentProcessInfo[GENERATE_DEM] = true;
        }
    }
    KRATOS_CATCH( "" )
}


/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    // Shape functions
    rThisKinematicVariables.N = row(GetGeometry().ShapeFunctionsValues(rIntegrationMethod), PointNumber);

    rThisKinematicVariables.detJ0 = this->CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    Matrix J;
    J = this->GetGeometry().Jacobian(J, PointNumber, rIntegrationMethod);

    GeometryUtils::DeformationGradient(J, rThisKinematicVariables.InvJ0,
                                        rThisKinematicVariables.F);
    this->CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.F,
                     rThisKinematicVariables.DN_DX);

    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void GenericTotalLagrangianFemDemElement<2,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<2,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<2,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<2,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<2,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<2,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<2,6>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rF, rDN_DX);
}
template<>
void GenericTotalLagrangianFemDemElement<3,6>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rF, rDN_DX);
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Calculate2DB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(2, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
    }

    KRATOS_CATCH( "" )
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Calculate3DB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateStress(Vector& rStrain,
                                      std::size_t IntegrationPoint,
                                      Vector& rStress,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS | ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetStrainVector(rStrain);
    cl_params.SetStressVector(rStress);
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateStress(Matrix const& rF,
                                      std::size_t IntegrationPoint,
                                      Vector& rStress,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    Vector strain(mConstitutiveLawVector[IntegrationPoint]->GetStrainSize());
    CalculateStrain(rF, *mConstitutiveLawVector[IntegrationPoint], strain, rCurrentProcessInfo);
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS | ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetStrainVector(strain);
    cl_params.SetStressVector(rStress);
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateStrain(Matrix const& rF,
                                      std::size_t IntegrationPoint,
                                      Vector& rStrain,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.SetDeformationGradientF(rF);
    cl_params.SetStrainVector(rStrain);
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateShapeSensitivity(
    ShapeParameter Deriv,
    Matrix& rDN_DX0,
    Matrix& rDN_DX0_Deriv,
    Matrix& rF_Deriv,
    double& rDetJ0_Deriv,
    std::size_t IntegrationPointIndex
    )
{
    KRATOS_TRY;
    const unsigned ws_dim = GetGeometry().WorkingSpaceDimension();
    const unsigned ls_dim = GetGeometry().LocalSpaceDimension();
    Matrix J0(ws_dim, ls_dim);
    GeometryUtils::JacobianOnInitialConfiguration(
        GetGeometry(), GetGeometry().IntegrationPoints(this->GetIntegrationMethod())[IntegrationPointIndex], J0);
    const Matrix& rDN_De = GetGeometry().ShapeFunctionsLocalGradients()[IntegrationPointIndex];
    auto sensitivity_utility =
        GeometricalSensitivityUtility(J0, rDN_De);
    sensitivity_utility.CalculateSensitivity(Deriv, rDetJ0_Deriv, rDN_DX0_Deriv);
    rF_Deriv.resize(ws_dim, ws_dim, false);
    rF_Deriv.clear();
    for (unsigned i = 0; i < ws_dim; ++i)
        for (unsigned j = 0; j < ws_dim; ++j) {
            for (unsigned k = 0; k < GetGeometry().PointsNumber(); ++k)
                rF_Deriv(i, j) += GetGeometry()[k].Coordinates()[i] * rDN_DX0_Deriv(k, j);
        }
    for (unsigned j = 0; j < ws_dim; ++j)
      rF_Deriv(Deriv.Direction, j) += rDN_DX0(Deriv.NodeIndex, j);
    KRATOS_CATCH("");
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateBSensitivity(
    Matrix const& rDN_DX,
    Matrix const& rF,
    Matrix const& rDN_DX_Deriv,
    Matrix const& rF_Deriv,
    Matrix& rB_Deriv
    )
{
    KRATOS_TRY;
    const unsigned dimension = GetGeometry().WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();
    rB_Deriv.resize(strain_size, dimension * GetGeometry().PointsNumber(), false);
    Matrix B_deriv1(rB_Deriv.size1(), rB_Deriv.size2());
    Matrix B_deriv2(rB_Deriv.size1(), rB_Deriv.size2());

    CalculateB(B_deriv1, rF_Deriv, rDN_DX);
    CalculateB(B_deriv2, rF, rDN_DX_Deriv);
    rB_Deriv = B_deriv1 + B_deriv2;
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
std::size_t GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::GetStrainSize() const
{
    return GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        const std::size_t ws_dim = r_geom.WorkingSpaceDimension();
        const std::size_t nnodes = r_geom.PointsNumber();
        const std::size_t mat_dim = nnodes * ws_dim;
        rOutput.resize(mat_dim, mat_dim, false);
        rOutput.clear();
        Matrix F, F_deriv, DN_DX0_deriv, strain_tensor_deriv, DN_DX0, B, B_deriv;
        Matrix M_deriv;
        const auto strain_size = GetStrainSize();
        B.resize(strain_size, ws_dim * nnodes, false);
        Vector strain_vector_deriv(strain_size);
        Vector stress_vector(strain_size), stress_vector_deriv(strain_size);
        Vector residual_deriv(ws_dim * nnodes);
        Vector body_force, acceleration;
        double detJ0_deriv;
        LargeDisplacementKinematics large_disp_kinematics(r_geom, this->GetIntegrationMethod());
        for (std::size_t g = 0; g < r_geom.IntegrationPointsNumber(); ++g) {
            large_disp_kinematics.DeformationGradient(g, F);
            large_disp_kinematics.ShapeFunctionsGradients(g, DN_DX0);
            CalculateB(B, F, DN_DX0);
            CalculateStress(F, g, stress_vector, rCurrentProcessInfo);
            double weight = GetIntegrationWeight(
                r_geom.IntegrationPoints(this->GetIntegrationMethod()), g, large_disp_kinematics.DetJ0(g));
            const Vector& rN = row(GetGeometry().ShapeFunctionsValues(), g);
            body_force = GetBodyForce(r_geom.IntegrationPoints(this->GetIntegrationMethod()), g);

            for (auto s = ShapeParameter::Sequence(nnodes, ws_dim); s; ++s) {
                const auto& deriv = s.CurrentValue();
                CalculateShapeSensitivity(deriv, DN_DX0, DN_DX0_deriv, F_deriv, detJ0_deriv, g);
                CalculateGreenLagrangeStrainSensitivity(F, F_deriv, strain_tensor_deriv);
                noalias(strain_vector_deriv) =
                    MathUtils<double>::StrainTensorToVector(strain_tensor_deriv);
                CalculateStress(strain_vector_deriv, g, stress_vector_deriv, rCurrentProcessInfo);
                CalculateBSensitivity(DN_DX0, F, DN_DX0_deriv, F_deriv, B_deriv);
                const double weight_deriv =
                    GetIntegrationWeight(r_geom.IntegrationPoints(this->GetIntegrationMethod()), g, detJ0_deriv);
                noalias(residual_deriv) = -weight_deriv * prod(trans(B), stress_vector);
                noalias(residual_deriv) -= weight * prod(trans(B_deriv), stress_vector);
                noalias(residual_deriv) -= weight * prod(trans(B), stress_vector_deriv);
                CalculateAndAddExtForceContribution(
                    rN, rCurrentProcessInfo, body_force, residual_deriv, weight_deriv);
                for (std::size_t k = 0; k < residual_deriv.size(); ++k)
                    rOutput(deriv.NodeIndex * ws_dim + deriv.Direction, k) +=
                        residual_deriv(k);
            }
        }

        for (auto s = ShapeParameter::Sequence(nnodes, ws_dim); s; ++s) {
            const auto& deriv = s.CurrentValue();
            CalculateShapeGradientOfMassMatrix(M_deriv, deriv);
            GetSecondDerivativesVector(acceleration);
            noalias(residual_deriv) = -prod(M_deriv, acceleration);
            for (std::size_t k = 0; k < residual_deriv.size(); ++k)
                rOutput(deriv.NodeIndex * ws_dim + deriv.Direction, k) += residual_deriv(k);
        }
    }
    else
        KRATOS_ERROR << "Unsupported variable: " << rDesignVariable << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::IntegrateStressDamageMechanics(
    double& rThreshold,
    double& rDamage,
    const Vector& rStrainVector,
    const Vector& rStressVector,
    const int Edge,
    const double CharacteristicLength,
    ConstitutiveLaw::Parameters& rValues,
    bool& rIsDamaging
    )
{
    double uniaxial_stress;
    this->CalculateEquivalentStress(rStressVector, rStrainVector, uniaxial_stress, rValues);

    double initial_threshold;
    this->GetInitialUniaxialThreshold(rValues, initial_threshold);

    if (rThreshold < tolerance) {
        rThreshold = initial_threshold; // 1st iteration sets threshold as c_max
    } else if (initial_threshold > rThreshold) { // remeshing stuff
        rThreshold = initial_threshold;
    }

    const double F = uniaxial_stress - rThreshold;
    if (F <= tolerance) { // Elastic region --> Damage is constant
        rDamage = mDamages[Edge];
    } else {
        double damage_parameter; // A parameter
        this->CalculateDamageParameter(rValues, damage_parameter, CharacteristicLength);
        this->CalculateExponentialDamage(rDamage, damage_parameter, uniaxial_stress, initial_threshold);            
        rThreshold = uniaxial_stress;
        rIsDamaging = true;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericTotalLagrangianFemDemElement<3,0>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<3,1>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<3,2>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<3,3>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<3,4>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<3,5>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<3,6>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,0>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,1>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,2>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,3>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,4>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,5>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericTotalLagrangianFemDemElement<2,6>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateAverageVariableOnEdge2D(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector, // must be assigned to the current element outside this method
    const int edge
    )
{
    auto& r_elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_elem_neigb.size() == 0) << " Neighbour Elements not calculated" << std::endl;
    rAverageVector += r_elem_neigb[edge].GetValue(ThisVariable);
    rAverageVector *= 0.5;
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateAverageVariableOnEdge3D(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector, // must be assigned to the current element outside this method
    const int edge
    )
{
    std::vector<Element*> p_edge_neighbours = this->GetEdgeNeighbourElements(edge);
    int counter = 0;

    for (unsigned int elem = 0; elem < p_edge_neighbours.size(); elem++) {
        rAverageVector += p_edge_neighbours[elem]->GetValue(ThisVariable);
        counter++;
    }
    rAverageVector /= (counter + 1);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double GenericTotalLagrangianFemDemElement<3,0>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<3,1>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<3,2>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<3,3>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<3,4>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<3,5>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<3,6>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,0>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,1>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,2>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,3>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,4>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,5>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericTotalLagrangianFemDemElement<2,6>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}

template<unsigned int TDim, unsigned int TyieldSurf>
double GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateElementalDamage3D(const Vector& rEdgeDamages)
{
    // 7 modes of fracture of the tetrahedron
    Vector damage_mode_fracture = ZeroVector(7);
    const double one_third = 1.0 / 3.0;
    damage_mode_fracture[0] = one_third * (rEdgeDamages[0] + rEdgeDamages[1] + rEdgeDamages[2]);
    damage_mode_fracture[1] = one_third * (rEdgeDamages[0] + rEdgeDamages[3] + rEdgeDamages[4]);
    damage_mode_fracture[2] = one_third * (rEdgeDamages[1] + rEdgeDamages[3] + rEdgeDamages[5]);
    damage_mode_fracture[3] = 0.25 * (rEdgeDamages[1] + rEdgeDamages[2] + rEdgeDamages[3] + rEdgeDamages[4]);
    damage_mode_fracture[4] = 0.25 * (rEdgeDamages[0] + rEdgeDamages[1] + rEdgeDamages[4] + rEdgeDamages[5]);
    damage_mode_fracture[5] = one_third * (rEdgeDamages[2] + rEdgeDamages[4] + rEdgeDamages[5]);
    damage_mode_fracture[6] = 0.25 * (rEdgeDamages[0] + rEdgeDamages[2] + rEdgeDamages[3] + rEdgeDamages[5]);
    return ConstitutiveLawUtilities<VoigtSize>::GetMaxValue(damage_mode_fracture);
}

template<unsigned int TDim, unsigned int TyieldSurf>
double GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateElementalDamage2D(const Vector& rEdgeDamages)
{
    Vector two_max_values;
    ConstitutiveLawUtilities<VoigtSize>::Get2MaxValues(two_max_values, rEdgeDamages[0], rEdgeDamages[1], rEdgeDamages[2]);
    return 0.5*(two_max_values[0] + two_max_values[1]);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::InitializeInternalVariablesAfterMapping()
{
    // After the mapping, the thresholds of the edges (are equal to 0.0) are imposed equal to the IP threshold
    const double element_threhsold = mThreshold;
    if (norm_2(mThresholds) < tolerance) {
        for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
            mThresholds[edge] = element_threhsold;
        }
    }

    // IDEM with the edge damages
    const double damage_element = mDamage;
    if (norm_2(mDamages) < tolerance) {
        for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
            mDamages[edge] = damage_element;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericTotalLagrangianFemDemElement<3,0>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericTotalLagrangianFemDemElement<3,1>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericTotalLagrangianFemDemElement<3,2>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericTotalLagrangianFemDemElement<3,3>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericTotalLagrangianFemDemElement<3,4>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericTotalLagrangianFemDemElement<3,5>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericTotalLagrangianFemDemElement<3,6>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericTotalLagrangianFemDemElement<2,0>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericTotalLagrangianFemDemElement<2,1>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericTotalLagrangianFemDemElement<2,2>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericTotalLagrangianFemDemElement<2,3>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericTotalLagrangianFemDemElement<2,4>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericTotalLagrangianFemDemElement<2,5>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericTotalLagrangianFemDemElement<2,6>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::AuxComputeEdgeNeighbours(
    ProcessInfo& rCurrentProcessInfo
    )
{
    std::vector<std::vector<Element*>> edge_neighbours_container;
    auto& r_nodes_current_element = this->GetGeometry();

    auto& pNode0 = r_nodes_current_element[0];
    auto& pNode1 = r_nodes_current_element[1];
    auto& pNode2 = r_nodes_current_element[2];
    auto& pNode3 = r_nodes_current_element[3];

    // Neighbour elements of each node of the current element
    GlobalPointersVector<Element>& r_neigh_node_0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
    GlobalPointersVector<Element>& r_neigh_node_1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
    GlobalPointersVector<Element>& r_neigh_node_2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
    GlobalPointersVector<Element>& r_neigh_node_3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);

    // Nodal neighbours container
    std::vector<GlobalPointersVector<Element>> nodal_neighbours;
    nodal_neighbours.push_back(r_neigh_node_0);
    nodal_neighbours.push_back(r_neigh_node_1);
    nodal_neighbours.push_back(r_neigh_node_2);
    nodal_neighbours.push_back(r_neigh_node_3);

    // Aux indexes
    Matrix nodes_indexes = ZeroMatrix(6, 2);
    this->SetNodeIndexes(nodes_indexes);

    // Loop over EDGES to assign the elements that share that edge -> Fill mEdgeNeighboursContainer
    for (unsigned int edge = 0; edge < 6; edge++) {
        const int node_index_1 = nodes_indexes(edge, 0);
        const int node_index_2 = nodes_indexes(edge, 1);

        // Neigh elements of local node 1 and 2  //
        auto& r_neigh_of_node_1 = nodal_neighbours[node_index_1];
        auto& r_neigh_of_node_2 = nodal_neighbours[node_index_2];

        const int node_id_1 = r_nodes_current_element[node_index_1].Id();
        const int node_id_2 = r_nodes_current_element[node_index_2].Id();

        std::vector<Element*> edge_shared_elements_node_1;
        // Loop over neigh elements of the node 1
        for (unsigned int neigh_elem = 0; neigh_elem < r_neigh_of_node_1.size(); neigh_elem++) {
            // Nodes of the neigh element
            auto& r_nodes_neigh_elem = r_neigh_of_node_1[neigh_elem].GetGeometry();

            // Loop over the nodes of the neigh element
            for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
                const int neigh_element_node_id = r_nodes_neigh_elem[neigh_elem_node].Id();

                if (neigh_element_node_id == node_id_2 && this->Id() != r_neigh_of_node_1[neigh_elem].Id()) {
                    edge_shared_elements_node_1.push_back(&r_neigh_of_node_1[neigh_elem]); // ( [] returns an Element object!!)
                }
            }
        }

        std::vector<Element *> edge_shared_elements_node_2;
        // Loop over neigh elements of the node 2
        for (unsigned int neigh_elem = 0; neigh_elem < r_neigh_of_node_2.size(); neigh_elem++) {
            // Nodes of the neigh element
            auto &r_nodes_neigh_elem = r_neigh_of_node_2[neigh_elem].GetGeometry();

            // Loop over the nodes of the neigh element
            for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
                const int neigh_element_node_id = r_nodes_neigh_elem[neigh_elem_node].Id();

                if (neigh_element_node_id == node_id_1 && this->Id() != r_neigh_of_node_2[neigh_elem].Id()) {
                    edge_shared_elements_node_2.push_back(&r_neigh_of_node_2[neigh_elem]);
                }
            }
        }
        // Let's create the vector of neighbour elements for this edge
        std::vector<Element *> edge_shared_elements = edge_shared_elements_node_1;
        // Add the neigh elements from the node 2
        for (unsigned int i = 0; i < edge_shared_elements_node_2.size(); i++) {
            int aux = 0;

            for (unsigned int j = 0; j < edge_shared_elements.size(); j++) {
                if (edge_shared_elements_node_2[i]->Id() == edge_shared_elements[j]->Id())
                    aux++;
            }
            if (aux == 0)
                edge_shared_elements.push_back(edge_shared_elements_node_2[i]);
        }
        edge_neighbours_container.push_back(edge_shared_elements);
    } // End loop edges

    // Storages the information inside the element
    this->SaveEdgeNeighboursContainer(edge_neighbours_container);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericTotalLagrangianFemDemElement<2,0>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressModifiedMohrCoulomb(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,0>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressModifiedMohrCoulomb(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<2,1>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressRankine(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,1>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressRankine(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<2,2>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressSimoJu(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,2>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressSimoJu(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<2,3>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressDruckerPrager(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,3>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressDruckerPrager(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<2,4>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressHuberVonMises(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,4>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressHuberVonMises(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<2,5>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressTresca(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,5>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressTresca(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<2,6>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressMohrCoulomb(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericTotalLagrangianFemDemElement<3,6>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressMohrCoulomb(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericTotalLagrangianFemDemElement<2,0>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdModifiedMohrCoulomb(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,0>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdModifiedMohrCoulomb(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<2,1>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdRankine(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,1>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdRankine(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<2,2>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdSimoJu(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,2>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdSimoJu(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<2,3>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdDruckerPrager(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,3>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdDruckerPrager(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<2,4>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdHuberVonMises(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,4>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdHuberVonMises(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<2,5>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdTresca(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,5>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdTresca(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<2,6>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdMohrCoulomb(rValues, rThreshold);
}
template<>
void GenericTotalLagrangianFemDemElement<3,6>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdMohrCoulomb(rValues, rThreshold);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericTotalLagrangianFemDemElement<2,0>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterModifiedMohrCoulomb(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,0>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterModifiedMohrCoulomb(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<2,1>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterRankine(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,1>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterRankine(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<2,2>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterSimoJu(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,2>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterSimoJu(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<2,3>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterDruckerPrager(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,3>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterDruckerPrager(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<2,4>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterHuberVonMises(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,4>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterHuberVonMises(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<2,5>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterTresca(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,5>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterTresca(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<2,6>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterMohrCoulomb(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericTotalLagrangianFemDemElement<3,6>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterMohrCoulomb(rValues, rAParameter, CharacteristicLength);
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
double GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateCharacteristicLength(
    GenericTotalLagrangianFemDemElement<TDim,TyieldSurf> *pCurrentElement
    )
{
    auto& r_geometry = pCurrentElement->GetGeometry();
    const auto& r_edges = r_geometry.Edges();

    double sum_of_lengths = 0.0;
    for (IndexType i = 0; i < NumberOfEdges; ++i) {
        auto& node_1 = r_edges[i][0];
        auto& node_2 = r_edges[i][1];
        auto coordinates_1 = node_1.GetInitialPosition();
        auto coordinates_2 = node_2.GetInitialPosition();
        sum_of_lengths += std::sqrt(std::pow(coordinates_1[0] - coordinates_2[0], 2) + 
            std::pow(coordinates_1[1] - coordinates_2[1], 2) + std::pow(coordinates_1[2] - coordinates_2[2], 2));
    }
    return sum_of_lengths / NumberOfEdges;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateExponentialDamage(
    double& rDamage,
    const double DamageParameter,
    const double UniaxialStress,
    const double InitialThrehsold
    )
{
    rDamage = 1.0 - (InitialThrehsold / UniaxialStress) * std::exp(DamageParameter *
             (1.0 - UniaxialStress / InitialThrehsold)); // Exponential softening law
    if (rDamage > 0.999) rDamage = 0.999;
}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void  GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateTangentTensor(
    Matrix& rTangentTensor,
    const Vector& rStrainVectorGP,
    const Vector& rStressVectorGP,
    const Matrix& rDeformationGradientGP,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double number_components = rStrainVectorGP.size();
    rTangentTensor.resize(number_components, number_components);
    Vector perturbed_stress, perturbed_strain;
    perturbed_strain.resize(number_components);
    perturbed_stress.resize(number_components);
    Matrix perturbed_deformation_gradient;
    perturbed_deformation_gradient.resize(number_components, number_components);
    const double size_1 = rDeformationGradientGP.size1();
    const double size_2 = rDeformationGradientGP.size2();
    
    for (unsigned int i_component = 0; i_component < size_1; i_component++) {
        for (unsigned int j_component = i_component; j_component < size_2; j_component++) {
            double perturbation;
            const int component_voigt_index = this->CalculateVoigtIndex(number_components, i_component, j_component);
            this->CalculatePerturbation(rStrainVectorGP, perturbation, component_voigt_index);
            this->PerturbateDeformationGradient(perturbed_deformation_gradient, rDeformationGradientGP, perturbation, i_component, j_component);
            this->CalculateGreenLagrangeStrainVector(perturbed_strain, perturbed_deformation_gradient);
            this->IntegratePerturbedStrain(perturbed_stress, perturbed_strain, rElasticMatrix, rValues);
            const Vector& r_delta_stress = perturbed_stress - rStressVectorGP;
            this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress, perturbation, component_voigt_index);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void  GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateTangentTensorSecondOrder(
    Matrix& rTangentTensor,
    const Vector& rStrainVectorGP,
    const Vector& rStressVectorGP,
    const Matrix& rDeformationGradientGP,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double number_components = rStrainVectorGP.size();
    rTangentTensor.resize(number_components, number_components);
    Vector perturbed_stress_plus, perturbed_strain_plus;
    Vector perturbed_stress_minus, perturbed_strain_minus;

    perturbed_stress_plus.resize(number_components);
    perturbed_strain_plus.resize(number_components);
    perturbed_stress_minus.resize(number_components);
    perturbed_strain_minus.resize(number_components);

    Matrix perturbed_deformation_gradient_plus, perturbed_deformation_gradient_minus;
    perturbed_deformation_gradient_plus.resize(number_components, number_components);
    perturbed_deformation_gradient_minus.resize(number_components, number_components);
    const double size_1 = rDeformationGradientGP.size1();
    const double size_2 = rDeformationGradientGP.size2();
    
    for (unsigned int i_component = 0; i_component < size_1; i_component++) {
        for (unsigned int j_component = i_component; j_component < size_2; j_component++) {
            double perturbation;
            const int component_voigt_index = this->CalculateVoigtIndex(number_components, i_component, j_component);
            this->CalculatePerturbation(rStrainVectorGP, perturbation, component_voigt_index);

            this->PerturbateDeformationGradient(perturbed_deformation_gradient_plus, rDeformationGradientGP, perturbation, i_component, j_component);
            this->CalculateGreenLagrangeStrainVector(perturbed_strain_plus, perturbed_deformation_gradient_plus);
            this->IntegratePerturbedStrain(perturbed_stress_plus, perturbed_strain_plus, rElasticMatrix, rValues);

            this->PerturbateDeformationGradient(perturbed_deformation_gradient_minus, rDeformationGradientGP, -perturbation, i_component, j_component);
            this->CalculateGreenLagrangeStrainVector(perturbed_strain_minus, perturbed_deformation_gradient_minus);
            this->IntegratePerturbedStrain(perturbed_stress_minus, perturbed_strain_minus, rElasticMatrix, rValues);

            const Vector& r_delta_stress = perturbed_stress_plus - perturbed_stress_minus;
            this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress,2.0*perturbation, component_voigt_index);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::PerturbateDeformationGradient(
    Matrix& rPerturbedDeformationGradient,
    const Matrix& rDeformationGradientGP,
    const double Perturbation,
    const int ComponentI,
    const int ComponentJ
    )
{
    rPerturbedDeformationGradient = rDeformationGradientGP;
    if (ComponentI == ComponentJ) {
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += Perturbation;
    } else {
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += 0.5 * Perturbation;
        rPerturbedDeformationGradient(ComponentJ, ComponentI) += 0.5 * Perturbation;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculatePerturbation(
    const Vector& rStrainVectorGP,
    double& rPerturbation,
    const int Component
    )
{
    double perturbation_1, perturbation_2;
    if (std::abs(rStrainVectorGP[Component]) > tolerance) {
        perturbation_1 = 1.0e-7 * rStrainVectorGP[Component];
    } else {
        double min_strain_component = ConstitutiveLawUtilities<VoigtSize>::GetMinAbsValue(rStrainVectorGP);
        perturbation_1 = 1.0e-7 * min_strain_component;
    }
    const double max_strain_component = ConstitutiveLawUtilities<VoigtSize>::GetMaxAbsValue(rStrainVectorGP);
    perturbation_2 = 1.0e-12 * max_strain_component;
    rPerturbation = std::max(perturbation_1, perturbation_2);
    if (rPerturbation < 1e-9) rPerturbation = 1e-9;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::PerturbateStrainVector(
    Vector& rPerturbedStrainVector,
    const Vector& rStrainVectorGP,
    const double Perturbation,
    const int Component
    )
{
    noalias(rPerturbedStrainVector) = rStrainVectorGP;
    rPerturbedStrainVector[Component] += Perturbation;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::IntegratePerturbedStrain(
    Vector& rPerturbedStressVector,
    const Vector& rPerturbedStrainVector,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Vector& r_perturbed_predictive_stress = prod(rElasticMatrix, rPerturbedStrainVector);
    Vector damages_edges = ZeroVector(NumberOfEdges);
    const double characteristic_length = this->CalculateCharacteristicLength(this);
    bool dummy = false;
    Vector average_stress_edge(VoigtSize), average_strain_edge(VoigtSize);

    for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
        average_stress_edge = r_perturbed_predictive_stress;
        average_strain_edge = rPerturbedStrainVector;
        this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
        this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

        double damage_edge = mDamages[edge];
        double threshold = mThresholds[edge];

        this->IntegrateStressDamageMechanics(threshold, damage_edge, average_strain_edge, 
                                             average_stress_edge, edge, characteristic_length, 
                                             rValues, dummy);

        damages_edges[edge] = damage_edge;
    } // Loop edges
    const double damage_element = this->CalculateElementalDamage(damages_edges);
    rPerturbedStressVector = (1.0 - damage_element) * r_perturbed_predictive_stress;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::AssignComponentsToTangentTensor(
    Matrix& rTangentTensor,
    const Vector& rDeltaStress,
    const double Perturbation,
    const int Component
    )
{
    const IndexType voigt_size = rDeltaStress.size();
    for (IndexType row = 0; row < voigt_size; ++row) {
        rTangentTensor(row, Component) = rDeltaStress[row] / Perturbation;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
int GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateVoigtIndex(
    const SizeType VoigtSize,
    const int ComponentI,
    const int ComponentJ
    )
{
    if (VoigtSize == 6) {
        switch(ComponentI) {
            case 0:
                switch(ComponentJ) {
                    case 0:
                        return 0;
                    case 1:
                        return 3;
                    case 2:
                        return 5;
                    default:
                        return 0;
                }
            case 1:
                switch(ComponentJ) {
                    case 0:
                        return 3;
                    case 1:
                        return 1;
                    case 2:
                        return 4;
                    default:
                        return 0;
                }
            case 2:
                switch(ComponentJ) {
                    case 0:
                        return 5;
                    case 1:
                        return 4;
                    case 2:
                        return 2;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    } else {
        switch(ComponentI) {
            case 0:
                switch(ComponentJ) {
                    case 0:
                        return 0;
                    case 1:
                        return 2;
                    default:
                        return 0;
                }
            case 1:
                switch(ComponentJ) {
                    case 0:
                        return 2;
                    case 1:
                        return 1;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateGreenLagrangeStrainVector(
    Vector& rStrainVector,
    const Matrix& rF
    )
{
    Matrix strain_tensor;
    if (strain_tensor.size1() != TDim)
        strain_tensor.resize(TDim, TDim);
    
    if (rStrainVector.size() != VoigtSize)
        rStrainVector.resize(VoigtSize);

    Matrix identity = identity_matrix<double>(TDim);
    noalias(strain_tensor) = 0.5 * (prod(trans(rF), rF) - identity);
    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(strain_tensor, rStrainVector.size());
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::GetValueOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DAMAGE_ELEMENT || 
    rVariable == IS_DAMAGED || 
    rVariable == STRESS_THRESHOLD || 
    rVariable == EQUIVALENT_STRESS_VM) {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DAMAGE_ELEMENT) {
        rOutput.resize(1);
        for (unsigned int point_number = 0; point_number < 1; point_number++) {
            rOutput[point_number] = mDamage;
        }
    } else if (rVariable == STRESS_THRESHOLD) {
        rOutput.resize(1);
        for (unsigned int point_number = 0; point_number < 1; point_number++) {
            rOutput[point_number] = mThreshold;
        }
    }
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    if (rVariable == STRESS_VECTOR_INTEGRATED) {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
        const SizeType number_of_integration_points = integration_points.size();
        if ( rOutput.size() != number_of_integration_points )
            rOutput.resize( number_of_integration_points );
        // Create and initialize element variables:
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& cl_options=cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            //call the constitutive law to update material variables
            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);

            if ( rOutput[point_number].size() != strain_size )
                rOutput[point_number].resize(strain_size, false);
            rOutput[point_number] = (1.0-mDamage)*this_constitutive_variables.StressVector;
        }
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    std::vector<Vector> stress_vector;

    if (rVariable == STRESS_TENSOR_INTEGRATED) {
        this->CalculateOnIntegrationPoints(STRESS_VECTOR_INTEGRATED, stress_vector, rCurrentProcessInfo);

        // Loop integration points
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            if (rOutput[point_number].size2() != dimension)
                rOutput[point_number].resize( dimension, dimension, false);

            rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
        }
    }

}


/***********************************************************************************/
/***********************************************************************************/
template class GenericTotalLagrangianFemDemElement<2,0>;
template class GenericTotalLagrangianFemDemElement<2,1>;
template class GenericTotalLagrangianFemDemElement<2,2>;
template class GenericTotalLagrangianFemDemElement<2,3>;
template class GenericTotalLagrangianFemDemElement<2,4>;
template class GenericTotalLagrangianFemDemElement<2,5>;
template class GenericTotalLagrangianFemDemElement<2,6>;
template class GenericTotalLagrangianFemDemElement<3,0>;
template class GenericTotalLagrangianFemDemElement<3,1>;
template class GenericTotalLagrangianFemDemElement<3,2>;
template class GenericTotalLagrangianFemDemElement<3,3>;
template class GenericTotalLagrangianFemDemElement<3,4>;
template class GenericTotalLagrangianFemDemElement<3,5>;
template class GenericTotalLagrangianFemDemElement<3,6>;
} // Namespace Kratos


