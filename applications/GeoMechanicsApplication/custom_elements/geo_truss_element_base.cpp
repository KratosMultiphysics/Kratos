// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter,
//                   Vahid Galavi
//

// System includes

// External includes

// Project includes
#include "custom_elements/geo_truss_element_base.hpp"
#include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementBase<TDim, TNumNodes>::GeoTrussElementBase(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementBase<TDim, TNumNodes>::GeoTrussElementBase(IndexType               NewId,
                                                          GeometryType::Pointer   pGeometry,
                                                          PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElementBase<TDim, TNumNodes>::Create(IndexType             NewId,
                                                              NodesArrayType const& rThisNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoTrussElementBase>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElementBase<TDim, TNumNodes>::Create(IndexType             NewId,
                                                              GeometryType::Pointer pGeom,
                                                              PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoTrussElementBase>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementBase<TDim, TNumNodes>::~GeoTrussElementBase()
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
            mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        } else {
            KRATOS_ERROR << "A constitutive law needs to be specified for the "
                            "element with ID "
                         << Id() << std::endl;
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CreateElementStiffnessMatrix(MatrixType& rLocalStiffnessMatrix,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    CalculateElasticStiffnessMatrix(rLocalStiffnessMatrix, rCurrentProcessInfo);

    FullDofMatrixType K_geo;
    CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

    noalias(rLocalStiffnessMatrix) += K_geo;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this, rDampingMatrix, rCurrentProcessInfo, TDim * TNumNodes);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Clear matrix
    if (rMassMatrix.size1() != TDim * TNumNodes || rMassMatrix.size2() != TDim * TNumNodes) {
        rMassMatrix.resize(TDim * TNumNodes, TDim * TNumNodes, false);
    }
    rMassMatrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);

    if (StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(GetProperties(), rCurrentProcessInfo)) {
        // Compute lumped mass matrix
        VectorType temp_vector(TDim * TNumNodes);
        CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);

        // Fill the matrix
        for (IndexType i = 0; i < TDim * TNumNodes; ++i) {
            rMassMatrix(i, i) = temp_vector[i];
        }
    } else {
        CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateConsistentMassMatrix(MatrixType& rMassMatrix,
                                                                         const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const double A   = GetProperties()[CROSS_AREA];
    const double L   = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = GetProperties()[DENSITY];

    const double factor = TDim * TNumNodes;
    const BoundedMatrix<double, TDim, TDim> fill_matrix = A * L * rho * (IdentityMatrix(TDim) / factor);

    project(rMassMatrix, range(0, TDim), range(0, TDim)) += fill_matrix * 2.0;
    project(rMassMatrix, range(0, TDim), range(TDim, NDof)) += fill_matrix;
    project(rMassMatrix, range(TDim, NDof), range(0, TDim)) += fill_matrix;
    project(rMassMatrix, range(TDim, NDof), range(TDim, NDof)) += fill_matrix * 2.0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateBodyForces(FullDofVectorType& rGlobalBodyForces)
{
    KRATOS_TRY

    // getting shapefunctionvalues
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

    // creating necessary values
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    double total_mass = A * L * rho;

    rGlobalBodyForces = ZeroVector(TDim * TNumNodes);

    // assemble global Vector
    BoundedVector<double, 3> body_forces_node = ZeroVector(3);
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        body_forces_node = total_mass * GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
                           Ncontainer(0, i);

        for (unsigned int j = 0; j < TDim; ++j) {
            rGlobalBodyForces[(i * TDim) + j] = body_forces_node[j];
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSolutionStepValues(GetDofs(), Step);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractFirstTimeDerivatives(GetDofs(), Step);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSecondTimeDerivatives(GetDofs(), Step);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                                VectorType& rRightHandSideVector,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rRightHandSideVector = ZeroVector(TDim * TNumNodes);

    FullDofVectorType internal_forces;
    UpdateInternalForces(internal_forces, rCurrentProcessInfo);

    noalias(rRightHandSideVector) -= internal_forces;

    FullDofVectorType GlobalBodyForces;
    CalculateBodyForces(GlobalBodyForces);

    noalias(rRightHandSideVector) += GlobalBodyForces;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // creating LHS
    CreateElementStiffnessMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::Calculate(const Variable<Matrix>& rVariable,
                                                     Matrix&                 rOutput,
                                                     const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == MEMBRANE_PRESTRESS) {
        std::vector<Vector> prestress_matrix;
        CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, prestress_matrix, rCurrentProcessInfo);
        const int manual_integraton_points_size(1);
        rOutput = ZeroMatrix(3, manual_integraton_points_size);

        // each column represents 1 GP
        for (SizeType i = 0; i < manual_integraton_points_size; ++i) {
            column(rOutput, i) = prestress_matrix[i];
        }
    }
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::Calculate(const Variable<double>& rVariable,
                                                     double&                 rOutput,
                                                     const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STRAIN_ENERGY) {
        const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double A  = GetProperties()[CROSS_AREA];

        // material strain energy
        double strain_energy(0.00);
        Vector strain_vector = ZeroVector(mpConstitutiveLaw->GetStrainSize());
        strain_vector[0]     = CalculateGreenLagrangeStrain();
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        Values.SetStrainVector(strain_vector);
        mpConstitutiveLaw->CalculateValue(Values, STRAIN_ENERGY, strain_energy);

        // constant pre_stress strain energy
        if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            const double pre_stress = GetProperties()[TRUSS_PRESTRESS_PK2];
            strain_energy += pre_stress * strain_vector[0];
        }

        // integration over reference domain
        strain_energy *= A * L0;

        rOutput = strain_energy;
    } else if (rVariable == KINETIC_ENERGY) {
        Matrix mass_matrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);
        CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
        Vector current_nodal_velocities = ZeroVector(TDim * TNumNodes);
        GetFirstDerivativesVector(current_nodal_velocities);
        rOutput = 0.50 * inner_prod(current_nodal_velocities, prod(mass_matrix, current_nodal_velocities));
    } else if (rVariable == ENERGY_DAMPING_DISSIPATION) {
        // Attention! this is only the current state and must be integrated over time (*dt)
        Matrix damping_matrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);
        CalculateDampingMatrix(damping_matrix, rCurrentProcessInfo);
        Vector current_nodal_velocities = ZeroVector(TDim * TNumNodes);
        GetFirstDerivativesVector(current_nodal_velocities);
        rOutput = inner_prod(current_nodal_velocities, prod(damping_matrix, current_nodal_velocities));
    } else if (rVariable == EXTERNAL_ENERGY) {
        // Dead Load contribution to external energy

        FullDofVectorType dead_load_rhs;
        CalculateBodyForces(dead_load_rhs);
        Vector current_nodal_displacements = ZeroVector(TDim * TNumNodes);
        GetValuesVector(current_nodal_displacements, 0);
        rOutput = inner_prod(dead_load_rhs, current_nodal_displacements);
    }
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                        std::vector<double>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == TRUSS_PRESTRESS_PK2) {
        rOutput[0] = 0.00;
        if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            rOutput[0] = GetProperties()[TRUSS_PRESTRESS_PK2];
        }
    }
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
        const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        rOutput[0]      = l / L0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                        std::vector<Vector>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        Vector strain = ZeroVector(TDim);
        strain[0]     = CalculateGreenLagrangeStrain();
        rOutput[0]    = strain;
    }

    if ((rVariable == CAUCHY_STRESS_VECTOR) || (rVariable == PK2_STRESS_VECTOR)) {
        const double prestress =
            (GetProperties().Has(TRUSS_PRESTRESS_PK2) ? GetProperties()[TRUSS_PRESTRESS_PK2] : 0.0);

        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        Vector                      temp_strain = ZeroVector(1);
        Vector                      temp_stress = ZeroVector(1);
        temp_strain[0]                          = CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);
        Values.SetStressVector(temp_stress);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
        const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

        temp_stress[0] += prestress;
        rOutput[0] = temp_stress;

        if (rVariable == CAUCHY_STRESS_VECTOR) {
            rOutput[0] *= l / L0;
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                        std::vector<array_1d<double, 3>>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == FORCE) {
        std::vector<Vector> array_output;
        CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, array_output, rCurrentProcessInfo);

        array_1d<double, 3> temp_internal_stresses = ZeroVector(3);
        temp_internal_stresses[0]                  = array_output[0][0];

        rOutput[0] = temp_internal_stresses * GetProperties()[CROSS_AREA];
    }
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int GeoTrussElementBase<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const double   numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    if (mpConstitutiveLaw != nullptr) {
        mpConstitutiveLaw->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);
    }

    if (dimension != TDim || number_of_nodes != TNumNodes) {
        KRATOS_ERROR << "Wrong dimension or number of nodes in the truss "
                        "element."
                     << std::endl
                     << "This element works with dimension: " << TDim << std::endl
                     << "and number of nodes: " << TNumNodes << std::endl;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    if constexpr (TDim > 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const NodeType& rnode = GetGeometry()[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode);

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode);
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode);
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode);
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const NodeType& rnode = GetGeometry()[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode);

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode);
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode);
        }
    }

    if (GetProperties().Has(CROSS_AREA) == false || GetProperties()[CROSS_AREA] <= numerical_limit) {
        KRATOS_ERROR << "CROSS_AREA not provided for this element" << Id() << std::endl;
    }

    if (GetProperties().Has(DENSITY) == false) {
        KRATOS_ERROR << "DENSITY not provided for this element" << Id() << std::endl;
    }

    KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this) < std::numeric_limits<double>::epsilon())
        << "Element #" << Id() << " has a length of zero!" << std::endl;

    if (TDim == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (GetGeometry()[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << GetGeometry()[i].Id()
                             << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoTrussElementBase<TDim, TNumNodes>::CalculateGreenLagrangeStrain() const
{
    KRATOS_TRY

    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double e = ((l * l - L * L) / (2.00 * L * L));
    return e;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::UpdateInternalForces(FullDofVectorType& rInternalForces,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    FullDofMatrixType transformation_matrix;
    CreateTransformationMatrix(transformation_matrix);

    const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double A  = GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    Vector                      temp_strain = ZeroVector(1);
    Vector                      temp_stress = ZeroVector(1);
    temp_strain[0]                          = CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    const double normal_force = ((temp_stress[0] + prestress) * l * A) / L0;

    // internal force vectors
    FullDofVectorType f_local = ZeroVector(TDim * TNumNodes);
    f_local[0]                = -1.00 * normal_force;
    f_local[TDim]             = 1.00 * normal_force;
    rInternalForces           = ZeroVector(TDim * TNumNodes);
    noalias(rInternalForces)  = prod(transformation_matrix, f_local);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<3, 2>::WriteTransformationCoordinates(FullDofVectorType& rReferenceCoordinates)
{
    KRATOS_TRY

    Vector current_displacement = ZeroVector(DIM * NUM_NODES);
    GetValuesVector(current_displacement, 0);
    rReferenceCoordinates[0] = GetGeometry()[0].X0() + current_displacement[0];
    rReferenceCoordinates[1] = GetGeometry()[0].Y0() + current_displacement[1];
    rReferenceCoordinates[2] = GetGeometry()[0].Z0() + current_displacement[2];
    rReferenceCoordinates[3] = GetGeometry()[1].X0() + current_displacement[3];
    rReferenceCoordinates[4] = GetGeometry()[1].Y0() + current_displacement[4];
    rReferenceCoordinates[5] = GetGeometry()[1].Z0() + current_displacement[5];

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<2, 2>::WriteTransformationCoordinates(FullDofVectorType& rReferenceCoordinates)
{
    KRATOS_TRY

    Vector current_displacement = ZeroVector(DIM * NUM_NODES);
    GetValuesVector(current_displacement, 0);
    rReferenceCoordinates[0] = GetGeometry()[0].X0() + current_displacement[0];
    rReferenceCoordinates[1] = GetGeometry()[0].Y0() + current_displacement[1];
    rReferenceCoordinates[2] = GetGeometry()[1].X0() + current_displacement[2];
    rReferenceCoordinates[3] = GetGeometry()[1].Y0() + current_displacement[3];

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<3, 2>::CreateTransformationMatrix(FullDofMatrixType& rRotationMatrix)
{
    KRATOS_TRY

    // 1st calculate transformation matrix
    using arraydim                          = BoundedVector<double, DIM>;
    arraydim          direction_vector_x    = ZeroVector(DIM);
    arraydim          direction_vector_y    = ZeroVector(DIM);
    arraydim          direction_vector_z    = ZeroVector(DIM);
    FullDofVectorType reference_coordinates = ZeroVector(DIM * NUM_NODES);
    arraydim          global_z_vector       = ZeroVector(DIM);
    global_z_vector[2]                      = 1.0;

    this->WriteTransformationCoordinates(reference_coordinates);

    for (unsigned int i = 0; i < DIM; ++i) {
        direction_vector_x[i] = (reference_coordinates[i + DIM] - reference_coordinates[i]);
    }

    // local x-axis (e1_local) is the beam axis  (in GID is e3_local)
    const double numeric_limit = std::numeric_limits<double>::epsilon();
    double       VectorNorm;
    VectorNorm = MathUtils<double>::Norm(direction_vector_x);
    if (VectorNorm > numeric_limit) {
        direction_vector_x /= VectorNorm;
    } else {
        KRATOS_ERROR << "length of element " << Id() << " is ~ zero: " << VectorNorm << std::endl;
    }

    if (std::abs(direction_vector_x[2] - 1.00) <= numeric_limit) {
        direction_vector_y[1] = 1.0;
        direction_vector_z[0] = -1.0;
    }

    else if (std::abs(direction_vector_x[2] + 1.00) <= numeric_limit) {
        direction_vector_y[1] = 1.0;
        direction_vector_z[0] = 1.0;
    }

    else {
        MathUtils<double>::UnitCrossProduct(direction_vector_y, direction_vector_x, global_z_vector);
        MathUtils<double>::UnitCrossProduct(direction_vector_z, direction_vector_y, direction_vector_x);
    }

    // 2nd fill big rotation matrix
    BoundedMatrix<double, DIM, DIM> current_coordinate_system = ZeroMatrix(DIM, DIM);
    for (unsigned int i = 0; i < DIM; ++i) {
        current_coordinate_system(i, 0) = direction_vector_x[i];
        current_coordinate_system(i, 1) = direction_vector_y[i];
        current_coordinate_system(i, 2) = direction_vector_z[i];
    }

    rRotationMatrix = ZeroMatrix(DIM * NUM_NODES, DIM * NUM_NODES);
    // Building the rotation matrix for the local element matrix
    for (unsigned int kk = 0; kk < DIM * NUM_NODES; kk += DIM) {
        for (unsigned int i = 0; i < DIM; ++i) {
            for (unsigned int j = 0; j < DIM; ++j) {
                rRotationMatrix(i + kk, j + kk) = current_coordinate_system(i, j);
            }
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<2, 2>::CreateTransformationMatrix(FullDofMatrixType& rRotationMatrix)
{
    KRATOS_TRY

    static constexpr unsigned int DIM_2D = 2;
    static constexpr unsigned int DIM_3D = 3;

    // 1st calculate transformation matrix

    BoundedVector<double, DIM_2D * NUM_NODES> reference_coordinates_2D = ZeroVector(DIM_2D * NUM_NODES);
    this->WriteTransformationCoordinates(reference_coordinates_2D);

    BoundedVector<double, DIM_3D> direction_vector_x_3D = ZeroVector(DIM_3D);
    for (unsigned int i = 0; i < DIM_2D; ++i) {
        direction_vector_x_3D[i] = (reference_coordinates_2D[i + DIM_2D] - reference_coordinates_2D[i]);
    }

    // local x-axis (e1_local) is the beam axis  (in GID is e3_local)
    double VectorNorm;
    VectorNorm = MathUtils<double>::Norm(direction_vector_x_3D);
    if (VectorNorm > std::numeric_limits<double>::epsilon()) {
        direction_vector_x_3D /= VectorNorm;
    } else {
        KRATOS_ERROR << "length of element " << Id() << " is ~ zero: " << VectorNorm << std::endl;
    }

    BoundedVector<double, DIM_3D> direction_vector_y_3D = ZeroVector(DIM_3D);
    BoundedVector<double, DIM_3D> global_z_vector_3D    = ZeroVector(DIM_3D);
    global_z_vector_3D[INDEX_Z]                         = 1.0;
    MathUtils<double>::UnitCrossProduct(direction_vector_y_3D, direction_vector_x_3D, global_z_vector_3D);

    // 2nd fill big rotation matrix
    BoundedMatrix<double, DIM_2D, DIM_2D> current_coordinate_system = ZeroMatrix(DIM_2D, DIM_2D);

    for (unsigned int i = 0; i < DIM_2D; ++i) {
        current_coordinate_system(i, 0) = direction_vector_x_3D[i];
        current_coordinate_system(i, 1) = direction_vector_y_3D[i];
    }

    rRotationMatrix = ZeroMatrix(DIM_2D * NUM_NODES, DIM_2D * NUM_NODES);
    // copying the 3D rotation matrix to 2D matrix
    for (unsigned int kk = 0; kk < DIM_2D * NUM_NODES; kk += DIM_2D) {
        for (unsigned int i = 0; i < DIM_2D; ++i) {
            for (unsigned int j = 0; j < DIM_2D; ++j) {
                rRotationMatrix(i + kk, j + kk) = current_coordinate_system(i, j);
            }
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::AddExplicitContribution(const VectorType& rRHSVector,
                                                                   const Variable<VectorType>& rRHSVariable,
                                                                   const Variable<double>& rDestinationVariable,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();

    if (rDestinationVariable == NODAL_MASS) {
        VectorType element_mass_vector(TDim * TNumNodes);
        CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);

        for (SizeType i = 0; i < TNumNodes; ++i) {
            double& r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int     index        = i * TDim;

            AtomicAdd(r_nodal_mass, element_mass_vector(index));
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::AddExplicitContribution(const VectorType& rRHSVector,
                                                                   const Variable<VectorType>& rRHSVariable,
                                                                   const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        FullDofVectorType damping_residual_contribution = ZeroVector(TDim * TNumNodes);
        Vector            current_nodal_velocities      = ZeroVector(TDim * TNumNodes);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        CalculateDampingMatrix(damping_matrix, rCurrentProcessInfo);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (size_t i = 0; i < TNumNodes; ++i) {
            size_t index = TDim * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < TDim; ++j) {
                AtomicAdd(r_force_residual[j],
                          (rRHSVector[index + j] - damping_residual_contribution[index + j]));
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {
        // Getting the vector mass
        VectorType mass_vector(TDim * TNumNodes);
        CalculateLumpedMassVector(mass_vector, rCurrentProcessInfo);

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            double& r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
            int     index        = i * TDim;

            AtomicAdd(r_nodal_mass, mass_vector[index]);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<3, 2>::CalculateGeometricStiffnessMatrix(FullDofMatrixType& rGeometricStiffnessMatrix,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    double       E = ReturnTangentModulus1D(rCurrentProcessInfo);
    const double A = GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    rGeometricStiffnessMatrix = ZeroMatrix(DIM * NUM_NODES, DIM * NUM_NODES);

    // du... delta displacement in x-direction
    // dv... delta displacement in y-direction
    // dw... delta displacement in z-direction
    // L... inital member length
    // l... deformed member length
    // e_gl... green_lagrange strain

    double du = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
                GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    double dv = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
                GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    double dw = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
                GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
    const double dx   = GetGeometry()[1].X0() - GetGeometry()[0].X0();
    const double dy   = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
    const double dz   = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
    const double L    = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l    = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    double       e_gL = (l * l - L * L) / (2.00 * L * L);
    const double L3   = L * L * L;

    const double K_sigma = ((E * A * e_gL) / L) + ((prestress * A) / L);
    const double K_uij   = (E * A) / L3;

    rGeometricStiffnessMatrix(0, 0) = K_sigma + K_uij * (2 * du * dx + du * du);
    rGeometricStiffnessMatrix(3, 3) = rGeometricStiffnessMatrix(0, 0);

    rGeometricStiffnessMatrix(1, 1) = K_sigma + K_uij * (2 * dv * dy + dv * dv);
    rGeometricStiffnessMatrix(4, 4) = rGeometricStiffnessMatrix(1, 1);

    rGeometricStiffnessMatrix(2, 2) = K_sigma + K_uij * (2 * dw * dz + dw * dw);
    rGeometricStiffnessMatrix(5, 5) = rGeometricStiffnessMatrix(2, 2);

    rGeometricStiffnessMatrix(0, 1) = K_uij * (dx * dv + dy * du + du * dv);
    rGeometricStiffnessMatrix(1, 0) = rGeometricStiffnessMatrix(0, 1);

    rGeometricStiffnessMatrix(0, 2) = K_uij * (dx * dw + dz * du + du * dw);
    rGeometricStiffnessMatrix(2, 0) = rGeometricStiffnessMatrix(0, 2);

    rGeometricStiffnessMatrix(0, 3) = -rGeometricStiffnessMatrix(0, 0);
    rGeometricStiffnessMatrix(3, 0) = rGeometricStiffnessMatrix(0, 3);

    rGeometricStiffnessMatrix(0, 4) = -rGeometricStiffnessMatrix(0, 1);
    rGeometricStiffnessMatrix(4, 0) = rGeometricStiffnessMatrix(0, 4);

    rGeometricStiffnessMatrix(0, 5) = -rGeometricStiffnessMatrix(0, 2);
    rGeometricStiffnessMatrix(5, 0) = rGeometricStiffnessMatrix(0, 5);

    rGeometricStiffnessMatrix(1, 2) = K_uij * (dy * dw + dz * dv + dv * dw);
    rGeometricStiffnessMatrix(2, 1) = rGeometricStiffnessMatrix(1, 2);

    rGeometricStiffnessMatrix(1, 3) = rGeometricStiffnessMatrix(0, 4);
    rGeometricStiffnessMatrix(3, 1) = rGeometricStiffnessMatrix(1, 3);

    rGeometricStiffnessMatrix(1, 4) = -rGeometricStiffnessMatrix(1, 1);
    rGeometricStiffnessMatrix(4, 1) = rGeometricStiffnessMatrix(1, 4);

    rGeometricStiffnessMatrix(1, 5) = -rGeometricStiffnessMatrix(1, 2);
    rGeometricStiffnessMatrix(5, 1) = rGeometricStiffnessMatrix(1, 5);

    rGeometricStiffnessMatrix(2, 3) = -rGeometricStiffnessMatrix(0, 2);
    rGeometricStiffnessMatrix(3, 2) = rGeometricStiffnessMatrix(2, 3);

    rGeometricStiffnessMatrix(2, 4) = -rGeometricStiffnessMatrix(1, 2);
    rGeometricStiffnessMatrix(4, 2) = rGeometricStiffnessMatrix(2, 4);

    rGeometricStiffnessMatrix(2, 5) = -rGeometricStiffnessMatrix(2, 2);
    rGeometricStiffnessMatrix(5, 2) = rGeometricStiffnessMatrix(2, 5);

    rGeometricStiffnessMatrix(3, 4) = rGeometricStiffnessMatrix(0, 1);
    rGeometricStiffnessMatrix(4, 3) = rGeometricStiffnessMatrix(3, 4);

    rGeometricStiffnessMatrix(3, 5) = rGeometricStiffnessMatrix(0, 2);
    rGeometricStiffnessMatrix(5, 3) = rGeometricStiffnessMatrix(3, 5);

    rGeometricStiffnessMatrix(4, 5) = rGeometricStiffnessMatrix(1, 2);
    rGeometricStiffnessMatrix(5, 4) = rGeometricStiffnessMatrix(4, 5);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<2, 2>::CalculateGeometricStiffnessMatrix(FullDofMatrixType& rGeometricStiffnessMatrix,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    double       E = ReturnTangentModulus1D(rCurrentProcessInfo);
    const double A = GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    rGeometricStiffnessMatrix = ZeroMatrix(DIM * NUM_NODES, DIM * NUM_NODES);

    // du... delta displacement in x-direction
    // dv... delta displacement in y-direction
    // dw... delta displacement in z-direction
    // L... inital member length
    // l... deformed member length
    // e_gl... green_lagrange strain

    double du = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
                GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    double dv = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
                GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double dx   = GetGeometry()[1].X0() - GetGeometry()[0].X0();
    const double dy   = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
    const double L    = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l    = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    double       e_gL = (l * l - L * L) / (2.00 * L * L);
    const double L3   = L * L * L;

    const double K_sigma = ((E * A * e_gL) / L) + ((prestress * A) / L);
    const double K_uij   = (E * A) / L3;

    rGeometricStiffnessMatrix(0, 0)     = K_sigma + K_uij * (2 * du * dx + du * du);
    rGeometricStiffnessMatrix(DIM, DIM) = rGeometricStiffnessMatrix(0, 0);

    rGeometricStiffnessMatrix(1, 1)             = K_sigma + K_uij * (2 * dv * dy + dv * dv);
    rGeometricStiffnessMatrix(DIM + 1, DIM + 1) = rGeometricStiffnessMatrix(1, 1);

    rGeometricStiffnessMatrix(0, 1) = K_uij * (dx * dv + dy * du + du * dv);
    rGeometricStiffnessMatrix(1, 0) = rGeometricStiffnessMatrix(0, 1);

    rGeometricStiffnessMatrix(0, DIM) = -rGeometricStiffnessMatrix(0, 0);
    rGeometricStiffnessMatrix(DIM, 0) = rGeometricStiffnessMatrix(0, DIM);

    rGeometricStiffnessMatrix(0, DIM + 1) = -rGeometricStiffnessMatrix(0, 1);
    rGeometricStiffnessMatrix(DIM + 1, 0) = rGeometricStiffnessMatrix(0, DIM + 1);

    rGeometricStiffnessMatrix(1, DIM) = rGeometricStiffnessMatrix(0, DIM + 1);
    rGeometricStiffnessMatrix(DIM, 1) = rGeometricStiffnessMatrix(1, DIM);

    rGeometricStiffnessMatrix(1, DIM + 1) = -rGeometricStiffnessMatrix(1, 1);
    rGeometricStiffnessMatrix(DIM + 1, 1) = rGeometricStiffnessMatrix(1, DIM + 1);

    rGeometricStiffnessMatrix(DIM, DIM + 1) = rGeometricStiffnessMatrix(0, 1);
    rGeometricStiffnessMatrix(DIM + 1, DIM) = rGeometricStiffnessMatrix(DIM, DIM + 1);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<3, 2>::CalculateElasticStiffnessMatrix(MatrixType& rElasticStiffnessMatrix,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    double E = ReturnTangentModulus1D(rCurrentProcessInfo);

    double A = GetProperties()[CROSS_AREA];

    rElasticStiffnessMatrix = ZeroMatrix(DIM * NUM_NODES, DIM * NUM_NODES);

    const double dx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
    const double dy = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
    const double dz = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
    const double L  = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double L3 = L * L * L;

    const double EA = E * A;

    rElasticStiffnessMatrix(0, 0) = (EA * dx * dx) / L3;
    rElasticStiffnessMatrix(3, 3) = rElasticStiffnessMatrix(0, 0);

    rElasticStiffnessMatrix(1, 1) = (EA * dy * dy) / L3;
    rElasticStiffnessMatrix(4, 4) = rElasticStiffnessMatrix(1, 1);

    rElasticStiffnessMatrix(2, 2) = (EA * dz * dz) / L3;
    rElasticStiffnessMatrix(5, 5) = rElasticStiffnessMatrix(2, 2);

    rElasticStiffnessMatrix(0, 1) = (EA * dx * dy) / L3;
    rElasticStiffnessMatrix(1, 0) = rElasticStiffnessMatrix(0, 1);

    rElasticStiffnessMatrix(0, 2) = (EA * dx * dz) / L3;
    rElasticStiffnessMatrix(2, 0) = rElasticStiffnessMatrix(0, 2);

    rElasticStiffnessMatrix(0, 3) = -rElasticStiffnessMatrix(0, 0);
    rElasticStiffnessMatrix(3, 0) = rElasticStiffnessMatrix(0, 3);

    rElasticStiffnessMatrix(0, 4) = -rElasticStiffnessMatrix(0, 1);
    rElasticStiffnessMatrix(4, 0) = rElasticStiffnessMatrix(0, 4);

    rElasticStiffnessMatrix(0, 5) = -rElasticStiffnessMatrix(0, 2);
    rElasticStiffnessMatrix(5, 0) = rElasticStiffnessMatrix(0, 5);

    rElasticStiffnessMatrix(1, 2) = (EA * dy * dz) / L3;
    rElasticStiffnessMatrix(2, 1) = rElasticStiffnessMatrix(1, 2);

    rElasticStiffnessMatrix(1, 3) = rElasticStiffnessMatrix(0, 4);
    rElasticStiffnessMatrix(3, 1) = rElasticStiffnessMatrix(1, 3);

    rElasticStiffnessMatrix(1, 4) = -rElasticStiffnessMatrix(1, 1);
    rElasticStiffnessMatrix(4, 1) = rElasticStiffnessMatrix(1, 4);

    rElasticStiffnessMatrix(1, 5) = -rElasticStiffnessMatrix(1, 2);
    rElasticStiffnessMatrix(5, 1) = rElasticStiffnessMatrix(1, 5);

    rElasticStiffnessMatrix(2, 3) = -rElasticStiffnessMatrix(0, 2);
    rElasticStiffnessMatrix(3, 2) = rElasticStiffnessMatrix(2, 3);

    rElasticStiffnessMatrix(2, 4) = -rElasticStiffnessMatrix(1, 2);
    rElasticStiffnessMatrix(4, 2) = rElasticStiffnessMatrix(2, 4);

    rElasticStiffnessMatrix(2, 5) = -rElasticStiffnessMatrix(2, 2);
    rElasticStiffnessMatrix(5, 2) = rElasticStiffnessMatrix(2, 5);

    rElasticStiffnessMatrix(3, 4) = rElasticStiffnessMatrix(0, 1);
    rElasticStiffnessMatrix(4, 3) = rElasticStiffnessMatrix(3, 4);

    rElasticStiffnessMatrix(3, 5) = rElasticStiffnessMatrix(0, 2);
    rElasticStiffnessMatrix(5, 3) = rElasticStiffnessMatrix(3, 5);

    rElasticStiffnessMatrix(4, 5) = rElasticStiffnessMatrix(1, 2);
    rElasticStiffnessMatrix(5, 4) = rElasticStiffnessMatrix(4, 5);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementBase<2, 2>::CalculateElasticStiffnessMatrix(MatrixType& rElasticStiffnessMatrix,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    double E = ReturnTangentModulus1D(rCurrentProcessInfo);
    double A = GetProperties()[CROSS_AREA];

    rElasticStiffnessMatrix = ZeroMatrix(DIM * NUM_NODES, DIM * NUM_NODES);

    const double dx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
    const double dy = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
    const double L  = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double L3 = L * L * L;

    const double EA = E * A;

    rElasticStiffnessMatrix(0, 0)     = (EA * dx * dx) / L3;
    rElasticStiffnessMatrix(DIM, DIM) = rElasticStiffnessMatrix(0, 0);

    rElasticStiffnessMatrix(1, 1)             = (EA * dy * dy) / L3;
    rElasticStiffnessMatrix(DIM + 1, DIM + 1) = rElasticStiffnessMatrix(1, 1);

    rElasticStiffnessMatrix(0, 1) = (EA * dx * dy) / L3;
    rElasticStiffnessMatrix(1, 0) = rElasticStiffnessMatrix(0, 1);

    rElasticStiffnessMatrix(0, DIM) = -rElasticStiffnessMatrix(0, 0);
    rElasticStiffnessMatrix(DIM, 0) = rElasticStiffnessMatrix(0, DIM);

    rElasticStiffnessMatrix(0, DIM + 1) = -rElasticStiffnessMatrix(0, 1);
    rElasticStiffnessMatrix(DIM + 1, 0) = rElasticStiffnessMatrix(0, DIM + 1);

    rElasticStiffnessMatrix(1, DIM) = rElasticStiffnessMatrix(0, DIM + 1);
    rElasticStiffnessMatrix(DIM, 1) = rElasticStiffnessMatrix(1, DIM);

    rElasticStiffnessMatrix(1, DIM + 1) = -rElasticStiffnessMatrix(1, 1);
    rElasticStiffnessMatrix(DIM + 1, 1) = rElasticStiffnessMatrix(1, DIM + 1);

    rElasticStiffnessMatrix(DIM, DIM + 1) = rElasticStiffnessMatrix(0, 1);
    rElasticStiffnessMatrix(DIM + 1, DIM) = rElasticStiffnessMatrix(DIM, DIM + 1);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    Vector                      temp_strain = ZeroVector(1);
    Vector                      temp_stress = ZeroVector(1);
    temp_strain[0]                          = CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("mpConstitutiveLaw", mpConstitutiveLaw);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mpConstitutiveLaw", mpConstitutiveLaw);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementBase<TDim, TNumNodes>::CalculateLumpedMassVector(VectorType& rLumpedMassVector,
                                                                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Clear matrix
    if (rLumpedMassVector.size() != TDim * TNumNodes) {
        rLumpedMassVector.resize(TDim * TNumNodes, false);
    }

    const double A   = GetProperties()[CROSS_AREA];
    const double L   = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = GetProperties()[DENSITY];

    const double total_mass = A * L * rho;

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        for (unsigned int j = 0; j < TDim; ++j) {
            unsigned int index = i * TDim + j;

            rLumpedMassVector[index] = total_mass * 0.50;
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoTrussElementBase<TDim, TNumNodes>::ReturnTangentModulus1D(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    double tangent_modulus(0.00);

    Vector strain_vector = ZeroVector(mpConstitutiveLaw->GetStrainSize());
    strain_vector[0]     = CalculateGreenLagrangeStrain();

    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    Values.SetStrainVector(strain_vector);

    mpConstitutiveLaw->CalculateValue(Values, TANGENT_MODULUS, tangent_modulus);

    return tangent_modulus;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType GeoTrussElementBase<TDim, TNumNodes>::GetDofs() const
{
    auto result = Element::DofsVectorType{};
    for (const auto& r_node : GetGeometry()) {
        result.push_back(r_node.pGetDof(DISPLACEMENT_X));
        result.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        if constexpr (TDim == 3) result.push_back(r_node.pGetDof(DISPLACEMENT_Z));
    }
    return result;
}

//--------------------------------------------------------------------------------------------
template class GeoTrussElementBase<2, 2>;
template class GeoTrussElementBase<3, 2>;

} // namespace Kratos.
