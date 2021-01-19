// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/beam_element_2D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
BeamElement2D2N::BeamElement2D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

BeamElement2D2N::BeamElement2D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
BeamElement2D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<BeamElement2D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
BeamElement2D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<BeamElement2D2N>(NewId, pGeom,
            pProperties);
}

BeamElement2D2N::~BeamElement2D2N() {}

void BeamElement2D2N::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
if (rResult.size() != msElementSize) {
        rResult.resize(msElementSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();

        rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }
}

void BeamElement2D2N::GetDofList(DofsVectorType& rElementalDofList,
                                   ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != msElementSize) {
        rElementalDofList.resize(msElementSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Y);

        rElementalDofList[index + 2] = GetGeometry()[i].pGetDof(ROTATION_Z);
    }
}

void BeamElement2D2N::Initialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("")
}

void BeamElement2D2N::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;

        rValues[index] =
            GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[index + 1] =
            GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[index + 2] = GetGeometry()[i].FastGetSolutionStepValue(
                                 ANGULAR_ACCELERATION_Z, Step);
    }
    KRATOS_CATCH("")
}


void BeamElement2D2N::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        rValues[index] =
            GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, Step);
        rValues[index + 1] =
            GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
        rValues[index + 2] = GetGeometry()[i].FastGetSolutionStepValue(
                                 ANGULAR_VELOCITY_Z, Step);
    }

    KRATOS_CATCH("")
}

void BeamElement2D2N::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        rValues[index] =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
        rValues[index + 1] =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
        rValues[index + 2] =
            GetGeometry()[i].FastGetSolutionStepValue(ROTATION_Z, Step);
    }
    KRATOS_CATCH("")
}


Matrix BeamElement2D2N::CalculateInitialLocalCS() const
{
    KRATOS_TRY
    Matrix temp_matrix = ZeroMatrix(msLocalSize);
    return temp_matrix;
    KRATOS_CATCH("")
}

void BeamElement2D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    ConstCalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BeamElement2D2N::ConstCalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    ConstCalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    ConstCalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    KRATOS_CATCH("");
}

void BeamElement2D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstCalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BeamElement2D2N::ConstCalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    rRightHandSideVector = ZeroVector(msElementSize);

    Vector f_int_global = ZeroVector(msElementSize);
    GetGlobalInternalForces(f_int_global);


    noalias(rRightHandSideVector) -= f_int_global;
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}


void BeamElement2D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    ConstCalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BeamElement2D2N::ConstCalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
    Matrix K_el = ZeroMatrix(msElementSize);
    MaterialTangentMatrix(K_el);
    Matrix K_geo = ZeroMatrix(msElementSize);
    GeometricTangentMatrix(K_geo);

    rLeftHandSideMatrix = K_el + K_geo;
    KRATOS_CATCH("")
}


void BeamElement2D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

void BeamElement2D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

int BeamElement2D2N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 2 || GetGeometry().size() != 2)
            << "The beam element works only in 2D and with 2 noded elements" << std::endl;

    // verify that the variables are correctly initialized
    if (VELOCITY.Key() == 0) {
        KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (DISPLACEMENT.Key() == 0) {
        KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (ACCELERATION.Key() == 0) {
        KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (DENSITY.Key() == 0) {
        KRATOS_ERROR << "DENSITY has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (CROSS_AREA.Key() == 0) {
        KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }

    // verify that the dofs exist
    for (unsigned int i = 0; i < GetGeometry().size(); ++i) {
        if (GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
            KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                         << GetGeometry()[i].Id() << std::endl;
        }
        if (GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false) {
            KRATOS_ERROR
                    << "missing one of the dofs for the variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
        }
    }

    KRATOS_ERROR_IF(!GetProperties().Has(CROSS_AREA) ||
                    GetProperties()[CROSS_AREA] <= numerical_limit)
            << "Please provide a reasonable value for \"CROSS_AREA\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(YOUNG_MODULUS) ||
                    GetProperties()[YOUNG_MODULUS] <= numerical_limit)
            << "Please provide a reasonable value for \"YOUNG_MODULUS\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(DENSITY) ||
                    GetProperties()[DENSITY] <= numerical_limit)
            << "Please provide a reasonable value for \"DENSITY\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(I33) ||
                    GetProperties()[I33] <= numerical_limit)
            << "Please provide a reasonable value for \"I33\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(POISSON_RATIO))
            << "\"POISSON_RATIO\" not provided for element #" << Id() << std::endl;

    KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this)
                    < std::numeric_limits<double>::epsilon())
            << "Element #" << Id() << " has a length of zero!" << std::endl;

    return 0;

    KRATOS_CATCH("")
}

double BeamElement2D2N::RigidBodyOrientation(const ConfigurationType& rConfiguration) const {
    KRATOS_TRY;
    double dx = GetGeometry()[1].X0()-GetGeometry()[0].X0();
    double dy = GetGeometry()[1].Y0()-GetGeometry()[0].Y0();
    if (rConfiguration == ConfigurationType::Current){
        dx += GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)-GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        dy += GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)-GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    }
    return std::atan2(dy,dx);
    KRATOS_CATCH("");
}

void BeamElement2D2N::CalculateBMatrix(Matrix& rBMatrix) const{
    KRATOS_TRY;
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double b = RigidBodyOrientation(ConfigurationType::Current);
    const double c = std::cos(b);
    const double s = std::sin(b);

    rBMatrix = ZeroMatrix(msLocalSize,msElementSize);
    rBMatrix(0,0) = -c;
    rBMatrix(0,1) = -s;
    rBMatrix(0,3) =  c;
    rBMatrix(0,4) =  s;

    rBMatrix(1,0) = -s/l;
    rBMatrix(1,1) =  c/l;
    rBMatrix(1,2) =  1.0;
    rBMatrix(1,3) =  s/l;
    rBMatrix(1,4) =  -c/l;

    rBMatrix(2,0) = -s/l;
    rBMatrix(2,1) =  c/l;
    rBMatrix(2,3) =  s/l;
    rBMatrix(2,4) =  -c/l;
    rBMatrix(2,5) =  1.0;
    KRATOS_CATCH("")
}

void BeamElement2D2N::LocalMaterialMatrix(Matrix& rMaterialMatrix) const{
    KRATOS_TRY;

    const double Iz = GetProperties()[I33];
    const double E = GetProperties()[YOUNG_MODULUS];
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double r_g2 = Iz/A;

    rMaterialMatrix = ZeroMatrix(msLocalSize);

    rMaterialMatrix(0,0) = 1.0;

    rMaterialMatrix(1,1) = 4.0*r_g2;
    rMaterialMatrix(1,2) = 2.0*r_g2;

    rMaterialMatrix(2,1) = 2.0*r_g2;
    rMaterialMatrix(2,2) = 4.0*r_g2;

    rMaterialMatrix *= E*A/L;

    KRATOS_CATCH("");
}

void BeamElement2D2N::MaterialTangentMatrix(Matrix& rMaterialMatrix) const{
    KRATOS_TRY;

    Matrix B_matrix = ZeroMatrix(msLocalSize,msElementSize);
    CalculateBMatrix(B_matrix);

    Matrix K_mat_local = ZeroMatrix(msLocalSize);
    LocalMaterialMatrix(K_mat_local);

    rMaterialMatrix = ZeroMatrix(msElementSize);
    rMaterialMatrix = prod(K_mat_local,B_matrix);
    rMaterialMatrix = prod(trans(B_matrix),rMaterialMatrix);

    KRATOS_CATCH("");
}

void BeamElement2D2N::GeometricTangentMatrix(Matrix& rGeometricMatrix) const{
    KRATOS_TRY;

    Vector f_int_local = ZeroVector(msLocalSize);
    GetLocalInternalForces(f_int_local);
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double b = RigidBodyOrientation(ConfigurationType::Current);
    const double c = std::cos(b);
    const double s = std::sin(b);

    Matrix r = ZeroMatrix(msElementSize,1);
    r(0,0) = -c;
    r(1,0) = -s;
    r(3,0) =  c;
    r(4,0) =  s;

    Matrix z = ZeroMatrix(msElementSize,1);
    r(0,0) =  s;
    r(1,0) = -c;
    r(3,0) = -s;
    r(4,0) =  c;

    rGeometricMatrix = ZeroMatrix(msElementSize);
    rGeometricMatrix += prod(z,trans(z)) * f_int_local[0] / l;
    rGeometricMatrix += (f_int_local[1] + f_int_local[2]) * (prod(r,trans(z)) + prod(z,trans(r)));

    KRATOS_CATCH("");
}


void BeamElement2D2N::GetLocalDeformations(Vector& rLocalDeformation) const{
    KRATOS_TRY;
    const double b_0 = RigidBodyOrientation(ConfigurationType::Reference);
    const double b_1 = b_0 + GetGeometry()[0].FastGetSolutionStepValue(ROTATION_Z);
    const double b_2 = b_0 + GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z);
    const double b = RigidBodyOrientation(ConfigurationType::Current);

    const double phi_1_local = std::atan2((std::cos(b)*std::sin(b_1)) - (std::sin(b)*std::cos(b_1)),(std::cos(b)*std::cos(b_1)) + (std::sin(b)*std::sin(b_1)));
    const double phi_2_local = std::atan2((std::cos(b)*std::sin(b_2)) - (std::sin(b)*std::cos(b_2)),(std::cos(b)*std::cos(b_2)) + (std::sin(b)*std::sin(b_2)));

    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);

    const double u_local = (std::pow(l,2.0)-std::pow(L,2.0)) / (l+L);  //Crisfield

    rLocalDeformation = ZeroVector(msLocalSize);
    rLocalDeformation[0] = u_local;
    rLocalDeformation[1] = phi_1_local;
    rLocalDeformation[2] = phi_2_local;

    KRATOS_CATCH("");
}

void BeamElement2D2N::GetLocalInternalForces(Vector& rLocalForces) const{
    KRATOS_TRY;

    Vector d_local = ZeroVector(msLocalSize);
    GetLocalDeformations(d_local);

    Matrix K_mat_local = ZeroMatrix(msLocalSize);
    LocalMaterialMatrix(K_mat_local);

    rLocalForces = ZeroVector(msLocalSize);
    rLocalForces = prod(K_mat_local,d_local);

    KRATOS_CATCH("");
}

void BeamElement2D2N::GetGlobalInternalForces(Vector& rGlobalForces) const{
    KRATOS_TRY;

    Vector f_int_local = ZeroVector(msLocalSize);
    GetLocalInternalForces(f_int_local);

    Matrix B_matrix = ZeroMatrix(msLocalSize,msElementSize);
    CalculateBMatrix(B_matrix);

    rGlobalForces = ZeroVector(msElementSize);
    rGlobalForces = prod(trans(B_matrix),f_int_local);

    KRATOS_CATCH("");
}



void BeamElement2D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double A = GetProperties()[CROSS_AREA];
    const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    bool use_consistent_mass_matrix = true;
    if (GetProperties().Has(USE_CONSISTENT_MASS_MATRIX)) {
        use_consistent_mass_matrix = GetProperties()[USE_CONSISTENT_MASS_MATRIX];
    }

    if (use_consistent_mass_matrix) {
        const double pre_beam = (rho * A * L) / 420.00;
        const double pre_bar = (rho * A * L) / 6.00;

        // bar part
        rMassMatrix(0, 0) = 2.00 * pre_bar;
        rMassMatrix(0, 3) = 1.00 * pre_bar;
        rMassMatrix(3, 0) = 1.00 * pre_bar;
        rMassMatrix(3, 3) = 2.00 * pre_bar;

        // beam part

        rMassMatrix(1, 1) = pre_beam * 156.00;
        rMassMatrix(1, 2) = pre_beam * 22.00 * L;
        rMassMatrix(1, 4) = pre_beam * 54.00;
        rMassMatrix(1, 5) = pre_beam * (-13.00) * L;

        rMassMatrix(2, 1) = pre_beam * 22.00 * L;
        rMassMatrix(2, 2) = pre_beam * 4.00 * L * L;
        rMassMatrix(2, 4) = pre_beam * 13.00 * L;
        rMassMatrix(2, 5) = pre_beam * (-3.00) * L * L;

        rMassMatrix(4, 1) = pre_beam * 54.00;
        rMassMatrix(4, 2) = pre_beam * 13.00 * L;
        rMassMatrix(4, 4) = pre_beam * 156.00;
        rMassMatrix(4, 5) = pre_beam * (-22.00) * L;

        rMassMatrix(5, 1) = pre_beam * (-13.00) * L;
        rMassMatrix(5, 2) = pre_beam * (-3.00) * L * L;
        rMassMatrix(5, 4) = pre_beam * (-22.00) * L;
        rMassMatrix(5, 5) = pre_beam * (4.00) * L * L;

        // Globalize Mass Matrix

        const double b_0 = RigidBodyOrientation(ConfigurationType::Reference);
        const double b_1 = b_0 + GetGeometry()[0].FastGetSolutionStepValue(ROTATION_Z);
        const double b_2 = b_0 + GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z);
        const double c_1 = std::cos(b_1);
        const double s_1 = std::sin(b_1);

        const double c_2 = std::cos(b_2);
        const double s_2 = std::sin(b_2);

        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix =
        ZeroMatrix(msElementSize, msElementSize);

        rotation_matrix(0, 0) = c_1;
        rotation_matrix(0, 1) = -s_1;
        rotation_matrix(1, 0) = s_1;
        rotation_matrix(1, 1) = c_1;
        rotation_matrix(2, 2) = 1.00;

        rotation_matrix(3, 3) = c_2;
        rotation_matrix(3, 4) = -s_2;
        rotation_matrix(4, 3) = s_2;
        rotation_matrix(4, 4) = c_2;
        rotation_matrix(5, 5) = 1.00;

        rMassMatrix = prod(rMassMatrix,trans(rotation_matrix));
        rMassMatrix = prod(rotation_matrix,rMassMatrix);

    }

    else {
        const double lumped_mass = A * L * rho;

        // w.r.t. Felippa - Chapter 31: LUMPED AND CONSISTENT MASS MATRICES - p.31–10
        double alpha = 0.00;
        if (GetProperties().Has(LUMPED_MASS_ROTATION_COEFFICIENT)) {
            alpha = GetProperties()[LUMPED_MASS_ROTATION_COEFFICIENT];
        }
        rMassMatrix(0, 0) = lumped_mass * 0.50;
        rMassMatrix(1, 1) = lumped_mass * 0.50;
        rMassMatrix(2, 2) = lumped_mass * L * L * alpha;
        rMassMatrix(3, 3) = lumped_mass * 0.50;
        rMassMatrix(4, 4) = lumped_mass * 0.50;
        rMassMatrix(5, 5) = lumped_mass * L * L * alpha;
    }



    KRATOS_CATCH("")
}


void BeamElement2D2N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        msElementSize);
}


void BeamElement2D2N::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    // FORCE- & Moment- Residual is 3D vector
    KRATOS_TRY;

    BoundedVector<double, msElementSize> damping_residual_contribution =
        ZeroVector(msElementSize);
    // calculate damping contribution to residual -->
    if ((GetProperties().Has(RAYLEIGH_ALPHA) ||
            GetProperties().Has(RAYLEIGH_BETA)) &&
            (rDestinationVariable != NODAL_INERTIA)) {
        Vector current_nodal_velocities = ZeroVector(msElementSize);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) =
            prod(damping_matrix, current_nodal_velocities);
    }

    if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == FORCE_RESIDUAL) {

        for (SizeType i = 0; i < msNumberOfNodes; ++i) {
            SizeType index = msLocalSize * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3>& r_force_residual =
                GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (SizeType j = 0; j < msDimension; ++j) {
                r_force_residual[j] +=
                    rRHSVector[index + j] - damping_residual_contribution[index + j];
            }

            r_force_residual[msDimension] = 0.00;
            GetGeometry()[i].UnSetLock();
        }
    }

    if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == MOMENT_RESIDUAL) {

        for (SizeType i = 0; i < msNumberOfNodes; ++i) {
            SizeType index = (msLocalSize * i) + msDimension;

            GetGeometry()[i].SetLock();

            array_1d<double, 3>& r_moment_residual =
                GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

            for (SizeType j = 0; j < msDimension; ++j) {
                r_moment_residual[j] = 0.00;
            }
            r_moment_residual[msDimension] +=
                rRHSVector[index] - damping_residual_contribution[index];
            GetGeometry()[i].UnSetLock();
        }
    }

    if (rDestinationVariable == NODAL_INERTIA) {
        Matrix element_mass_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_info; // Dummy
        CalculateMassMatrix(element_mass_matrix, temp_info);

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            double aux_nodal_mass = 0.0;
            double aux_nodal_inertia = 0.0;

            const SizeType index = i * msLocalSize;

            for (IndexType j = 0; j < msElementSize; ++j) {
                aux_nodal_mass += element_mass_matrix(index, j);
                aux_nodal_inertia += element_mass_matrix(index + msDimension, j);
            }

            #pragma omp atomic
            GetGeometry()[i].GetValue(NODAL_MASS) += aux_nodal_mass;

            array_1d<double, 3>& r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
            #pragma omp atomic
            r_nodal_inertia[msDimension] += std::abs(aux_nodal_inertia);
        }
    }

    KRATOS_CATCH("")
}

BoundedVector<double, BeamElement2D2N::msElementSize>
    BeamElement2D2N::CalculateBodyForces() const
{
    KRATOS_TRY
    // getting shapefunctionvalues for linear SF
    const Matrix& Ncontainer =
        GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    BoundedVector<double, 3> equivalent_line_load = ZeroVector(3);
    BoundedVector<double, msElementSize> body_forces_global =
        ZeroVector(msElementSize);

    const double A = GetProperties()[CROSS_AREA];
    const double l = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double rho = GetProperties()[DENSITY];

    // calculating equivalent line load
    for (int i = 0; i < msNumberOfNodes; ++i) {
        equivalent_line_load +=
            A * rho *
            GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
            Ncontainer(0, i);
    }

    // adding the nodal forces
    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        for (int j = 0; j < msDimension; ++j) {
            body_forces_global[j + index] =
                equivalent_line_load[j] * Ncontainer(0, i) * l;
        }
    }

    // adding the nodal moments
    CalculateAndAddWorkEquivalentNodalForcesLineLoad(equivalent_line_load,
            body_forces_global, l);
    // return the total ForceVector
    return body_forces_global;
    KRATOS_CATCH("")
}

void BeamElement2D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
    const BoundedVector<double, 3> ForceInput,
    BoundedVector<double, BeamElement2D2N::msElementSize>
    & rRightHandSideVector,
    const double GeometryLength) const
{
    KRATOS_TRY;
    // calculate orthogonal load vector
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    Vector geometric_orientation = ZeroVector(3);
    geometric_orientation[0] =
        GetGeometry()[1].X() - GetGeometry()[0].X();
    geometric_orientation[1] =
        GetGeometry()[1].Y() - GetGeometry()[0].Y();
    geometric_orientation[2] = 0.000;

    const double vector_norm_a = MathUtils<double>::Norm(geometric_orientation);
    if (vector_norm_a > numerical_limit) {
        geometric_orientation /= vector_norm_a;
    }

    Vector line_load_direction = ZeroVector(3);
    for (int i = 0; i < msDimension; ++i) {
        line_load_direction[i] = ForceInput[i];
    }

    const double vector_norm_b = MathUtils<double>::Norm(line_load_direction);
    if (vector_norm_b > numerical_limit) {
        line_load_direction /= vector_norm_b;
    }

    double cos_angle = 0.00;
    for (int i = 0; i < msDimension; ++i) {
        cos_angle += line_load_direction[i] * geometric_orientation[i];
    }

    const double sin_angle = std::sqrt(1.00 - (cos_angle * cos_angle));
    const double norm_force_vector_orthogonal = sin_angle * vector_norm_b;

    Vector node_a = ZeroVector(3);
    node_a[0] = GetGeometry()[0].X();
    node_a[1] = GetGeometry()[0].Y();
    node_a[2] = 0.00;

    Vector node_b = ZeroVector(3);
    node_b = node_a + line_load_direction;

    Vector node_c = ZeroVector(3);
    node_c = node_a + (geometric_orientation * cos_angle);

    Vector load_orthogonal_direction = ZeroVector(3);
    load_orthogonal_direction = node_b - node_c;
    const double vector_norm_c =
        MathUtils<double>::Norm(load_orthogonal_direction);
    if (vector_norm_c > numerical_limit) {
        load_orthogonal_direction /= vector_norm_c;
    }

    // now caluclate respective work equivilent nodal moments

    const double custom_moment =
        norm_force_vector_orthogonal * GeometryLength * GeometryLength / 12.00;

    Vector moment_node_a = ZeroVector(3);
    moment_node_a = MathUtils<double>::CrossProduct(geometric_orientation,
                    load_orthogonal_direction);
    moment_node_a *= custom_moment;

    rRightHandSideVector[msDimension] += moment_node_a[2];
    rRightHandSideVector[(2 * msDimension) + 1] -= moment_node_a[2];

    KRATOS_CATCH("")
}


} // namespace Kratos.