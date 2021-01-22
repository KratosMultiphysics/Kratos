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
#include "custom_elements/beam_element_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
BeamElement3D2N::BeamElement3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

BeamElement3D2N::BeamElement3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
BeamElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<BeamElement3D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
BeamElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<BeamElement3D2N>(NewId, pGeom,
            pProperties);
}

BeamElement3D2N::~BeamElement3D2N() {}

void BeamElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != msElementSize) {
        rResult.resize(msElementSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msNumberOfNodes * msDimension;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }
}

void BeamElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
                                   ProcessInfo& rCurrentProcessInfo)
{

    if (rElementalDofList.size() != msElementSize) {
        rElementalDofList.resize(msElementSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msNumberOfNodes * msDimension;
        rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Z);

        rElementalDofList[index + 3] = GetGeometry()[i].pGetDof(ROTATION_X);
        rElementalDofList[index + 4] = GetGeometry()[i].pGetDof(ROTATION_Y);
        rElementalDofList[index + 5] = GetGeometry()[i].pGetDof(ROTATION_Z);
    }
}

void BeamElement3D2N::Initialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("")
}

void BeamElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step) const
{

    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& acc =
            GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const auto& ang_acc = GetGeometry()[i].FastGetSolutionStepValue(
                                  ANGULAR_ACCELERATION, Step);

        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];

        rValues[index + 3] = ang_acc[0];
        rValues[index + 4] = ang_acc[1];
        rValues[index + 5] = ang_acc[2];
    }
    KRATOS_CATCH("")
}

void BeamElement3D2N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mDeformationPreviousIteration = mDeformationCurrentIteration;
    GetValuesVector(mDeformationCurrentIteration, 0);

    //maybe we dont need the members "mDeformationPreviousIteration" anymore ....
    UpdateGlobalNodalRotations();
    KRATOS_CATCH("")
}

void BeamElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step) const
{

    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& vel =
            GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const auto& ang_vel =
            GetGeometry()[i].FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);

        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];

        rValues[index + 3] = ang_vel[0];
        rValues[index + 4] = ang_vel[1];
        rValues[index + 5] = ang_vel[2];
    }
    KRATOS_CATCH("")
}

void BeamElement3D2N::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& disp =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const auto& rot =
            GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step);

        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
    KRATOS_CATCH("")
}

Matrix BeamElement3D2N::CreateElementStiffnessMatrix_Material() const
{

    KRATOS_TRY;
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = CalculateShearModulus();
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

    const double It = GetProperties()[TORSIONAL_INERTIA];
    const double Iy = GetProperties()[I22];
    const double Iz = GetProperties()[I33];

    const int local_deformation_possibilities(7);


    Matrix local_stiffness_matrix =  ZeroMatrix(local_deformation_possibilities);


    local_stiffness_matrix(0, 0) = E * A / L;

    local_stiffness_matrix(1, 1) = G*It/L;
    local_stiffness_matrix(4, 4) = G*It/L;
    local_stiffness_matrix(1, 4) = -G*It/L;
    local_stiffness_matrix(4, 1) = -G*It/L;

    local_stiffness_matrix(2, 2) = 4.0*E*Iz/L;
    local_stiffness_matrix(5, 5) = 4.0*E*Iz/L;
    local_stiffness_matrix(2, 5) = 2.0*E*Iz/L;
    local_stiffness_matrix(5, 2) = 2.0*E*Iz/L;

    local_stiffness_matrix(3, 3) = 4.0*E*Iy/L;
    local_stiffness_matrix(6, 6) = 4.0*E*Iy/L;
    local_stiffness_matrix(3, 6) = 2.0*E*Iy/L;
    local_stiffness_matrix(6, 3) = 2.0*E*Iy/L;

    return local_stiffness_matrix;
    KRATOS_CATCH("")
}

Matrix BeamElement3D2N::CreateElementStiffnessMatrixIntermediate() const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    Matrix b_a = BMatrixIntermediateA();
    Matrix local_stiffness_matrix = prod(CreateElementStiffnessMatrix_Material(),b_a);
    local_stiffness_matrix = prod(trans(b_a),local_stiffness_matrix);

    Vector local_deformation = LocalDeformations();
    Vector local_internal_forces = LocalInternalForces();

    for (int node_nr=0;node_nr<2;++node_nr)
    {
        Vector m_i = ZeroVector(msDimension);
        Vector phi_i = ZeroVector(msDimension);

        if (node_nr==0)
        {
            for (int i=0;i<msDimension;++i)
            {
                m_i[i] = local_internal_forces[1+i];
                phi_i[i] = local_deformation[1+i];
            }
        }
        else if (node_nr==1)
        {
            for (int i=0;i<msDimension;++i)
            {
                m_i[i] = local_internal_forces[4+i];
                phi_i[i] = local_deformation[4+i];
            }
        }
        else KRATOS_ERROR << "wrong node number in CreateElementStiffnessMatrixIntermediate" << std::endl;


        const double alpha = MathUtils<double>::Norm(phi_i);
        if (alpha>numerical_limit)
        {
            const double eta = ( (2.0*std::sin(alpha)) - (alpha*(1.0+std::cos(alpha))) ) / (2.0*std::pow(alpha,2.0)*std::sin(alpha));
            const double nu  = ( (alpha*(alpha+std::sin(alpha))) - (8.0*std::pow(std::sin(alpha/2.0),2.0)) ) / ( 4.0*std::pow(alpha,4.0)*std::pow(std::sin(alpha/2.0),2.0));

            Matrix phi_skew = SkewSymmetricMatrix(phi_i);
            Matrix m_skew = SkewSymmetricMatrix(m_i);

            Matrix k_h_i = outer_prod(phi_i,m_i);
            k_h_i -= 2.0*outer_prod(m_i,phi_i);
            k_h_i += inner_prod(phi_i,m_i) * IdentityMatrix(msDimension);
            k_h_i *= eta;

            Matrix temp_mat_1 = prod(phi_skew,phi_skew);
            k_h_i += nu * prod(temp_mat_1,outer_prod(m_i,phi_i));
            k_h_i -= 0.50 * m_skew;

            k_h_i = prod(k_h_i,InverseLocalRotation(node_nr+1));


            if (node_nr==0) project(local_stiffness_matrix, range(1,4), range(1,4)) += k_h_i;
            else if (node_nr==1) project(local_stiffness_matrix, range(4,7), range(4,7)) += k_h_i;
            else KRATOS_ERROR << "wrong node number in CreateElementStiffnessMatrixIntermediate" << std::endl;

        }


    }
    return local_stiffness_matrix;
}

Matrix BeamElement3D2N::GlobalTangentStiffnessMatrix() const
{
    Matrix b_g = BMatrixGlobal();
    Matrix tangent_stiffness_matrix = prod(CreateElementStiffnessMatrixIntermediate(),b_g);
    tangent_stiffness_matrix = prod(trans(b_g),tangent_stiffness_matrix);

    Vector intermediate_internal_forces = LocalInternalIntermediateForces();

    ///////
    tangent_stiffness_matrix += intermediate_internal_forces[0]*DMatrix();

    Matrix e = EMatrix();
    Matrix g = GMatrix();

    Matrix temp_1 = prod(trans(g),trans(e));
    temp_1 = prod(QMatrix(),temp_1);
    temp_1 = prod(e,temp_1);


    ///////
    tangent_stiffness_matrix -= temp_1;


    temp_1 = outer_prod(AVector(),RVector());
    temp_1 = prod(g,temp_1);
    temp_1 = prod(e,temp_1);

    ///////
    tangent_stiffness_matrix += temp_1;


    return tangent_stiffness_matrix;
}

Matrix BeamElement3D2N::CalculateInitialLocalCS() const
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    array_1d<double, msDimension> direction_vector_x = ZeroVector(msDimension);
    array_1d<double, msDimension> direction_vector_y = ZeroVector(msDimension);
    array_1d<double, msDimension> direction_vector_z = ZeroVector(msDimension);
    array_1d<double, msLocalSize> reference_coordinates = ZeroVector(msLocalSize);

    reference_coordinates[0] = GetGeometry()[0].X0();
    reference_coordinates[1] = GetGeometry()[0].Y0();
    reference_coordinates[2] = GetGeometry()[0].Z0();
    reference_coordinates[3] = GetGeometry()[1].X0();
    reference_coordinates[4] = GetGeometry()[1].Y0();
    reference_coordinates[5] = GetGeometry()[1].Z0();

    for (unsigned int i = 0; i < msDimension; ++i) {
        direction_vector_x[i] =
            (reference_coordinates[i + msDimension] - reference_coordinates[i]);
    }
    Matrix temp_matrix = ZeroMatrix(msDimension);

    // take user defined local axis 2 from GID input
    if (Has(LOCAL_AXIS_2)) {
        double vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        direction_vector_y = GetValue(LOCAL_AXIS_2);

        direction_vector_z[0] = direction_vector_x[1] * direction_vector_y[2] -
                                direction_vector_x[2] * direction_vector_y[1];
        direction_vector_z[1] = direction_vector_x[2] * direction_vector_y[0] -
                                direction_vector_x[0] * direction_vector_y[2];
        direction_vector_z[2] = direction_vector_x[0] * direction_vector_y[1] -
                                direction_vector_x[1] * direction_vector_y[0];

        vector_norm = MathUtils<double>::Norm(direction_vector_z);
        if (vector_norm > numerical_limit) {
            direction_vector_z /= vector_norm;
        } else
            KRATOS_ERROR << "LOCAL_AXIS_3 has length 0 for element " << Id()
                         << std::endl;

        for (int i = 0; i < msDimension; ++i) {
            temp_matrix(i, 0) = direction_vector_x[i];
            temp_matrix(i, 1) = direction_vector_y[i];
            temp_matrix(i, 2) = direction_vector_z[i];
        }
    }

    // if no user defined local axis 2 input available use this
    else {
        // use orientation class 1st constructor
        double theta_custom = 0.00;
        if (GetProperties().Has(ANG_ROT)) {
            theta_custom = GetProperties()[ANG_ROT];
        }

        typedef array_1d<double, msDimension> arraydim;
        arraydim global_z = ZeroVector(msDimension);
        global_z[2] = 1.0;

        arraydim v2 = ZeroVector(msDimension);
        arraydim v3 = ZeroVector(msDimension);

        double vector_norm;
        vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        if (std::abs(direction_vector_x[2] - 1.00) < numerical_limit) {
            v2[1] = 1.0;
            v3[0] = -1.0;
        }

        else if (std::abs(direction_vector_x[2] + 1.00) < numerical_limit) {
            v2[1] = 1.0;
            v3[0] = 1.0;
        }

        else {
            MathUtils<double>::UnitCrossProduct(v2, global_z, direction_vector_x);
            MathUtils<double>::UnitCrossProduct(v3, direction_vector_x, v2);
        }

        // manual rotation around the beam axis
        if (std::abs(theta_custom) > numerical_limit) {
            const Vector nz_temp = v3;
            const Vector ny_temp = v2;
            const double cos_theta = std::cos(theta_custom);
            const double sin_theta = std::sin(theta_custom);

            v2 = ny_temp * cos_theta + nz_temp * sin_theta;
            vector_norm = MathUtils<double>::Norm(v2);
            if (vector_norm > numerical_limit) {
                v2 /= vector_norm;
            }

            v3 = nz_temp * cos_theta - ny_temp * sin_theta;
            vector_norm = MathUtils<double>::Norm(v3);
            if (vector_norm > numerical_limit) {
                v3 /= vector_norm;
            }
        }

        for (int i = 0; i < msDimension; ++i) {
            temp_matrix(i, 0) = direction_vector_x[i];
            temp_matrix(i, 1) = v2[i];
            temp_matrix(i, 2) = v3[i];
        }
    }


    return temp_matrix;
    KRATOS_CATCH("")
}

Vector BeamElement3D2N::CurrentLocalAxis1() const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    BoundedVector<double, msLocalSize> current_nodal_position = GetCurrentNodalPosition();
    Vector r_1 = ZeroVector(msDimension);

    for (unsigned int i = 0; i < msDimension; ++i) {
        r_1[i] = (current_nodal_position[i + msDimension] - current_nodal_position[i]);
    }

    const double vector_norm = MathUtils<double>::Norm(r_1);
    if (vector_norm > numerical_limit) r_1 /= vector_norm;
    return r_1;
}

Matrix BeamElement3D2N::CoRotatingCS() const
{
    Matrix reference_rotation = CalculateInitialLocalCS();
    Vector r_1 = CurrentLocalAxis1();
    Vector temp_vec = ZeroVector(msDimension);
    temp_vec[1] = 1.0;

    Vector q_1 = prod(reference_rotation,trans(temp_vec));
    q_1 = prod(mGlobalRotationNode1,q_1);

    Vector q_2 = prod(reference_rotation,trans(temp_vec));
    q_2 = prod(mGlobalRotationNode2,q_2);

    Vector q_mean = 0.5 * (q_1+q_2);

    Vector r_3 = ZeroVector(msDimension);
    MathUtils<double>::UnitCrossProduct(r_3, r_1, q_mean);

    Vector r_2 = ZeroVector(msDimension);
    MathUtils<double>::UnitCrossProduct(r_2, r_3, r_1);

    Matrix co_rotating_rotation = ZeroMatrix(msDimension);

    for (int i = 0; i < msDimension; ++i) {
        co_rotating_rotation(i, 0) = r_1[i];
        co_rotating_rotation(i, 1) = r_2[i];
        co_rotating_rotation(i, 2) = r_3[i];
    }
    return co_rotating_rotation;
}

Matrix BeamElement3D2N::SkewSymmetricMatrix(const Vector& rinput_vec) const
{
    Matrix skew = ZeroMatrix(msDimension);
    skew(0, 1) = -rinput_vec[2];
    skew(0, 2) =  rinput_vec[1];
    skew(1, 0) =  rinput_vec[2];
    skew(1, 2) = -rinput_vec[0];
    skew(2, 0) = -rinput_vec[1];
    skew(2, 1) =  rinput_vec[0];
    return skew;
}

Vector BeamElement3D2N::LocalDeformations() const
{
    const int local_deformation_possibilities(7);
    Vector local_deformations = ZeroVector(local_deformation_possibilities);

    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

    //local_deformations[0] = l-L;
    local_deformations[0] = (std::pow(l,2.0)-std::pow(L,2.0)) / (l+L);  //Crisfield

    Matrix reference_co_rot_matrix = CalculateInitialLocalCS();
    Matrix current_co_rot_matrix = CoRotatingCS();


    Matrix local_rot_node_1 = prod(mGlobalRotationNode1,reference_co_rot_matrix);
    local_rot_node_1 = prod(trans(current_co_rot_matrix),local_rot_node_1);

    Matrix local_rot_node_2 = prod(mGlobalRotationNode2,reference_co_rot_matrix);
    local_rot_node_2 = prod(trans(current_co_rot_matrix),local_rot_node_2);

    Matrix log_rotation_node_1 = LogRotationMatrix(local_rot_node_1);
    Matrix log_rotation_node_2 = LogRotationMatrix(local_rot_node_2);


    local_deformations[1] = log_rotation_node_1(2, 1);
    local_deformations[2] = log_rotation_node_1(0, 2);
    local_deformations[3] = log_rotation_node_1(1, 0);

    local_deformations[4] = log_rotation_node_2(2, 1);
    local_deformations[5] = log_rotation_node_2(0, 2);
    local_deformations[6] = log_rotation_node_2(1, 0);

    return local_deformations;
}

Matrix BeamElement3D2N::LogRotationMatrix(const Matrix& rRotationMatrix) const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double trace = rRotationMatrix(0,0)+rRotationMatrix(1,1)+rRotationMatrix(2,2);
    double phi = std::acos((trace - 1.0)/2.0);

    if (std::abs(std::abs(rRotationMatrix(0,0))-1.0) < numerical_limit) phi = std::asin(rRotationMatrix(0,1));

    Matrix log_matrix = ZeroMatrix(msDimension);

    if (std::abs(phi)>numerical_limit){
        log_matrix = rRotationMatrix-trans(rRotationMatrix);
        log_matrix *= phi / (2.0*std::sin(phi));
    }

    return log_matrix;
}

void BeamElement3D2N::UpdateGlobalNodalRotations()
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const Vector deformation_increment = GetIncrementDeformation();
    Vector d_phi_node_1 = ZeroVector(msDimension);
    Vector d_phi_node_2 = ZeroVector(msDimension);

    for (int i = 0; i < msDimension; ++i) {
        d_phi_node_1[i] = deformation_increment[msDimension+i];
        d_phi_node_2[i] = deformation_increment[(msDimension*3)+i];
    }

    // update node 1
    const Matrix d_phi_1_skew = SkewSymmetricMatrix(d_phi_node_1);
    const double d_phi_1_norm = MathUtils<double>::Norm(d_phi_node_1);

    if (d_phi_1_norm>numerical_limit){
        Matrix update_rotation = IdentityMatrix(msDimension);
        update_rotation += (std::sin(d_phi_1_norm)/d_phi_1_norm) * d_phi_1_skew;
        update_rotation += 0.50 * std::pow((std::sin(d_phi_1_norm/2.0)/(d_phi_1_norm/2.0)),2.0) * prod(d_phi_1_skew,d_phi_1_skew);

        mGlobalRotationNode1 = prod(update_rotation,mGlobalRotationNode1);
    }

    // update node 2
    const Matrix d_phi_2_skew = SkewSymmetricMatrix(d_phi_node_2);
    const double d_phi_2_norm = MathUtils<double>::Norm(d_phi_node_2);

    if (d_phi_2_norm>numerical_limit){
        Matrix update_rotation = IdentityMatrix(msDimension);
        update_rotation += (std::sin(d_phi_2_norm)/d_phi_2_norm) * d_phi_2_skew;
        update_rotation += 0.50 * std::pow((std::sin(d_phi_2_norm/2.0)/(d_phi_2_norm/2.0)),2.0) * prod(d_phi_2_skew,d_phi_2_skew);

        mGlobalRotationNode2 = prod(update_rotation,mGlobalRotationNode2);
    }

}


Vector BeamElement3D2N::GetIncrementDeformation() const
{
    KRATOS_TRY;
    return mDeformationCurrentIteration - mDeformationPreviousIteration;
    KRATOS_CATCH("")
}

Vector BeamElement3D2N::LocalInternalForces() const
{
    return prod(CreateElementStiffnessMatrix_Material(),LocalDeformations());
}

Vector BeamElement3D2N::LocalInternalIntermediateForces() const
{
    return prod(trans(BMatrixIntermediateA()),LocalInternalForces());
}

Vector BeamElement3D2N::CalculateGlobalNodalForces() const
{
    return prod(trans(BMatrixGlobal()),LocalInternalIntermediateForces());
}

Matrix BeamElement3D2N::DMatrix() const
{
    Vector r_1 = CurrentLocalAxis1();
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

    Matrix d_3 = (IdentityMatrix(msDimension)-outer_prod(r_1,r_1)) / l;
    Matrix d_total = ZeroMatrix(msElementSize);

    project(d_total, range(0,3),range(0,3)) += d_3;
    project(d_total, range(0,3),range(6,9)) -= d_3;
    project(d_total, range(6,9),range(0,3)) -= d_3;
    project(d_total, range(6,9),range(6,9)) += d_3;

    return d_total;
}

Matrix BeamElement3D2N::GMatrix() const
{
    Matrix reference_rotation = CalculateInitialLocalCS();
    Matrix current_cor_rotation = CoRotatingCS();

    Vector temp_vec = ZeroVector(msDimension);
    temp_vec[1] = 1.0;

    Vector q_1 = prod(reference_rotation,trans(temp_vec));
    q_1 = prod(mGlobalRotationNode1,q_1);

    Vector q_2 = prod(reference_rotation,trans(temp_vec));
    q_2 = prod(mGlobalRotationNode2,q_2);

    Vector q_mean = 0.5 * (q_1+q_2);


    Vector q_i = prod(trans(current_cor_rotation),q_mean);
    Vector q_1i = prod(trans(current_cor_rotation),q_1);
    Vector q_2i = prod(trans(current_cor_rotation),q_2);

    const double eta = q_i[0]/q_i[1];
    const double eta_11 = q_1i[0]/q_i[1];
    const double eta_12 = q_1i[1]/q_i[1];
    const double eta_21 = q_2i[0]/q_i[1];
    const double eta_22 = q_2i[1]/q_i[1];

    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

    Matrix G_transposed = ZeroMatrix(msDimension,msElementSize);
    G_transposed(0, 2) = eta/l;
    G_transposed(0, 3) = eta_12/2.0;
    G_transposed(0, 4) = -eta_11/2.0;
    G_transposed(0, 8) = -eta/l;
    G_transposed(0, 9) = eta_22/2.0;
    G_transposed(0, 10) = -eta_21/2.0;
    G_transposed(1, 2) = 1.0/l;
    G_transposed(1, 8) = -1.0/l;
    G_transposed(2, 1) = -1.0/l;
    G_transposed(2, 7) = 1.0/l;


    return trans(G_transposed);
}

Vector BeamElement3D2N::RVector() const
{
    Vector r_1 = CurrentLocalAxis1();
    Vector r = ZeroVector(msElementSize);

    for (int i=0;i<3;++i) r[i] = -r_1[i];
    for (int i=0;i<3;++i) r[msLocalSize+i] = r_1[i];
    return r;
}

Matrix BeamElement3D2N::PMatrix() const
{
    Matrix p = ZeroMatrix(msLocalSize,msElementSize);
    project(p, range(0,3),range(3,6)) += IdentityMatrix(msDimension);
    project(p, range(3,6),range(9,12)) += IdentityMatrix(msDimension);

    Matrix g = GMatrix();

    project(p, range(0,3),range(0,12)) -= trans(g);
    project(p, range(3,6),range(0,12)) -= trans(g);

    return p;
}

Matrix BeamElement3D2N::EMatrix() const
{
    Matrix e = ZeroMatrix(msElementSize);
    Matrix current_co_rotation_cs = CoRotatingCS();

    project(e, range(0,3),range(0,3)) += current_co_rotation_cs;
    project(e, range(3,6),range(3,6)) += current_co_rotation_cs;
    project(e, range(6,9),range(6,9)) += current_co_rotation_cs;
    project(e, range(9,12),range(9,12)) += current_co_rotation_cs;

    return e;
}

Matrix BeamElement3D2N::BMatrixGlobal() const
{
    const int local_deformation_possibilities(7);
    Matrix b_global = ZeroMatrix(local_deformation_possibilities,msElementSize);

    Vector r = RVector();
    for (SizeType i=0;i<msElementSize;++i) b_global(0,i) = r[i];

    Matrix P_E_t = prod(PMatrix(),trans(EMatrix()));

    project(b_global, range(1,7), range(0,12)) += P_E_t;

    return b_global;
}

Matrix BeamElement3D2N::QMatrix() const
{
    Vector internal_forces_a = LocalInternalIntermediateForces();
    Vector internal_moments = ZeroVector(msLocalSize);
    for (SizeType i=0;i<msLocalSize;++i) internal_moments[i] = internal_forces_a[i+1];

    Vector q_vec = prod(trans(PMatrix()),internal_moments);

    Vector n_1 = ZeroVector(msDimension);
    Vector n_2 = ZeroVector(msDimension);
    Vector m_1 = ZeroVector(msDimension);
    Vector m_2 = ZeroVector(msDimension);

    for (int i=0;i<msDimension;++i)
    {
        n_1[i] = q_vec[i];
        m_1[i] = q_vec[3+i];
        n_2[i] = q_vec[6+i];
        m_2[i] = q_vec[9+i];
    }

    Matrix q_total = ZeroMatrix(msElementSize,msDimension);
    project(q_total, range(0,3), range(0,3)) += SkewSymmetricMatrix(n_1);
    project(q_total, range(3,6), range(0,3)) += SkewSymmetricMatrix(m_1);
    project(q_total, range(6,9), range(0,3)) += SkewSymmetricMatrix(n_2);
    project(q_total, range(9,12), range(0,3)) += SkewSymmetricMatrix(m_2);

    return q_total;
}

Matrix BeamElement3D2N::InverseLocalRotation(const int node_nr) const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    Vector local_deformation = LocalDeformations();
    Vector nodal_rot = ZeroVector(msDimension);

    if (node_nr==1) for (int i=0;i<msDimension;++i) nodal_rot[i] = local_deformation[1+i];
    else if (node_nr==2) for (int i=0;i<msDimension;++i) nodal_rot[i] = local_deformation[4+i];
    else KRATOS_ERROR << "wrong node number in InverseLocalRotation" << std::endl;

    const double phi_norm = MathUtils<double>::Norm(nodal_rot);
    Matrix rot_inverse = IdentityMatrix(msDimension);

    if (phi_norm>numerical_limit){
        Vector u = nodal_rot/phi_norm;

        rot_inverse = IdentityMatrix(msDimension) * (phi_norm/2.0)/(std::tan(phi_norm/2.0));
        rot_inverse += outer_prod(u,u) * (1.0 - ((phi_norm/2.0)/std::tan(phi_norm/2.0)));
        rot_inverse -= 0.50 * SkewSymmetricMatrix(nodal_rot);
    }

    return rot_inverse;
}

Matrix BeamElement3D2N::BMatrixIntermediateA() const
{
    const int local_deformation_possibilities(7);
    Matrix b_a = ZeroMatrix(local_deformation_possibilities);
    b_a(0, 0) = 1.0;

    project(b_a, range(1,4), range(1,4)) += InverseLocalRotation(1);
    project(b_a, range(4,7), range(4,7)) += InverseLocalRotation(2);

    return b_a;
}

Vector BeamElement3D2N::AVector() const
{
    Vector a = ZeroVector(msDimension);
    Vector internal_forces_a = LocalInternalIntermediateForces();

    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double eta = GMatrix()(2,0)*l;

    a[1] = ((eta/l) * (internal_forces_a[1]+internal_forces_a[4])) - ((1.0/l) * (internal_forces_a[2]+internal_forces_a[5]));
    a[2] = (1.0 / l) * (internal_forces_a[3]+internal_forces_a[6]);

    return a;
}

void BeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    ConstCalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BeamElement3D2N::ConstCalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    ConstCalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    ConstCalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    KRATOS_CATCH("");
}

void BeamElement3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstCalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BeamElement3D2N::ConstCalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -= CalculateGlobalNodalForces();
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}


void BeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    ConstCalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BeamElement3D2N::ConstCalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
    rLeftHandSideMatrix += GlobalTangentStiffnessMatrix();
    KRATOS_CATCH("")
}

BoundedVector<double, BeamElement3D2N::msLocalSize>
BeamElement3D2N::GetCurrentNodalPosition() const
{
    BoundedVector<double, msLocalSize> current_nodal_position =
        ZeroVector(msLocalSize);
    for (unsigned int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension;
        current_nodal_position[index] =
            GetGeometry()[i].X0() +
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, 0);
        current_nodal_position[index + 1] =
            GetGeometry()[i].Y0() +
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, 0);
        current_nodal_position[index + 2] =
            GetGeometry()[i].Z0() +
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z, 0);
    }

    return current_nodal_position;
}

double BeamElement3D2N::CalculateShearModulus() const
{
    KRATOS_TRY;
    const double nu = GetProperties()[POISSON_RATIO];
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = E / (2.0 * (1.0 + nu));
    return G;
    KRATOS_CATCH("")
}

void BeamElement3D2N::CalculateConsistentMassMatrix(
    MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double A = GetProperties()[CROSS_AREA];
    const double Iy = GetProperties()[I22];
    const double Iz = GetProperties()[I33];


    double I_t = 0.00;
    if (GetProperties().Has(MASS_MOMENT_OF_INERTIA)){
        I_t = GetProperties()[MASS_MOMENT_OF_INERTIA];
    }
    else {
        //This is an approximation for a circular cross section
        I_t = Iy + Iz;
    }

    rMassMatrix(0,0) = (1.0/3.0)*A*L*rho;
    rMassMatrix(0,6) = 0.16666666666666669*A*L*rho;
    rMassMatrix(1,1) = 0.371428571428571*A*L*rho;
    rMassMatrix(1,5) = 0.052380952380952528*A*std::pow(L, 2)*rho;
    rMassMatrix(1,7) = 0.12857142857142856*A*L*rho;
    rMassMatrix(1,11) = -0.030952380952380898*A*std::pow(L, 2)*rho;
    rMassMatrix(2,2) = 0.371428571428571*A*L*rho;
    rMassMatrix(2,4) = -0.052380952380952639*A*std::pow(L, 2)*rho;
    rMassMatrix(2,8) = 0.12857142857142856*A*L*rho;
    rMassMatrix(2,10) = 0.030952380952380731*A*std::pow(L, 2)*rho;
    rMassMatrix(3,3) = (1.0/3.0)*I_t*L*rho;
    rMassMatrix(3,9) = 0.16666666666666669*I_t*L*rho;
    rMassMatrix(4,2) = -0.052380952380952639*A*std::pow(L, 2)*rho;
    rMassMatrix(4,4) = 0.009523809523809601*A*std::pow(L, 3)*rho;
    rMassMatrix(4,8) = -0.030952380952380842*A*std::pow(L, 2)*rho;
    rMassMatrix(4,10) = -0.0071428571428571175*A*std::pow(L, 3)*rho;
    rMassMatrix(5,1) = 0.052380952380952528*A*std::pow(L, 2)*rho;
    rMassMatrix(5,5) = 0.009523809523809601*A*std::pow(L, 3)*rho;
    rMassMatrix(5,7) = 0.030952380952380842*A*std::pow(L, 2)*rho;
    rMassMatrix(5,11) = -0.0071428571428571175*A*std::pow(L, 3)*rho;
    rMassMatrix(6,0) = 0.16666666666666669*A*L*rho;
    rMassMatrix(6,6) = 0.33333333333333326*A*L*rho;
    rMassMatrix(7,1) = 0.12857142857142856*A*L*rho;
    rMassMatrix(7,5) = 0.030952380952380842*A*std::pow(L, 2)*rho;
    rMassMatrix(7,7) = 0.37142857142857144*A*L*rho;
    rMassMatrix(7,11) = -0.052380952380952528*A*std::pow(L, 2)*rho;
    rMassMatrix(8,2) = 0.12857142857142856*A*L*rho;
    rMassMatrix(8,4) = -0.030952380952380842*A*std::pow(L, 2)*rho;
    rMassMatrix(8,8) = 0.37142857142857144*A*L*rho;
    rMassMatrix(8,10) = 0.052380952380952528*A*std::pow(L, 2)*rho;
    rMassMatrix(9,3) = 0.16666666666666669*I_t*L*rho;
    rMassMatrix(9,9) = 0.33333333333333326*I_t*L*rho;
    rMassMatrix(10,2) = 0.030952380952380731*A*std::pow(L, 2)*rho;
    rMassMatrix(10,4) = -0.0071428571428571175*A*std::pow(L, 3)*rho;
    rMassMatrix(10,8) = 0.052380952380952528*A*std::pow(L, 2)*rho;
    rMassMatrix(10,10) = 0.0095238095238095455*A*std::pow(L, 3)*rho;
    rMassMatrix(11,1) = -0.030952380952380898*A*std::pow(L, 2)*rho;
    rMassMatrix(11,5) = -0.0071428571428571175*A*std::pow(L, 3)*rho;
    rMassMatrix(11,7) = -0.052380952380952528*A*std::pow(L, 2)*rho;
    rMassMatrix(11,11) = 0.0095238095238095455*A*std::pow(L, 3)*rho;




    Matrix reference_co_rot_matrix = CalculateInitialLocalCS();
    Matrix local_rot_node_1 = prod(mGlobalRotationNode1,reference_co_rot_matrix);
    Matrix local_rot_node_2 = prod(mGlobalRotationNode2,reference_co_rot_matrix);

    Matrix current_co_rotating_rotation_matrix = ZeroMatrix(msElementSize);

    project(current_co_rotating_rotation_matrix, range(0,3),range(0,3)) += local_rot_node_1;
    project(current_co_rotating_rotation_matrix, range(3,6),range(3,6)) += local_rot_node_1;
    project(current_co_rotating_rotation_matrix, range(6,9),range(6,9)) += local_rot_node_2;
    project(current_co_rotating_rotation_matrix, range(9,12),range(9,12)) += local_rot_node_2;


    // rotate the consistent mass matrix
    /* Matrix current_co_rotating_rotation_matrix = EMatrix(); */



    rMassMatrix = prod(rMassMatrix,trans(current_co_rotating_rotation_matrix));
    rMassMatrix = prod(current_co_rotating_rotation_matrix,rMassMatrix);
    KRATOS_CATCH("")
}

void BeamElement3D2N::CalculateLumpedMassMatrix(
    MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    const double total_mass = A * L * rho;
    const double temp = 0.50 * total_mass;

    // w.r.t. Felippa - Chapter 31: LUMPED AND CONSISTENT MASS MATRICES - p.31–10
    const double rotational_inertia_lumped = total_mass * L * L * GetProperties()[LUMPED_MASS_ROTATION_COEFFICIENT];

    // translatonal mass
    for (int i = 0; i < msNumberOfNodes; ++i) {
        for (int j = 0; j < msDimension; ++j) {
            int index = i * (msDimension * 2) + j;
            rMassMatrix(index, index) = temp;

            // add rotational inertia
            rMassMatrix(index+msDimension, index+msDimension) = rotational_inertia_lumped;
        }
    }

    KRATOS_CATCH("")
}

void BeamElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    bool use_consistent_mass_matrix = false;

    if (GetProperties().Has(USE_CONSISTENT_MASS_MATRIX)) {
        use_consistent_mass_matrix = GetProperties()[USE_CONSISTENT_MASS_MATRIX];
    }

    if (use_consistent_mass_matrix) CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);
    else CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void BeamElement3D2N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        msElementSize);
}

BoundedVector<double, BeamElement3D2N::msElementSize> BeamElement3D2N::CalculateBodyForces() const
{
    KRATOS_TRY
    // getting shapefunctionvalues for linear SF
    const Matrix& Ncontainer =
        GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    BoundedVector<double, msDimension> equivalent_line_load =
        ZeroVector(msDimension);
    BoundedVector<double, msElementSize> body_forces_global =
        ZeroVector(msElementSize);

    const double A = GetProperties()[CROSS_AREA];
    const double l = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    // calculating equivalent line load
    for (int i = 0; i < msNumberOfNodes; ++i) {
        noalias(equivalent_line_load) +=
            (A * rho * Ncontainer(0, i)) *
            GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    // adding the nodal forces
    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        for (int j = 0; j < msDimension; ++j) {
            body_forces_global[j + index] =
                equivalent_line_load[j] * Ncontainer(0, i) * l;
        }
    }

    /* // adding the nodal moments
    CalculateAndAddWorkEquivalentNodalForcesLineLoad(equivalent_line_load,
            body_forces_global, l); */

    // return the total ForceVector
    return body_forces_global;
    KRATOS_CATCH("")
}

void BeamElement3D2N::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    BoundedVector<double, msElementSize> damping_residual_contribution = ZeroVector(msElementSize);
    // Calculate damping contribution to residual -->
    if ((GetProperties().Has(RAYLEIGH_ALPHA) ||
            GetProperties().Has(RAYLEIGH_BETA)) &&
            (rDestinationVariable != NODAL_INERTIA)) {
        Vector current_nodal_velocities = ZeroVector(msElementSize);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == FORCE_RESIDUAL) {

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            const SizeType index = msLocalSize * i;

            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == MOMENT_RESIDUAL) {

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            const SizeType index = (msLocalSize * i) + msDimension;

            array_1d<double, 3>& r_moment_residual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

            for (IndexType j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_moment_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    if (rDestinationVariable == NODAL_INERTIA) {
        Matrix element_mass_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_info; // Dummy
        CalculateMassMatrix(element_mass_matrix, temp_info);

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            double aux_nodal_mass = 0.0;
            array_1d<double, 3> aux_nodal_inertia(3, 0.0);

            const SizeType index = i * msLocalSize;

            for (IndexType j = 0; j < msElementSize; ++j) {
                aux_nodal_mass += element_mass_matrix(index, j);
                for (IndexType k = 0; k < msDimension; ++k) {
                    aux_nodal_inertia[k] += element_mass_matrix(index + msDimension + k, j);
                }
            }

            #pragma omp atomic
            GetGeometry()[i].GetValue(NODAL_MASS) += aux_nodal_mass;

            array_1d<double, 3>& r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
            for (IndexType k = 0; k < msDimension; ++k) {
                #pragma omp atomic
                r_nodal_inertia[k] += std::abs(aux_nodal_inertia[k]);
            }
        }
    }

    KRATOS_CATCH("")
}

void BeamElement3D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
    const BoundedVector<double, BeamElement3D2N::msDimension> ForceInput,
    BoundedVector<double, BeamElement3D2N::msElementSize>
    & rRightHandSideVector,
    const double GeometryLength) const
{
    KRATOS_TRY;
    // calculate orthogonal load vector
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    Vector geometric_orientation = ZeroVector(msDimension);
    geometric_orientation[0] =
        GetGeometry()[1].X() - GetGeometry()[0].X();
    geometric_orientation[1] =
        GetGeometry()[1].Y() - GetGeometry()[0].Y();
    if (msDimension == 3) {
        geometric_orientation[2] =
            GetGeometry()[1].Z() - GetGeometry()[0].Z();
    }

    const double vector_norm_a = MathUtils<double>::Norm(geometric_orientation);
    if (vector_norm_a > numerical_limit) {
        geometric_orientation /= vector_norm_a;
    }

    Vector line_load_direction = ZeroVector(msDimension);
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

    Vector node_a = ZeroVector(msDimension);
    node_a[0] = GetGeometry()[0].X();
    node_a[1] = GetGeometry()[0].Y();
    if (msDimension == 3) {
        node_a[2] = GetGeometry()[0].Z();
    }

    Vector node_b = ZeroVector(msDimension);
    node_b = node_a + line_load_direction;

    Vector node_c = ZeroVector(msDimension);
    node_c = node_a + (geometric_orientation * cos_angle);

    Vector load_orthogonal_direction = ZeroVector(msDimension);
    load_orthogonal_direction = node_b - node_c;
    const double vector_norm_c =
        MathUtils<double>::Norm(load_orthogonal_direction);
    if (vector_norm_c > numerical_limit) {
        load_orthogonal_direction /= vector_norm_c;
    }

    // now caluclate respective work equivilent nodal moments

    const double custom_moment =
        norm_force_vector_orthogonal * GeometryLength * GeometryLength / 12.00;

    Vector moment_a = ZeroVector(msDimension);
    moment_a = MathUtils<double>::CrossProduct(geometric_orientation,
               load_orthogonal_direction);
    moment_a *= custom_moment;

    for (int i = 0; i < msDimension; ++i) {
        rRightHandSideVector[(1 * msDimension) + i] += moment_a[i];
        rRightHandSideVector[(3 * msDimension) + i] -= moment_a[i];
    }

    KRATOS_CATCH("")
}

BeamElement3D2N::IntegrationMethod BeamElement3D2N::GetIntegrationMethod() const
{
    // do this to have 3GP as an output in GID
    return Kratos::GeometryData::GI_GAUSS_3;
}

void BeamElement3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }


    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        Vector nodal_forces_local_qe = LocalInternalForces();
        rOutput[0][0] = -1.0 * nodal_forces_local_qe[1] * 0.75 + nodal_forces_local_qe[4] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[1] * 0.50 + nodal_forces_local_qe[4] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[1] * 0.25 + nodal_forces_local_qe[4] * 0.75;

        rOutput[0][1] = -1.0 * (nodal_forces_local_qe[2] - (nodal_forces_local_qe[2]+nodal_forces_local_qe[5]) * 0.25);
        rOutput[1][1] = -1.0 * (nodal_forces_local_qe[2] - (nodal_forces_local_qe[2]+nodal_forces_local_qe[5]) * 0.50);
        rOutput[2][1] = -1.0 * (nodal_forces_local_qe[2] - (nodal_forces_local_qe[2]+nodal_forces_local_qe[5]) * 0.75);

        rOutput[0][2] = (nodal_forces_local_qe[3] - (nodal_forces_local_qe[3]+nodal_forces_local_qe[6]) * 0.25);
        rOutput[1][2] = (nodal_forces_local_qe[3] - (nodal_forces_local_qe[3]+nodal_forces_local_qe[6]) * 0.50);
        rOutput[2][2] = (nodal_forces_local_qe[3] - (nodal_forces_local_qe[3]+nodal_forces_local_qe[6]) * 0.75);
    } else if (rVariable == FORCE) {
        Vector nodal_forces_local_qe = LocalInternalForces();
        const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

        rOutput[0][0] = nodal_forces_local_qe[0];
        rOutput[1][0] = nodal_forces_local_qe[0];
        rOutput[2][0] = nodal_forces_local_qe[0];

        rOutput[0][1] = -1.0 * (nodal_forces_local_qe[3]+nodal_forces_local_qe[6]) / l;
        rOutput[1][1] = -1.0 * (nodal_forces_local_qe[3]+nodal_forces_local_qe[6]) / l;
        rOutput[2][1] = -1.0 * (nodal_forces_local_qe[3]+nodal_forces_local_qe[6]) / l;

        rOutput[0][2] = (nodal_forces_local_qe[2]+nodal_forces_local_qe[5]) / l;
        rOutput[1][2] = (nodal_forces_local_qe[2]+nodal_forces_local_qe[5]) / l;
        rOutput[2][2] = (nodal_forces_local_qe[2]+nodal_forces_local_qe[5]) / l;
    }

    else if (rVariable == LOCAL_AXIS_1) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = CoRotatingCS();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 0)[i];
        }
    } else if (rVariable == LOCAL_AXIS_2) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = CoRotatingCS();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 1)[i];
        }
    } else if (rVariable == LOCAL_AXIS_3) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = CoRotatingCS();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 2)[i];
        }
    }


    KRATOS_CATCH("")
}

void BeamElement3D2N::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    KRATOS_CATCH("")
}


int BeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
            << "The beam element works only in 3D and with 2 noded elements" << std::endl;

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
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
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

    KRATOS_ERROR_IF(!GetProperties().Has(I22) ||
                    GetProperties()[I22] <= numerical_limit)
            << "Please provide a reasonable value for \"I22\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(I33) ||
                    GetProperties()[I33] <= numerical_limit)
            << "Please provide a reasonable value for \"I33\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(TORSIONAL_INERTIA) ||
                    GetProperties()[TORSIONAL_INERTIA] <= numerical_limit)
            << "Please provide a reasonable value for \"TORSIONAL_INERTIA\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(POISSON_RATIO))
            << "\"POISSON_RATIO\" not provided for element #" << Id() << std::endl;

    if (Has(LOCAL_AXIS_2)) {
        array_1d<double, msDimension> direction_vector_x = ZeroVector(msDimension);
        array_1d<double, msDimension> direction_vector_y = ZeroVector(msDimension);
        array_1d<double, msLocalSize> reference_coordinates = ZeroVector(msLocalSize);

        reference_coordinates[0] = GetGeometry()[0].X0();
        reference_coordinates[1] = GetGeometry()[0].Y0();
        reference_coordinates[2] = GetGeometry()[0].Z0();
        reference_coordinates[3] = GetGeometry()[1].X0();
        reference_coordinates[4] = GetGeometry()[1].Y0();
        reference_coordinates[5] = GetGeometry()[1].Z0();

        for (unsigned int i = 0; i < msDimension; ++i) {
            direction_vector_x[i] = (reference_coordinates[i + msDimension] - reference_coordinates[i]);
        }

        const double vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        direction_vector_y = GetValue(LOCAL_AXIS_2);

        KRATOS_ERROR_IF(MathUtils<double>::Norm(direction_vector_y)<numerical_limit)
                << "Given LOCAL_AXIS_2 has length 0 for element " << Id() << std::endl;

        // a tollerance of 1e-3 allows for a rough deviation of 0.06 degrees from 90.0 degrees
        KRATOS_ERROR_IF(std::abs(MathUtils<double>::Dot(direction_vector_x,direction_vector_y))>1e-3)
                << "LOCAL_AXIS_1 is not perpendicular to LOCAL_AXIS_2 for element " << Id() << std::endl;
    }

    KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this)
                    < std::numeric_limits<double>::epsilon())
            << "Element #" << Id() << " has a length of zero!" << std::endl;

    return 0;

    KRATOS_CATCH("")
}

void BeamElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("NodalDeformationCurrent", mDeformationCurrentIteration);
    rSerializer.save("NodalDeformationPrevious", mDeformationPreviousIteration);
    rSerializer.save("GlobalRotationNode1", mGlobalRotationNode1);
    rSerializer.save("GlobalRotationNode2", mGlobalRotationNode2);
}

void BeamElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("NodalDeformationCurrent", mDeformationCurrentIteration);
    rSerializer.load("NodalDeformationPrevious", mDeformationPreviousIteration);
    rSerializer.load("GlobalRotationNode1", mGlobalRotationNode1);
    rSerializer.load("GlobalRotationNode2", mGlobalRotationNode2);
}

} // namespace Kratos.