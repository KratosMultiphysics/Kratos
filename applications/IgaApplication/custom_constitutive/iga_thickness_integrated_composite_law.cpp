// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: IgaApplication/license.txt
//
//  Main authors:    Max Friedrichs-Dachale, based on implementation by Malte Woidt and Kai Sommerwerk

#include "includes/mat_variables.h"
#include "utilities/math_utils.h"

#include "custom_constitutive/iga_thickness_integrated_composite_law.h"

#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{


// CLT and energy-method for transverse-shear correction after Whitney/ Chow 
IgaThicknessIntegratedCompositeLaw::IgaThicknessIntegratedCompositeLaw()
    : BaseType()
{
}

IgaThicknessIntegratedCompositeLaw::IgaThicknessIntegratedCompositeLaw(
    const std::vector<double>& rZCoordinates,
    const std::vector<double>& rEulerAngles,
    const std::vector<double>& rThicknesses)
    : BaseType(rZCoordinates, rEulerAngles, rThicknesses)
{
}

IgaThicknessIntegratedCompositeLaw::IgaThicknessIntegratedCompositeLaw(
    const IgaThicknessIntegratedCompositeLaw& rOther)
    : BaseType(rOther)
{
}

IgaThicknessIntegratedCompositeLaw::~IgaThicknessIntegratedCompositeLaw() = default;

ConstitutiveLaw::Pointer IgaThicknessIntegratedCompositeLaw::Clone() const
{
    return Kratos::make_shared<IgaThicknessIntegratedCompositeLaw>(*this);
}

ConstitutiveLaw::Pointer IgaThicknessIntegratedCompositeLaw::Create(
    Kratos::Parameters NewParameters) const
{

    KRATOS_ERROR_IF_NOT(NewParameters.Has("z_layer_coordinate_vector"))
        << "IgaThicknessIntegratedCompositeLaw: missing z_layer_coordinate_vector"
        << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters.Has("Euler_angle_layer_vector"))
        << "IgaThicknessIntegratedCompositeLaw: missing Euler_angle_layer_vector"
        << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters.Has("thickness_layer_vector"))
        << "IgaThicknessIntegratedCompositeLaw: missing thickness_layer_vector"
        << std::endl;

    const SizeType n = NewParameters["thickness_layer_vector"].size();
    KRATOS_ERROR_IF_NOT(NewParameters["z_layer_coordinate_vector"].size() == n)
        << "IgaThicknessIntegratedCompositeLaw: z/thickness layer vectors size mismatch"
        << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters["Euler_angle_layer_vector"].size() == n)
        << "IgaThicknessIntegratedCompositeLaw: Euler/thickness layer vectors size mismatch"
        << std::endl;

    std::vector<double> z(n), euler(n), thicknesses(n);
    for (IndexType i = 0; i < n; ++i) {
        z[i]           = NewParameters["z_layer_coordinate_vector"][i].GetDouble();
        euler[i]       = NewParameters["Euler_angle_layer_vector"][i].GetDouble();
        thicknesses[i] = NewParameters["thickness_layer_vector"][i].GetDouble();
    }

    return Kratos::make_shared<IgaThicknessIntegratedCompositeLaw>(z, euler, thicknesses);
}

void IgaThicknessIntegratedCompositeLaw::InitializeShearReductionFactors(
    const Properties& rMaterialProperties)
{
    const IndexType number_of_laws = mConstitutiveLaws.size();
    const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize();
    const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension();

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(rMaterialProperties);

    Vector null_strain_vector(subprop_strain_size);
    null_strain_vector.clear();
    parameters.SetStrainVector(null_strain_vector);

    Vector null_stress_vector(subprop_strain_size);
    null_stress_vector.clear();
    parameters.SetStressVector(null_stress_vector);

    Matrix constitutive_matrix(subprop_strain_size, subprop_strain_size);
    constitutive_matrix.clear();
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    Matrix generalized_constitutive_matrix(BaseType::VoigtSize, BaseType::VoigtSize);
    generalized_constitutive_matrix.clear();

    Flags& r_flags = parameters.GetOptions();
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    Matrix F(subprop_dimension, subprop_dimension);
    double aux_weight, detF, Euler_angle;

    double Gyz = 0.0;
    double Gxz = 0.0;

    BoundedMatrix<double, 2, 2> T;
    BoundedMatrix<double, 3, 3> Tvoigt;

    std::vector<Matrix> Qbar_layer(number_of_laws);
    std::vector<BoundedMatrix<double, 2, 2>> S_compliance(number_of_laws);

    const auto it_prop_begin = rMaterialProperties.GetSubProperties().begin();

    for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {

        const double half_thickness = 0.5 * mThicknesses[i_layer];
        const double z_inf_layer = mZCoordinates[i_layer] - half_thickness;
        const double z_sup_layer = mZCoordinates[i_layer] + half_thickness;
        const double weight = mThicknesses[i_layer];
        aux_weight = weight * mZCoordinates[i_layer];
        Euler_angle = mEulerAngles[i_layer];

        Properties &r_subprop = *(it_prop_begin + i_layer);
        parameters.SetMaterialProperties(r_subprop);
        CalculateShearModulus(Gyz, Gxz, parameters);

        //rotate a ply stiffness from the fibre axes to the laminate frame (MA Gl. 2.59 /  2.66)
        AdvancedConstitutiveLawUtilities<3>::CalculateRotationOperatorEuler1(Euler_angle, T);
        ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(T, Tvoigt);

        null_strain_vector.clear();
        noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(null_strain_vector);
        detF = MathUtils<double>::Det2(F);
        parameters.SetDeterminantF(detF);
        parameters.SetDeformationGradientF(F);

        // single UD-ply in-plane stiffness C_OA in the fibre frame (MA Gl. 2.57/2.58)
        mConstitutiveLaws[i_layer]->CalculateMaterialResponsePK2(parameters);

        // rotate the ply stiffness into the laminate frame, Qbar = T^T C_OA T (MA Gl. 2.59)
        Matrix Qbar = prod(trans(Tvoigt), Matrix(prod(constitutive_matrix, Tvoigt)));
        Qbar_layer[i_layer] = Qbar;

        // CLT integration into the ABD blocks (MA Gl. 2.60 - 2.62)
        noalias(project(generalized_constitutive_matrix, range(0, 3), range(0, 3))) += weight * Qbar;
        noalias(project(generalized_constitutive_matrix, range(3, 6), range(3, 6))) += (std::pow(z_sup_layer, 3) - std::pow(z_inf_layer, 3)) / 3.0 * Qbar;
        noalias(project(generalized_constitutive_matrix, range(0, 3), range(3, 6))) += aux_weight * Qbar;

        const double cos_e = std::cos(Euler_angle * Globals::Pi / 180.0);
        const double sin_e = std::sin(Euler_angle * Globals::Pi / 180.0);
        const double c2 = cos_e * cos_e;
        const double s2 = sin_e * sin_e;
        const double cs = cos_e * sin_e;
        //single-ply transverse-shear stiffness S rotated by the fibre angle (MA Gl. 2.65/2.66) then integrated over the thickness into G = integral S dz = sum S_k t_k (MA Gl. 2.67)
        const double S_xz_xz_shell = Gxz * c2 + Gyz * s2;
        const double S_yz_yz_shell = Gxz * s2 + Gyz * c2;
        const double S_yz_xz_shell = (Gxz - Gyz) * cs;
        generalized_constitutive_matrix(6, 6) += weight * S_yz_yz_shell;
        generalized_constitutive_matrix(7, 7) += weight * S_xz_xz_shell;
        generalized_constitutive_matrix(6, 7) += weight * S_yz_xz_shell;
        generalized_constitutive_matrix(7, 6) += weight * S_yz_xz_shell;

        BoundedMatrix<double, 2, 2> S_local;
        S_local(0, 0) = S_xz_xz_shell; S_local(0, 1) = S_yz_xz_shell;
        S_local(1, 0) = S_yz_xz_shell; S_local(1, 1) = S_yz_yz_shell;
        //per-layer shear compliance (MA Gl. 2.84/2.86)
        BoundedMatrix<double, 2, 2> S_inv;
        double det_S;
        MathUtils<double>::InvertMatrix2(S_local, S_inv, det_S);
        S_compliance[i_layer] = S_inv;
    }

    generalized_constitutive_matrix(3, 0) = generalized_constitutive_matrix(0, 3);
    generalized_constitutive_matrix(4, 0) = generalized_constitutive_matrix(0, 4);
    generalized_constitutive_matrix(5, 0) = generalized_constitutive_matrix(0, 5);
    generalized_constitutive_matrix(3, 1) = generalized_constitutive_matrix(1, 3);
    generalized_constitutive_matrix(4, 1) = generalized_constitutive_matrix(1, 4);
    generalized_constitutive_matrix(5, 1) = generalized_constitutive_matrix(1, 5);
    generalized_constitutive_matrix(3, 2) = generalized_constitutive_matrix(2, 3);
    generalized_constitutive_matrix(4, 2) = generalized_constitutive_matrix(2, 4);
    generalized_constitutive_matrix(5, 2) = generalized_constitutive_matrix(2, 5);

    Matrix ABD(6, 6);
    for (IndexType i = 0; i < 6; ++i)
        for (IndexType j = 0; j < 6; ++j)
            ABD(i, j) = generalized_constitutive_matrix(i, j);
    //MA Gl. 2.78)
    Matrix ABD_inv(6, 6);
    double det_ABD;
    MathUtils<double>::InvertMatrix(ABD, ABD_inv, det_ABD);

    array_1d<double, 3> Bstar0, Dstar0, Bstar1, Dstar1;
    for (IndexType i = 0; i < 3; ++i) {
        Bstar0[i] = ABD_inv(i, 3);   Dstar0[i] = ABD_inv(i + 3, 3);
        Bstar1[i] = ABD_inv(i, 4);   Dstar1[i] = ABD_inv(i + 3, 4);
    }

    const std::vector<double> gauss_w({5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0});
    const std::vector<double> gauss_xi({-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)});

    double g11_off = 0.0, g31_off = 0.0, g32_off = 0.0, g22_off = 0.0;
    double I1 = 0.0, I2 = 0.0;

    for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {
        const Matrix& Qbar = Qbar_layer[i_layer];
        const double half_thickness = 0.5 * mThicknesses[i_layer];
        const double z_lo = mZCoordinates[i_layer] - half_thickness;
        const double z_hi = mZCoordinates[i_layer] + half_thickness;

        auto dot_row = [&Qbar](IndexType r, const array_1d<double, 3>& v) {
            return Qbar(r, 0) * v[0] + Qbar(r, 1) * v[1] + Qbar(r, 2) * v[2];
        };

        const double p11 = dot_row(0, Bstar0), q11 = dot_row(0, Dstar0);//->g_11: sigma_xx from M_xx
        const double p31 = dot_row(2, Bstar0), q31 = dot_row(2, Dstar0);//->g_31: sigma_xy from M_xx
        const double p32 = dot_row(2, Bstar1), q32 = dot_row(2, Dstar1);//->g_32: sigma_xy from M_yy
        const double p22 = dot_row(1, Bstar1), q22 = dot_row(1, Dstar1);  //->g_22: sigma_yy from M_yy

        const BoundedMatrix<double, 2, 2>& Sc = S_compliance[i_layer];
        const double s44 = Sc(0, 0);
        const double s45 = Sc(0, 1);
        const double s55 = Sc(1, 1);

        //reconstruct the through-thickness transverse-shear-stress distribution from 3D  equilibrium: g_ik(z) = -integral_{-h/2}^z (p_ik + q_ik z') dz' (MA Gl. 2.74/2.82/2.83);
        for (IndexType gp = 0; gp < gauss_w.size(); ++gp) {
            const double z = mZCoordinates[i_layer] + half_thickness * gauss_xi[gp];
            const double dz = z - z_lo;
            const double dz2 = 0.5 * (z * z - z_lo * z_lo);
            const double g11 = g11_off - (p11 * dz + q11 * dz2);
            const double g31 = g31_off - (p31 * dz + q31 * dz2);
            const double g32 = g32_off - (p32 * dz + q32 * dz2);
            const double g22 = g22_off - (p22 * dz + q22 * dz2);
            const double w = gauss_w[gp] * half_thickness;
            //transverse-shear strain energy of the reconstructed distribution I = integral g^T S^-1 g dz (MA Gl. 2.84/2.86): I1 for shear direction xz (M_xx case), I2 for yz (M_yy case).
            I1 += w * (s44 * g11 * g11 + 2.0 * s45 * g11 * g31 + s55 * g31 * g31);
            I2 += w * (s44 * g32 * g32 + 2.0 * s45 * g22 * g32 + s55 * g22 * g22);
        }

        // advance the continuity offset by the integral of (p + q z) over the full layer, so the next layer starts from the accumulated shear stress at this interface (MA Gl. 2.82/2.83).
        const double dz_full = z_hi - z_lo;
        const double dz2_full = 0.5 * (z_hi * z_hi - z_lo * z_lo);
        g11_off -= (p11 * dz_full + q11 * dz2_full);
        g31_off -= (p31 * dz_full + q31 * dz2_full);
        g32_off -= (p32 * dz_full + q32 * dz2_full);
        g22_off -= (p22 * dz_full + q22 * dz2_full);
    }

    //energy equivalence of the FSDT (constant-shear) and the reconstructed-distribution energies (MA Gl. 2.86/2.87)
    mShearReductionFactors[1] = 1.0 / (generalized_constitutive_matrix(7, 7) * I1);
    mShearReductionFactors[0] = 1.0 / (generalized_constitutive_matrix(6, 6) * I2);
}

void IgaThicknessIntegratedCompositeLaw::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateLayeredResponse(rValues, PlyStressMeasure::Cauchy);
}

void IgaThicknessIntegratedCompositeLaw::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateLayeredResponse(rValues, PlyStressMeasure::PK2);
}

void IgaThicknessIntegratedCompositeLaw::CalculateLayeredResponse(
    ConstitutiveLaw::Parameters& rValues,
    PlyStressMeasure PlyMeasure)
{
    KRATOS_TRY

    Flags& r_flags = rValues.GetOptions();
    const auto& r_material_properties = rValues.GetMaterialProperties();
    const IndexType number_of_laws = mConstitutiveLaws.size();
    const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize();
    const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension();

    const bool flag_compute_constitutive_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_compute_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Vector generalized_strain_vector = rValues.GetStrainVector();
    Vector generalized_stress_vector(BaseType::VoigtSize);
    Matrix generalized_constitutive_matrix(BaseType::VoigtSize, BaseType::VoigtSize);
    generalized_constitutive_matrix.clear();
    generalized_stress_vector.clear();

    if (flag_compute_stress || flag_compute_constitutive_tensor) {

        Vector& r_stress_vector = rValues.GetStressVector();
        Vector& r_strain_vector = rValues.GetStrainVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        r_strain_vector.resize(subprop_strain_size, false);
        r_stress_vector.resize(subprop_strain_size, false);
        r_constitutive_matrix.resize(subprop_strain_size, subprop_strain_size, false);
        r_strain_vector.clear();
        r_stress_vector.clear();
        r_constitutive_matrix.clear();

        Matrix F(subprop_dimension, subprop_dimension);
        double weight, z_coord, aux_weight, detF, Euler_angle;

        const double stenberg_stabilization = 1.0; //TBD maybe add later(?)

        double Gyz = 0.0;
        double Gxz = 0.0;

        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();

        BoundedMatrix<double, 2, 2> T;
        BoundedMatrix<double, 3, 3> Tvoigt;

        for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {

            const double z_inf_layer = mZCoordinates[i_layer] - 0.5 * mThicknesses[i_layer];
            const double z_sup_layer = mZCoordinates[i_layer] + 0.5 * mThicknesses[i_layer];

            Properties &r_subprop = *(it_prop_begin + i_layer);
            rValues.SetMaterialProperties(r_subprop);

            CalculateShearModulus(Gyz, Gxz, rValues);

            weight = mThicknesses[i_layer];
            z_coord = mZCoordinates[i_layer];
            Euler_angle = mEulerAngles[i_layer];

            aux_weight = weight * z_coord;

            //in-plane strain at layer mid-plane height z: epsilon(z) = gamma + z chi (MA Gl. 2.45/2.76)
            r_strain_vector[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3];
            r_strain_vector[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4];
            r_strain_vector[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5];

            // rotate the laminate strain into the ply frame
            AdvancedConstitutiveLawUtilities<3>::CalculateRotationOperatorEuler1(Euler_angle, T);
            ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(T, Tvoigt);
            r_strain_vector = prod(Tvoigt, r_strain_vector);

            noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(r_strain_vector);
            detF = MathUtils<double>::Det2(F);
            rValues.SetDeterminantF(detF);
            rValues.SetDeformationGradientF(F);

            switch (PlyMeasure) {
                case PlyStressMeasure::Cauchy:
                    mConstitutiveLaws[i_layer]->CalculateMaterialResponseCauchy(rValues);
                    break;
                case PlyStressMeasure::PK2:
                    mConstitutiveLaws[i_layer]->CalculateMaterialResponsePK2(rValues);
                    break;
            }

            const double cos_e = std::cos(Euler_angle * Globals::Pi / 180.0);
            const double sin_e = std::sin(Euler_angle * Globals::Pi / 180.0);
            const double c2 = cos_e * cos_e;
            const double s2 = sin_e * sin_e;
            const double cs = cos_e * sin_e;
            const double S_xz_xz_shell = Gxz * c2 + Gyz * s2;
            const double S_yz_yz_shell = Gxz * s2 + Gyz * c2;
            const double S_yz_xz_shell = (Gxz - Gyz) * cs;
            // layer transverse-shear weight = thickness * energy shear-correction factor k^2
            const double w_red0 = weight * stenberg_stabilization * mShearReductionFactors[0];
            const double w_red1 = weight * stenberg_stabilization * mShearReductionFactors[1];

            if (flag_compute_stress) {
                r_stress_vector = prod(trans(Tvoigt), r_stress_vector);

                generalized_stress_vector[0] += r_stress_vector[0] * weight;
                generalized_stress_vector[1] += r_stress_vector[1] * weight;
                generalized_stress_vector[2] += r_stress_vector[2] * weight;

                generalized_stress_vector[3] += r_stress_vector[0] * aux_weight;
                generalized_stress_vector[4] += r_stress_vector[1] * aux_weight;
                generalized_stress_vector[5] += r_stress_vector[2] * aux_weight;

                // (MA Gl. 2.72/2.85)
                generalized_stress_vector[6] +=
                    w_red0 * ( S_yz_yz_shell * generalized_strain_vector[6]
                             + S_yz_xz_shell * generalized_strain_vector[7] );
                generalized_stress_vector[7] +=
                    w_red1 * ( S_yz_xz_shell * generalized_strain_vector[6]
                             + S_xz_xz_shell * generalized_strain_vector[7] );
            }

            if (flag_compute_constitutive_tensor) {
                r_constitutive_matrix = prod(trans(Tvoigt), Matrix(prod(r_constitutive_matrix, Tvoigt)));

                noalias(project(generalized_constitutive_matrix, range(0, 3), range(0, 3))) += weight * r_constitutive_matrix;
                noalias(project(generalized_constitutive_matrix, range(3, 6), range(3, 6))) += (std::pow(z_sup_layer, 3) - std::pow(z_inf_layer, 3)) / 3.0 * r_constitutive_matrix;
                noalias(project(generalized_constitutive_matrix, range(0, 3), range(3, 6))) += aux_weight * r_constitutive_matrix;

                generalized_constitutive_matrix(3, 0) = generalized_constitutive_matrix(0, 3);
                generalized_constitutive_matrix(4, 0) = generalized_constitutive_matrix(0, 4);
                generalized_constitutive_matrix(5, 0) = generalized_constitutive_matrix(0, 5);
                generalized_constitutive_matrix(3, 1) = generalized_constitutive_matrix(1, 3);
                generalized_constitutive_matrix(4, 1) = generalized_constitutive_matrix(1, 4);
                generalized_constitutive_matrix(5, 1) = generalized_constitutive_matrix(1, 5);
                generalized_constitutive_matrix(3, 2) = generalized_constitutive_matrix(2, 3);
                generalized_constitutive_matrix(4, 2) = generalized_constitutive_matrix(2, 4);
                generalized_constitutive_matrix(5, 2) = generalized_constitutive_matrix(2, 5);

                // transverse-shear stiffness block, k^2-corrected (MA Gl. 2.72/2.85)
                generalized_constitutive_matrix(6, 6) += S_yz_yz_shell * w_red0;
                generalized_constitutive_matrix(7, 7) += S_xz_xz_shell * w_red1;
                generalized_constitutive_matrix(6, 7) += S_yz_xz_shell * w_red0;
                generalized_constitutive_matrix(7, 6) += S_yz_xz_shell * w_red1;
            }
        }

        rValues.SetMaterialProperties(r_material_properties);
        r_strain_vector.resize(BaseType::VoigtSize, false);
        noalias(r_strain_vector) = generalized_strain_vector;

        if (flag_compute_stress) {
            r_stress_vector.resize(BaseType::VoigtSize, false);
            noalias(r_stress_vector) = generalized_stress_vector;
        }

        if (flag_compute_constitutive_tensor) {
            r_constitutive_matrix.resize(BaseType::VoigtSize, BaseType::VoigtSize, false);
            noalias(r_constitutive_matrix) = generalized_constitutive_matrix;
        }
    }
    KRATOS_CATCH("IgaThicknessIntegratedCompositeLaw::CalculateLayeredResponse")
}

void IgaThicknessIntegratedCompositeLaw::CalculateLayerGeneralizedStiffness(
    IndexType LayerIndex,
    ConstitutiveLaw::Parameters& rValues,
    Matrix& rQbar,
    Matrix& rShear,
    double& rZLo,
    double& rZHi)
{
    KRATOS_TRY

    const auto& r_material_properties = rValues.GetMaterialProperties();

    rZLo = mZCoordinates[LayerIndex] - 0.5 * mThicknesses[LayerIndex];
    rZHi = mZCoordinates[LayerIndex] + 0.5 * mThicknesses[LayerIndex];
    const double Euler_angle = mEulerAngles[LayerIndex];

    Properties& r_subprop = *(r_material_properties.GetSubProperties().begin() + LayerIndex);
    rValues.SetMaterialProperties(r_subprop);

    double Gyz = 0.0, Gxz = 0.0;
    CalculateShearModulus(Gyz, Gxz, rValues);

    BoundedMatrix<double, 2, 2> T;
    BoundedMatrix<double, 3, 3> Tvoigt;
    AdvancedConstitutiveLawUtilities<3>::CalculateRotationOperatorEuler1(Euler_angle, T);
    ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(T, Tvoigt);

    Vector strain3 = ZeroVector(3);
    Vector stress3 = ZeroVector(3);
    Matrix C3 = ZeroMatrix(3, 3);
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    rValues.SetStrainVector(strain3);
    rValues.SetStressVector(stress3);
    rValues.SetConstitutiveMatrix(C3);
    Matrix F = IdentityMatrix(2);
    rValues.SetDeterminantF(1.0);
    rValues.SetDeformationGradientF(F);
    mConstitutiveLaws[LayerIndex]->CalculateMaterialResponsePK2(rValues);

    //per-layer in-plane stiffness rotated into the laminate frame (MA Gl. 2.57; rotation: MA Gl. 2.59), used for stress recovery
    if (rQbar.size1() != 3 || rQbar.size2() != 3) rQbar.resize(3, 3, false);
    noalias(rQbar) = prod(trans(Tvoigt), Matrix(prod(rValues.GetConstitutiveMatrix(), Tvoigt)));

    const double cos_e = std::cos(Euler_angle * Globals::Pi / 180.0);
    const double sin_e = std::sin(Euler_angle * Globals::Pi / 180.0);
    const double c2 = cos_e * cos_e, s2 = sin_e * sin_e, cs = cos_e * sin_e;
    const double S_xz = Gxz * c2 + Gyz * s2;
    const double S_yz = Gxz * s2 + Gyz * c2;
    const double S_cross = (Gxz - Gyz) * cs;
    //per-layer transverse-shear stiffness rotated by the fibre angle (MA Gl. 2.65/2.66) weighted by the energy shear-correction factors k^2 (MA Gl. 2.85-2.87).
    if (rShear.size1() != 2 || rShear.size2() != 2) rShear.resize(2, 2, false);
    rShear(0, 0) = S_yz * mShearReductionFactors[0];
    rShear(1, 1) = S_xz * mShearReductionFactors[1];
    rShear(0, 1) = S_cross * mShearReductionFactors[0];
    rShear(1, 0) = S_cross * mShearReductionFactors[1];

    rValues.SetMaterialProperties(r_material_properties);

    KRATOS_CATCH("IgaThicknessIntegratedCompositeLaw::CalculateLayerGeneralizedStiffness")
}

int IgaThicknessIntegratedCompositeLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY


    int aux_out = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const auto it_prop_begin = rMaterialProperties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        const auto& r_subprop = *(it_prop_begin + i_layer);

        if (r_subprop.Has(SHEAR_MODULUS_XZ)) {
            KRATOS_ERROR_IF(r_subprop[SHEAR_MODULUS_XZ] <= 0.0)
                << "IgaThicknessIntegratedCompositeLaw: layer " << i_layer
                << " SHEAR_MODULUS_XZ (Gxz) must be positive, got "
                << r_subprop[SHEAR_MODULUS_XZ] << "." << std::endl;
        }
        if (r_subprop.Has(SHEAR_MODULUS_YZ)) {
            KRATOS_ERROR_IF(r_subprop[SHEAR_MODULUS_YZ] <= 0.0)
                << "IgaThicknessIntegratedCompositeLaw: layer " << i_layer
                << " SHEAR_MODULUS_YZ (Gyz) must be positive, got "
                << r_subprop[SHEAR_MODULUS_YZ] << "." << std::endl;
        }
    }

    return aux_out;

    KRATOS_CATCH("IgaThicknessIntegratedCompositeLaw::Check")
}

}
