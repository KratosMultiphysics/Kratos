//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alessandro Franci
//                   Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_cut_fem_element.h"
#include "includes/cfd_variables.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include <cmath>

namespace Kratos
{

    template <unsigned int TDim>
    Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::Clone(
        IndexType NewId,
        NodesArrayType const &rThisNodes) const
    {

        TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement new_element(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

        new_element.SetData(this->GetData());
        new_element.SetFlags(this->GetFlags());

        return Kratos::make_intrusive<TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement>(new_element);
    }

    template <unsigned int TDim>
    int TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        int ierr = BaseType::Check(rCurrentProcessInfo);

        return ierr;

        KRATOS_CATCH("");
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateLocalMomentumEquations(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {
        // Volume Navier-Stokes contribution
        // Note that this uses the CalculateGeometryData below, meaning that if it is cut, it already does the subintegration
        BaseType::CalculateLocalMomentumEquations(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

        // If intersected, add the boundary contribution
        if (IsCut())
        {
            // Calculate intersection Gauss points geometry data
            Matrix interface_N;
            ShapeFunctionDerivativesArrayType interface_DN_DX;
            Vector interface_gauss_weights;
            ModifiedShapeFunctions::AreaNormalsContainerType interface_unit_normals;
            CalculateIntersectionGeometryData(interface_DN_DX, interface_N, interface_gauss_weights, interface_unit_normals);

            // Create an auxiliary elemental data container for the boundary terms integration
            ElementalVariables elemental_variables;
            this->InitializeElementalVariables(elemental_variables);

            // Get other data
            const double kappa = rCurrentProcessInfo[PENALTY_COEFFICIENT];
            KRATOS_ERROR_IF(kappa < 1.0e-12) << "'PENALTY_COEFFICIENT' is zero." << std::endl;
            const double h = this->ElementSize();

            // Check if the element is "slip"
            // Note that if one node is flagged as SLIP the entire element is considered slip
            bool is_slip = false;
            const auto &r_geom = this->GetGeometry();
            for (const auto& r_node : r_geom) {
                if (r_node.Is(SLIP)) {
                    is_slip = true;
                    break;
                }
            }

            // Interface Gauss points loop
            array_1d<double,3> vel_gauss;
            const double rho = this->mMaterialDensity;
            const double dt = rCurrentProcessInfo[DELTA_TIME];
            const double theta = this->GetThetaMomentum();
            array_1d<double, TDim> proj_dev_stress;
            const std::size_t n_nodes = r_geom.PointsNumber();
            const unsigned int n_int_gauss_pts = interface_gauss_weights.size();
            for (unsigned int g = 0; g < n_int_gauss_pts; g++)
            {
                // Get interface Gauss point data
                const double g_weight = interface_gauss_weights[g];
                const auto g_DN_DX = interface_DN_DX[g];
                const auto g_shape_functions = row(interface_N, g);
                const auto &r_g_unit_normal = interface_unit_normals[g];

                // Calculate the mechanical response at the interface Gauss point to get the viscous stress
                this->CalcMechanicsUpdated(elemental_variables, rCurrentProcessInfo, g_DN_DX);

                auto &r_strain_vector = elemental_variables.SpatialDefRate;
                auto &r_dev_stress_vector = elemental_variables.UpdatedDeviatoricCauchyStress;
                auto &r_constitutive_matrix = elemental_variables.ConstitutiveMatrix;

                auto constitutive_law_values = ConstitutiveLaw::Parameters(
                    r_geom,
                    this->GetProperties(),
                    rCurrentProcessInfo);

                auto &r_constitutive_law_options = constitutive_law_values.GetOptions();
                r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
                r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

                constitutive_law_values.SetShapeFunctionsValues(g_shape_functions);
                constitutive_law_values.SetStrainVector(r_strain_vector);
                constitutive_law_values.SetStressVector(r_dev_stress_vector);
                constitutive_law_values.SetConstitutiveMatrix(r_constitutive_matrix);

                this->mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);
                this->UpdateStressTensor(elemental_variables);

                // Take dynamic viscosity from the bottom right corner of the constitutive matrix
                double mu = r_constitutive_matrix(StrainSize-1, StrainSize-1);
                const double max_mu_value = 10000; // TODO: check it
                if (mu > max_mu_value)
                {
                    mu = max_mu_value;
                }
                // Interpolate the pressure and velocity at the interface Gauss point to calculate the isochoric stress
                double pres_gauss = 0.0;
                vel_gauss.clear();
                for (std::size_t j = 0; j < n_nodes; ++j)
                {
                    const double p_0 = r_geom[j].FastGetSolutionStepValue(PRESSURE);
                    const double p_1 = r_geom[j].FastGetSolutionStepValue(PRESSURE,1);
                    pres_gauss += g_shape_functions[j] * (theta * p_0 + (1 - theta) * p_1);
                    const auto &r_vel_j_0 = r_geom[j].FastGetSolutionStepValue(VELOCITY);
                    const auto &r_vel_j_1 = r_geom[j].FastGetSolutionStepValue(VELOCITY,1);
                    noalias(vel_gauss) += g_shape_functions[j] * (theta * r_vel_j_0 + (1.0 - theta) * r_vel_j_1) ;
                }

                // Navier-Stokes traction boundary term
                VoigtStressNormalProjection(r_dev_stress_vector, r_g_unit_normal, proj_dev_stress);
                for (std::size_t i = 0; i < n_nodes; ++i)
                {
                    const double aux = g_weight * g_shape_functions[i];
                    for (std::size_t d = 0; d < TDim; ++d)
                    {
                        // Add Right Hand Side contribution
                        // TODO: Add the LHS terms
                        rRightHandSideVector(i * TDim + d) += aux * (pres_gauss * r_g_unit_normal[d] - proj_dev_stress[d]);
                    }
                }

                // Allocate and calculate auxiliary arrays
                array_1d<double,3> vel_j;
                array_1d<double,3> bc_vel;
                array_1d<double,3> wall_vel;
                BoundedMatrix<double, 3, 3> norm_proj_mat;
                BoundedMatrix<double, StrainSize, TDim * NumNodes> B;
                BoundedMatrix<double, TDim, StrainSize> voigt_normal;
                CalculateBMatrix(g_DN_DX, B);
                VoigtTransformForProduct(r_g_unit_normal, voigt_normal);
                BoundedMatrix<double, StrainSize, TDim *NumNodes> aux_BC = prod(trans(B), trans(r_constitutive_matrix));
                BoundedMatrix<double, TDim * NumNodes, TDim> aux_BC_proj = prod(aux_BC, trans(voigt_normal));

                // Cut-FEM boundary condition Nitsche imposition
                const double penalty_parameter = kappa * (mu / h + rho * (norm_2(vel_gauss) + h / dt));
                for (IndexType i = 0; i < n_nodes; ++i)
                {
                    for (IndexType j = 0; j < n_nodes; ++j)
                    {
                        // j-node data
                        const auto &r_vel_j_0 = r_geom[j].FastGetSolutionStepValue(VELOCITY);
                        const auto &r_vel_j_1 = r_geom[j].FastGetSolutionStepValue(VELOCITY,1);
                        vel_j = theta * r_vel_j_0 + (1.0 - theta) * r_vel_j_1;
                        noalias(wall_vel) = ZeroVector(3); // TODO: This should be the interpolation of the "structure" velocity in the future

                        // Check the boundary condition to be imposed (no-slip or pure slip)
                        if (is_slip) {
                            noalias(norm_proj_mat) = outer_prod(r_g_unit_normal, r_g_unit_normal);
                            noalias(bc_vel) = prod(norm_proj_mat, vel_j - wall_vel);
                        } else {
                            noalias(bc_vel) = vel_j - wall_vel;
                        }

                        // Assemble boundary condition RHS and LHS contributions
                        const double aux_1 = g_weight * penalty_parameter * g_shape_functions[i] * g_shape_functions[j];
                        const double aux_2 = g_weight * g_shape_functions[j];
                        for (IndexType d1 = 0; d1 < TDim; ++d1)
                        {
                            // Penalty term
                            rLeftHandSideMatrix(i * TDim + d1, j * TDim + d1) += aux_1;
                            rRightHandSideVector(i * TDim + d1) -= aux_1 * bc_vel[d1];
                            // Nitsche term (only viscous component)
                            for (IndexType d2 = 0; d2 < TDim; ++d2)
                            {
                                rLeftHandSideMatrix(i * TDim + d1, j * TDim + d2) -= aux_BC_proj(i * TDim + d1,  d2) * aux_2;
                                rRightHandSideVector(i * TDim + d1) += aux_BC_proj(i * TDim + d1, d2) * aux_2 * bc_vel[d1];
                            }
                        }
                    }
                }
            }
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateGeometryData(
        ShapeFunctionDerivativesArrayType &rDN_DX,
        Matrix &NContainer,
        Vector &rGaussWeights)
    {
        if (IsCut())
        {
            // Calculate cut element Gauss point values
            CalculateCutGeometryData(rDN_DX, NContainer, rGaussWeights);
        }
        else
        {
            // If not cut, we use the standard shape functions data calculator from the parent
            BaseType::CalculateGeometryData(rDN_DX, NContainer, rGaussWeights);
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateGeometryData(Vector &rGaussWeights)
    {
        if (IsCut())
        {
            CalculateCutGeometryData(rGaussWeights);
        }
        else
        {
            BaseType::CalculateGeometryData(rGaussWeights);
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateIntersectionGeometryData(
        ShapeFunctionDerivativesArrayType &rInterfaceDNDX,
        Matrix &rInterfaceN,
        Vector &rInterfaceGaussWeights,
        ModifiedShapeFunctions::AreaNormalsContainerType &rInterfaceUnitNormals)
    {
        const auto &r_geom = this->GetGeometry();

        // Auxiliary distance vector for the element subdivision utility
        Vector distances_vector(NumNodes);
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
            distances_vector[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Get the subintegration utility
        ModifiedShapeFunctions::Pointer p_mod_sh_func = nullptr;
        if constexpr (TDim == 2)
        {
            p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }
        else
        {
            p_mod_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }

        // Fluid side interface
        p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            rInterfaceN,
            rInterfaceDNDX,
            rInterfaceGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_1);

        // Fluid side interface normals
        p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(
            rInterfaceUnitNormals,
            GeometryData::IntegrationMethod::GI_GAUSS_1);

        for (unsigned int i = 0; i < rInterfaceUnitNormals.size(); ++i)
        {
            const double norm = norm_2(rInterfaceUnitNormals[i]);
            KRATOS_WARNING_IF("CalculateIntersectionGeometryData", norm < 1.0e-12) << "Normal is close to zero in element " << this->Id() << " cut interface." << std::endl;
            rInterfaceUnitNormals[i] /= norm;
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateCutGeometryData(
        ShapeFunctionDerivativesArrayType &rDNDX,
        Matrix &rN,
        Vector &rGaussWeights)
    {
        const auto &r_geom = this->GetGeometry();

        // Auxiliary distance vector for the element subdivision utility
        Vector distances_vector(NumNodes);
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
            distances_vector[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Get the subintegration utility
        ModifiedShapeFunctions::Pointer p_mod_sh_func = nullptr;
        if constexpr (TDim == 2)
        {
            p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }
        else
        {
            p_mod_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }

        // Fluid side
        p_mod_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            rN,
            rDNDX,
            rGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_1);
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateCutGeometryData(Vector &rGaussWeights)
    {
        const auto &r_geom = this->GetGeometry();

        // Auxiliary distance vector for the element subdivision utility
        Vector distances_vector(NumNodes);
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
            distances_vector[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Get the subintegration utility
        ModifiedShapeFunctions::Pointer p_mod_sh_func = nullptr;
        if constexpr (TDim == 2)
        {
            p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }
        else
        {
            p_mod_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }

        // Fluid side
        Matrix aux_N_container;
        p_mod_sh_func->ComputePositiveSideShapeFunctionsAndWeights(
            aux_N_container,
            rGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_1);
    }

    template <unsigned int TDim>
    bool TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::IsCut() const
    {
        const auto &r_geom = this->GetGeometry();
        SizeType n_pos = 0;
        SizeType n_neg = 0;
        for (const auto &r_node : r_geom)
        {
            if (r_node.FastGetSolutionStepValue(DISTANCE) > 0.0)
            {
                n_pos++;
            }
            else
            {
                n_neg++;
            }
        }

        return n_pos != 0 && n_neg != 0 ? true : false;
    }

    template <unsigned int TDim>
    bool TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::IsPositive() const
    {
        const auto &r_geom = this->GetGeometry();
        SizeType n_pos = 0;
        for (const auto &r_node : r_geom)
        {
            if (r_node.FastGetSolutionStepValue(DISTANCE) > 0.0)
            {
                n_pos++;
            }
        }

        return n_pos == NumNodes ? true : false;
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::VoigtStressNormalProjection(
        const Vector &rVoigtStress,
        const array_1d<double, 3> &rUnitNormal,
        array_1d<double, 2> &rProjectedStress)
    {
        rProjectedStress[0] = rVoigtStress[0] * rUnitNormal[0] + rVoigtStress[2] * rUnitNormal[1];
        rProjectedStress[1] = rVoigtStress[2] * rUnitNormal[0] + rVoigtStress[1] * rUnitNormal[1];
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::VoigtStressNormalProjection(
        const Vector &rVoigtStress,
        const array_1d<double, 3> &rUnitNormal,
        array_1d<double, 3> &rProjectedStress)
    {
        rProjectedStress[0] = rVoigtStress[0] * rUnitNormal[0] + rVoigtStress[3] * rUnitNormal[1] + rVoigtStress[5] * rUnitNormal[2];
        rProjectedStress[1] = rVoigtStress[3] * rUnitNormal[0] + rVoigtStress[1] * rUnitNormal[1] + rVoigtStress[4] * rUnitNormal[2];
        rProjectedStress[2] = rVoigtStress[5] * rUnitNormal[0] + rVoigtStress[4] * rUnitNormal[1] + rVoigtStress[2] * rUnitNormal[2];
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::CalculateBMatrix(
        const Matrix &rDNDX,
        BoundedMatrix<double, StrainSize, 2 * NumNodes> &rB)
    {
        rB.clear();
        IndexType index;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            index = 2 * i;
            rB(0, index + 0) = rDNDX(i, 0);
            rB(1, index + 1) = rDNDX(i, 1);
            rB(2, index + 0) = rDNDX(i, 1);
            rB(2, index + 1) = rDNDX(i, 0);
        }
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::CalculateBMatrix(
        const Matrix &rDNDX,
        BoundedMatrix<double, StrainSize, 3 * NumNodes> &rB)
    {
        rB.clear();
        IndexType index;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            index = 3 * i;
            rB(0, index + 0) = rDNDX(i, 0);
            rB(1, index + 1) = rDNDX(i, 1);
            rB(2, index + 2) = rDNDX(i, 2);
            rB(3, index + 0) = rDNDX(i, 1);
            rB(3, index + 1) = rDNDX(i, 0);
            rB(4, index + 1) = rDNDX(i, 2);
            rB(4, index + 2) = rDNDX(i, 1);
            rB(5, index + 0) = rDNDX(i, 2);
            rB(5, index + 2) = rDNDX(i, 0);
        }
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::VoigtTransformForProduct(
        const array_1d<double, 3> &rVector,
        BoundedMatrix<double, 2, StrainSize> &rVoigtMatrix)
    {

        rVoigtMatrix.clear();

        rVoigtMatrix(0, 0) = rVector(0);
        rVoigtMatrix(0, 2) = rVector(1);
        rVoigtMatrix(1, 1) = rVector(1);
        rVoigtMatrix(1, 2) = rVector(0);
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::VoigtTransformForProduct(
        const array_1d<double, 3> &rVector,
        BoundedMatrix<double, 3, StrainSize> &rVoigtMatrix)
    {

        rVoigtMatrix.clear();

        rVoigtMatrix(0, 0) = rVector(0);
        rVoigtMatrix(0, 3) = rVector(1);
        rVoigtMatrix(0, 5) = rVector(2);
        rVoigtMatrix(1, 1) = rVector(1);
        rVoigtMatrix(1, 3) = rVector(0);
        rVoigtMatrix(1, 4) = rVector(2);
        rVoigtMatrix(2, 2) = rVector(2);
        rVoigtMatrix(2, 4) = rVector(1);
        rVoigtMatrix(2, 5) = rVector(0);
    }

    template class TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>;
    template class TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>;

} // namespace Kratos
