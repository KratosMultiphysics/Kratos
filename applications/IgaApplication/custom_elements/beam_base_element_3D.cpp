//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti,
//                   Ricky Aristio, Alicia Knauer
//

// System includes
#include <cstddef>
#include <cmath>

// External includes

// Project includes
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/beam_base_element_3D.h"

namespace Kratos
{

    void BeamBaseElement3D::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 4 * number_of_control_points) {
            rResult.resize(4 * number_of_control_points, false);
        }

        const IndexType dof_position =
            r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 4;

            rResult[index] =
                r_geometry[i]
                    .GetDof(DISPLACEMENT_X, dof_position)
                    .EquationId();

            rResult[index + 1] =
                r_geometry[i]
                    .GetDof(DISPLACEMENT_Y, dof_position + 1)
                    .EquationId();

            rResult[index + 2] =
                r_geometry[i]
                    .GetDof(DISPLACEMENT_Z, dof_position + 2)
                    .EquationId();

            rResult[index + 3] =
                r_geometry[i]
                    .GetDof(ROTATION_X, dof_position + 3)
                    .EquationId();
        }

        KRATOS_CATCH("")
    }


    void BeamBaseElement3D::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(4 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(
                r_geometry[i].pGetDof(DISPLACEMENT_X));

            rElementalDofList.push_back(
                r_geometry[i].pGetDof(DISPLACEMENT_Y));

            rElementalDofList.push_back(
                r_geometry[i].pGetDof(DISPLACEMENT_Z));

            rElementalDofList.push_back(
                r_geometry[i].pGetDof(ROTATION_X));
        }

        KRATOS_CATCH("")
    }


    void BeamBaseElement3D::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_shape_functions =
            r_geometry.ShapeFunctionsValues();

        const SizeType number_of_integration_points =
            r_geometry.IntegrationPointsNumber();

        if (mConstitutiveLawVector.size() != number_of_integration_points) {
            mConstitutiveLawVector.resize(number_of_integration_points);
        }

        for (
            IndexType point_number = 0;
            point_number < number_of_integration_points;
            ++point_number)
        {
            mConstitutiveLawVector[point_number] =
                GetProperties()[CONSTITUTIVE_LAW]->Clone();

            mConstitutiveLawVector[point_number]->InitializeMaterial(
                r_properties,
                r_geometry,
                row(r_shape_functions, point_number));
        }

        KRATOS_CATCH("")
    }


    void BeamBaseElement3D::GetValuesVector(
        Vector& rValues,
        const int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        if (rValues.size() != number_of_nodes * 4) {
            rValues.resize(number_of_nodes * 4, false);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 4;

            const auto& r_displacement =
                r_geometry[i].FastGetSolutionStepValue(
                    DISPLACEMENT,
                    Step);

            const double rotation =
                r_geometry[i].FastGetSolutionStepValue(
                    ROTATION_X,
                    Step);

            rValues[index] = r_displacement[0];
            rValues[index + 1] = r_displacement[1];
            rValues[index + 2] = r_displacement[2];
            rValues[index + 3] = rotation;
        }
    }

    int BeamBaseElement3D::Check(
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        const double numerical_limit =
            std::numeric_limits<double>::epsilon();

        const auto& r_geometry = GetGeometry();

        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            if (r_geometry[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
                KRATOS_ERROR
                    << "missing variable DISPLACEMENT on node "
                    << r_geometry[i].Id() << std::endl;
            }

            if (r_geometry[i].SolutionStepsDataHas(ROTATION_X) == false) {
                KRATOS_ERROR
                    << "missing variable ROTATION_X on node "
                    << r_geometry[i].Id() << std::endl;
            }

            if (
                r_geometry[i].HasDofFor(DISPLACEMENT_X) == false ||
                r_geometry[i].HasDofFor(DISPLACEMENT_Y) == false ||
                r_geometry[i].HasDofFor(DISPLACEMENT_Z) == false)
            {
                KRATOS_ERROR
                    << "missing one of the dofs for the variable "
                    << "DISPLACEMENT on node "
                    << r_geometry[i].Id() << std::endl;
            }

            if (r_geometry[i].HasDofFor(ROTATION_X) == false) {
                KRATOS_ERROR
                    << "missing dof for the variable ROTATION_X on node "
                    << r_geometry[i].Id() << std::endl;
            }
        }

        KRATOS_ERROR_IF(
            !GetProperties().Has(CROSS_AREA) ||
            GetProperties()[CROSS_AREA] <= numerical_limit)
            << "Please provide a reasonable value for \"CROSS_AREA\" "
            << "for element #" << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(YOUNG_MODULUS) ||
            GetProperties()[YOUNG_MODULUS] <= numerical_limit)
            << "Please provide a reasonable value for \"YOUNG_MODULUS\" "
            << "for element #" << Id() << std::endl;

        KRATOS_ERROR_IF(
            !GetProperties().Has(DENSITY) ||
            GetProperties()[DENSITY] <= numerical_limit)
            << "Please provide a reasonable value for \"DENSITY\" "
            << "for element #" << Id()
            << ". Provided density: "
            << GetProperties()[DENSITY] << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(POISSON_RATIO))
            << "\"POISSON_RATIO\" not provided for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW))
            << "\"CONSTITUTIVE_LAW\" not provided for element #"
            << Id() << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    // ---------------------------------------------------------------------
    // Beam kinematics
    // ---------------------------------------------------------------------

    void BeamBaseElement3D::ComputeTangentVariation(
            const Vector& rShapeFunctionDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        Vector& rTangentVariation) const
        {
            rTangentVariation.resize(3 * mNumberOfDofs);
            rTangentVariation.clear();

            const double first_base_vector_norm =
                norm_2(rFirstBaseVector);

            const double first_base_vector_norm_cubed =
                pow(first_base_vector_norm, 3);

            Vector first_base_vector_variation;
            first_base_vector_variation.resize(
                3 * mNumberOfDofs);
            first_base_vector_variation.clear();

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0;
                    dof_index < mNumberOfDofs;
                    ++dof_index)
                {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    const size_t control_point_index =
                        dof_index / mDofsPerNode;

                    if (component == dof_component) {
                        first_base_vector_variation(
                            component * mNumberOfDofs +
                            dof_index) +=
                            rShapeFunctionDerivatives[
                                control_point_index];
                    }
                }
            }

            Vector base_vector_product_variation;
            base_vector_product_variation.resize(
                mNumberOfDofs);
            base_vector_product_variation.clear();

            for (size_t dof_index = 0;
                dof_index < mNumberOfDofs;
                ++dof_index)
            {
                for (size_t component = 0;
                    component < 3;
                    ++component)
                {
                    base_vector_product_variation(dof_index) +=
                        first_base_vector_variation[
                            component * mNumberOfDofs +
                            dof_index] *
                        rFirstBaseVector[component];
                }
            }

            for (size_t component = 0;
                component < 3;
                ++component)
            {
                for (size_t dof_index = 0;
                    dof_index < mNumberOfDofs;
                    ++dof_index)
                {
                    rTangentVariation(
                        component * mNumberOfDofs +
                        dof_index) +=
                        first_base_vector_variation[
                            component * mNumberOfDofs +
                            dof_index] /
                        first_base_vector_norm
                        -
                        rFirstBaseVector[component] *
                        base_vector_product_variation[dof_index] /
                        first_base_vector_norm_cubed;
                }
            }
        }


    void BeamBaseElement3D::ComputeTangentSecondVariation(
            const Vector& rShapeFunctionDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        Matrix& rTangentSecondVariation) const
        {
            rTangentSecondVariation.resize(
                3 * mNumberOfDofs,
                mNumberOfDofs);

            rTangentSecondVariation.clear();

            const double first_base_vector_norm =
                norm_2(rFirstBaseVector);

            const double first_base_vector_norm_cubed =
                pow(first_base_vector_norm, 3);

            const double first_base_vector_norm_fifth =
                pow(first_base_vector_norm, 5);

            Vector first_base_vector_variation;
            first_base_vector_variation.resize(
                3 * mNumberOfDofs);
            first_base_vector_variation.clear();

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0;
                    dof_index < mNumberOfDofs;
                    ++dof_index)
                {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    const size_t control_point_index =
                        dof_index / mDofsPerNode;

                    if (
                        component == dof_component &&
                        dof_component < 3)
                    {
                        first_base_vector_variation(
                            component * mNumberOfDofs +
                            dof_index) =
                            rShapeFunctionDerivatives[
                                control_point_index];
                    }
                }
            }

            Vector base_vector_product_variation;
            base_vector_product_variation.resize(
                mNumberOfDofs);
            base_vector_product_variation.clear();

            for (size_t dof_index = 0;
                dof_index < mNumberOfDofs;
                ++dof_index)
            {
                const size_t dof_component =
                    dof_index % mDofsPerNode;

                if (dof_component > 2) {
                    base_vector_product_variation(dof_index) = 0;
                } else {
                    for (size_t component = 0;
                        component < 3;
                        ++component)
                    {
                        base_vector_product_variation(dof_index) +=
                            first_base_vector_variation(
                                component * mNumberOfDofs +
                                dof_index) *
                            rFirstBaseVector[component];
                    }
                }
            }

            Matrix variation_product;
            variation_product.resize(
                mNumberOfDofs,
                mNumberOfDofs);

            variation_product.clear();

            for (size_t component = 0;
                component < 3;
                ++component)
            {
                for (size_t dof_index = 0;
                    dof_index < mNumberOfDofs;
                    ++dof_index)
                {
                    for (size_t second_dof_index = 0;
                        second_dof_index < mNumberOfDofs;
                        ++second_dof_index)
                    {
                        variation_product(
                            dof_index,
                            second_dof_index) +=
                            first_base_vector_variation[
                                component * mNumberOfDofs +
                                dof_index] *
                            first_base_vector_variation[
                                component * mNumberOfDofs +
                                second_dof_index];
                    }
                }
            }

            for (size_t component = 0;
                component < 3;
                ++component)
            {
                for (size_t dof_index = 0;
                    dof_index < mNumberOfDofs;
                    ++dof_index)
                {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    for (size_t second_dof_index = 0;
                        second_dof_index < mNumberOfDofs;
                        ++second_dof_index)
                    {
                        const size_t second_dof_component =
                            second_dof_index % mDofsPerNode;

                        if (
                            dof_component > 2 ||
                            second_dof_component > 2)
                        {
                            continue;
                        }

                        rTangentSecondVariation(
                            component * mNumberOfDofs +
                            dof_index,
                            second_dof_index) +=
                            3.0 *
                            base_vector_product_variation[dof_index] *
                            base_vector_product_variation[
                                second_dof_index] *
                            rFirstBaseVector[component] /
                            first_base_vector_norm_fifth;

                        rTangentSecondVariation(
                            component * mNumberOfDofs +
                            dof_index,
                            second_dof_index) +=
                            -(
                                first_base_vector_variation[
                                    component * mNumberOfDofs +
                                    dof_index] *
                                base_vector_product_variation[
                                    second_dof_index]
                                +
                                first_base_vector_variation[
                                    component * mNumberOfDofs +
                                    second_dof_index] *
                                base_vector_product_variation[
                                    dof_index])
                            /
                            first_base_vector_norm_cubed;

                        rTangentSecondVariation(
                            component * mNumberOfDofs +
                            dof_index,
                            second_dof_index) +=
                            -variation_product(
                                dof_index,
                                second_dof_index) *
                            rFirstBaseVector[component] /
                            first_base_vector_norm_cubed;
                    }
                }
            }
        }


    void BeamBaseElement3D::ComputeTangentDerivativeVariation(
            const Vector& rShapeFunctionDerivatives,
            const Vector& rShapeFunctionSecondDerivatives,
            const array_1d<double, 3>& rFirstBaseVector,
        const array_1d<double, 3>& rSecondBaseVector,
        Vector& rTangentDerivativeVariation) const
        {
            rTangentDerivativeVariation.resize(3 * mNumberOfDofs);
            rTangentDerivativeVariation.clear();

            const double first_second_base_vector_product =
                inner_prod(rFirstBaseVector, rSecondBaseVector);

            const double first_base_vector_norm =
                norm_2(rFirstBaseVector);

            Vector first_base_vector_variation;
            first_base_vector_variation.resize(3 * mNumberOfDofs);
            first_base_vector_variation.clear();

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    const size_t control_point_index =
                        dof_index / mDofsPerNode;

                    if (component == dof_component) {
                        first_base_vector_variation(
                            component * mNumberOfDofs + dof_index) =
                            rShapeFunctionDerivatives[control_point_index];
                    }
                }
            }

            Vector second_base_vector_variation;
            second_base_vector_variation.resize(3 * mNumberOfDofs);
            second_base_vector_variation.clear();

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    const size_t control_point_index =
                        dof_index / mDofsPerNode;

                    if (component == dof_component) {
                        second_base_vector_variation(
                            component * mNumberOfDofs + dof_index) =
                            rShapeFunctionSecondDerivatives[control_point_index];
                    }
                }
            }

            Vector first_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            Vector second_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            Vector first_base_vector_dot_second_variation =
                ZeroVector(mNumberOfDofs);

            for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                const size_t dof_component =
                    dof_index % mDofsPerNode;

                if (dof_component > 2) {
                    continue;
                }

                for (size_t component = 0; component < 3; ++component) {
                    first_base_vector_dot_first_variation[dof_index] +=
                        first_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rFirstBaseVector[component];

                    second_base_vector_dot_first_variation[dof_index] +=
                        first_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rSecondBaseVector[component];

                    first_base_vector_dot_second_variation[dof_index] +=
                        second_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rFirstBaseVector[component];
                }
            }

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    if (dof_component > 2) {
                        rTangentDerivativeVariation[
                            component * mNumberOfDofs + dof_index] = 0.0;

                        continue;
                    }

                    rTangentDerivativeVariation[
                        component * mNumberOfDofs + dof_index] +=
                        second_base_vector_variation[
                            component * mNumberOfDofs + dof_index] /
                            first_base_vector_norm
                        -
                        rSecondBaseVector[component] *
                            first_base_vector_dot_first_variation[dof_index] /
                            pow(first_base_vector_norm, 3);

                    rTangentDerivativeVariation[
                        component * mNumberOfDofs + dof_index] +=
                        3.0 *
                            first_base_vector_dot_first_variation[dof_index] *
                            first_second_base_vector_product *
                            rFirstBaseVector[component] /
                            pow(first_base_vector_norm, 5)
                        -
                        (
                            first_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                                first_second_base_vector_product
                            +
                            (
                                first_base_vector_dot_second_variation[dof_index]
                                +
                                second_base_vector_dot_first_variation[dof_index]
                            ) *
                            rFirstBaseVector[component]
                        ) /
                        pow(first_base_vector_norm, 3);
                }
            }
        }


    void BeamBaseElement3D::ComputeTangentDerivativeSecondVariation(
            const Vector& rShapeFunctionDerivatives,
            const Vector& rShapeFunctionSecondDerivatives,
            const array_1d<double, 3>& rFirstBaseVector,
        const array_1d<double, 3>& rSecondBaseVector,
        Matrix& rTangentDerivativeSecondVariation) const
        {
            rTangentDerivativeSecondVariation.resize(
                3 * mNumberOfDofs,
                mNumberOfDofs);

            rTangentDerivativeSecondVariation.clear();

            const double first_second_base_vector_product =
                inner_prod(rFirstBaseVector, rSecondBaseVector);

            const double first_base_vector_norm =
                norm_2(rFirstBaseVector);

            const double first_base_vector_norm_cubed =
                pow(first_base_vector_norm, 3);

            const double first_base_vector_norm_fifth =
                pow(first_base_vector_norm, 5);

            const double first_base_vector_norm_seventh =
                pow(first_base_vector_norm, 7);

            Vector first_base_vector_variation =
                ZeroVector(3 * mNumberOfDofs);

            Vector second_base_vector_variation =
                ZeroVector(3 * mNumberOfDofs);

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    const size_t control_point_index =
                        dof_index / mDofsPerNode;

                    if (component == dof_component) {
                        first_base_vector_variation(
                            component * mNumberOfDofs + dof_index) +=
                            rShapeFunctionDerivatives[control_point_index];

                        second_base_vector_variation(
                            component * mNumberOfDofs + dof_index) +=
                            rShapeFunctionSecondDerivatives[control_point_index];
                    }
                }
            }

            Vector first_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            Vector second_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            Vector first_base_vector_dot_second_variation =
                ZeroVector(mNumberOfDofs);

            for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                const size_t dof_component =
                    dof_index % mDofsPerNode;

                if (dof_component > 2) {
                    continue;
                }

                for (size_t component = 0; component < 3; ++component) {
                    first_base_vector_dot_first_variation[dof_index] +=
                        first_base_vector_variation(
                            component * mNumberOfDofs + dof_index) *
                        rFirstBaseVector[component];

                    second_base_vector_dot_first_variation[dof_index] +=
                        rSecondBaseVector[component] *
                        first_base_vector_variation(
                            component * mNumberOfDofs + dof_index);

                    first_base_vector_dot_second_variation[dof_index] +=
                        rFirstBaseVector[component] *
                        second_base_vector_variation(
                            component * mNumberOfDofs + dof_index);
                }
            }

            Matrix first_variation_product(
                mNumberOfDofs,
                mNumberOfDofs,
                0.0);

            Matrix first_second_variation_product(
                mNumberOfDofs,
                mNumberOfDofs,
                0.0);

            Matrix second_first_variation_product(
                mNumberOfDofs,
                mNumberOfDofs,
                0.0);

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    for (
                        size_t second_dof_index = 0;
                        second_dof_index < mNumberOfDofs;
                        ++second_dof_index)
                    {
                        first_variation_product(
                            dof_index,
                            second_dof_index) +=
                            first_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                            first_base_vector_variation[
                                component * mNumberOfDofs + second_dof_index];

                        first_second_variation_product(
                            dof_index,
                            second_dof_index) +=
                            first_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                            second_base_vector_variation[
                                component * mNumberOfDofs + second_dof_index];

                        second_first_variation_product(
                            dof_index,
                            second_dof_index) +=
                            second_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                            first_base_vector_variation[
                                component * mNumberOfDofs + second_dof_index];
                    }
                }
            }

            for (size_t component = 0; component < 3; ++component) {
                for (
                    size_t second_dof_index = 0;
                    second_dof_index < mNumberOfDofs;
                    ++second_dof_index)
                {
                    for (
                        size_t dof_index = 0;
                        dof_index < mNumberOfDofs;
                        ++dof_index)
                    {
                        const size_t dof_component =
                            dof_index % mNumberOfDofs;

                        const size_t second_dof_component =
                            second_dof_index % mNumberOfDofs;

                        if (
                            dof_component > 2 ||
                            second_dof_component > 2)
                        {
                            continue;
                        }

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            -(
                                second_base_vector_variation(
                                    component * mNumberOfDofs + dof_index) *
                                first_base_vector_dot_first_variation[
                                    second_dof_index]
                                +
                                second_base_vector_variation(
                                    component * mNumberOfDofs + second_dof_index) *
                                first_base_vector_dot_first_variation[
                                    dof_index]
                                +
                                rSecondBaseVector[component] *
                                first_variation_product(
                                    dof_index,
                                    second_dof_index)
                            ) /
                            first_base_vector_norm_cubed;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            3.0 *
                            first_base_vector_dot_first_variation[dof_index] *
                            rSecondBaseVector[component] *
                            first_base_vector_dot_first_variation[
                                second_dof_index] /
                            first_base_vector_norm_fifth;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            3.0 *
                            first_variation_product(
                                dof_index,
                                second_dof_index) *
                            first_second_base_vector_product *
                            rFirstBaseVector[component] /
                            first_base_vector_norm_fifth;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            3.0 *
                            first_base_vector_dot_first_variation[dof_index] *
                            (
                                second_base_vector_dot_first_variation[
                                    second_dof_index]
                                +
                                first_base_vector_dot_second_variation[
                                    second_dof_index]
                            ) *
                            rFirstBaseVector[component] /
                            first_base_vector_norm_fifth;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            3.0 *
                            first_base_vector_dot_first_variation[dof_index] *
                            first_base_vector_variation[
                                component * mNumberOfDofs +
                                second_dof_index] *
                            first_second_base_vector_product /
                            first_base_vector_norm_fifth;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            3.0 *
                            first_base_vector_dot_first_variation[
                                second_dof_index] *
                            first_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                            first_second_base_vector_product /
                            first_base_vector_norm_fifth;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            -15.0 *
                            first_base_vector_dot_first_variation[dof_index] *
                            first_base_vector_dot_first_variation[
                                second_dof_index] *
                            first_second_base_vector_product *
                            rFirstBaseVector[component] /
                            first_base_vector_norm_seventh;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            -(
                                (
                                    first_second_variation_product(
                                        dof_index,
                                        second_dof_index)
                                    +
                                    second_first_variation_product(
                                        dof_index,
                                        second_dof_index)
                                ) *
                                rFirstBaseVector[component]
                                +
                                (
                                    first_base_vector_dot_second_variation[
                                        dof_index]
                                    +
                                    second_base_vector_dot_first_variation[
                                        dof_index]
                                ) *
                                first_base_vector_variation[
                                    component * mNumberOfDofs +
                                    second_dof_index]
                                +
                                (
                                    first_base_vector_dot_second_variation[
                                        second_dof_index]
                                    +
                                    second_base_vector_dot_first_variation[
                                        second_dof_index]
                                ) *
                                first_base_vector_variation[
                                    component * mNumberOfDofs + dof_index]
                            ) /
                            first_base_vector_norm_cubed;

                        rTangentDerivativeSecondVariation(
                            component * mNumberOfDofs + dof_index,
                            second_dof_index) +=
                            3.0 *
                            first_base_vector_dot_first_variation[
                                second_dof_index] *
                            (
                                second_base_vector_dot_first_variation[
                                    dof_index]
                                +
                                first_base_vector_dot_second_variation[
                                    dof_index]
                            ) *
                            rFirstBaseVector[component] /
                            first_base_vector_norm_fifth;
                    }
                }
            }
        }

    
    void BeamBaseElement3D::ComputeTangentSecondDerivativeVariation(
            const Vector& rShapeFunctionDerivatives,
            const Vector& rShapeFunctionSecondDerivatives,
            const Vector& rShapeFunctionThirdDerivatives,
            const array_1d<double, 3>& rFirstBaseVector,
            const array_1d<double, 3>& rSecondBaseVector,
        const array_1d<double, 3>& rThirdBaseVector,
        Vector& rTangentSecondDerivativeVariation) const
        {
            rTangentSecondDerivativeVariation.resize(
                3 * mNumberOfDofs);

            rTangentSecondDerivativeVariation.clear();

            const double first_second_base_vector_product =
                inner_prod(
                    rFirstBaseVector,
                    rSecondBaseVector);

            const double first_second_base_vector_product_derivative =
                inner_prod(
                    rSecondBaseVector,
                    rSecondBaseVector)
                +
                inner_prod(
                    rFirstBaseVector,
                    rThirdBaseVector);

            const double first_base_vector_norm =
                norm_2(rFirstBaseVector);

            Vector first_base_vector_variation =
                ZeroVector(3 * mNumberOfDofs);

            Vector second_base_vector_variation =
                ZeroVector(3 * mNumberOfDofs);

            Vector third_base_vector_variation =
                ZeroVector(3 * mNumberOfDofs);

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    const size_t control_point_index =
                        dof_index / mDofsPerNode;

                    if (component == dof_component) {
                        first_base_vector_variation(
                            component * mNumberOfDofs + dof_index) =
                            rShapeFunctionDerivatives[
                                control_point_index];

                        second_base_vector_variation(
                            component * mNumberOfDofs + dof_index) =
                            rShapeFunctionSecondDerivatives[
                                control_point_index];

                        third_base_vector_variation(
                            component * mNumberOfDofs + dof_index) =
                            rShapeFunctionThirdDerivatives[
                                control_point_index];
                    }
                }
            }

            Vector first_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            Vector second_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            Vector first_base_vector_dot_second_variation =
                ZeroVector(mNumberOfDofs);

            Vector first_base_vector_dot_third_variation =
                ZeroVector(mNumberOfDofs);

            Vector second_base_vector_dot_second_variation =
                ZeroVector(mNumberOfDofs);

            Vector third_base_vector_dot_first_variation =
                ZeroVector(mNumberOfDofs);

            for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                const size_t dof_component =
                    dof_index % mDofsPerNode;

                if (dof_component > 2) {
                    continue;
                }

                for (size_t component = 0; component < 3; ++component) {
                    first_base_vector_dot_first_variation[dof_index] +=
                        first_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rFirstBaseVector[component];

                    second_base_vector_dot_first_variation[dof_index] +=
                        first_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rSecondBaseVector[component];

                    first_base_vector_dot_second_variation[dof_index] +=
                        second_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rFirstBaseVector[component];

                    first_base_vector_dot_third_variation[dof_index] +=
                        third_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rFirstBaseVector[component];

                    second_base_vector_dot_second_variation[dof_index] +=
                        second_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rSecondBaseVector[component];

                    third_base_vector_dot_first_variation[dof_index] +=
                        first_base_vector_variation[
                            component * mNumberOfDofs + dof_index] *
                        rThirdBaseVector[component];
                }
            }

            const Vector first_second_product_variation =
                first_base_vector_dot_second_variation +
                second_base_vector_dot_first_variation;

            const Vector first_second_product_derivative_variation =
                2.0 * second_base_vector_dot_second_variation
                +
                third_base_vector_dot_first_variation
                +
                first_base_vector_dot_third_variation;

            for (size_t component = 0; component < 3; ++component) {
                for (size_t dof_index = 0; dof_index < mNumberOfDofs; ++dof_index) {
                    const size_t dof_component =
                        dof_index % mDofsPerNode;

                    if (dof_component > 2) {
                        rTangentSecondDerivativeVariation[
                            component * mNumberOfDofs + dof_index] = 0.0;

                        continue;
                    }

                    rTangentSecondDerivativeVariation[
                        component * mNumberOfDofs + dof_index] +=
                        third_base_vector_variation[
                            component * mNumberOfDofs + dof_index] /
                            first_base_vector_norm
                        -
                        rThirdBaseVector[component] *
                            first_base_vector_dot_first_variation[dof_index] /
                            pow(first_base_vector_norm, 3);

                    rTangentSecondDerivativeVariation[
                        component * mNumberOfDofs + dof_index] +=
                        -(
                            second_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                                first_second_base_vector_product
                            +
                            rSecondBaseVector[component] *
                                first_second_product_variation[dof_index]
                        ) /
                        pow(first_base_vector_norm, 3)
                        +
                        3.0 *
                            rSecondBaseVector[component] *
                            first_second_base_vector_product *
                            first_base_vector_dot_first_variation[dof_index] /
                            pow(first_base_vector_norm, 5);

                    rTangentSecondDerivativeVariation[
                        component * mNumberOfDofs + dof_index] +=
                        -(
                            second_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                                first_second_base_vector_product
                            +
                            rSecondBaseVector[component] *
                                first_second_product_variation[dof_index]
                            +
                            first_base_vector_variation[
                                component * mNumberOfDofs + dof_index] *
                                first_second_base_vector_product_derivative
                            +
                            rFirstBaseVector[component] *
                                first_second_product_derivative_variation[dof_index]
                        ) /
                        pow(first_base_vector_norm, 3)
                        +
                        3.0 *
                            (
                                rSecondBaseVector[component] *
                                    first_second_base_vector_product
                                +
                                rFirstBaseVector[component] *
                                    first_second_base_vector_product_derivative
                            ) *
                            first_base_vector_dot_first_variation[dof_index] /
                            pow(first_base_vector_norm, 5);

                    rTangentSecondDerivativeVariation[
                        component * mNumberOfDofs + dof_index] +=
                        3.0 *
                            (
                                first_base_vector_variation[
                                    component * mNumberOfDofs + dof_index] *
                                    pow(first_second_base_vector_product, 2)
                                +
                                rFirstBaseVector[component] *
                                    2.0 *
                                    first_second_base_vector_product *
                                    first_second_product_variation[dof_index]
                            ) /
                            pow(first_base_vector_norm, 5)
                        -
                        15.0 *
                            rFirstBaseVector[component] *
                            pow(first_second_base_vector_product, 2) *
                            first_base_vector_dot_first_variation[dof_index] /
                            pow(first_base_vector_norm, 7);
                }
            }
        }
    // ---------------------------------------------------------------------
    // Beam rotation helpers
    // ---------------------------------------------------------------------
    // These helpers are local to this translation unit and are not part of
    // the BeamBaseElement3D interface.

    namespace
    {

        using Matrix3d = BoundedMatrix<double, 3, 3>;

        constexpr double tolerance = 1.0e-8;

            template <typename... Args>
            auto CrossProduct(Args&&... args)
            {
                return MathUtils<double>::CrossProduct(
                    std::forward<Args>(args)...);
            }

            Matrix3d CrossProductVectorMatrix(
                const array_1d<double, 3>& rVector,
                const Matrix3d& rMatrix)
            {
                Matrix3d matrix_vector;
                matrix_vector.clear();

                int levi_civita[3][3][3] = {};

                levi_civita[0][1][2] = 1;
                levi_civita[2][0][1] = 1;
                levi_civita[1][2][0] = 1;

                levi_civita[0][2][1] = -1;
                levi_civita[1][0][2] = -1;
                levi_civita[2][1][0] = -1;

                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t j = 0; j < 3; ++j) {
                        for (std::size_t k = 0; k < 3; ++k) {
                            for (std::size_t l = 0; l < 3; ++l) {
                                matrix_vector(i, j) +=
                                    levi_civita[j][k][l] *
                                    rMatrix(i, k) *
                                    rVector(l);
                            }
                        }
                    }
                }

                return matrix_vector;
            }

        }
      

    void BeamBaseElement3D::ComputeRodriguesMatrix(const array_1d<double, 3>& rAxis, 
            const double Phi, Matrix3d& rRodriguesMatrix)
        {
            rRodriguesMatrix.clear();

            Matrix3d identity_matrix;
            identity_matrix.resize(3, 3, false);
            identity_matrix.clear();

            for (SizeType i = 0; i < 3; ++i) {
                identity_matrix(i, i) = 1.0;
            }

            for (SizeType i = 0; i < 3; ++i) {
                rRodriguesMatrix(i, i) = cos(Phi);
            }

            rRodriguesMatrix +=
                CrossProductVectorMatrix(
                    rAxis,
                    identity_matrix) *
                sin(Phi);
        }

    void BeamBaseElement3D::ComputeRodriguesMatrixDerivative(
            const array_1d<double, 3>& rAxis,
            const array_1d<double, 3>& rAxisDerivative,
            const double Phi,
            const double PhiDerivative,
            Matrix3d& rRodriguesMatrixDerivative)
        {
            rRodriguesMatrixDerivative.clear();

            Matrix3d identity_matrix;
            identity_matrix.resize(3, 3, false);
            identity_matrix.clear();

            for (SizeType i = 0; i < 3; ++i) {
                identity_matrix(i, i) = 1.0;
            }

            for (SizeType i = 0; i < 3; ++i) {
                rRodriguesMatrixDerivative(i, i) =
                    -PhiDerivative * sin(Phi);
            }

            rRodriguesMatrixDerivative +=
                CrossProductVectorMatrix(
                    rAxis,
                    identity_matrix) *
                cos(Phi) *
                PhiDerivative;

            rRodriguesMatrixDerivative +=
                CrossProductVectorMatrix(
                    rAxisDerivative,
                    identity_matrix) *
                sin(Phi);
        }

    void BeamBaseElement3D::ComputeRodriguesMatrixVariation(
            const array_1d<double, 3>& rAxis,
            const Vector& rAxisVariation,
            const Vector& rShapeFunctions,
            const double Phi,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rRodriguesMatrixVariation)
        {
            rRodriguesMatrixVariation.clear();

            int levi_civita[3][3][3] = {};

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            for (size_t row_component = 0; row_component < 3; ++row_component) {
                for (size_t column_component = 0; column_component < 3; ++column_component) {
                    for (size_t dof_index = 0; dof_index < NumberOfDofs; ++dof_index) {
                        const size_t dof_component =
                            dof_index % DofsPerNode;

                        const size_t control_point_index =
                            dof_index / DofsPerNode;

                        if (row_component == column_component) {
                            if (dof_component > 2) {
                                rRodriguesMatrixVariation(
                                    row_component * NumberOfDofs + dof_index,
                                    column_component) +=
                                    -sin(Phi) *
                                    rShapeFunctions[control_point_index];
                            }
                        } else if (dof_component > 2) {
                            for (SizeType k = 0; k < 3; ++k) {
                                rRodriguesMatrixVariation(
                                    row_component * NumberOfDofs + dof_index,
                                    column_component) +=
                                    cos(Phi) *
                                    rShapeFunctions[control_point_index] *
                                    levi_civita[row_component][k][column_component] *
                                    rAxis[k];
                            }
                        }

                        for (SizeType k = 0; k < 3; ++k) {
                            rRodriguesMatrixVariation(
                                row_component * NumberOfDofs + dof_index,
                                column_component) +=
                                sin(Phi) *
                                levi_civita[row_component][k][column_component] *
                                rAxisVariation[
                                    k * NumberOfDofs + dof_index];
                        }
                    }
                }
            }
        }


    void BeamBaseElement3D::ComputeRodriguesMatrixSecondVariation(
            const array_1d<double, 3>& rAxis,
            const Vector& rAxisVariation,
            const Matrix& rAxisSecondVariation,
            const Vector& rShapeFunctions,
            const double Phi,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rRodriguesMatrixSecondVariation)
        {
            rRodriguesMatrixSecondVariation.clear();

            int levi_civita[3][3][3] = {};

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector phi_variation =
                ZeroVector(NumberOfDofs);

            for (size_t dof_index = 0; dof_index < NumberOfDofs; ++dof_index) {
                const size_t dof_component =
                    dof_index % DofsPerNode;

                const size_t control_point_index =
                    dof_index / DofsPerNode;

                if (dof_component > 2) {
                    phi_variation(dof_index) =
                        rShapeFunctions[control_point_index];
                }
            }

            for (size_t row_component = 0; row_component < 3; ++row_component) {
                for (size_t column_component = 0; column_component < 3; ++column_component) {
                    for (size_t second_dof_index = 0;
                        second_dof_index < NumberOfDofs;
                        ++second_dof_index)
                    {
                        const size_t second_dof_component =
                            second_dof_index % DofsPerNode;

                        for (size_t dof_index = 0;
                            dof_index < NumberOfDofs;
                            ++dof_index)
                        {
                            const size_t dof_component =
                                dof_index % DofsPerNode;

                            if (
                                row_component == column_component &&
                                (dof_component > 2 ||
                                second_dof_component > 2))
                            {
                                rRodriguesMatrixSecondVariation(
                                    row_component * NumberOfDofs + dof_index,
                                    column_component * NumberOfDofs +
                                        second_dof_index) +=
                                    -cos(Phi) *
                                    phi_variation[dof_index] *
                                    phi_variation[second_dof_index];
                            }

                            for (SizeType k = 0; k < 3; ++k) {
                                if (
                                    dof_component > 2 ||
                                    second_dof_component > 2)
                                {
                                    rRodriguesMatrixSecondVariation(
                                        row_component * NumberOfDofs + dof_index,
                                        column_component * NumberOfDofs +
                                            second_dof_index) +=
                                        -sin(Phi) *
                                            phi_variation[dof_index] *
                                            phi_variation[second_dof_index] *
                                            levi_civita[row_component][k][column_component] *
                                            rAxis[k]
                                        +
                                        phi_variation[dof_index] *
                                            cos(Phi) *
                                            levi_civita[row_component][k][column_component] *
                                            rAxisVariation[
                                                k * NumberOfDofs +
                                                second_dof_index]
                                        +
                                        phi_variation[second_dof_index] *
                                            cos(Phi) *
                                            levi_civita[row_component][k][column_component] *
                                            rAxisVariation[
                                                k * NumberOfDofs +
                                                dof_index];
                                } else {
                                    rRodriguesMatrixSecondVariation(
                                        row_component * NumberOfDofs + dof_index,
                                        column_component * NumberOfDofs +
                                            second_dof_index) +=
                                        sin(Phi) *
                                        levi_civita[row_component][k][column_component] *
                                        rAxisSecondVariation(
                                            k * NumberOfDofs + dof_index,
                                            second_dof_index);
                                }
                            }
                        }
                    }
                }
            }
        }


    void BeamBaseElement3D::ComputeRodriguesMatrixDerivativeVariation(
            const array_1d<double, 3>& rAxis,
            const Vector& rAxisVariation,
            const array_1d<double, 3>& rAxisDerivative,
            const Vector& rAxisDerivativeVariation,
            const Vector& rShapeFunctions,
            const Vector& rShapeFunctionDerivatives,
            const double Phi,
            const double PhiDerivative,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rRodriguesMatrixDerivativeVariation)
        {
            rRodriguesMatrixDerivativeVariation.clear();

            int levi_civita[3][3][3] = {};

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector phi_variation =
                ZeroVector(NumberOfDofs);

            for (
                size_t control_point_index = 0;
                control_point_index < NumberOfDofs / DofsPerNode;
                ++control_point_index)
            {
                phi_variation(
                    control_point_index * DofsPerNode + 3) =
                    rShapeFunctions(control_point_index);
            }

            Vector phi_derivative_variation =
                ZeroVector(NumberOfDofs);

            for (
                size_t control_point_index = 0;
                control_point_index < NumberOfDofs / DofsPerNode;
                ++control_point_index)
            {
                phi_derivative_variation(
                    control_point_index * DofsPerNode + 3) =
                    rShapeFunctionDerivatives(control_point_index);
            }

            for (size_t row_component = 0; row_component < 3; ++row_component) {
                for (size_t column_component = 0; column_component < 3; ++column_component) {
                    for (size_t dof_index = 0; dof_index < NumberOfDofs; ++dof_index) {
                        const size_t dof_component =
                            dof_index % DofsPerNode;

                        const size_t control_point_index =
                            dof_index / DofsPerNode;

                        if (
                            row_component == column_component &&
                            dof_component > 2)
                        {
                            rRodriguesMatrixDerivativeVariation(
                                row_component * NumberOfDofs + dof_index,
                                column_component) +=
                                -phi_derivative_variation(dof_index) *
                                    sin(Phi)
                                -
                                cos(Phi) *
                                    PhiDerivative *
                                    rShapeFunctions[control_point_index];
                        }

                        for (SizeType k = 0; k < 3; ++k) {
                            if (dof_component > 2) {
                                rRodriguesMatrixDerivativeVariation(
                                    row_component * NumberOfDofs + dof_index,
                                    column_component) +=
                                    (
                                        phi_derivative_variation(dof_index) *
                                            cos(Phi)
                                        -
                                        PhiDerivative *
                                            rShapeFunctions[control_point_index] *
                                            sin(Phi)
                                    ) *
                                    levi_civita[row_component][k][column_component] *
                                    rAxis[k]
                                    +
                                    cos(Phi) *
                                    phi_variation(dof_index) *
                                    levi_civita[row_component][k][column_component] *
                                    rAxisDerivative[k];
                            } else {
                                rRodriguesMatrixDerivativeVariation(
                                    row_component * NumberOfDofs + dof_index,
                                    column_component) +=
                                    cos(Phi) *
                                        PhiDerivative *
                                        levi_civita[row_component][k][column_component] *
                                        rAxisVariation[
                                            k * NumberOfDofs + dof_index]
                                    +
                                    sin(Phi) *
                                        levi_civita[row_component][k][column_component] *
                                        rAxisDerivativeVariation[
                                            k * NumberOfDofs + dof_index];
                            }
                        }
                    }
                }
            }
        }


    void BeamBaseElement3D::ComputeRodriguesMatrixDerivativeSecondVariation(
            const array_1d<double, 3>& rAxis,
            const Vector& rAxisVariation,
            const array_1d<double, 3>& rAxisDerivative,
            const Vector& rAxisDerivativeVariation,
            const Matrix& rAxisSecondVariation,
            const Matrix& rAxisDerivativeSecondVariation,
            const Vector& rShapeFunctions,
            const Vector& rShapeFunctionDerivatives,
            const double Phi,
            const double PhiDerivative,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rRodriguesMatrixDerivativeSecondVariation)
        {
                rRodriguesMatrixDerivativeSecondVariation.clear();

                int levi_civita[3][3][3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            levi_civita[i][j][k] = 0;
                        }
                    }
                }

                levi_civita[0][1][2] = 1;
                levi_civita[2][0][1] = 1;
                levi_civita[1][2][0] = 1;

                levi_civita[0][2][1] = -1;
                levi_civita[1][0][2] = -1;
                levi_civita[2][1][0] = -1;


                Vector phi_variation;
                phi_variation.resize(NumberOfDofs);
                phi_variation.clear();
                Vector phi_derivative_variation;
                phi_derivative_variation.resize(NumberOfDofs);
                phi_derivative_variation.clear();

                for (size_t  r  = 0;r < NumberOfDofs / DofsPerNode;r++)
                {
                    phi_variation(r * DofsPerNode + 3) = rShapeFunctions[r];
                    phi_derivative_variation(r * DofsPerNode + 3) = rShapeFunctionDerivatives(r);
                }

                double cos_phi;
                cos_phi = cos(Phi);
                double sin_phi;
                sin_phi = sin(Phi);

                for (size_t t = 0;t < 3;t++) 
                {
                    for (size_t u = 0;u < 3;u++)
                    {
                        for (size_t  r  = 0;r < NumberOfDofs;r++)
                        {
                            
                            for (size_t  s  = 0;s < NumberOfDofs;s++)
                            {
                                
                                
                                if (t == u)
                                {
                                    rRodriguesMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -phi_derivative_variation(r) * phi_variation(s) * cos_phi - phi_derivative_variation(s) * phi_variation(r) * cos_phi + sin_phi * PhiDerivative * phi_variation[r] * phi_variation[s];
                                }
                                
                                {
                                    for (int k = 0; k < 3; k++)
                                    {
                                        rRodriguesMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -phi_derivative_variation(r) * phi_variation(s) * sin_phi * levi_civita[t][k][u] * rAxis[k]
                                            - phi_derivative_variation(s) * phi_variation(r) * sin_phi * levi_civita[t][k][u] * rAxis[k]
                                                + phi_derivative_variation(r) * cos_phi * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + s]
                                                    + phi_derivative_variation(s) * cos_phi * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r]
                                                    - cos_phi * PhiDerivative * phi_variation[r] * phi_variation[s] * levi_civita[t][k][u] * rAxis[k]
                                                    - sin_phi * PhiDerivative * phi_variation[r] * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + s]
                                                        - sin_phi * PhiDerivative * phi_variation[s] * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r]
                                                        + cos_phi * PhiDerivative * levi_civita[t][k][u] * rAxisSecondVariation(k * NumberOfDofs + r, s)
                                                        - phi_variation[r] * phi_variation[s] * sin_phi * levi_civita[t][k][u] * rAxisDerivative[k]
                                                        + phi_variation[r] * cos_phi * levi_civita[t][k][u] * rAxisDerivativeVariation[k * NumberOfDofs + s]
                                                            + phi_variation[s] * cos_phi * levi_civita[t][k][u] * rAxisDerivativeVariation[k * NumberOfDofs + r]
                                                        ;

                                                        
                                                        rRodriguesMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += sin_phi * levi_civita[t][k][u] * rAxisDerivativeSecondVariation(k * NumberOfDofs + r, s); 

                                    }
                                }
                            }
                        }
                    }
                }
            }

        void BeamBaseElement3D::ComputeRodriguesMatrixSecondDerivativeVariation(
            const array_1d<double, 3>& rAxis,
            const Vector& rAxisVariation,
            const array_1d<double, 3>& rAxisDerivative,
            const Vector& rAxisDerivativeVariation,
            const array_1d<double, 3>& rAxisSecondDerivative,
            const Vector& rAxisSecondDerivativeVariation,
            const Vector& rShapeFunctions,
            const Vector& rShapeFunctionDerivatives,
            const Vector& rShapeFunctionSecondDerivatives,
            const double Phi,
            const double PhiDerivative,
            const double PhiSecondDerivative,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rRodriguesMatrixSecondDerivativeVariation)
        {
                
                rRodriguesMatrixSecondDerivativeVariation.resize(3 * NumberOfDofs, 3);
                rRodriguesMatrixSecondDerivativeVariation.clear();

                int levi_civita[3][3][3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            levi_civita[i][j][k] = 0;
                        }
                    }
                }

                levi_civita[0][1][2] = 1;
                levi_civita[2][0][1] = 1;
                levi_civita[1][2][0] = 1;

                levi_civita[0][2][1] = -1;
                levi_civita[1][0][2] = -1;
                levi_civita[2][1][0] = -1;

                for (size_t t = 0; t < 3; t++) 
                {
                    for (size_t u = 0; u < 3; u++)
                    {
                        for (size_t  r  = 0; r < NumberOfDofs; r++)
                        {
                            size_t xyz = r % DofsPerNode; 
                            size_t i= r / DofsPerNode;     
                            if (t == u)
                            {
                                if (xyz > 2)
                                    rRodriguesMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += -rShapeFunctionSecondDerivatives[i] * sin(Phi) - cos(Phi) * PhiSecondDerivative * rShapeFunctions[i] - 2 * rShapeFunctionDerivatives[i] * PhiDerivative * cos(Phi) + sin(Phi) * pow(PhiDerivative, 2) * rShapeFunctions[i];
                                else
                                    rRodriguesMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += 0;
                            }

                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    if (xyz > 2)
                                        rRodriguesMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += (rShapeFunctionSecondDerivatives[i] * cos(Phi) - PhiSecondDerivative * rShapeFunctions[i] * sin(Phi) - 2 * rShapeFunctionDerivatives[i] * PhiDerivative * sin(Phi) - pow(PhiDerivative, 2) * rShapeFunctions[i] * cos(Phi)) * levi_civita[t][k][u] * rAxis[k] + 2 * (cos(Phi) * rShapeFunctionDerivatives[i] - sin(Phi) * rShapeFunctions[i] * PhiDerivative) * levi_civita[t][k][u] * rAxisDerivative[k] + cos(Phi) * rShapeFunctions[i] * levi_civita[t][k][u] * rAxisSecondDerivative[k];
                                    else
                                        rRodriguesMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += (cos(Phi) * PhiSecondDerivative - pow(PhiDerivative, 2) * sin(Phi)) * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r] + 2 * PhiDerivative * cos(Phi) * levi_civita[t][k][u] * rAxisDerivativeVariation[k * NumberOfDofs + r] + sin(Phi) * levi_civita[t][k][u] * rAxisSecondDerivativeVariation[k * NumberOfDofs + r]; 
                                }
                            }
                        }
                    }
                }

            }

    void BeamBaseElement3D::ComputeRodriguesMatrixVariations(
        const array_1d<double, 3>& rAxis,
        const array_1d<double, 3>& rAxisVariation,
        const array_1d<double, 3>& rAxisDerivative,
        const Vector& rAxisDerivativeVariation,
        const Matrix& rAxisSecondVariation,
        const Matrix& rAxisDerivativeSecondVariation,
        const Vector& rShapeFunctions,
        const Vector& rShapeFunctionDerivatives,
        const double Phi,
        const double PhiDerivative,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixVariation,
        Matrix& rRodriguesMatrixDerivativeVariation,
        Matrix& rRodriguesMatrixSecondVariation,
        Matrix& rRodriguesMatrixDerivativeSecondVariation)
    {
    
            rRodriguesMatrixVariation.clear();
            
            rRodriguesMatrixDerivativeVariation.clear();
            
            rRodriguesMatrixSecondVariation.clear();
            
            rRodriguesMatrixDerivativeSecondVariation.clear();

            int levi_civita[3][3][3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        levi_civita[i][j][k] = 0;
                    }
                }
            }

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector phi_variation;
            phi_variation.resize(NumberOfDofs);
            phi_variation.clear();
            Vector phi_derivative_variation;
            phi_derivative_variation.resize(NumberOfDofs);
            phi_derivative_variation.clear();

            for (size_t  r  = 0;r < NumberOfDofs / DofsPerNode;r++)
            {
                phi_variation(r * DofsPerNode + 3) = rShapeFunctions[r];
                phi_derivative_variation(r * DofsPerNode + 3) = rShapeFunctionDerivatives(r);
            }

            double cos_phi;
            cos_phi = cos(Phi);
            double sin_phi;
            sin_phi = sin(Phi);

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz_r = r % DofsPerNode; 
                        size_t i= r / DofsPerNode;     
                        if (t == u)
                        {
                            if (xyz_r > 2)
                            {
                                rRodriguesMatrixVariation(t * NumberOfDofs + r, u) += -sin(Phi) * rShapeFunctions[i];
                                rRodriguesMatrixDerivativeVariation(t * NumberOfDofs + r, u) += -phi_derivative_variation(r) * sin(Phi) - cos(Phi) * PhiDerivative * rShapeFunctions[i];
                            }
                        }
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz_r > 2)
                            {
                                rRodriguesMatrixVariation(t * NumberOfDofs + r, u) += cos(Phi) * rShapeFunctions[i] * levi_civita[t][k][u] * rAxis[k];
                                rRodriguesMatrixDerivativeVariation(t * NumberOfDofs + r, u) += (phi_derivative_variation(r) * cos(Phi) - PhiDerivative * rShapeFunctions[i] * sin(Phi)) * levi_civita[t][k][u] * rAxis[k] + cos(Phi) * phi_variation(r) * levi_civita[t][k][u] * rAxisDerivative[k];
                            }
                            else
                            {
                                rRodriguesMatrixVariation(t * NumberOfDofs + r, u) += sin(Phi) * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r];
                                rRodriguesMatrixDerivativeVariation(t * NumberOfDofs + r, u) += cos(Phi) * PhiDerivative * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r] + sin(Phi) * levi_civita[t][k][u] * rAxisDerivativeVariation[k * NumberOfDofs + r]; 
                            }
                        }
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyz_s = s % DofsPerNode; 
                            
                            if (t == u)
                            {
                                rRodriguesMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -cos(Phi) * phi_variation[r] * phi_variation[s];
                                rRodriguesMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -phi_derivative_variation(r) * phi_variation(s) * cos_phi - phi_derivative_variation(s) * phi_variation(r) * cos_phi + sin_phi * PhiDerivative * phi_variation[r] * phi_variation[s];
                            }
                            
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    if (xyz_r > 2 || xyz_s > 2)
                                        rRodriguesMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -sin(Phi) * phi_variation[r] * phi_variation[s] * levi_civita[t][k][u] * rAxis[k] + phi_variation[r] * cos(Phi) * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + s] + phi_variation[s] * cos(Phi) * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r];
                                    else
                                        rRodriguesMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += sin(Phi) * levi_civita[t][k][u] * rAxisSecondVariation(k * NumberOfDofs + r, s); 

                                    rRodriguesMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -phi_derivative_variation(r) * phi_variation(s) * sin_phi * levi_civita[t][k][u] * rAxis[k]
                                        - phi_derivative_variation(s) * phi_variation(r) * sin_phi * levi_civita[t][k][u] * rAxis[k]
                                            + phi_derivative_variation(r) * cos_phi * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + s]
                                                + phi_derivative_variation(s) * cos_phi * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r]
                                                - cos_phi * PhiDerivative * phi_variation[r] * phi_variation[s] * levi_civita[t][k][u] * rAxis[k]
                                                - sin_phi * PhiDerivative * phi_variation[r] * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + s]
                                                    - sin_phi * PhiDerivative * phi_variation[s] * levi_civita[t][k][u] * rAxisVariation[k * NumberOfDofs + r]
                                                    + cos_phi * PhiDerivative * levi_civita[t][k][u] * rAxisSecondVariation(k * NumberOfDofs + r, s)
                                                    - phi_variation[r] * phi_variation[s] * sin_phi * levi_civita[t][k][u] * rAxisDerivative[k]
                                                    + phi_variation[r] * cos_phi * levi_civita[t][k][u] * rAxisDerivativeVariation[k * NumberOfDofs + s]
                                                        + phi_variation[s] * cos_phi * levi_civita[t][k][u] * rAxisDerivativeVariation[k * NumberOfDofs + r]
                                                    ;

                                                    
                                                    rRodriguesMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += sin_phi * levi_civita[t][k][u] * rAxisDerivativeSecondVariation(k * NumberOfDofs + r, s); 
                                }
                            }
                        }
                    }
                }
            }
        }



    void BeamBaseElement3D::ComputeRodriguesMatrixSecondDerivative(
            const array_1d<double, 3>& rAxis,
            const array_1d<double, 3>& rAxisDerivative,
            const array_1d<double, 3>& rAxisSecondDerivative,
            const double Phi,
            const double PhiDerivative,
            const double PhiSecondDerivative,
            Matrix3d& rRodriguesMatrixSecondDerivative)
        {
            rRodriguesMatrixSecondDerivative.clear();

            Matrix3d identity_matrix;
            identity_matrix.clear();  
            for (int i = 0; i < 3; i++) { identity_matrix(i, i) = 1; }

            for (int i = 0; i < 3; i++) { rRodriguesMatrixSecondDerivative(i, i) = -PhiSecondDerivative * sin(Phi) - pow(PhiDerivative, 2) * cos(Phi); }
            rRodriguesMatrixSecondDerivative += CrossProductVectorMatrix(rAxis, identity_matrix) * (cos(Phi) * PhiSecondDerivative - pow(PhiDerivative, 2) * sin(Phi));
            rRodriguesMatrixSecondDerivative += CrossProductVectorMatrix(rAxisDerivative, identity_matrix) * 2 * PhiDerivative * cos(Phi);
            rRodriguesMatrixSecondDerivative += CrossProductVectorMatrix(rAxisSecondDerivative, identity_matrix) * sin(Phi);
        }

    void BeamBaseElement3D::ComputeLambdaMatrix(const array_1d<double, 3>& rReferenceTangent, const array_1d<double, 3>& rCurrentTangent, Matrix3d& rLambdaMatrix)
        {
            rLambdaMatrix.clear();  
            Matrix3d lambda_matrix_auxiliary;
            lambda_matrix_auxiliary.clear();  
            double scaling_factor;

            Matrix3d identity_matrix;
            identity_matrix.resize(3, 3, false);
            identity_matrix.clear();  
            for (int i = 0; i < 3; i++) { identity_matrix(i, i) = 1.; }
            array_1d<double, 3> reference_cross_current;
            CrossProduct(reference_cross_current, rReferenceTangent, rCurrentTangent);
            double cross_product_norm = norm_2(reference_cross_current);
            array_1d<double, 3>& rotation_axis = reference_cross_current;
            if (cross_product_norm > 0.000000000001) rotation_axis = rotation_axis / cross_product_norm;

            if ((inner_prod(rReferenceTangent, rCurrentTangent) + 1.0)> tolerance)
            {
                for (int i = 0; i < 3; i++) { rLambdaMatrix(i, i) = inner_prod(rReferenceTangent, rCurrentTangent); }
                rLambdaMatrix += CrossProductVectorMatrix(CrossProduct(rReferenceTangent, rCurrentTangent), identity_matrix);

                array_1d<double, 3> cross_product = CrossProduct(rReferenceTangent, rCurrentTangent);
                lambda_matrix_auxiliary = outer_prod(cross_product, cross_product);
                scaling_factor = 1.0 / (1.0 + inner_prod(rReferenceTangent, rCurrentTangent));
                lambda_matrix_auxiliary = lambda_matrix_auxiliary * scaling_factor;
                rLambdaMatrix += lambda_matrix_auxiliary;
            }
            else
            {
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        rLambdaMatrix(i, j) = (rotation_axis(i) * rotation_axis(j)) * (1 - inner_prod(rReferenceTangent, rCurrentTangent)) + inner_prod(rReferenceTangent, rCurrentTangent) * identity_matrix(i, j);
                    }
                }
                rLambdaMatrix += CrossProductVectorMatrix(reference_cross_current, identity_matrix);
            }

        }


    void BeamBaseElement3D::ComputeLambdaMatrixDerivative(
            const array_1d<double, 3>& rReferenceTangent,
            const array_1d<double, 3>& rCurrentTangent,
            const array_1d<double, 3>& rReferenceTangentDerivative,
            const array_1d<double, 3>& rCurrentTangentDerivative,
            Matrix3d& rLambdaMatrixDerivative)
        {
            rLambdaMatrixDerivative.clear();  

            double scaling_factor;

            Matrix3d identity_matrix;
            identity_matrix.clear();  
            for (int i = 0; i < 3; i++) { identity_matrix(i, i) = 1.; }

            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
            double reference_dot_current_derivative = inner_prod(rReferenceTangent, rCurrentTangentDerivative);
            double reference_derivative_dot_current = inner_prod(rReferenceTangentDerivative, rCurrentTangent);
            array_1d<double, 3> reference_cross_current = CrossProduct(rReferenceTangent, rCurrentTangent);
            array_1d<double, 3> reference_cross_current_derivative = CrossProduct(rReferenceTangent, rCurrentTangentDerivative);
            array_1d<double, 3> reference_derivative_cross_current = CrossProduct(rReferenceTangentDerivative, rCurrentTangent);
            array_1d<double, 3> tangent_cross_product_derivative = reference_cross_current_derivative + reference_derivative_cross_current;

            for (int i = 0; i < 3; i++) { rLambdaMatrixDerivative(i, i) = reference_dot_current_derivative + reference_derivative_dot_current; }
            rLambdaMatrixDerivative += CrossProductVectorMatrix(tangent_cross_product_derivative, identity_matrix);
            scaling_factor = -(reference_dot_current_derivative + reference_derivative_dot_current) / pow((1.0 + tangent_dot_product), 2);
            rLambdaMatrixDerivative += outer_prod(reference_cross_current, reference_cross_current) * scaling_factor;
            scaling_factor = 1.0 / (1.0 + tangent_dot_product);
            rLambdaMatrixDerivative += outer_prod(tangent_cross_product_derivative, reference_cross_current) * scaling_factor;
            rLambdaMatrixDerivative += outer_prod(reference_cross_current, tangent_cross_product_derivative) * scaling_factor;

        }


    void BeamBaseElement3D::ComputeLambdaMatrixSecondDerivative(
            const array_1d<double, 3>& rReferenceTangent,
            const array_1d<double, 3>& rCurrentTangent,
            const array_1d<double, 3>& rReferenceTangentDerivative,
            const array_1d<double, 3>& rCurrentTangentDerivative,
            const array_1d<double, 3>& rReferenceTangentSecondDerivative,
            const array_1d<double, 3>& rCurrentTangentSecondDerivative,
            Matrix3d& rLambdaMatrixSecondDerivative)
        {
            rLambdaMatrixSecondDerivative.clear();  

            double scaling_factor;

            Matrix3d identity_matrix;
            identity_matrix.clear();  
            for (int i = 0; i < 3; i++) { identity_matrix(i, i) = 1.; }

            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
            double reference_dot_current_derivative = inner_prod(rReferenceTangent, rCurrentTangentDerivative);
            double reference_derivative_dot_current = inner_prod(rReferenceTangentDerivative, rCurrentTangent);
            double tangent_derivative_dot_product = inner_prod(rReferenceTangentDerivative, rCurrentTangentDerivative);
            double reference_second_derivative_dot_current = inner_prod(rReferenceTangentSecondDerivative, rCurrentTangent);
            double reference_dot_current_second_derivative = inner_prod(rReferenceTangent, rCurrentTangentSecondDerivative);
            double tangent_dot_product_derivative = reference_derivative_dot_current + reference_dot_current_derivative;
            double tangent_dot_product_second_derivative = reference_second_derivative_dot_current + 2 * tangent_derivative_dot_product + reference_dot_current_second_derivative;
            array_1d<double, 3> reference_cross_current = CrossProduct(rReferenceTangent, rCurrentTangent);
            array_1d<double, 3> reference_cross_current_derivative = CrossProduct(rReferenceTangent, rCurrentTangentDerivative);
            array_1d<double, 3> reference_derivative_cross_current = CrossProduct(rReferenceTangentDerivative, rCurrentTangent);
            array_1d<double, 3> reference_derivative_cross_current_derivative = CrossProduct(rReferenceTangentDerivative, rCurrentTangentDerivative);
            array_1d<double, 3> reference_second_derivative_cross_current = CrossProduct(rReferenceTangentSecondDerivative, rCurrentTangent);
            array_1d<double, 3> reference_cross_current_second_derivative = CrossProduct(rReferenceTangent, rCurrentTangentSecondDerivative);
            array_1d<double, 3> tangent_cross_product_derivative = reference_cross_current_derivative + reference_derivative_cross_current;
            array_1d<double, 3> tangent_cross_product_second_derivative = reference_second_derivative_cross_current + 2 * reference_derivative_cross_current_derivative + reference_cross_current_second_derivative;

            for (int i = 0; i < 3; i++) { rLambdaMatrixSecondDerivative(i, i) = tangent_dot_product_second_derivative; }

            rLambdaMatrixSecondDerivative += CrossProductVectorMatrix(tangent_cross_product_second_derivative, identity_matrix);
            scaling_factor = 2 * pow(tangent_dot_product_derivative, 2) / pow((1.0 + tangent_dot_product), 3) - tangent_dot_product_second_derivative / pow((1.0 + tangent_dot_product), 2);
            rLambdaMatrixSecondDerivative += outer_prod(reference_cross_current, reference_cross_current) * scaling_factor;
            scaling_factor = -tangent_dot_product_derivative / pow((1.0 + tangent_dot_product), 2) * 2;
            rLambdaMatrixSecondDerivative += outer_prod(tangent_cross_product_derivative, reference_cross_current) * scaling_factor;
            rLambdaMatrixSecondDerivative += outer_prod(reference_cross_current, tangent_cross_product_derivative) * scaling_factor;
            scaling_factor = 1.0 / (1.0 + tangent_dot_product);
            rLambdaMatrixSecondDerivative += outer_prod(tangent_cross_product_second_derivative, reference_cross_current) * scaling_factor;
            rLambdaMatrixSecondDerivative += outer_prod(tangent_cross_product_derivative, tangent_cross_product_derivative) * 2 * scaling_factor;
            rLambdaMatrixSecondDerivative += outer_prod(reference_cross_current, tangent_cross_product_second_derivative) * scaling_factor;

        }


    void BeamBaseElement3D::ComputeLambdaMatrixVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const Vector& rCurrentTangentVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixVariation)
        {
            

            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
            
            rLambdaMatrixVariation.clear();

            Matrix3d identity_matrix;
            identity_matrix.clear();  
            for (int i = 0; i < 3; i++) { identity_matrix(i, i) = 1.; }

            int levi_civita[3][3][3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        levi_civita[i][j][k] = 0;
                    }
                }
            }

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector reference_cross_current_variation;
            Vector reference_cross_current;

            reference_cross_current_variation.resize(NumberOfDofs * 3);
            reference_cross_current.resize(3);

            reference_cross_current_variation.clear();
            reference_cross_current.clear();

            reference_cross_current = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangent);

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++)
                {
                    size_t xyz = r % DofsPerNode; 
                    for (size_t u = 0;u < 3;u++)
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            if (xyz > 2)
                                reference_cross_current_variation[t * NumberOfDofs + r] += 0;
                            else
                                reference_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                        }
                    }
                }
            }

            Vector tangent_dot_product_variation;
            tangent_dot_product_variation.resize(NumberOfDofs);
            tangent_dot_product_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyz = r % DofsPerNode; 
                    

                    if (xyz > 2)
                        tangent_dot_product_variation(r) = 0;
                    else
                        tangent_dot_product_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz = r % DofsPerNode; 
                        
                        if (t == u)
                        {
                            if (xyz > 2)
                                rLambdaMatrixVariation(t * NumberOfDofs + r, u) += 0;
                            else
                                rLambdaMatrixVariation(t * NumberOfDofs + r, u) += tangent_dot_product_variation(r);
                        }
                        for (int k = 0; k < 3; k++)
                        {
                            rLambdaMatrixVariation(t * NumberOfDofs + r, u) += levi_civita[t][k][u] * reference_cross_current_variation[r + k * NumberOfDofs]; 
                        }
                        rLambdaMatrixVariation(t * NumberOfDofs + r, u) += -tangent_dot_product_variation[r] / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current[t] * reference_cross_current[u]);
                        rLambdaMatrixVariation(t * NumberOfDofs + r, u) += +1.0 / (1.0 + tangent_dot_product) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + r]);
                    }
                }
            }
        }

    void BeamBaseElement3D::ComputeLambdaMatrixSecondVariation(
            const array_1d<double, 3>& rReferenceTangent,
            const array_1d<double, 3>& rCurrentTangent,
            const Vector& rCurrentTangentVariation,
            const Matrix& rCurrentTangentSecondVariation,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rLambdaMatrixSecondVariation)
        {   
            
            
            rLambdaMatrixSecondVariation.clear();  

            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);

            int levi_civita[3][3][3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        levi_civita[i][j][k] = 0;
                    }
                }
            }

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector reference_cross_current_variation;
            Matrix reference_cross_current_second_variation;
            Vector reference_cross_current;

            reference_cross_current_variation.resize(NumberOfDofs * 3);
            reference_cross_current_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_cross_current.resize(3);

            reference_cross_current_variation.clear();
            reference_cross_current_second_variation.clear();
            reference_cross_current.clear();

            reference_cross_current = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangent);

            Vector tangent_dot_product_variation;
            tangent_dot_product_variation.resize(NumberOfDofs);
            tangent_dot_product_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyz = r % DofsPerNode; 

                    if (xyz > 2)
                        tangent_dot_product_variation(r) = 0;
                    else
                        tangent_dot_product_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                }
            }

            Matrix tangent_dot_product_second_variation;
            tangent_dot_product_second_variation.resize(NumberOfDofs, NumberOfDofs);
            tangent_dot_product_second_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    for (size_t  s  = 0;s < NumberOfDofs;s++) 
                    {
                        size_t xyzr = r % DofsPerNode; 
                        size_t xyzs = s % DofsPerNode; 

                        if (xyzr > 2 || xyzs > 2)
                            tangent_dot_product_second_variation(r, s) += 0;
                        else
                            tangent_dot_product_second_variation(r, s) += rCurrentTangentSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangent[t];
                    }
                }
            }



            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++)
                {
                    size_t xyz_r = r % DofsPerNode; 
                    for (size_t u = 0;u < 3;u++)
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            if (xyz_r > 2)
                                reference_cross_current_variation[t * NumberOfDofs + r] += 0;
                            else
                                reference_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++) 
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz_r = r % DofsPerNode; 
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyz_s = s % DofsPerNode; 
                            for (size_t k = 0;k < 3;k++)
                            {
                                if (xyz_r > 2 || xyz_s > 2)
                                    reference_cross_current_second_variation(t * NumberOfDofs + r, s) += 0;
                                else
                                    reference_cross_current_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentSecondVariation(u * NumberOfDofs + r, s);
                            }
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyzr = r % DofsPerNode; 
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyzs = s % DofsPerNode; 
                            if (xyzr > 2 || xyzs > 2) rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) = 0;
                            else
                            {
                                if (t == u)
                                    rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += tangent_dot_product_second_variation(r, s);
                                else
                                {
                                    for (int k = 0; k < 3; k++)
                                        rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += levi_civita[t][k][u] * reference_cross_current_second_variation(k * NumberOfDofs + r, s); 
                                }
                                rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (2 * tangent_dot_product_variation(r) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 3) - tangent_dot_product_second_variation(r, s) / pow(1.0 + tangent_dot_product, 2)) * reference_cross_current[t] * reference_cross_current[u];
                                rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -tangent_dot_product_variation(r) / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current_variation[t * NumberOfDofs + s] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + s])
                                    - tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + r]);
                                rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += 1.0 / (1.0 + tangent_dot_product) * (reference_cross_current_second_variation(t * NumberOfDofs + r, s) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + r) * reference_cross_current_variation(u * NumberOfDofs + s) +
                                    reference_cross_current_variation(t * NumberOfDofs + s) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current[t] * reference_cross_current_second_variation(u * NumberOfDofs + s, r));
                            }
                        }
                    }
                }
            }
        }

    void BeamBaseElement3D::ComputeLambdaMatrixDerivativeVariation(
            const array_1d<double, 3>& rReferenceTangent,
            const array_1d<double, 3>& rCurrentTangent,
            const array_1d<double, 3>& rReferenceTangentDerivative,
            const Vector& rCurrentTangentVariation,
            const array_1d<double, 3>& rCurrentTangentDerivative,
            const Vector& rCurrentTangentDerivativeVariation,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rLambdaMatrixDerivativeVariation)
        {
            
            rLambdaMatrixDerivativeVariation.clear();  


            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
            double reference_dot_current_derivative = inner_prod(rReferenceTangent, rCurrentTangentDerivative);
            double reference_derivative_dot_current = inner_prod(rReferenceTangentDerivative, rCurrentTangent);

            int levi_civita[3][3][3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        levi_civita[i][j][k] = 0;
                    }
                }
            }

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector reference_cross_current_variation;
            Vector reference_cross_current_derivative_variation;
            Vector reference_derivative_cross_current_variation;
            Vector reference_cross_current;
            Vector tangent_cross_product_derivative;
            Vector reference_derivative_cross_current;

            reference_cross_current_variation.resize(NumberOfDofs * 3);
            reference_cross_current_derivative_variation.resize(NumberOfDofs * 3);
            reference_derivative_cross_current_variation.resize(NumberOfDofs * 3);
            reference_cross_current.resize(3);
            tangent_cross_product_derivative.resize(3);
            reference_derivative_cross_current.resize(3);

            reference_cross_current_variation.clear();
            reference_cross_current_derivative_variation.clear();
            reference_derivative_cross_current_variation.clear();
            reference_cross_current.clear();
            tangent_cross_product_derivative.clear();
            reference_derivative_cross_current.clear();

            reference_cross_current = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangent);
            tangent_cross_product_derivative = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangentDerivative);
            reference_derivative_cross_current = MathUtils<double>::CrossProduct(rReferenceTangentDerivative, rCurrentTangent);

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++)
                {
                    for (size_t u = 0;u < 3;u++)
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            size_t xyz = r % DofsPerNode;
                            if (xyz > 2)
                            {
                                reference_cross_current_variation[t * NumberOfDofs + r] = 0;
                                reference_cross_current_derivative_variation[t * NumberOfDofs + r] = 0;
                                reference_derivative_cross_current_variation[t * NumberOfDofs + r] = 0;
                            }
                            else
                            {
                                reference_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                                reference_cross_current_derivative_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentDerivativeVariation[u * NumberOfDofs + r];
                                reference_derivative_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                            }
                        }
                    }
                }
            }

            Vector tangent_dot_product_variation;
            tangent_dot_product_variation.resize(NumberOfDofs);
            tangent_dot_product_variation.clear();
            Vector tangent_dot_product_derivative_variation;
            tangent_dot_product_derivative_variation.resize(NumberOfDofs);
            tangent_dot_product_derivative_variation.clear();
            Vector reference_derivative_dot_current_variation;
            reference_derivative_dot_current_variation.resize(NumberOfDofs);
            reference_derivative_dot_current_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyz = r % DofsPerNode; 

                    if (xyz > 2)
                    {
                        tangent_dot_product_variation(r) += 0;
                        tangent_dot_product_derivative_variation(r) += 0;
                        reference_derivative_dot_current_variation(r) += 0;
                    }
                    else
                    {
                        
                        tangent_dot_product_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                        tangent_dot_product_derivative_variation(r) += rCurrentTangentDerivativeVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                        reference_derivative_dot_current_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangentDerivative[t];
                    }
                }
            }


            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz = r % DofsPerNode; 

                        if (xyz > 2)
                                rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += 0;
                        else
                        {
                            if (t == u)
                                rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r);

                            else
                            {
                                {
                                    for (int k = 0; k < 3; k++)
                                        rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += levi_civita[t][k][u] * reference_cross_current_derivative_variation[r + k * NumberOfDofs] + levi_civita[t][k][u] * reference_derivative_cross_current_variation[r + k * NumberOfDofs]; 
                                }

                            }
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += (2 * (tangent_dot_product_variation(r)) * (reference_dot_current_derivative + reference_derivative_dot_current) / pow(1.0 + tangent_dot_product, 3) - (tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r)) / pow(1.0 + tangent_dot_product, 2)) * reference_cross_current[t] * reference_cross_current[u];
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += -(reference_dot_current_derivative + reference_derivative_dot_current) / pow(1.0 + tangent_dot_product, 2) * ((reference_cross_current_variation[t * NumberOfDofs + r]) * reference_cross_current[u] + reference_cross_current[t] * (reference_cross_current_variation[u * NumberOfDofs + r]));
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += -(tangent_dot_product_variation(r)) / pow(1.0 + tangent_dot_product, 2) * ((tangent_cross_product_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current[u] + reference_cross_current[t] * (tangent_cross_product_derivative[u] + reference_derivative_cross_current[u]));
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += 1.0 / (1.0 + tangent_dot_product) * ((reference_cross_current_derivative_variation[t * NumberOfDofs + r] + reference_derivative_cross_current_variation[t * NumberOfDofs + r]) * reference_cross_current[u] + (reference_cross_current_variation[t * NumberOfDofs + r]) * (tangent_cross_product_derivative[u] + reference_derivative_cross_current[u]) + (tangent_cross_product_derivative[t] + reference_derivative_cross_current[t]) * (reference_cross_current_variation[u * NumberOfDofs + r]) + reference_cross_current[t] * (reference_cross_current_derivative_variation[u * NumberOfDofs + r] + reference_derivative_cross_current_variation[u * NumberOfDofs + r]));
                        }

                    }
                }
            }

        }


        void BeamBaseElement3D::ComputeLambdaMatrixSecondDerivativeVariation(
            const array_1d<double, 3>& rReferenceTangent,
            const array_1d<double, 3>& rCurrentTangent,
            const array_1d<double, 3>& rReferenceTangentDerivative,
            const array_1d<double, 3>& rReferenceTangentSecondDerivative,
            const Vector& rCurrentTangentVariation,
            const array_1d<double, 3>& rCurrentTangentDerivative,
            const array_1d<double, 3>& rCurrentTangentSecondDerivative,
            const Vector& rCurrentTangentDerivativeVariation,
            const Vector& rCurrentTangentSecondDerivativeVariation,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rLambdaMatrixSecondDerivativeVariation)
        {
                    
                    rLambdaMatrixSecondDerivativeVariation.resize(3 * NumberOfDofs, 3);
                    rLambdaMatrixSecondDerivativeVariation.clear();  
    
                    double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
                    double reference_dot_current_derivative = inner_prod(rReferenceTangent, rCurrentTangentDerivative);
                    double reference_derivative_dot_current = inner_prod(rReferenceTangentDerivative, rCurrentTangent);
                    double reference_second_derivative_dot_current = inner_prod(rReferenceTangentSecondDerivative, rCurrentTangent);
                    double tangent_derivative_dot_product = inner_prod(rReferenceTangentDerivative, rCurrentTangentDerivative);
                    double reference_dot_current_second_derivative = inner_prod(rReferenceTangent, rCurrentTangentSecondDerivative);
                    double tangent_dot_product_derivative = reference_derivative_dot_current + reference_dot_current_derivative;
                    double tangent_dot_product_second_derivative = reference_second_derivative_dot_current + 2 * tangent_derivative_dot_product + reference_dot_current_second_derivative;
    
                    int levi_civita[3][3][3];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                levi_civita[i][j][k] = 0;
                            }
                        }
                    }
    
                    levi_civita[0][1][2] = 1;
                    levi_civita[2][0][1] = 1;
                    levi_civita[1][2][0] = 1;
    
                    levi_civita[0][2][1] = -1;
                    levi_civita[1][0][2] = -1;
                    levi_civita[2][1][0] = -1;
    
                    Vector reference_cross_current_variation;           
                    Vector reference_cross_current_derivative_variation;        
                    Vector reference_derivative_cross_current_variation;        
                    Vector reference_second_derivative_cross_current_variation;     
                    Vector reference_derivative_cross_current_derivative_variation;     
                    Vector reference_cross_current_second_derivative_variation;     
                    Vector reference_cross_current;              
                    Vector reference_cross_current_derivative;           
                    Vector reference_derivative_cross_current;           
                    Vector reference_cross_current_second_derivative;        
                    Vector reference_second_derivative_cross_current;        
                    Vector reference_derivative_cross_current_derivative;        
                    Vector tangent_cross_product_derivative;          
                    Vector tangent_cross_product_second_derivative;       
                    Vector tangent_cross_product_derivative_variation;       
                    Vector tangent_cross_product_second_derivative_variation;    
    
                    reference_cross_current_variation.resize(NumberOfDofs * 3);
                    reference_cross_current_derivative_variation.resize(NumberOfDofs * 3);
                    reference_derivative_cross_current_variation.resize(NumberOfDofs * 3);
                    reference_second_derivative_cross_current_variation.resize(NumberOfDofs * 3);
                    reference_derivative_cross_current_derivative_variation.resize(NumberOfDofs * 3);
                    reference_cross_current_second_derivative_variation.resize(NumberOfDofs * 3);
                    tangent_cross_product_derivative_variation.resize(NumberOfDofs * 3);
                    tangent_cross_product_second_derivative_variation.resize(NumberOfDofs * 3);
                    reference_cross_current.resize(3);
                    reference_cross_current_derivative.resize(3);
                    reference_derivative_cross_current.resize(3);
                    reference_second_derivative_cross_current.resize(3);
                    reference_derivative_cross_current_derivative.resize(3);
                    reference_cross_current_second_derivative.resize(3);
                    tangent_cross_product_derivative.resize(3);
                    tangent_cross_product_second_derivative.resize(3);
    
                    reference_cross_current_variation.clear();
                    reference_cross_current_derivative_variation.clear();
                    reference_derivative_cross_current_variation.clear();
                    reference_second_derivative_cross_current_variation.clear();
                    reference_derivative_cross_current_derivative_variation.clear();
                    reference_cross_current_second_derivative_variation.clear();
                    reference_cross_current.clear();
                    reference_cross_current_derivative.clear();
                    reference_derivative_cross_current.clear();
                    reference_second_derivative_cross_current.clear();
                    reference_derivative_cross_current_derivative.clear();
                    reference_cross_current_second_derivative.clear();
                    tangent_cross_product_derivative.clear();
                    tangent_cross_product_second_derivative.clear();
                    tangent_cross_product_derivative_variation.clear();
                    tangent_cross_product_second_derivative_variation.clear();
    
                    reference_cross_current = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangent);
                    reference_cross_current_derivative = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangentDerivative);
                    reference_derivative_cross_current = MathUtils<double>::CrossProduct(rReferenceTangentDerivative, rCurrentTangent);
                    reference_second_derivative_cross_current = MathUtils<double>::CrossProduct(rReferenceTangentSecondDerivative, rCurrentTangent);
                    reference_derivative_cross_current_derivative = MathUtils<double>::CrossProduct(rReferenceTangentDerivative, rCurrentTangentDerivative);
                    reference_cross_current_second_derivative = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangentSecondDerivative);
                    tangent_cross_product_derivative = reference_cross_current_derivative + reference_derivative_cross_current;
                    tangent_cross_product_second_derivative = 2 * reference_derivative_cross_current_derivative + reference_second_derivative_cross_current + reference_cross_current_second_derivative;
    
                    for (size_t t = 0; t < 3; t++) 
                    {
                        for (size_t  r  = 0; r < NumberOfDofs; r++)
                        {
                            for (size_t u = 0; u < 3; u++)
                            {
                                for (size_t k = 0; k < 3; k++)
                                {
                                    size_t xyz = r % DofsPerNode;
                                    if (xyz > 2)
                                    {
                                        reference_cross_current_variation[t * NumberOfDofs + r] = 0;
                                        reference_cross_current_derivative_variation[t * NumberOfDofs + r] = 0;
                                        reference_derivative_cross_current_variation[t * NumberOfDofs + r] = 0;
                                        reference_second_derivative_cross_current_variation[t * NumberOfDofs + r] = 0;
                                        reference_derivative_cross_current_derivative_variation[t * NumberOfDofs + r] = 0;
                                        reference_cross_current_second_derivative_variation[t * NumberOfDofs + r] = 0;
                                    }
                                    else
                                    {
                                        reference_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                                        reference_cross_current_derivative_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentDerivativeVariation[u * NumberOfDofs + r];
                                        reference_derivative_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                                        reference_second_derivative_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangentSecondDerivative[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                                        reference_derivative_cross_current_derivative_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentDerivativeVariation[u * NumberOfDofs + r];
                                        reference_cross_current_second_derivative_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentSecondDerivativeVariation[u * NumberOfDofs + r];
                                    }
                                }
                            }
                        }
                    }
                    tangent_cross_product_derivative_variation = reference_cross_current_derivative_variation + reference_derivative_cross_current_variation;
                    tangent_cross_product_second_derivative_variation = 2 * reference_derivative_cross_current_derivative_variation + reference_second_derivative_cross_current_variation + reference_cross_current_second_derivative_variation;
    
                    Vector tangent_dot_product_variation;
                    tangent_dot_product_variation.resize(NumberOfDofs);
                    tangent_dot_product_variation.clear();
                    Vector reference_dot_current_derivative_variation;
                    reference_dot_current_derivative_variation.resize(NumberOfDofs);
                    reference_dot_current_derivative_variation.clear();
                    Vector reference_derivative_dot_current_variation;
                    reference_derivative_dot_current_variation.resize(NumberOfDofs);
                    reference_derivative_dot_current_variation.clear();
                    Vector tangent_derivative_dot_product_variation;
                    tangent_derivative_dot_product_variation.resize(NumberOfDofs);
                    tangent_derivative_dot_product_variation.clear();
                    Vector reference_second_derivative_dot_current_variation;
                    reference_second_derivative_dot_current_variation.resize(NumberOfDofs);
                    reference_second_derivative_dot_current_variation.clear();
                    Vector reference_dot_current_second_derivative_variation;
                    reference_dot_current_second_derivative_variation.resize(NumberOfDofs);
                    reference_dot_current_second_derivative_variation.clear();
    
                    for (size_t t = 0; t < 3; t++)
                    {
                        for (size_t  r  = 0; r < NumberOfDofs; r++) 
                        {
                            size_t xyz = r % DofsPerNode; 
    
                            if (xyz > 2)
                            {
                                tangent_dot_product_variation(r) += 0;
                                reference_dot_current_derivative_variation(r) += 0;
                                reference_derivative_dot_current_variation(r) += 0;
                                tangent_derivative_dot_product_variation(r) += 0;
                                reference_second_derivative_dot_current_variation(r) += 0;
                                reference_dot_current_second_derivative_variation(r) += 0;
                            }
                            else
                            {
                                
                                tangent_dot_product_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                                reference_dot_current_derivative_variation(r) += rCurrentTangentDerivativeVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                                reference_derivative_dot_current_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangentDerivative[t];
                                reference_second_derivative_dot_current_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangentSecondDerivative[t];
                                tangent_derivative_dot_product_variation(r) += rCurrentTangentDerivativeVariation[t * NumberOfDofs + r] * rReferenceTangentDerivative[t];
                                reference_dot_current_second_derivative_variation(r) += rCurrentTangentSecondDerivativeVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                            }
                        }
                    }
                    Vector tangent_dot_product_derivative_variation;
                    tangent_dot_product_derivative_variation.resize(NumberOfDofs);
                    tangent_dot_product_derivative_variation.clear();
                    Vector tangent_dot_product_second_derivative_variation;
                    tangent_dot_product_second_derivative_variation.resize(NumberOfDofs);
                    tangent_dot_product_second_derivative_variation.clear();
    
                    tangent_dot_product_derivative_variation = reference_derivative_dot_current_variation + reference_dot_current_derivative_variation;
                    tangent_dot_product_second_derivative_variation = reference_second_derivative_dot_current_variation + 2 * tangent_derivative_dot_product_variation + reference_dot_current_second_derivative_variation;
    
                    double inverse_one_plus_tangent_dot_squared = 1.0 / pow(1.0 + tangent_dot_product, 2);
                    double inverse_one_plus_tangent_dot_cubed = 1.0 / pow(1.0 + tangent_dot_product, 3);
                    double inverse_one_plus_tangent_dot_fourth = 1.0 / pow(1.0 + tangent_dot_product, 4);
    
                    for (size_t t = 0; t < 3; t++) 
                    {
                        for (size_t u = 0; u < 3; u++)
                        {
                            for (size_t  r  = 0; r < NumberOfDofs; r++)
                            {
                                size_t xyz = r % DofsPerNode; 
    
                                if (xyz > 2)
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += 0;
                                else
                                {
                                    if (t == u)
                                        rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += tangent_dot_product_second_derivative_variation(r);
    
                                    else
                                    {
                                        {
                                            for (int k = 0; k < 3; k++)
                                                rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += levi_civita[t][k][u] * tangent_cross_product_second_derivative_variation[r + k * NumberOfDofs]; 
                                        }
    
                                    }
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += (-tangent_dot_product_second_derivative_variation(r) * inverse_one_plus_tangent_dot_squared + 2 * (tangent_dot_product_variation(r) * tangent_dot_product_second_derivative + 2 * tangent_dot_product_derivative_variation[r] * tangent_dot_product_derivative) * inverse_one_plus_tangent_dot_cubed - 6 * (pow(tangent_dot_product_derivative, 2) * tangent_dot_product_variation[r]) * inverse_one_plus_tangent_dot_fourth) * reference_cross_current[t] * reference_cross_current[u];
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += (-tangent_dot_product_second_derivative * inverse_one_plus_tangent_dot_squared + 2 * pow(tangent_dot_product_derivative, 2) * inverse_one_plus_tangent_dot_cubed) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * (reference_cross_current_variation[u * NumberOfDofs + r]));
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += (4 * tangent_dot_product_derivative * tangent_dot_product_variation[r] * inverse_one_plus_tangent_dot_cubed - 2 * tangent_dot_product_derivative_variation[r] * inverse_one_plus_tangent_dot_squared) * (tangent_cross_product_derivative[t] * reference_cross_current[u] + reference_cross_current[t] * tangent_cross_product_derivative[u]);
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += -2 * tangent_dot_product_derivative * inverse_one_plus_tangent_dot_squared * (tangent_cross_product_derivative_variation[t * NumberOfDofs + r] * reference_cross_current[u] + tangent_cross_product_derivative[t] * reference_cross_current_variation[u * NumberOfDofs + r] + reference_cross_current_variation[t * NumberOfDofs + r] * tangent_cross_product_derivative[u] + reference_cross_current[t] * tangent_cross_product_derivative_variation[u * NumberOfDofs + r]);
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += -tangent_dot_product_variation(r) * inverse_one_plus_tangent_dot_squared * (tangent_cross_product_second_derivative[t] * reference_cross_current[u] + 2 * tangent_cross_product_derivative[t] * tangent_cross_product_derivative[u] + reference_cross_current[t] * tangent_cross_product_second_derivative[u]);
                                    rLambdaMatrixSecondDerivativeVariation(t * NumberOfDofs + r, u) += 1.0 / (1.0 + tangent_dot_product) * (tangent_cross_product_second_derivative_variation[t * NumberOfDofs + r] * reference_cross_current[u] + tangent_cross_product_second_derivative[t] * reference_cross_current_variation[u * NumberOfDofs + r] + 2 * tangent_cross_product_derivative_variation[t * NumberOfDofs + r] * tangent_cross_product_derivative[u] + 2 * tangent_cross_product_derivative[t] * tangent_cross_product_derivative_variation[u * NumberOfDofs + r] + reference_cross_current_variation[t * NumberOfDofs + r] * tangent_cross_product_second_derivative[u] + reference_cross_current[t] * tangent_cross_product_second_derivative_variation[u * NumberOfDofs + r]);
                                }
    
                            }
                        }
                    }
    
                }

void BeamBaseElement3D::ComputeLambdaMatrixDerivativeSecondVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const Vector& rCurrentTangentVariation,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        const Vector& rCurrentTangentDerivativeVariation,
        const Matrix& rCurrentTangentSecondVariation,
        const Matrix& rCurrentTangentDerivativeSecondVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixDerivativeSecondVariation)
        {
            
            rLambdaMatrixDerivativeSecondVariation.clear();  

            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
            double reference_dot_current_derivative = inner_prod(rReferenceTangent, rCurrentTangentDerivative);
            double reference_derivative_dot_current = inner_prod(rReferenceTangentDerivative, rCurrentTangent);
            double tangent_dot_product_derivative = reference_dot_current_derivative + reference_derivative_dot_current;

            int levi_civita[3][3][3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        levi_civita[i][j][k] = 0;
                    }
                }
            }

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector reference_cross_current_variation;
            Vector reference_cross_current_derivative_variation;
            Vector reference_cross_current;
            Vector reference_cross_current_derivative;
            Vector reference_derivative_cross_current;
            Vector reference_derivative_cross_current_variation;

            reference_cross_current_variation.resize(NumberOfDofs * 3);
            reference_cross_current_derivative_variation.resize(NumberOfDofs * 3);
            reference_cross_current.resize(3);
            reference_cross_current_derivative.resize(3);
            reference_derivative_cross_current.resize(3);
            reference_derivative_cross_current_variation.resize(NumberOfDofs * 3);

            reference_cross_current_variation.clear();
            reference_cross_current_derivative_variation.clear();
            reference_cross_current.clear();
            reference_cross_current_derivative.clear();
            reference_derivative_cross_current.clear();
            reference_derivative_cross_current_variation.clear();

            reference_cross_current = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangent);
            reference_cross_current_derivative = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangentDerivative);
            reference_derivative_cross_current = MathUtils<double>::CrossProduct(rReferenceTangentDerivative, rCurrentTangent);

            Vector tangent_dot_product_variation;
            tangent_dot_product_variation.resize(NumberOfDofs);
            tangent_dot_product_variation.clear();
            Vector reference_derivative_dot_current_variation;
            reference_derivative_dot_current_variation.resize(NumberOfDofs);
            reference_derivative_dot_current_variation.clear();
            Vector tangent_dot_product_derivative_variation;
            tangent_dot_product_derivative_variation.resize(NumberOfDofs);
            tangent_dot_product_derivative_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyz = r % DofsPerNode; 

                    if (xyz > 2)
                    {
                        tangent_dot_product_variation(r) += 0;
                        reference_derivative_dot_current_variation(r) += 0;
                        tangent_dot_product_derivative_variation(r) += 0;
                    }
                    else
                    {
                        tangent_dot_product_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                        reference_derivative_dot_current_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangentDerivative[t];
                        tangent_dot_product_derivative_variation(r) += rCurrentTangentDerivativeVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                    }
                }
            }

            Matrix tangent_dot_product_second_variation;
            tangent_dot_product_second_variation.resize(NumberOfDofs, NumberOfDofs);
            tangent_dot_product_second_variation.clear();
            Matrix tangent_dot_product_derivative_second_variation;
            tangent_dot_product_derivative_second_variation.resize(NumberOfDofs, NumberOfDofs);
            tangent_dot_product_derivative_second_variation.clear();
            Matrix reference_derivative_dot_current_second_variation;
            reference_derivative_dot_current_second_variation.resize(NumberOfDofs, NumberOfDofs);
            reference_derivative_dot_current_second_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyzr = r % DofsPerNode; 
                    
                    for (size_t  s  = 0;s < NumberOfDofs;s++) 
                    {
                        size_t xyzs = s % DofsPerNode; 
                        

                        if (xyzr > 2 || xyzs > 2)
                        {
                            tangent_dot_product_second_variation(r, s) = 0;
                        }
                        else
                        {
                            tangent_dot_product_second_variation(r, s) += rCurrentTangentSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangent[t];
                            tangent_dot_product_derivative_second_variation(r, s) += rCurrentTangentDerivativeSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangent[t];
                            reference_derivative_dot_current_second_variation(r, s) += rCurrentTangentSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangentDerivative[t];
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++)
                {
                    size_t xyz_r = r % DofsPerNode; 
                    for (size_t u = 0;u < 3;u++)
                    {
                        if (xyz_r > 2)
                        {
                            reference_cross_current_variation[t * NumberOfDofs + r] += 0;
                        }
                        else
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                reference_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                                reference_cross_current_derivative_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentDerivativeVariation[u * NumberOfDofs + r];
                                reference_derivative_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                            }
                        }
                    }
                }
            }

            Matrix reference_cross_current_second_variation;
            reference_cross_current_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_cross_current_second_variation.clear();
            Matrix reference_cross_current_derivative_second_variation;
            reference_cross_current_derivative_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_cross_current_derivative_second_variation.clear();
            Matrix reference_derivative_cross_current_second_variation;
            reference_derivative_cross_current_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_derivative_cross_current_second_variation.clear();

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++) 
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz_r = r % DofsPerNode; 
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyz_s = s % DofsPerNode; 
                            if (xyz_r < 3 || xyz_s < 3)
                            {
                                for (size_t k = 0;k < 3;k++)
                                {
                                    reference_cross_current_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentSecondVariation(u * NumberOfDofs + r, s);
                                    reference_cross_current_derivative_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentDerivativeSecondVariation(u * NumberOfDofs + r, s);
                                    reference_derivative_cross_current_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentSecondVariation(u * NumberOfDofs + r, s);
                                }
                            }
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++) 
                {
                    for (size_t  s  = 0;s < NumberOfDofs;s++)
                    {
                        size_t xyz_s = s % DofsPerNode; 

                        for (size_t  r  = 0;r < NumberOfDofs;r++)
                        {
                            size_t xyz_r = r % DofsPerNode; 
                            if (t == u)
                            {
                                if (xyz_r > 2 || xyz_s > 2)
                                    rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * DofsPerNode + s) += 0;
                                else
                                    rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += tangent_dot_product_derivative_second_variation(r, s) + reference_derivative_dot_current_second_variation(r, s);
                            }
                            else
                            {
                                if (xyz_r > 2 || xyz_s > 2)
                                    for (int k = 0; k < 3; k++)
                                        rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, s) += 0;
                                else
                                {
                                    for (int k = 0; k < 3; k++)
                                        rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += levi_civita[t][k][u] * reference_cross_current_derivative_second_variation(k * NumberOfDofs + r, s) + levi_civita[t][k][u] * reference_derivative_cross_current_second_variation(k * NumberOfDofs + r, s); 
                                }
                            }
                        }
                    }
                }
            }
            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyzr = r % DofsPerNode; 
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyzs = s % DofsPerNode; 
                            if (xyzr > 2 || xyzs > 2)
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += 0;
                            else
                            {
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (2 * (tangent_dot_product_variation(r) * (tangent_dot_product_derivative_variation(s) + reference_derivative_dot_current_variation(s)) + tangent_dot_product_derivative * tangent_dot_product_second_variation(r, s)) / pow(1.0 + tangent_dot_product, 3) - 6 * tangent_dot_product_derivative * tangent_dot_product_variation(r) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 4) - (tangent_dot_product_derivative_second_variation(r, s) + reference_derivative_dot_current_second_variation(r, s)) / pow(1.0 + tangent_dot_product, 2) + 2 * (tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r)) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 3)) * reference_cross_current[t] * reference_cross_current[u];
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (2 * tangent_dot_product_variation(r) * tangent_dot_product_derivative / pow(1.0 + tangent_dot_product, 3) - (tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r)) / pow(1.0 + tangent_dot_product, 2)) * (reference_cross_current_variation[t * NumberOfDofs + s] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + s])
                                    + (2 * tangent_dot_product_variation(s) * tangent_dot_product_derivative / pow(1.0 + tangent_dot_product, 3) - (tangent_dot_product_derivative_variation(s) + reference_derivative_dot_current_variation(s)) / pow(1.0 + tangent_dot_product, 2)) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + r]);
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -tangent_dot_product_derivative / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current_second_variation(t * NumberOfDofs + r, s) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + r) * reference_cross_current_variation(u * NumberOfDofs + s) +
                                    reference_cross_current_variation(t * NumberOfDofs + s) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current[t] * reference_cross_current_second_variation(u * NumberOfDofs + r, s));   
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (-tangent_dot_product_second_variation(r, s) / pow(1.0 + tangent_dot_product, 2) + 2 * tangent_dot_product_variation(r) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 3)) * ((reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current[u] + reference_cross_current[t] * (reference_cross_current_derivative[u] + reference_derivative_cross_current[u]));
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (-tangent_dot_product_variation(r) / pow(1.0 + tangent_dot_product, 2)) * ((reference_cross_current_derivative_variation(t * NumberOfDofs + s) + reference_derivative_cross_current_variation(t * NumberOfDofs + s)) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + s) * (reference_cross_current_derivative(u) + reference_derivative_cross_current(u)) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current_variation(u * NumberOfDofs + s) + reference_cross_current[t] * (reference_cross_current_derivative_variation(u * NumberOfDofs + s) + reference_derivative_cross_current_variation(u * NumberOfDofs + s)));
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (-tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 2)) * ((reference_cross_current_derivative_variation(t * NumberOfDofs + r) + reference_derivative_cross_current_variation(t * NumberOfDofs + r)) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + r) * (reference_cross_current_derivative(u) + reference_derivative_cross_current(u)) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current[t] * (reference_cross_current_derivative_variation(u * NumberOfDofs + r) + reference_derivative_cross_current_variation(u * NumberOfDofs + r)));
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += 1.0 / (1.0 + tangent_dot_product) * ((reference_cross_current_derivative_second_variation(t * NumberOfDofs + r, s) + reference_derivative_cross_current_second_variation(t * NumberOfDofs + r, s)) * reference_cross_current[u] + reference_cross_current_second_variation(t * NumberOfDofs + r, s) * (reference_cross_current_derivative[u] + reference_derivative_cross_current[u]) + (reference_cross_current_derivative_variation(t * NumberOfDofs + s) + reference_derivative_cross_current_variation(t * NumberOfDofs + s)) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current_variation(t * NumberOfDofs + s) * (reference_cross_current_derivative_variation(u * NumberOfDofs + r) + reference_derivative_cross_current_variation(u * NumberOfDofs + r)) + (reference_cross_current_derivative_variation(t * NumberOfDofs + r) + reference_derivative_cross_current_variation(t * NumberOfDofs + r)) * reference_cross_current_variation(u * NumberOfDofs + s) + reference_cross_current_variation(t * NumberOfDofs + r) * (reference_cross_current_derivative_variation(u * NumberOfDofs + s) + reference_derivative_cross_current_variation(u * NumberOfDofs + s)) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current_second_variation(u * NumberOfDofs + r, s) + reference_cross_current[t] * (reference_cross_current_derivative_second_variation(u * NumberOfDofs + r, s) + reference_derivative_cross_current_second_variation(u * NumberOfDofs + r, s)));  
                            }
                        }
                    }
                }
            }
        }

    void BeamBaseElement3D::ComputeLambdaMatrixVariations(
            const array_1d<double, 3>& rReferenceTangent,
            const array_1d<double, 3>& rCurrentTangent,
            const array_1d<double, 3>& rReferenceTangentDerivative,
            const Vector& rCurrentTangentVariation,
            const array_1d<double, 3>& rCurrentTangentDerivative,
            const Vector& rCurrentTangentDerivativeVariation,
            const Matrix& rCurrentTangentSecondVariation,
            const Matrix& rCurrentTangentDerivativeSecondVariation,
            const SizeType NumberOfDofs,
            const SizeType DofsPerNode,
            Matrix& rLambdaMatrixVariation,
            Matrix& rLambdaMatrixDerivativeVariation,
            Matrix& rLambdaMatrixSecondVariation,
            Matrix& rLambdaMatrixDerivativeSecondVariation)
        {
            

            
            rLambdaMatrixVariation.clear();  
            
            rLambdaMatrixDerivativeVariation.clear();  
            
            rLambdaMatrixSecondVariation.clear();  
            
            rLambdaMatrixDerivativeSecondVariation.clear();  

            double tangent_dot_product = inner_prod(rReferenceTangent, rCurrentTangent);
            double reference_dot_current_derivative = inner_prod(rReferenceTangent, rCurrentTangentDerivative);
            double reference_derivative_dot_current = inner_prod(rReferenceTangentDerivative, rCurrentTangent);
            double tangent_dot_product_derivative = reference_dot_current_derivative + reference_derivative_dot_current;

            int levi_civita[3][3][3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        levi_civita[i][j][k] = 0;
                    }
                }
            }

            levi_civita[0][1][2] = 1;
            levi_civita[2][0][1] = 1;
            levi_civita[1][2][0] = 1;

            levi_civita[0][2][1] = -1;
            levi_civita[1][0][2] = -1;
            levi_civita[2][1][0] = -1;

            Vector reference_cross_current_variation;
            Vector reference_cross_current_derivative_variation;
            Vector reference_cross_current;
            Vector reference_cross_current_derivative;
            Vector reference_derivative_cross_current;
            Vector reference_derivative_cross_current_variation;

            reference_cross_current_variation.resize(NumberOfDofs * 3);
            reference_cross_current_derivative_variation.resize(NumberOfDofs * 3);
            reference_cross_current.resize(3);
            reference_cross_current_derivative.resize(3);
            reference_derivative_cross_current.resize(3);
            reference_derivative_cross_current_variation.resize(NumberOfDofs * 3);

            reference_cross_current_variation.clear();
            reference_cross_current_derivative_variation.clear();
            reference_cross_current.clear();
            reference_cross_current_derivative.clear();
            reference_derivative_cross_current.clear();
            reference_derivative_cross_current_variation.clear();

            reference_cross_current = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangent);
            reference_cross_current_derivative = MathUtils<double>::CrossProduct(rReferenceTangent, rCurrentTangentDerivative);
            reference_derivative_cross_current = MathUtils<double>::CrossProduct(rReferenceTangentDerivative, rCurrentTangent);

            Vector tangent_dot_product_variation;
            tangent_dot_product_variation.resize(NumberOfDofs);
            tangent_dot_product_variation.clear();
            Vector reference_derivative_dot_current_variation;
            reference_derivative_dot_current_variation.resize(NumberOfDofs);
            reference_derivative_dot_current_variation.clear();
            Vector tangent_dot_product_derivative_variation;
            tangent_dot_product_derivative_variation.resize(NumberOfDofs);
            tangent_dot_product_derivative_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyz = r % DofsPerNode; 

                    if (xyz > 2)
                    {
                        tangent_dot_product_variation(r) += 0;
                        reference_derivative_dot_current_variation(r) += 0;
                        tangent_dot_product_derivative_variation(r) += 0;
                    }
                    else
                    {
                        tangent_dot_product_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                        reference_derivative_dot_current_variation(r) += rCurrentTangentVariation[t * NumberOfDofs + r] * rReferenceTangentDerivative[t];
                        tangent_dot_product_derivative_variation(r) += rCurrentTangentDerivativeVariation[t * NumberOfDofs + r] * rReferenceTangent[t];
                    }
                }
            }

            Matrix tangent_dot_product_second_variation;
            tangent_dot_product_second_variation.resize(NumberOfDofs, NumberOfDofs);
            tangent_dot_product_second_variation.clear();
            Matrix tangent_dot_product_derivative_second_variation;
            tangent_dot_product_derivative_second_variation.resize(NumberOfDofs, NumberOfDofs);
            tangent_dot_product_derivative_second_variation.clear();
            Matrix reference_derivative_dot_current_second_variation;
            reference_derivative_dot_current_second_variation.resize(NumberOfDofs, NumberOfDofs);
            reference_derivative_dot_current_second_variation.clear();

            for (size_t t = 0;t < 3;t++)
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++) 
                {
                    size_t xyzr = r % DofsPerNode; 
                    
                    for (size_t  s  = 0;s < NumberOfDofs;s++) 
                    {
                        size_t xyzs = s % DofsPerNode; 
                        

                        if (xyzr > 2 || xyzs > 2)
                        {
                            tangent_dot_product_second_variation(r, s) = 0;
                        }
                        else
                        {
                            tangent_dot_product_second_variation(r, s) += rCurrentTangentSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangent[t];
                            tangent_dot_product_derivative_second_variation(r, s) += rCurrentTangentDerivativeSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangent[t];
                            reference_derivative_dot_current_second_variation(r, s) += rCurrentTangentSecondVariation(t * NumberOfDofs + r, s) * rReferenceTangentDerivative[t];
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t  r  = 0;r < NumberOfDofs;r++)
                {
                    size_t xyz_r = r % DofsPerNode; 
                    for (size_t u = 0;u < 3;u++)
                    {
                        if (xyz_r > 2)
                        {
                            reference_cross_current_variation[t * NumberOfDofs + r] += 0;
                        }
                        else
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                reference_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                                reference_cross_current_derivative_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentDerivativeVariation[u * NumberOfDofs + r];
                                reference_derivative_cross_current_variation[t * NumberOfDofs + r] += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentVariation[u * NumberOfDofs + r];
                            }
                        }
                    }
                }
            }

            Matrix reference_cross_current_second_variation;
            reference_cross_current_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_cross_current_second_variation.clear();
            Matrix reference_cross_current_derivative_second_variation;
            reference_cross_current_derivative_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_cross_current_derivative_second_variation.clear();
            Matrix reference_derivative_cross_current_second_variation;
            reference_derivative_cross_current_second_variation.resize(3 * NumberOfDofs, NumberOfDofs);
            reference_derivative_cross_current_second_variation.clear();

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++) 
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz_r = r % DofsPerNode; 
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyz_s = s % DofsPerNode; 
                            if (xyz_r < 3 || xyz_s < 3)
                            {
                                for (size_t k = 0;k < 3;k++)
                                {
                                    reference_cross_current_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentSecondVariation(u * NumberOfDofs + r, s);
                                    reference_cross_current_derivative_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangent[k] * rCurrentTangentDerivativeSecondVariation(u * NumberOfDofs + r, s);
                                    reference_derivative_cross_current_second_variation(t * NumberOfDofs + r, s) += levi_civita[t][k][u] * rReferenceTangentDerivative[k] * rCurrentTangentSecondVariation(u * NumberOfDofs + r, s);
                                }
                            }
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyz = r % DofsPerNode; 

                        if (xyz > 2)
                        {
                            
                                        
                        }
                        else
                        {
                            if (t == u)
                            {
                                rLambdaMatrixVariation(t * NumberOfDofs + r, u) += tangent_dot_product_variation(r);
                                rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r);
                            }
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    rLambdaMatrixVariation(t * NumberOfDofs + r, u) += levi_civita[t][k][u] * reference_cross_current_variation[r + k * NumberOfDofs];
                                    rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += levi_civita[t][k][u] * reference_cross_current_derivative_variation[r + k * NumberOfDofs] + levi_civita[t][k][u] * reference_derivative_cross_current_variation[r + k * NumberOfDofs]; 
                                }
                            }
                            
                            rLambdaMatrixVariation(t * NumberOfDofs + r, u) += -tangent_dot_product_variation[r] / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current[t] * reference_cross_current[u]);
                            rLambdaMatrixVariation(t * NumberOfDofs + r, u) += +1.0 / (1.0 + tangent_dot_product) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + r]);
                            
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += (2 * (tangent_dot_product_variation(r)) * (reference_dot_current_derivative + reference_derivative_dot_current) / pow(1.0 + tangent_dot_product, 3) - (tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r)) / pow(1.0 + tangent_dot_product, 2)) * reference_cross_current[t] * reference_cross_current[u];
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += -(reference_dot_current_derivative + reference_derivative_dot_current) / pow(1.0 + tangent_dot_product, 2) * ((reference_cross_current_variation[t * NumberOfDofs + r]) * reference_cross_current[u] + reference_cross_current[t] * (reference_cross_current_variation[u * NumberOfDofs + r]));
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += -(tangent_dot_product_variation(r)) / pow(1.0 + tangent_dot_product, 2) * ((reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current[u] + reference_cross_current[t] * (reference_cross_current_derivative[u] + reference_derivative_cross_current[u]));
                            rLambdaMatrixDerivativeVariation(t * NumberOfDofs + r, u) += 1.0 / (1.0 + tangent_dot_product) * ((reference_cross_current_derivative_variation[t * NumberOfDofs + r] + reference_derivative_cross_current_variation[t * NumberOfDofs + r]) * reference_cross_current[u] + (reference_cross_current_variation[t * NumberOfDofs + r]) * (reference_cross_current_derivative[u] + reference_derivative_cross_current[u]) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * (reference_cross_current_variation[u * NumberOfDofs + r]) + reference_cross_current[t] * (reference_cross_current_derivative_variation[u * NumberOfDofs + r] + reference_derivative_cross_current_variation[u * NumberOfDofs + r]));
                        }
                    }
                }
            }

            for (size_t t = 0;t < 3;t++) 
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t  r  = 0;r < NumberOfDofs;r++)
                    {
                        size_t xyzr = r % DofsPerNode; 
                        for (size_t  s  = 0;s < NumberOfDofs;s++)
                        {
                            size_t xyzs = s % DofsPerNode; 
                            if (xyzr > 2 || xyzs > 2)
                            {
                                
                                
                            }
                            else
                            {
                                if (t == u)
                                {
                                    rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += tangent_dot_product_second_variation(r, s);
                                    rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += tangent_dot_product_derivative_second_variation(r, s) + reference_derivative_dot_current_second_variation(r, s);
                                }
                                else
                                {
                                    for (int k = 0; k < 3; k++)
                                    {
                                        rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += levi_civita[t][k][u] * reference_cross_current_second_variation(k * NumberOfDofs + r, s); 
                                        rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += levi_civita[t][k][u] * reference_cross_current_derivative_second_variation(k * NumberOfDofs + r, s) + levi_civita[t][k][u] * reference_derivative_cross_current_second_variation(k * NumberOfDofs + r, s); 
                                    }
                                }
                                
                                rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (2 * tangent_dot_product_variation(r) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 3) - tangent_dot_product_second_variation(r, s) / pow(1.0 + tangent_dot_product, 2)) * reference_cross_current[t] * reference_cross_current[u];
                                rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -tangent_dot_product_variation(r) / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current_variation[t * NumberOfDofs + s] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + s])
                                    - tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + r]);
                                rLambdaMatrixSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += 1.0 / (1.0 + tangent_dot_product) * (reference_cross_current_second_variation(t * NumberOfDofs + r, s) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + r) * reference_cross_current_variation(u * NumberOfDofs + s) +
                                    reference_cross_current_variation(t * NumberOfDofs + s) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current[t] * reference_cross_current_second_variation(u * NumberOfDofs + s, r));
                                
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (2 * (tangent_dot_product_variation(r) * (tangent_dot_product_derivative_variation(s) + reference_derivative_dot_current_variation(s)) + tangent_dot_product_derivative * tangent_dot_product_second_variation(r, s)) / pow(1.0 + tangent_dot_product, 3) - 6 * tangent_dot_product_derivative * tangent_dot_product_variation(r) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 4) - (tangent_dot_product_derivative_second_variation(r, s) + reference_derivative_dot_current_second_variation(r, s)) / pow(1.0 + tangent_dot_product, 2) + 2 * (tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r)) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 3)) * reference_cross_current[t] * reference_cross_current[u];
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (2 * tangent_dot_product_variation(r) * tangent_dot_product_derivative / pow(1.0 + tangent_dot_product, 3) - (tangent_dot_product_derivative_variation(r) + reference_derivative_dot_current_variation(r)) / pow(1.0 + tangent_dot_product, 2)) * (reference_cross_current_variation[t * NumberOfDofs + s] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + s])
                                    + (2 * tangent_dot_product_variation(s) * tangent_dot_product_derivative / pow(1.0 + tangent_dot_product, 3) - (tangent_dot_product_derivative_variation(s) + reference_derivative_dot_current_variation(s)) / pow(1.0 + tangent_dot_product, 2)) * (reference_cross_current_variation[t * NumberOfDofs + r] * reference_cross_current[u] + reference_cross_current[t] * reference_cross_current_variation[u * NumberOfDofs + r]);
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += -tangent_dot_product_derivative / pow(1.0 + tangent_dot_product, 2) * (reference_cross_current_second_variation(t * NumberOfDofs + r, s) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + r) * reference_cross_current_variation(u * NumberOfDofs + s) +
                                    reference_cross_current_variation(t * NumberOfDofs + s) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current[t] * reference_cross_current_second_variation(u * NumberOfDofs + r, s));   
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (-tangent_dot_product_second_variation(r, s) / pow(1.0 + tangent_dot_product, 2) + 2 * tangent_dot_product_variation(r) * tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 3)) * ((reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current[u] + reference_cross_current[t] * (reference_cross_current_derivative[u] + reference_derivative_cross_current[u]));
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (-tangent_dot_product_variation(r) / pow(1.0 + tangent_dot_product, 2)) * ((reference_cross_current_derivative_variation(t * NumberOfDofs + s) + reference_derivative_cross_current_variation(t * NumberOfDofs + s)) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + s) * (reference_cross_current_derivative(u) + reference_derivative_cross_current(u)) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current_variation(u * NumberOfDofs + s) + reference_cross_current[t] * (reference_cross_current_derivative_variation(u * NumberOfDofs + s) + reference_derivative_cross_current_variation(u * NumberOfDofs + s)));
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += (-tangent_dot_product_variation(s) / pow(1.0 + tangent_dot_product, 2)) * ((reference_cross_current_derivative_variation(t * NumberOfDofs + r) + reference_derivative_cross_current_variation(t * NumberOfDofs + r)) * reference_cross_current[u] + reference_cross_current_variation(t * NumberOfDofs + r) * (reference_cross_current_derivative(u) + reference_derivative_cross_current(u)) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current[t] * (reference_cross_current_derivative_variation(u * NumberOfDofs + r) + reference_derivative_cross_current_variation(u * NumberOfDofs + r)));
                                rLambdaMatrixDerivativeSecondVariation(t * NumberOfDofs + r, u * NumberOfDofs + s) += 1.0 / (1.0 + tangent_dot_product) * ((reference_cross_current_derivative_second_variation(t * NumberOfDofs + r, s) + reference_derivative_cross_current_second_variation(t * NumberOfDofs + r, s)) * reference_cross_current[u] + reference_cross_current_second_variation(t * NumberOfDofs + r, s) * (reference_cross_current_derivative[u] + reference_derivative_cross_current[u]) + (reference_cross_current_derivative_variation(t * NumberOfDofs + s) + reference_derivative_cross_current_variation(t * NumberOfDofs + s)) * reference_cross_current_variation(u * NumberOfDofs + r) + reference_cross_current_variation(t * NumberOfDofs + s) * (reference_cross_current_derivative_variation(u * NumberOfDofs + r) + reference_derivative_cross_current_variation(u * NumberOfDofs + r)) + (reference_cross_current_derivative_variation(t * NumberOfDofs + r) + reference_derivative_cross_current_variation(t * NumberOfDofs + r)) * reference_cross_current_variation(u * NumberOfDofs + s) + reference_cross_current_variation(t * NumberOfDofs + r) * (reference_cross_current_derivative_variation(u * NumberOfDofs + s) + reference_derivative_cross_current_variation(u * NumberOfDofs + s)) + (reference_cross_current_derivative[t] + reference_derivative_cross_current[t]) * reference_cross_current_second_variation(u * NumberOfDofs + r, s) + reference_cross_current[t] * (reference_cross_current_derivative_second_variation(u * NumberOfDofs + r, s) + reference_derivative_cross_current_second_variation(u * NumberOfDofs + r, s)));  

                            }
                        }
                    }
                }
            }

        }



} // namespace Kratos