//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala Pascual
//
//

#ifndef KRATOS_BLURRINESS_CALCULATOR
#define KRATOS_BLURRINESS_CALCULATOR

// System includes
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/openmp_utils.h"

// Application includes

namespace Kratos
{
    ///@addtogroup ConvectionDiffusionApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// @brief Blurriness of a stratified scalar field at cross-sections of a channel.
    ///
    /// For a section S (a surface model part of 2D conditions immersed in the 3D main
    /// model part) the blurriness is
    ///
    ///     B^2 = int_S (psi_h - f)^2 (u_h . n) dS  /  int_S (psi_inf - f)^2 (u_h . n) dS ,
    ///
    /// where psi_h is the scalar (interpolated from the 3D elements), f the inlet step
    /// profile (layer values separated by the chords z = z_i), u_h . n the normal
    /// velocity (the flow weight, replacing the analytical 1 - r^2) and psi_inf the far
    /// field. The integral runs over the WHOLE section: every layer contributes its
    /// flow-weighted squared departure from the inlet datum, so B runs from 0 (intact
    /// stratification) to 1 (fully mixed). This replaces the earlier per-interface
    /// normalisation (integral over the two layers adjacent to each chord), which is
    /// deprecated.
    ///
    /// By default psi_inf is the flow-weighted mean of the step profile evaluated with
    /// the same surface quadrature (the far field of the axial-diffusion-free problem);
    /// SetFarFieldValue overrides it (e.g. with the finite-Peclet far field or an
    /// outlet average).
    ///
    /// The Gauss points of the surface conditions are located in the 3D mesh with a
    /// bin-based locator, so the surface mesh is independent of the volume mesh; the
    /// step profile is evaluated at the Gauss point (not at the condition centre), so a
    /// surface mesh conforming to the chords integrates the discontinuity exactly.
    class BlurrinessCalculator
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of BlurrinessCalculator
        KRATOS_CLASS_POINTER_DEFINITION(BlurrinessCalculator);

        typedef BinBasedFastPointLocator<3> LocatorType;
        typedef LocatorType::ResultContainerType ResultContainerType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        /// @param rMainModelPart      3D model part carrying the scalar and the velocity.
        /// @param rSurfaceModelParts  cross-sections (2D conditions), one blurriness each.
        /// @param rInterfaces         chords z_i of the inlet profile, increasing, strictly inside (-R, R); N-1 values.
        /// @param rLayerValues        inlet value of each layer, from the bottom (z < z_0) up; N values.
        /// @param rScalarVariable     nodal scalar psi_h on the main model part.
        /// @param rVelocityVariable   nodal velocity on the main model part (flow weight u . n).
        BlurrinessCalculator(
            ModelPart &rMainModelPart,
            std::vector<ModelPart *> rSurfaceModelParts,
            const std::vector<double> &rInterfaces,
            const std::vector<double> &rLayerValues,
            const Variable<double> &rScalarVariable,
            const Variable<array_1d<double, 3>> &rVelocityVariable = VELOCITY)
            : mrMainModelPart(rMainModelPart), mSurfaceModelParts(rSurfaceModelParts),
              mInterfaces(rInterfaces), mLayerValues(rLayerValues),
              mrScalarVariable(rScalarVariable), mrVelocityVariable(rVelocityVariable),
              mFarFieldOverride(0.0), mHasFarFieldOverride(false), mIsComputed(false)
        {
            KRATOS_ERROR_IF(mLayerValues.empty()) << "BlurrinessCalculator: at least one layer value is required." << std::endl;
            KRATOS_ERROR_IF(mInterfaces.size() + 1 != mLayerValues.size())
                << "BlurrinessCalculator: " << mLayerValues.size() << " layer values need "
                << mLayerValues.size() - 1 << " interfaces (chords strictly inside the section), got "
                << mInterfaces.size() << "." << std::endl;
            for (std::size_t i = 1; i < mInterfaces.size(); ++i)
                KRATOS_ERROR_IF(mInterfaces[i] <= mInterfaces[i - 1])
                    << "BlurrinessCalculator: the interfaces must be strictly increasing." << std::endl;
            KRATOS_ERROR_IF_NOT(mrMainModelPart.HasNodalSolutionStepVariable(mrScalarVariable))
                << "BlurrinessCalculator: " << mrMainModelPart.Name() << " has no nodal variable "
                << mrScalarVariable.Name() << std::endl;
            KRATOS_ERROR_IF_NOT(mrMainModelPart.HasNodalSolutionStepVariable(mrVelocityVariable))
                << "BlurrinessCalculator: " << mrMainModelPart.Name() << " has no nodal variable "
                << mrVelocityVariable.Name() << std::endl;
            const std::size_t n = mSurfaceModelParts.size();
            mBlurriness.assign(n, 0.0);
            mFarField.assign(n, 0.0);
            mNumerator.assign(n, 0.0);
            mDenominator.assign(n, 0.0);
            mFlowRate.assign(n, 0.0);
            mArea.assign(n, 0.0);
            mExtrapolatedPoints.assign(n, 0);
        }

        /// Turn back information as a string.
        std::string Info() const
        {
            std::stringstream buffer;
            buffer << "BlurrinessCalculator";
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const { rOStream << "BlurrinessCalculator"; }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const {}

        ///@}
        ///@name Operations
        ///@{

        /// @brief Impose the far-field value used in the denominator (all sections).
        void SetFarFieldValue(const double FarFieldValue)
        {
            mFarFieldOverride = FarFieldValue;
            mHasFarFieldOverride = true;
            mIsComputed = false;
        }

        /// @brief Use again the flow-weighted inlet mean of each section as far field.
        void UnsetFarFieldValue()
        {
            mHasFarFieldOverride = false;
            mIsComputed = false;
        }

        /// @brief Compute the blurriness of every surface model part.
        void ComputeBlurriness()
        {
            // Bins of the point locator. Kratos' default takes N^(1/3) cells per
            // direction, which for a fine mesh gives cells smaller than the elements
            // and registers every element in many cells (gigabytes for 10^6 elements).
            // A cell of the order of the element size keeps the structure small; the
            // default here is twice the largest element diameter of the main model part.
            LocatorType locator(mrMainModelPart);
            double cell_size = mSearchCellSize;
            if (!(cell_size > 0.0))
            {
                double max_diameter = 0.0;
#pragma omp parallel for reduction(max : max_diameter)
                for (int e = 0; e < static_cast<int>(mrMainModelPart.NumberOfElements()); ++e)
                {
                    const Geometry<Node> &r_geom = (mrMainModelPart.ElementsBegin() + e)->GetGeometry();
                    const double d = r_geom.Length();
                    if (d > max_diameter) max_diameter = d;
                }
                cell_size = 2.0 * max_diameter;
            }
            locator.UpdateSearchDatabaseAssignedSize(cell_size);
            const unsigned max_results = 1000;
            // A section of the exact circle overhangs the polygonal boundary of the
            // volume mesh by the sagitta of its boundary edges (~h^2 / 8R). Gauss
            // points in that sliver are not inside any element: they are assigned
            // to the boundary element found with a relaxed tolerance on the local
            // coordinates, i.e. the finite-element field is extrapolated linearly
            // over a distance that vanishes with the mesh size.
            const double locator_tolerance = 1.0e-5;
            const double relaxed_tolerance = mRelaxedTolerance;
            const std::size_t n_layers = mLayerValues.size();

            for (std::size_t m = 0; m < mSurfaceModelParts.size(); ++m)
            {
                ModelPart &r_surface = *(mSurfaceModelParts[m]);
                const int number_of_conditions = static_cast<int>(r_surface.NumberOfConditions());
                KRATOS_ERROR_IF(number_of_conditions == 0)
                    << "BlurrinessCalculator: surface model part " << r_surface.Name() << " has no conditions." << std::endl;

                // Thread-local accumulators: flow rate q, area A, and the three moments
                // int (psi - f)^2 dq, int f dq, int f^2 dq, plus the flow of every layer.
                const int n_threads = ParallelUtilities::GetNumThreads();
                std::vector<double> q_t(n_threads, 0.0), a_t(n_threads, 0.0), e2_t(n_threads, 0.0),
                    f1_t(n_threads, 0.0), f2_t(n_threads, 0.0);
                std::vector<std::vector<double>> q_layer_t(n_threads, std::vector<double>(n_layers, 0.0));
                std::vector<int> missing_t(n_threads, 0), extrapolated_t(n_threads, 0);

#pragma omp parallel
                {
                    const int tid = OpenMPUtils::ThisThread();
                    ResultContainerType results(max_results);
                    Vector shape_functions;
                    Element::Pointer p_element;

#pragma omp for schedule(dynamic)
                    for (int c = 0; c < number_of_conditions; ++c)
                    {
                        auto it_cond = r_surface.ConditionsBegin() + c;
                        Geometry<Node> &r_geometry = it_cond->GetGeometry();
                        const unsigned num_nodes = r_geometry.size();
                        // 3-point rule on triangles (6 on quads): exact for the square of a
                        // linear interpolant, which is what the volume mesh provides.
                        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
                        const auto &r_integration_points = r_geometry.IntegrationPoints(integration_method);
                        const unsigned n_gauss = r_geometry.IntegrationPointsNumber(integration_method);
                        Vector det_j(n_gauss);
                        r_geometry.DeterminantOfJacobian(det_j, integration_method);
                        const Matrix &N_container = r_geometry.ShapeFunctionsValues(integration_method);
                        array_1d<double, 3> normal = r_geometry.UnitNormal(r_geometry.Center());

                        for (unsigned g = 0; g < n_gauss; ++g)
                        {
                            array_1d<double, 3> gauss_point = ZeroVector(3);
                            for (unsigned n = 0; n < num_nodes; ++n)
                                noalias(gauss_point) += N_container(g, n) * r_geometry[n].Coordinates();

                            bool found = locator.FindPointOnMesh(gauss_point, shape_functions, p_element,
                                                                 results.begin(), max_results, locator_tolerance);
                            if (!found)
                            {
                                found = locator.FindPointOnMesh(gauss_point, shape_functions, p_element,
                                                                results.begin(), max_results, relaxed_tolerance);
                                if (found)
                                    extrapolated_t[tid] += 1;
                            }
                            if (!found)
                            {
                                missing_t[tid] += 1;
                                continue;
                            }

                            // Interpolate the scalar and the normal velocity in the 3D element.
                            Geometry<Node> &r_elem_geometry = p_element->GetGeometry();
                            double scalar = 0.0, normal_velocity = 0.0;
                            for (unsigned n = 0; n < r_elem_geometry.size(); ++n)
                            {
                                const double N = shape_functions[n];
                                scalar += N * r_elem_geometry[n].FastGetSolutionStepValue(mrScalarVariable);
                                const array_1d<double, 3> &vel = r_elem_geometry[n].FastGetSolutionStepValue(mrVelocityVariable);
                                normal_velocity += N * inner_prod(vel, normal);
                            }

                            const std::size_t layer = LayerIndex(gauss_point[2]);
                            const double f = mLayerValues[layer];
                            const double weight = r_integration_points[g].Weight() * det_j[g];
                            const double dq = weight * normal_velocity;
                            q_t[tid] += dq;
                            a_t[tid] += weight;
                            q_layer_t[tid][layer] += dq;
                            e2_t[tid] += dq * (scalar - f) * (scalar - f);
                            f1_t[tid] += dq * f;
                            f2_t[tid] += dq * f * f;
                        }
                    }
                }

                double q = 0.0, area = 0.0, e2 = 0.0, f1 = 0.0, f2 = 0.0;
                std::vector<double> q_layer(n_layers, 0.0);
                int missing = 0, extrapolated = 0;
                for (int t = 0; t < n_threads; ++t)
                {
                    q += q_t[t]; area += a_t[t]; e2 += e2_t[t]; f1 += f1_t[t]; f2 += f2_t[t];
                    missing += missing_t[t]; extrapolated += extrapolated_t[t];
                    for (std::size_t i = 0; i < n_layers; ++i) q_layer[i] += q_layer_t[t][i];
                }
                KRATOS_ERROR_IF(missing > 0)
                    << "BlurrinessCalculator: " << missing << " Gauss points of surface model part "
                    << r_surface.Name() << " were not found inside " << mrMainModelPart.Name() << "." << std::endl;
                KRATOS_ERROR_IF(std::abs(q) <= 0.0)
                    << "BlurrinessCalculator: zero flow rate through " << r_surface.Name()
                    << "; the velocity field or the surface normal is not set." << std::endl;
                // Orient the flow weight with the flow: the sign of the surface normal
                // is a meshing convention that must not enter the ratio.
                if (q < 0.0)
                {
                    q = -q; e2 = -e2; f1 = -f1; f2 = -f2;
                    for (double &v : q_layer) v = -v;
                }

                // Far field: flow-weighted mean of the step profile through this section
                // unless a value was imposed.
                const double psi_inf = mHasFarFieldOverride ? mFarFieldOverride : f1 / q;
                const double denominator = psi_inf * psi_inf * q - 2.0 * psi_inf * f1 + f2;
                KRATOS_ERROR_IF(denominator <= 0.0)
                    << "BlurrinessCalculator: the inlet profile equals its far field on " << r_surface.Name()
                    << "; the blurriness is undefined." << std::endl;

                mExtrapolatedPoints[m] = extrapolated;
                mFlowRate[m] = q;
                mArea[m] = area;
                mFarField[m] = psi_inf;
                mNumerator[m] = e2;
                mDenominator[m] = denominator;
                mBlurriness[m] = std::sqrt(std::max(0.0, e2 / denominator));
                mLayerFlowRates.resize(mSurfaceModelParts.size());
                mLayerFlowRates[m] = q_layer;
            }
            mIsComputed = true;
        }

        /// @brief Blurriness of every surface model part (same order as the constructor).
        std::vector<double> GetBlurriness() const
        {
            CheckComputed("GetBlurriness");
            return mBlurriness;
        }

        /// @brief Far-field value used for every surface model part.
        std::vector<double> GetFarFieldValues() const
        {
            CheckComputed("GetFarFieldValues");
            return mFarField;
        }

        /// @brief Numerator int (psi_h - f)^2 dq of every surface.
        std::vector<double> GetNumerators() const
        {
            CheckComputed("GetNumerators");
            return mNumerator;
        }

        /// @brief Denominator int (psi_inf - f)^2 dq of every surface.
        std::vector<double> GetDenominators() const
        {
            CheckComputed("GetDenominators");
            return mDenominator;
        }

        /// @brief Flow rate int u . n dS through every surface (oriented with the flow).
        std::vector<double> GetFlowRates() const
        {
            CheckComputed("GetFlowRates");
            return mFlowRate;
        }

        /// @brief Area of every surface (a check of the surface mesh).
        std::vector<double> GetAreas() const
        {
            CheckComputed("GetAreas");
            return mArea;
        }

        /// @brief Number of Gauss points per surface assigned by extrapolation (outside every element).
        std::vector<int> GetExtrapolatedPointCounts() const
        {
            CheckComputed("GetExtrapolatedPointCounts");
            return mExtrapolatedPoints;
        }

        /// @brief Size of the cells of the bin-based point locator (default: twice the largest element diameter).
        void SetSearchCellSize(const double CellSize)
        {
            mSearchCellSize = CellSize;
            mIsComputed = false;
        }

        /// @brief Tolerance on the local coordinates accepted for points outside every element (default 0.05).
        void SetRelaxedTolerance(const double Tolerance)
        {
            mRelaxedTolerance = Tolerance;
            mIsComputed = false;
        }

        /// @brief Flow rate of every layer of every surface.
        std::vector<std::vector<double>> GetLayerFlowRates() const
        {
            CheckComputed("GetLayerFlowRates");
            return mLayerFlowRates;
        }

        ///@}

    private:
        ///@name Member Variables
        ///@{

        ModelPart &mrMainModelPart;
        std::vector<ModelPart *> mSurfaceModelParts;
        std::vector<double> mInterfaces;
        std::vector<double> mLayerValues;
        const Variable<double> &mrScalarVariable;
        const Variable<array_1d<double, 3>> &mrVelocityVariable;
        double mFarFieldOverride;
        bool mHasFarFieldOverride;
        bool mIsComputed;

        std::vector<double> mBlurriness, mFarField, mNumerator, mDenominator, mFlowRate, mArea;
        std::vector<std::vector<double>> mLayerFlowRates;
        std::vector<int> mExtrapolatedPoints;
        double mRelaxedTolerance = 0.05;
        double mSearchCellSize = 0.0;

        ///@}
        ///@name Deleted special members
        ///@{

        BlurrinessCalculator() = delete;
        BlurrinessCalculator &operator=(BlurrinessCalculator const &rOther) = delete;
        BlurrinessCalculator(BlurrinessCalculator const &rOther) = delete;

        ///@}
        ///@name Private Operations
        ///@{

        /// @brief Layer containing the height z: the first interface above z.
        std::size_t LayerIndex(const double z) const
        {
            std::size_t layer = 0;
            while (layer < mInterfaces.size() && z > mInterfaces[layer])
                ++layer;
            return layer;
        }

        void CheckComputed(const char *caller) const
        {
            KRATOS_ERROR_IF_NOT(mIsComputed)
                << "BlurrinessCalculator::" << caller << ": call ComputeBlurriness first." << std::endl;
        }

        ///@}
    }; // Class BlurrinessCalculator

    ///@}

} // namespace Kratos.

#endif // KRATOS_BLURRINESS_CALCULATOR
