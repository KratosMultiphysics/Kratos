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

#ifndef KRATOS_SURFACE_INTERPOLATED_L2_NORM
#define KRATOS_SURFACE_INTERPOLATED_L2_NORM

// System includes
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "includes/node.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
    ///@addtogroup ConvectionDiffusionApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// @brief Computes the standard L2 norm of a variable over a set of "surface" model parts.
    ///
    /// The surface model parts only contain 2D conditions whose nodes do not carry nodal
    /// solution values (empty geometries). For each Gauss point of a 2D condition, the 3D
    /// element of the main model part containing that point is located, and the variable is
    /// interpolated there. The L2 norm over each surface model part is then
    /// L2 = sqrt( integral_S |var|^2 dS ).
    ///
    /// @tparam TValueType Either `double` (scalar variable) or `array_1d<double,3>` (vector
    /// variable). For a vector, the magnitude is used: |var|^2 = var . var.
    template <class TValueType>
    class SurfaceInterpolatedL2Norm
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of SurfaceInterpolatedL2Norm
        KRATOS_CLASS_POINTER_DEFINITION(SurfaceInterpolatedL2Norm);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        SurfaceInterpolatedL2Norm(
            ModelPart &model_part,
            std::vector<ModelPart *> model_part_vector,
            const Variable<TValueType> &integration_variable)
            : mMainModelPart(model_part), mModelPartVector(model_part_vector), mVariable(integration_variable), mIsComputed(false)
        {
            mL2Norms.resize(mModelPartVector.size(), 0.0);
        }

        ///@}
        ///@name Operations
        ///@{

        /// Turn back information as a string.
        std::string Info() const
        {
            std::stringstream buffer;
            buffer << "SurfaceInterpolatedL2Norm";
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const { rOStream << "SurfaceInterpolatedL2Norm"; }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const {}

        /// @brief Compute the L2 norm of the variable in each surface model part
        void ComputeL2Norm()
        {
            // The variable is interpolated from the main model part nodes, so it must be part
            // of their solution-step data to be accessed through FastGetSolutionStepValue.
            KRATOS_ERROR_IF_NOT(mMainModelPart.HasNodalSolutionStepVariable(mVariable))
                << "The main model part " << mMainModelPart.Name() << " does not have the variable "
                << mVariable.Name() << " in its nodal solution-step data." << std::endl;

            for (unsigned m = 0; m < mModelPartVector.size(); m++)
            {
                ModelPart &r_model_part = *(mModelPartVector[m]);
                const unsigned number_of_conditions = r_model_part.NumberOfConditions();

                double integral = 0.0; // \int{ |var|^2 dS }
                double area = 0.0;     // Total surface area, for reporting

                for (unsigned c = 0; c < number_of_conditions; c++)
                {
                    ModelPart::ConditionsContainerType::iterator it_cond = r_model_part.ConditionsBegin() + c;

                    // Usual chunk of code necessary to integrate on a condition
                    Geometry<Node> &r_geometry = it_cond->GetGeometry();
                    unsigned int NumNodes = r_geometry.size();
                    GeometryData::IntegrationMethod integration_method = it_cond->GetIntegrationMethod();
                    const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints(integration_method);
                    unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
                    Vector detJ_vector(r_number_integration_points);
                    r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
                    Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

                    // Perform the surface integral
                    for (unsigned g = 0; g < r_number_integration_points; g++)
                    {
                        // Global position of this gauss point
                        array_1d<double, 3> gauss_point_global = ZeroVector(3);
                        for (unsigned n = 0; n < NumNodes; n++)
                        {
                            Point node_global_pos = r_geometry.GetPoint(n);
                            for (unsigned d = 0; d < 3; d++)
                            {
                                gauss_point_global[d] += node_global_pos[d] * NContainer(g, n);
                            }
                        }

                        // Find the 3D element containing this gauss point
                        array_1d<double, 3> gauss_point_local;
                        ModelPart::ElementType::Pointer p_elem;
                        FindElementContainingPoint(gauss_point_global, gauss_point_local, p_elem);
                        KRATOS_ERROR_IF(p_elem == nullptr)
                            << "Gauss point " << gauss_point_global << " in condition with ID = " << it_cond->Id() << " not found in model part " << mMainModelPart.Name() << std::endl;

                        // Interpolate the variable at the gauss point and accumulate
                        TValueType value = InterpolateValue(p_elem, gauss_point_local);
                        double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
                        integral += SquaredValue(value) * Weight;
                        area += Weight;
                    }
                }

                mL2Norms[m] = std::sqrt(integral / area);

                // Print info
                std::cout << "\nModel part " << m + 1 << " (" << r_model_part.Name() << "):" << std::endl;
                std::cout << "    - surface area = " << area << std::endl;
                std::cout << "    - L2 norm = " << mL2Norms[m] << std::endl;
            }
            mIsComputed = true;
        }

        /// @brief Return the L2 norm of the variable in each surface model part
        /// @return Vector with one L2 norm per surface model part
        std::vector<double> GetL2Norm()
        {
            if (!mIsComputed)
            {
                KRATOS_ERROR << "The L2 norm is not defined. Method `ComputeL2Norm` must be called before getting the L2 norm." << std::endl;
            }
            return mL2Norms;
        }

        ///@}

    private:
        ///@name Member Variables
        ///@{

        ModelPart &mMainModelPart;
        std::vector<ModelPart *> mModelPartVector;
        const Variable<TValueType> &mVariable;
        std::vector<double> mL2Norms;
        bool mIsComputed;

        ///@}
        ///@name Deleted special members
        ///@{

        /// Default constructor.
        SurfaceInterpolatedL2Norm() = delete;

        /// Assignment operator.
        SurfaceInterpolatedL2Norm &operator=(SurfaceInterpolatedL2Norm const &rOther) = delete;

        /// Copy constructor.
        SurfaceInterpolatedL2Norm(SurfaceInterpolatedL2Norm const &rOther) = delete;

        ///@}
        ///@name Private Operations
        ///@{

        /// @brief Squared magnitude of a scalar value
        static double SquaredValue(double v) { return v * v; }

        /// @brief Squared magnitude of a vector value
        static double SquaredValue(const array_1d<double, 3> &v) { return inner_prod(v, v); }

        /// @brief Find the element that contains the point
        /// @param point Point inside the element we want to find
        /// @param p_pos_local Local coordinates of the point inside the found element
        /// @param p_elem Pointer to the element containing the point
        void FindElementContainingPoint(const array_1d<double, 3> &point, array_1d<double, 3> &p_pos_local, ModelPart::ElementType::Pointer &p_elem)
        {
            p_elem = nullptr;
            const unsigned number_of_elements = mMainModelPart.NumberOfElements();
            for (unsigned int e = 0; e < number_of_elements; e++)
            {
                ModelPart::ElementsContainerType::iterator it_elem = mMainModelPart.ElementsBegin() + e;
                Geometry<Node> &r_geometry = it_elem->GetGeometry();
                if (r_geometry.IsInside(point, p_pos_local, 1e-10))
                {
                    p_elem = mMainModelPart.pGetElement(it_elem->Id());
                    break;
                }
            }
        }

        /// @brief Interpolate the variable at a local position inside a 3D element
        /// @param p_elem The element of the main model part
        /// @param p_pos_local The local coordinates at which the variable is interpolated
        /// @return Interpolated value of the variable
        TValueType InterpolateValue(Element::Pointer p_elem, const array_1d<double, 3> &p_pos_local)
        {
            TValueType value = mVariable.Zero();

            Geometry<Node> &r_geometry = p_elem->GetGeometry();
            unsigned int NumNodes = r_geometry.size();

            for (unsigned n = 0; n < NumNodes; n++)
            {
                const TValueType &nodal_value = r_geometry[n].FastGetSolutionStepValue(mVariable);
                double shape_function_value = r_geometry.ShapeFunctionValue(n, p_pos_local);
                value += nodal_value * shape_function_value;
            }
            return value;
        }

        ///@}
    }; // Class SurfaceInterpolatedL2Norm

    ///@}

}; // namespace Kratos.

#endif // KRATOS_SURFACE_INTERPOLATED_L2_NORM
