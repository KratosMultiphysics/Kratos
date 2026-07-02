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

#ifndef KRATOS_SURFACE_L2_NORM
#define KRATOS_SURFACE_L2_NORM

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
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
    ///@addtogroup ConvectionDiffusionApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// @brief Computes the L2 norm of a nodal variable over a set of 2D surface model parts.
    ///
    /// Unlike SurfaceInterpolatedL2Norm, no main 3D model part is required: the variable
    /// must already be set on the nodes of the surface model parts themselves. For each
    /// Gauss point of a 2D condition, the value is interpolated from the condition's own
    /// nodes via the shape functions. The area-averaged L2 norm over each surface is then
    /// L2 = sqrt( (1/A) * integral_S |var|^2 dS ).
    ///
    /// @tparam TValueType Either `double` (scalar variable) or `array_1d<double,3>` (vector
    /// variable). For a vector, the magnitude is used: |var|^2 = var . var.
    template <class TValueType>
    class SurfaceL2Norm
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of SurfaceL2Norm
        KRATOS_CLASS_POINTER_DEFINITION(SurfaceL2Norm);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        SurfaceL2Norm(
            std::vector<ModelPart *> model_part_vector,
            const Variable<TValueType> &variable)
            : mModelPartVector(model_part_vector), mVariable(variable), mIsComputed(false)
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
            buffer << "SurfaceL2Norm";
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const { rOStream << "SurfaceL2Norm"; }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const {}

        /// @brief Compute the area-averaged L2 norm of the variable in each surface model part.
        ///
        /// The variable must be registered as a nodal solution-step variable on every surface
        /// model part and its values must be set on all surface nodes before calling this method.
        void ComputeL2Norm()
        {
            for (unsigned m = 0; m < mModelPartVector.size(); m++)
            {
                ModelPart &r_model_part = *(mModelPartVector[m]);

                KRATOS_ERROR_IF_NOT(r_model_part.HasNodalSolutionStepVariable(mVariable))
                    << "Surface model part " << r_model_part.Name() << " does not have the variable "
                    << mVariable.Name() << " in its nodal solution-step data." << std::endl;

                using TwoReduction = CombinedReduction<SumReduction<double>, SumReduction<double>>;
                double integral, area;
                std::tie(integral, area) = block_for_each<TwoReduction>(
                    r_model_part.Conditions(),
                    [&](Condition &rCond) -> std::tuple<double, double> {
                        Geometry<Node> &r_geometry = rCond.GetGeometry();
                        unsigned int NumNodes = r_geometry.size();
                        GeometryData::IntegrationMethod integration_method = rCond.GetIntegrationMethod();
                        const std::vector<IntegrationPoint<3>> r_integration_points = r_geometry.IntegrationPoints(integration_method);
                        unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
                        Vector detJ_vector(r_number_integration_points);
                        r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
                        Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

                        double local_integral = 0.0;
                        double local_area = 0.0;
                        for (unsigned g = 0; g < r_number_integration_points; g++)
                        {
                            TValueType value = mVariable.Zero();
                            for (unsigned n = 0; n < NumNodes; n++)
                            {
                                value += r_geometry[n].FastGetSolutionStepValue(mVariable) * NContainer(g, n);
                            }
                            double Weight = r_integration_points[g].Weight() * detJ_vector[g];
                            local_integral += SquaredValue(value) * Weight;
                            local_area += Weight;
                        }
                        return {local_integral, local_area};
                    });

                mL2Norms[m] = std::sqrt(integral / area);

                // Print info
                std::cout << "\nModel part " << m + 1 << " (" << r_model_part.Name() << "):" << std::endl;
                std::cout << "    - surface area = " << area << std::endl;
                std::cout << "    - L2 norm = " << mL2Norms[m] << std::endl;
            }
            mIsComputed = true;
        }

        /// @brief Return the L2 norm of the variable in each surface model part.
        /// @return Vector with one L2 norm per surface model part.
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

        std::vector<ModelPart *> mModelPartVector;
        const Variable<TValueType> &mVariable;
        std::vector<double> mL2Norms;
        bool mIsComputed;

        ///@}
        ///@name Deleted special members
        ///@{

        /// Default constructor.
        SurfaceL2Norm() = delete;

        /// Assignment operator.
        SurfaceL2Norm &operator=(SurfaceL2Norm const &rOther) = delete;

        /// Copy constructor.
        SurfaceL2Norm(SurfaceL2Norm const &rOther) = delete;

        ///@}
        ///@name Private Operations
        ///@{

        /// @brief Squared magnitude of a scalar value.
        static double SquaredValue(double v) { return v * v; }

        /// @brief Squared magnitude of a vector value.
        static double SquaredValue(const array_1d<double, 3> &v) { return inner_prod(v, v); }

        ///@}
    }; // Class SurfaceL2Norm

    ///@}

}; // namespace Kratos.

#endif // KRATOS_SURFACE_L2_NORM
