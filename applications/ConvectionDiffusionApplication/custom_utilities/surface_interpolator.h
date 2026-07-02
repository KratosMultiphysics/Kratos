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

#ifndef KRATOS_SURFACE_INTERPOLATOR
#define KRATOS_SURFACE_INTERPOLATOR

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "includes/node.h"
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/binbased_fast_point_locator.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
    ///@addtogroup ConvectionDiffusionApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// @brief Interpolates a variable from a general 3D model part onto the nodes of a set of
    /// "surface" model parts.
    ///
    /// The surface model parts only contain 2D conditions whose nodes do not carry nodal
    /// solution values (empty geometries). For each node of a surface model part, the 3D
    /// element of the main model part containing that node is located, the variable is
    /// interpolated there, and the interpolated value is assigned to the surface node.
    ///
    /// @tparam TValueType Either `double` (scalar variable) or `array_1d<double,3>` (vector
    /// variable).
    template <class TValueType>
    class SurfaceInterpolator
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of SurfaceInterpolator
        KRATOS_CLASS_POINTER_DEFINITION(SurfaceInterpolator);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        SurfaceInterpolator(
            ModelPart &model_part,
            std::vector<ModelPart *> model_part_vector,
            const Variable<TValueType> &interpolation_variable)
            : mMainModelPart(model_part), mModelPartVector(model_part_vector), mVariable(interpolation_variable)
        {
        }

        ///@}
        ///@name Operations
        ///@{

        /// Turn back information as a string.
        std::string Info() const
        {
            std::stringstream buffer;
            buffer << "SurfaceInterpolator";
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const { rOStream << "SurfaceInterpolator"; }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const {}

        /// @brief Interpolate the variable from the main model part onto the nodes of each
        /// surface model part and assign the values to those nodes.
        void Interpolate()
        {
            // The variable is read from the main model part nodes, so it must be part of their
            // solution-step data to be accessed through FastGetSolutionStepValue.
            KRATOS_ERROR_IF_NOT(mMainModelPart.HasNodalSolutionStepVariable(mVariable))
                << "The main model part " << mMainModelPart.Name() << " does not have the variable "
                << mVariable.Name() << " in its nodal solution-step data." << std::endl;

            // The interpolated value is assigned to the surface nodes, so they must also have the
            // variable in their solution-step data.
            for (unsigned m = 0; m < mModelPartVector.size(); m++)
            {
                KRATOS_ERROR_IF_NOT(mModelPartVector[m]->HasNodalSolutionStepVariable(mVariable))
                    << "The surface model part " << mModelPartVector[m]->Name() << " does not have the variable "
                    << mVariable.Name() << " in its nodal solution-step data, so it cannot be assigned "
                    << "through FastGetSolutionStepValue." << std::endl;
            }

            BinBasedFastPointLocator<3> point_locator(mMainModelPart);
            point_locator.UpdateSearchDatabase();

            using ResultContainerType = BinBasedFastPointLocator<3>::ResultContainerType;

            for (unsigned m = 0; m < mModelPartVector.size(); m++)
            {
                ModelPart &r_model_part = *(mModelPartVector[m]);

                block_for_each(r_model_part.Nodes(), ResultContainerType(1000),
                    [&](Node &rNode, ResultContainerType &rResults) {
                        Vector shape_functions;
                        Element::Pointer p_elem;
                        const bool is_found = point_locator.FindPointOnMesh(
                            rNode.Coordinates(), shape_functions, p_elem, rResults.begin(), rResults.size(), 1e-6);
                        KRATOS_ERROR_IF_NOT(is_found)
                            << "Node " << rNode.Id() << " at " << rNode.Coordinates()
                            << " not found in model part " << mMainModelPart.Name() << std::endl;

                        Geometry<Node> &r_geometry = p_elem->GetGeometry();
                        TValueType value = mVariable.Zero();
                        for (unsigned n = 0; n < r_geometry.size(); n++)
                        {
                            value += r_geometry[n].FastGetSolutionStepValue(mVariable) * shape_functions[n];
                        }
                        rNode.FastGetSolutionStepValue(mVariable) = value;
                    });
            }
        }

        ///@}

    private:
        ///@name Member Variables
        ///@{

        ModelPart &mMainModelPart;
        std::vector<ModelPart *> mModelPartVector;
        const Variable<TValueType> &mVariable;

        ///@}
        ///@name Deleted special members
        ///@{

        /// Default constructor.
        SurfaceInterpolator() = delete;

        /// Assignment operator.
        SurfaceInterpolator &operator=(SurfaceInterpolator const &rOther) = delete;

        /// Copy constructor.
        SurfaceInterpolator(SurfaceInterpolator const &rOther) = delete;

        ///@}
        ///@}
    }; // Class SurfaceInterpolator

    ///@}

}; // namespace Kratos.

#endif // KRATOS_SURFACE_INTERPOLATOR
