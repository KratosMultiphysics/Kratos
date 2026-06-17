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
#include "spatial_containers/spatial_containers.h"
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

            for (unsigned m = 0; m < mModelPartVector.size(); m++)
            {
                ModelPart &r_model_part = *(mModelPartVector[m]);
                const unsigned number_of_nodes = r_model_part.NumberOfNodes();

                for (unsigned i = 0; i < number_of_nodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it_node = r_model_part.NodesBegin() + i;

                    // Global position of this surface node
                    const array_1d<double, 3> &node_global = it_node->Coordinates();

                    // Find the 3D element containing this node
                    array_1d<double, 3> node_local;
                    ModelPart::ElementType::Pointer p_elem;
                    FindElementContainingPoint(node_global, node_local, p_elem);
                    KRATOS_ERROR_IF(p_elem == nullptr)
                        << "Node " << it_node->Id() << " at " << node_global << " not found in model part " << mMainModelPart.Name() << std::endl;

                    // Interpolate the variable and assign it to the surface node
                    it_node->FastGetSolutionStepValue(mVariable) = InterpolateValue(p_elem, node_local);
                }
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
        ///@name Private Operations
        ///@{

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
                if (r_geometry.IsInside(point, p_pos_local, 1e-6))
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
    }; // Class SurfaceInterpolator

    ///@}

}; // namespace Kratos.

#endif // KRATOS_SURFACE_INTERPOLATOR
