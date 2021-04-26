// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef GEOMETRY_UTILITIES_H
#define GEOMETRY_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <functional>

#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/element_size_calculator.h"

#include "spatial_containers/spatial_containers.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

class GeometryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // For better reading
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    using GeometryType = ModelPart::ElementType::GeometryType;

    /// Pointer definition of GeometryUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GeometryUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryUtilities( ModelPart& modelPart )
        : mrModelPart( modelPart )
    {
    }

    /// Destructor.
    virtual ~GeometryUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void ComputeUnitSurfaceNormals()
    {
        KRATOS_TRY;

        const unsigned int domain_size = mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        KRATOS_ERROR_IF(mrModelPart.NumberOfConditions() == 0) <<
            "> Normal calculation requires surface or line conditions to be defined!" << std::endl;
        KRATOS_ERROR_IF((domain_size == 3 && mrModelPart.ConditionsBegin()->GetGeometry().size() == 2)) <<
            "> Normal calculation of 2-noded conditions in 3D domains is not possible!" << std::endl;
        CalculateAreaNormals(mrModelPart.Conditions(),domain_size);
        CalculateUnitNormals();

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ProjectNodalVariableOnUnitSurfaceNormals( const Variable<array_3d> &rNodalVariable )
    {
        KRATOS_TRY;

        ProjectNodalVariableOnDirection(rNodalVariable, NORMALIZED_SURFACE_NORMAL);

        KRATOS_CATCH("");
    }

    void ProjectNodalVariableOnDirection( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rDirectionVariable)
    {
        KRATOS_TRY;

        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            array_3d &node_normal = node_i->FastGetSolutionStepValue(rDirectionVariable);

            const double magnitude = inner_prod(nodal_variable, node_normal);
            noalias(nodal_variable) = magnitude * node_normal;
        }

        KRATOS_CATCH("");
    }

    void ProjectNodalVariableOnTangentPlane( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rPlaneNormalVariable)
    {
        KRATOS_TRY;

        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            array_3d &node_normal = node_i->FastGetSolutionStepValue(rPlaneNormalVariable);

            const double magnitude = inner_prod(nodal_variable, node_normal);
            nodal_variable -= magnitude * node_normal;
        }

        KRATOS_CATCH("");
    }
    // --------------------------------------------------------------------------
    void ExtractBoundaryNodes( std::string const& rBoundarySubModelPartName )
    {
    	KRATOS_TRY;

        ModelPart& r_boundary_model_part = mrModelPart.GetSubModelPart(rBoundarySubModelPartName);

        KRATOS_ERROR_IF(r_boundary_model_part.Nodes().size() != 0) << "ExtractBoundaryNodes: The boundary model part already has nodes!" << std::endl;

    	// Some type-definitions
        typedef std::unordered_map<vector<unsigned int>, unsigned int, KeyHasherRange<vector<unsigned int>>, KeyComparorRange<vector<unsigned int>> > hashmap;

    	// Create map to ask for number of faces for the given set of node ids representing one face in the model part
    	hashmap n_boundaries_map;

        unsigned int domain_size = static_cast<unsigned int>(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);

    	// Fill map that counts number of faces for given set of nodes
    	for (auto& elem_i : mrModelPart.Elements())
    	{
            KRATOS_ERROR_IF(elem_i.GetGeometry().Dimension() < domain_size) << "ExtractBoundaryNodes: This function does only work"
                <<" for solid elements in 3D and surface elements in 2D!" << std::endl;

            Element::GeometryType::GeometriesArrayType boundaries = elem_i.GetGeometry().GenerateBoundariesEntities();

            for(unsigned int boundary=0; boundary<boundaries.size(); boundary++)
            {
                // Create vector that stores all node is of current face
                DenseVector<unsigned int> ids(boundaries[boundary].size());

                // Store node ids
                for(unsigned int i=0; i<boundaries[boundary].size(); i++)
                    ids[i] = boundaries[boundary][i].Id();

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(ids.begin(), ids.end());

                // Fill the map
                n_boundaries_map[ids] += 1;
            }
    	}

    	// Vector to store all nodes on surface. Node ids may be listed several times
    	std::vector<std::size_t> temp_boundary_node_ids;

    	// Add surface nodes to sub-model part
    	for(auto it=n_boundaries_map.begin(); it!=n_boundaries_map.end(); it++)
    	{
    		// If given node set represents face that is not overlapping with a face of another element, add it as skin element
    		if(it->second == 1)
    		{
    			for(unsigned int i=0; i<it->first.size(); i++)
    				temp_boundary_node_ids.push_back(it->first[i]);
    		}
    	}

    	// Add nodes and remove double entries
    	r_boundary_model_part.AddNodes(temp_boundary_node_ids);

    	KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ComputeDistancesToBoundingModelPart(
        ModelPart& rBoundingModelPart,
        pybind11::list& rSignedDistances,
        pybind11::list& rDirections )
    {
        KRATOS_TRY;

        typedef Node < 3 > NodeType;
        typedef NodeType::Pointer NodeTypePointer;
        typedef std::vector<NodeType::Pointer> NodeVector;
        typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
        typedef std::vector<double>::iterator DoubleVectorIterator;
        typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > KDTree;

        KRATOS_ERROR_IF(rBoundingModelPart.NumberOfElements() != 0) <<
            "ComputeDistancesToBoundingModelPart: Model part must only contain conditions!" << std::endl;

        NodeVector all_bounding_nodes;
        all_bounding_nodes.reserve(rBoundingModelPart.Nodes().size());
        for (ModelPart::NodesContainerType::iterator node_it = rBoundingModelPart.NodesBegin(); node_it != rBoundingModelPart.NodesEnd(); ++node_it)
        {
            all_bounding_nodes.push_back(*(node_it.base()));
        }
        const size_t bucket_size = 100;
        KDTree search_tree(all_bounding_nodes.begin(), all_bounding_nodes.end(), bucket_size);

        GeometryUtilities(rBoundingModelPart).ComputeUnitSurfaceNormals();

        for (auto& r_node : mrModelPart.Nodes()){

            double distance;
            NodeTypePointer p_neighbor = search_tree.SearchNearestPoint(r_node, distance);

            const array_3d delta = r_node.Coordinates() - p_neighbor->Coordinates();
            const array_3d& bounding_normal = p_neighbor->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
            const double projected_length = inner_prod(delta, bounding_normal);

            rSignedDistances.append(projected_length);

            rDirections.append(bounding_normal[0]);
            rDirections.append(bounding_normal[1]);
            rDirections.append(bounding_normal[2]);
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    double ComputeVolume()
    {
        KRATOS_TRY

        std::function<double(const GeometryType&)> volume_method;

        const int domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        int number_of_nodes = -1;
        if (mrModelPart.NumberOfElements() > 0) {
            number_of_nodes = mrModelPart.Elements().front().GetGeometry().PointsNumber();

            block_for_each(mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
                KRATOS_ERROR_IF(static_cast<int>(rElement.GetGeometry().PointsNumber()) != number_of_nodes)
                    << "Mismatching geometries with different number of nodes "
                       "found in "
                    << mrModelPart.Name() << ". All geometries should be of same type. [ First element's geometry's nodes = "
                    << number_of_nodes << ", current_element_id = " << rElement.Id()
                    << ", current element's geometry's number of nodes = "
                    << rElement.GetGeometry().PointsNumber() << " ].\n";
            });
        }
        number_of_nodes = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(number_of_nodes);

        KRATOS_ERROR_IF(number_of_nodes == -1)
            << "No elements were found in " << mrModelPart.Name()
            << ". Please use a model part with elements.\n";

        if (domain_size == 2) {
            if (number_of_nodes == 3) {
                volume_method = [](const GeometryType& rGeometry) {
                    return std::pow(
                        ElementSizeCalculator<2, 3>::AverageElementSize(rGeometry), 2);
                };
            } else if (number_of_nodes == 4) {
                volume_method = [](const GeometryType& rGeometry) {
                    return std::pow(
                        ElementSizeCalculator<2, 4>::AverageElementSize(rGeometry), 2);
                };
            } else {
                KRATOS_ERROR
                    << "Unsupported geometry type. Only geometries with 3 and "
                       "4 nodes are supported in 2D. [ number_of_nodes = "
                    << number_of_nodes << " ].\n";
            }
        } else if (domain_size == 3) {
            if (number_of_nodes == 4) {
                volume_method = [](const GeometryType& rGeometry) {
                    return std::pow(
                        ElementSizeCalculator<3, 4>::AverageElementSize(rGeometry), 3);
                };
            } else if (number_of_nodes == 6) {
                volume_method = [](const GeometryType& rGeometry) {
                    return std::pow(
                        ElementSizeCalculator<3, 6>::AverageElementSize(rGeometry), 3);
                };
            } else if (number_of_nodes == 8) {
                volume_method = [](const GeometryType& rGeometry) {
                    return std::pow(
                        ElementSizeCalculator<3, 8>::AverageElementSize(rGeometry), 3);
                };
            } else {
                KRATOS_ERROR
                    << "Unsupported geometry type. Only geometries with 4, 6 and "
                       "8 nodes are supported in 3D. [ number_of_nodes = "
                    << number_of_nodes << " ].\n";
            }
        } else {
            KRATOS_ERROR << "Unsupported domain size. Only 2D and 3D is "
                            "supported [ domain_size = "
                         << domain_size << " ].\n";
        }

        const double volume = block_for_each<SumReduction<double>>(
            mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
                return volume_method(rElement.GetGeometry());
            });

        return mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(volume);

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ComputeVolumeShapeDerivatives(
        const Variable<array_3d>& rDerivativeVariable)
    {
        KRATOS_TRY

        using CUInt = const unsigned int;

        VariableUtils().SetNonHistoricalVariableToZero(rDerivativeVariable, mrModelPart.Nodes());

        const int domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        int number_of_nodes = -1;
        if (mrModelPart.NumberOfElements() > 0) {
            number_of_nodes = mrModelPart.Elements().front().GetGeometry().PointsNumber();

            block_for_each(mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
                KRATOS_ERROR_IF(static_cast<int>(rElement.GetGeometry().PointsNumber()) != number_of_nodes)
                    << "Mismatching geometries with different number of nodes "
                       "found in "
                    << mrModelPart.Name() << ". All geometries should be of same type. [ First element's geometry's nodes = "
                    << number_of_nodes << ", current_element_id = " << rElement.Id()
                    << ", current element's geometry's number of nodes = "
                    << rElement.GetGeometry().PointsNumber() << " ].\n";
            });
        }
        number_of_nodes = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(number_of_nodes);

        KRATOS_ERROR_IF(number_of_nodes == -1)
            << "No elements were found in " << mrModelPart.Name()
            << ". Please use a model part with elements.\n";

        std::function<double(CUInt, CUInt, const GeometryType&)> volume_derivative_method;

        if (domain_size == 2) {
            if (number_of_nodes == 3) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 2.0 * ElementSizeCalculator<2, 3>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 3>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else if (number_of_nodes == 4) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 2.0 * ElementSizeCalculator<2, 4>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else {
                KRATOS_ERROR
                    << "Unsupported geometry type. Only geometries with 3 and "
                       "4 nodes are supported in 2D. [ number_of_nodes = "
                    << number_of_nodes << " ].\n";
            }
        } else if (domain_size == 3) {
            if (number_of_nodes == 4) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 4>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else if (number_of_nodes == 6) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 6>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 6>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else if (number_of_nodes == 8) {
                volume_derivative_method = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 8>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 8>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
            } else {
                KRATOS_ERROR
                    << "Unsupported geometry type. Only geometries with 4, 6 and "
                       "8 nodes are supported in 3D. [ number_of_nodes = "
                    << number_of_nodes << " ].\n";
            }
        } else {
            KRATOS_ERROR << "Unsupported domain size. Only 2D and 3D is "
                            "supported [ domain_size = "
                         << domain_size << " ].\n";
        }

        block_for_each(mrModelPart.Elements(), [&](ModelPart::ElementType& rElement){
            auto& r_geometry = rElement.GetGeometry();
            for (int c = 0; c < static_cast<int>(r_geometry.PointsNumber()); ++c) {
                auto& r_node = r_geometry[c];

                for (int k = 0; k < domain_size; ++k) {
                    const double derivative_value = volume_derivative_method(c, k, r_geometry);

                    r_node.SetLock();
                    auto& r_derivative_value = r_node.GetValue(rDerivativeVariable);
                    r_derivative_value[k] += derivative_value;
                    r_node.UnSetLock();
                }
            }
        });

        mrModelPart.GetCommunicator().AssembleNonHistoricalData(rDerivativeVariable);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "GeometryUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CalculateAreaNormals(ConditionsArrayType& rConditions, int dimension)
    {
        KRATOS_TRY

        //resetting the normals
        array_1d<double,3> zero = Vector(3);
        noalias(zero) = ZeroVector(3);

        for(auto & cond_i : rConditions)
        {
            Element::GeometryType& rNodes = cond_i.GetGeometry();
            for(unsigned int in = 0; in<rNodes.size(); in++)
                noalias((rNodes[in]).GetSolutionStepValue(NORMAL)) = zero;
        }

        //calculating the normals and storing on the conditions
        array_1d<double,3> An;
        if(dimension == 2)
        {
            for(auto & cond_i : rConditions)
            {
                if (cond_i.GetGeometry().PointsNumber() == 2)
                    CalculateNormal2D(cond_i,An);
            }
        }
        else if(dimension == 3)
        {
            array_1d<double,3> v1;
            array_1d<double,3> v2;
            for(auto & cond_i : rConditions)
            {
                //calculate the normal on the given condition
                if (cond_i.GetGeometry().PointsNumber() == 3)
                    CalculateNormal3DTriangle(cond_i,An,v1,v2);
                else if (cond_i.GetGeometry().PointsNumber() == 4)
                    CalculateNormal3DQuad(cond_i,An,v1,v2);
                else
                    KRATOS_ERROR << "Calculation of surface normal not implemented for the given surface conditions!";
            }
        }

        //adding the normals to the nodes
        for(auto & cond_i : rConditions)
        {
            Geometry<Node<3> >& pGeometry = cond_i.GetGeometry();
            double coeff = 1.00/pGeometry.size();
	        const array_1d<double,3>& normal = cond_i.GetValue(NORMAL);
            for(unsigned int i = 0; i<pGeometry.size(); i++)
                noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * normal;
        }

        KRATOS_CATCH("")
    }

    // --------------------------------------------------------------------------
    static void CalculateNormal2D(Condition& cond, array_1d<double,3>& An)
    {
        Geometry<Node<3> >& pGeometry = cond.GetGeometry();

        An[0] =    pGeometry[1].Y() - pGeometry[0].Y();
        An[1] = - (pGeometry[1].X() - pGeometry[0].X());
        An[2] =    0.00;

        array_1d<double,3>& normal = cond.GetValue(NORMAL);
        noalias(normal) = An;
    }

    // --------------------------------------------------------------------------
    static void CalculateNormal3DTriangle(Condition& cond, array_1d<double,3>& An, array_1d<double,3>& v1,array_1d<double,3>& v2 )
    {
        Geometry<Node<3> >& pGeometry = cond.GetGeometry();

        v1[0] = pGeometry[1].X() - pGeometry[0].X();
        v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
        v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

        v2[0] = pGeometry[2].X() - pGeometry[0].X();
        v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
        v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

        MathUtils<double>::CrossProduct(An,v1,v2);
        An *= 0.5;

        array_1d<double,3>& normal = cond.GetValue(NORMAL);
        noalias(normal) = An;
    }

    // --------------------------------------------------------------------------
    static void CalculateNormal3DQuad(Condition& cond, array_1d<double,3>& An, array_1d<double,3>& v1,array_1d<double,3>& v2 )
    {
        Geometry<Node<3> >& pGeometry = cond.GetGeometry();

        v1[0] = pGeometry[2].X() - pGeometry[0].X();
        v1[1] = pGeometry[2].Y() - pGeometry[0].Y();
        v1[2] = pGeometry[2].Z() - pGeometry[0].Z();

        v2[0] = pGeometry[3].X() - pGeometry[1].X();
        v2[1] = pGeometry[3].Y() - pGeometry[1].Y();
        v2[2] = pGeometry[3].Z() - pGeometry[1].Z();

        MathUtils<double>::CrossProduct(An,v1,v2);
        An *= 0.5;

        array_1d<double,3>& normal = cond.GetValue(NORMAL);
        noalias(normal) = An;
    }

    // --------------------------------------------------------------------------
    void CalculateUnitNormals()
    {
        for (auto& node_i : mrModelPart.Nodes())
        {
            const array_1d<double,3>& area_normal = node_i.FastGetSolutionStepValue(NORMAL);
            array_3d& normalized_normal = node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

            const double norm2 = norm_2(area_normal);
            KRATOS_ERROR_IF(norm2<1e-10) << "CalculateUnitNormals: Norm2 of normal for node "
                << node_i.Id() << " is < 1e-10!" << std::endl;

            noalias(normalized_normal) = area_normal/norm2;
        }
    }

    // --------------------------------------------------------------------------


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      GeometryUtilities& operator=(GeometryUtilities const& rOther);

    /// Copy constructor.
//      GeometryUtilities(GeometryUtilities const& rOther);


    ///@}

}; // Class GeometryUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // GEOMETRY_UTILITIES_H
