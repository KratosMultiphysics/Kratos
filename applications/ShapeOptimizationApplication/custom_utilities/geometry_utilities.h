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

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "shape_optimization_application.h"

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

        // We loop over all nodes and compute the part of the sensitivity which is in direction to the surface normal
        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            array_3d &node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

            // We compute dFdX_n = (dFdX \cdot n) * n
            double surface_sens = inner_prod(nodal_variable,node_normal);
            nodal_variable = surface_sens * node_normal;
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

            Element::GeometryType::GeometriesArrayType boundaries;
            if (domain_size==3)
                boundaries = elem_i.GetGeometry().Faces();
            else if (domain_size == 2)
                boundaries = elem_i.GetGeometry().Edges();

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
                    KRATOS_ERROR << "> Calculation of surface normal not implemented for the given surface conditions!";
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
            noalias(normalized_normal) = area_normal/norm_2(area_normal);
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
