#ifndef BREP_MODEL_GEOMETRY_READER
#define BREP_MODEL_GEOMETRY_READER

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "BrepModel.h"
#include "ControlPoint.h"
#include "TrimmingCurve.h"
#include "BrepModel.h"
#include "Face.h"
#include "Edge.h"
#include "FaceTrim.h"
//#include "brep_gauss_point.h"

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

class BrepModelGeometryReader
{
public:
	///@name Type Definitions
	///@{
	typedef std::vector<BrepModel> BrepModelVector;
    typedef boost::python::extract<double> extractDouble;
    typedef boost::python::extract<int> extractInt;
	typedef boost::python::extract<std::string> extractString;

	typedef std::vector<double> DoubleVector;
	typedef std::vector<int> IntVector;

	//Face:
	typedef std::vector<std::vector<int>> TrimmingLoopVector;
	typedef std::vector<TrimmingCurve> TrimmingCurveVector;

	//BrepModel:
	typedef std::vector<Face> FacesVector;
	typedef std::vector<Edge> EdgesVector;

	typedef std::vector<ControlPoint> ControlPointVector;
	///@}
	

	/// Pointer definition of BrepModelGeometryReader
	//    KRATOS_CLASS_POINTER_DEFINITION[BrepModelGeometryReader];

	// Default constructor.
	BrepModelGeometryReader() 
	{	
	}

	// Constructor
	BrepModelGeometryReader(boost::python::dict cad_geometry_in_json) 
	: m_cad_geometry_in_json(cad_geometry_in_json)
	{	
	}

	/// Destructor.
	virtual ~BrepModelGeometryReader()
	{
	}

	// --------------------------------------------------------------------------
	void ReadGeometry(BrepModelVector& r_brep_model_vector, ModelPart& model_part)
	{
		std::cout << "\n> Start reading CAD geometry..." << std::endl;

		FacesVector faces_vector;
		EdgesVector edges_vector;

		// loop over faces
		for (int i = 0; i < boost::python::len(m_cad_geometry_in_json["faces"]); i++)
		{
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// 1. Step: faces
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			unsigned int face_id = extractInt(m_cad_geometry_in_json["faces"][i]["brep_id"]);

			std::cout << "> Reading face " << face_id << "..." << std::endl;

			// Variables needed later
			DoubleVector knot_vector_u;
			DoubleVector knot_vector_v;
			int p;
			int q;
			ControlPointVector control_points;

			// read and store knot_vector_u
			for (int u_idx = 0; u_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][0]); u_idx++)
			{
				double knot = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][0][u_idx]);
				knot_vector_u.push_back(knot);
			}

			// read and store knot_vector_v
			for (int v_idx = 0; v_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][1]); v_idx++)
			{
				double knot = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][1][v_idx]);
				knot_vector_v.push_back(knot);
			}

			// read and store polynamial degree p and q
			p = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["degrees"][0]);
			q = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["degrees"][1]);

			// read and store control_points
			// Control points in each patch get a global as well as a mapping matrix id
			// brep Id: Unique Id for each control point (given by json-file)
			// mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
			for (int cp_idx = 0; cp_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"]); cp_idx++)
			{
				unsigned int cp_id = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][0]);
				double x = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][0]);
				double y = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][1]);
				double z = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][2]);
				double w = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][3]);

				ControlPoint new_control_point(x, y, z, w, cp_id);
				//model_part.CreateNewNode(cp_id, x, y, z);
				control_points.push_back(new_control_point);
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// 2. step: boundary loops
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			TrimmingLoopVector trimming_loop;
			TrimmingCurveVector trimming_curves;

			// For better reading
			boost::python::list boundary_dict(m_cad_geometry_in_json["faces"][i]["boundary_loops"]);

			// Loop over all boundary loops
			for (int loop_idx = 0; loop_idx < boost::python::len(boundary_dict); loop_idx++)
			{
				IntVector loop;
				

				// Loop over all curves
				for (int edge_idx = 0; edge_idx < boost::python::len(boundary_dict[loop_idx]["trimming_curves"]); edge_idx++)
				{
					unsigned curve_id = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"]);

					loop.push_back(curve_id);

					// Variables needed later
					DoubleVector boundary_knot_vector_u;
					unsigned int boundary_p;
					ControlPointVector boundary_control_points;

					// read and store knot_vector_u
					for (int u_idx = 0; u_idx < boost::python::len(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"]); u_idx++)
					{
						double knot = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"][u_idx]);
						boundary_knot_vector_u.push_back(knot);
					}

					// read and store polynamial degree p and q
					boundary_p = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degrees"]);	

					// read and store control_points
					for (int cp_idx = 0; cp_idx < boost::python::len(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"]); cp_idx++)
					{
						unsigned int cp_id = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
						double x = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][0]);
						double y = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][1]);
						double z = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][2]);
						double w = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][3]);			

						ControlPoint new_control_point(x, y, z, w, cp_id);
						boundary_control_points.push_back(new_control_point);
					}							

					// Create and store edge
					TrimmingCurve new_boundary_curve(curve_id, boundary_knot_vector_u, boundary_p, boundary_control_points);
					trimming_curves.push_back(new_boundary_curve);
				}

				// Read loop type
				//std::string loop_type = extractString(boundary_dict[loop_idx]["loop_type"]);
				//std::string Inner("Inner");
				//bool is_inner_loop = false;
				//if(loop_type.compare(Inner) == 0)
				//	is_inner_loop = true;

				// Create and store boundary loop
				trimming_loop.push_back(loop);
			}		

			// create face
			Face face(face_id, trimming_curves, trimming_loop, knot_vector_u, knot_vector_v, p, q, control_points);
			faces_vector.push_back(face);
		}	

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 3. Step: Create BrepModel
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		BrepModel brep(faces_vector, edges_vector);
		r_brep_model_vector.push_back(brep);

		std::cout << "\n> Finished reading CAD geometry..." << std::endl;	
	}

	//// --------------------------------------------------------------------------
	//void ReadIntegrationData(BREPElementVector& r_brep_elements)
	//{
	//	std::cout << "\n> Starting reading CAD integration data of given json-file..." << std::endl;

	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// 1. Step: Loop over all 2D elements defined in the integration file to assign every element_id the 
	//	// 			corresponding patch_id. This is needed to store the patch id with each Gauss point that is read later.
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	std::map<unsigned int, unsigned int> corresponding_patch_id;

	//	// Loop over the patches
	//	for (unsigned int patch_itr = 0; patch_itr < boost::python::len(mr_cad_integration_data_in_json["2d_elements"]); patch_itr++)
	//	{		
	//		unsigned int patch_id = extractInt(mr_cad_integration_data_in_json["2d_elements"][patch_itr][0]);	

	//		// Loop over all 2D_elements
	//		for (unsigned int elem_itr = 0; elem_itr < boost::python::len(mr_cad_integration_data_in_json["2d_elements"][patch_itr][1]); elem_itr++)
	//		{
	//			// Read element Id specified in integration file
	//			unsigned int elem_id = extractInt(mr_cad_integration_data_in_json["2d_elements"][patch_itr][1][elem_itr][0]);

	//			// Create map between elem_id and corresponding patch_id
	//			corresponding_patch_id[elem_id] = patch_id;
	//		}
	//	}

	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// 2. Step: Read and store brep elements with corresponding Gauss points
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////		

	//	// Loop over all boundary edges
	//	for (unsigned int edge_itr = 0; edge_itr < boost::python::len(mr_cad_integration_data_in_json["brep_elements"]); edge_itr++)
	//	{
	//		// Extract Id of brep edge
	//		unsigned int brep_edge_id = extractInt(mr_cad_integration_data_in_json["brep_elements"][edge_itr][0]);

	//		// loop over all brep elements on boundary edge
	//		for (unsigned int elem_itr = 0; elem_itr < boost::python::len(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1]); elem_itr++)
	//		{
	//			BREPGaussPointVector gauss_points;

	//			// Variables needed to distinguish different types of boundary elements
	//			bool has_coupling_condition = false;
	//			bool has_dirichlet_condition = false;

	//			// Declare variables for Gauss point;
	//			unsigned int master_patch_id;
	//			unsigned int gauss_point_id;				
	//			unsigned int master_elem_id;
	//			double weight;
	//			Vector location = ZeroVector(2);
	//			Vector tangent = ZeroVector(2);

	//			// Check if current edge only carries information about a master element or also about a corresponding slave element (as only needed for coupling)
	//			// To this end we check if the first gauss point of this element has master-slave information 
	//			if(boost::python::len(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][0][0])==2)
	//			{
	//				// if above statement is true, than we are given a coupling condition

	//				// Additional variables in case of a coupling edge
	//				unsigned int slave_patch_id;
	//				unsigned int slave_elem_id;
	//				Vector location_slave = ZeroVector(2);
	//				Vector tangent_slave = ZeroVector(2);

	//				// Loop over all Gauss points on the current brep element
	//				for (unsigned int gp_itr = 0; gp_itr < boost::python::len(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1]); gp_itr++)
	//				{
	//					// Read integration data from file
	//					master_elem_id = extractInt(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][0][0]);
	//					slave_elem_id = extractInt(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][0][1]);
	//					gauss_point_id = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][0]);
	//					weight = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][1]);
	//					location(0) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][2][0]);
	//					location(1) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][2][1]);
	//					tangent(0) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][3][0]);
	//					tangent(1) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][3][1]);	
	//					location_slave(0) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][4][0]);
	//					location_slave(1) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][4][1]);					
	//					tangent_slave(0) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][5][0]);
	//					tangent_slave(1) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][5][1]);	

	//					// identify patch Id corresponding to above read element id (each for master and slave element)
	//					master_patch_id = corresponding_patch_id[master_elem_id];
	//					slave_patch_id = corresponding_patch_id[slave_elem_id];			

	//					// Create new Gauss point on brep element
	//					BREPGaussPoint new_brep_gp(master_patch_id, slave_patch_id, gauss_point_id, weight, location, tangent, location_slave, tangent_slave);

	//					// Add Gauss point to list of all brep Gauss points
	//					gauss_points.push_back(new_brep_gp);
	//				}	

	//				// Identify edge element as part of a coupling edge
	//				has_coupling_condition = true;
	//			}
	//			else // we are given a Dirichlet condition
	//			{
	//				// Loop over all Gauss points on the current brep element
	//				for (unsigned int gp_itr = 0; gp_itr < boost::python::len(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1]); gp_itr++)
	//				{
	//					// Read integration data from file
	//					master_elem_id = extractInt(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][0][0]);
	//					gauss_point_id = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][0]);
	//					weight = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][1]);
	//					location(0) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][2][0]);
	//					location(1) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][2][1]);
	//					tangent(0) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][3][0]);
	//					tangent(1) = extractDouble(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][1][gp_itr][1][3][1]);	

	//					// identify patch Id corresponding to above read element id (each for master and slave element)
	//					master_patch_id = corresponding_patch_id[master_elem_id];	

	//					// Create new Gauss point on brep element
	//					BREPGaussPoint new_brep_gp(master_patch_id, gauss_point_id, weight, location, tangent);

	//					// Add Gauss point to list of all brep Gauss points
	//					gauss_points.push_back(new_brep_gp);
	//				}

	//				// Identify edge element as part of a Dirichlet edge
	//				has_dirichlet_condition = true;
	//			}

	//			// Create new brep element and add to list of all brep elements
	//			unsigned int brep_elem_id = extractInt(mr_cad_integration_data_in_json["brep_elements"][edge_itr][1][elem_itr][0]);
	//			BREPElement new_brep_element(brep_elem_id, brep_edge_id, gauss_points,has_coupling_condition, has_dirichlet_condition);
	//			r_brep_elements.push_back(new_brep_element);
	//		}
	//	}
	//	std::cout << "\n> Finished reading CAD integration data of given json-file..." << std::endl;
	//}

	// --------------------------------------------------------------------------
	//void UpdateControlPoints(PatchVector& patches, Vector& ds)
	//{
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// 1. Step: Update C++ data base
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	// Map to identify control point for given global id (needed for python update later)
	//	std::map<unsigned int, ControlPoint*> control_point_corresponding_to_global_id;

	//	for (PatchVector::iterator patch_i = patches.begin(); patch_i != patches.end(); ++patch_i)
	//	{
	//		for (ControlPointVector::iterator cp_i = patch_i->GetSurface().GetControlPoints().begin(); cp_i != patch_i->GetSurface().GetControlPoints().end(); ++cp_i)
	//		{
	//			if(cp_i->IsRelevantForMapping())
	//			{
	//				// Updating c++ data base
	//				unsigned int cp_mapping_matrix_id = cp_i->GetMappingMatrixId();
	//				cp_i->setdX( ds[3*cp_mapping_matrix_id+0] );
	//				cp_i->setdY( ds[3*cp_mapping_matrix_id+1] );
	//				cp_i->setdZ( ds[3*cp_mapping_matrix_id+2] );

	//			}
	//			// Filling map to be used later
	//			unsigned int cp_global_id = cp_i->getGlobalId();
	//			control_point_corresponding_to_global_id[cp_global_id] = &(*cp_i);
	//		}
	//	}

	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// 2. Step: Update pyhon data base
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	
	//	// loop over patches / faces in cad geometry file
	//	for (int i = 0; i < boost::python::len(m_cad_geometry_in_json["faces"]); i++)
	//	{
	//		for (int cp_idx = 0; cp_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"]); cp_idx++)
	//		{
	//			unsigned int global_id = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][0]);

	//			ControlPoint* cp_j = control_point_corresponding_to_global_id[global_id];

	//			double sx = cp_j->getX();
	//			double sy = cp_j->getY();
	//			double sz = cp_j->getZ();

	//			// Update python data base which is store as a reference (so it is also updated on python level)
	//			m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][0] = sx;
	//			m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][1] = sy;
	//			m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][2] = sz;
	//		}
	//	}
	//}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BrepModelGeometryReader";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BrepModelGeometryReader";
	}

	// ==============================================================================
	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}


private:
	// ==============================================================================
	// Initialized by class constructor
	// ==============================================================================
	boost::python::dict m_cad_geometry_in_json;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BrepModelGeometryReader& operator=[BrepModelGeometryReader const& rOther];

	/// Copy constructor.
	//      BrepModelGeometryReader[BrepModelGeometryReader const& rOther];

}; // Class BrepModelGeometryReader

} // namespace Kratos.

#endif // BREP_MODEL_GEOMETRY_READER
