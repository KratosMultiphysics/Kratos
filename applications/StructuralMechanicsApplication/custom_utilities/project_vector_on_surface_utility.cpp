//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus Sautter
//                   Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "project_vector_on_surface_utility.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

void ProjectVectorOnSurfaceUtility::Execute(ModelPart& rModelPart, const Parameters MethodParameters)
{
	const auto specific_parameters = MethodParameters["method_specific_settings"];
	const std::string& r_projection_type = specific_parameters["projection_type"].GetString();
	const std::string& r_local_variable_name = specific_parameters["local_variable_name"].GetString();

	std::cout << std::endl
			  << "Assigning " << r_local_variable_name << " orientation to elements using method: " << r_projection_type << std::endl;

	if (r_projection_type=="planar") PlanarProjection(rModelPart,MethodParameters);
	else if (r_projection_type=="radial") RadialProjection(rModelPart,MethodParameters);
	else if (r_projection_type=="spherical") SphericalProjection(rModelPart,MethodParameters);
	else
	{
		KRATOS_ERROR << "projection type: " << r_projection_type << " not available, please use planar,radial,spherical" << std::endl;
	}

	std::cout << std::endl
			  << ".......done assigning direction for all elements......." << std::endl;
}

void ProjectVectorOnSurfaceUtility::PlanarProjection(ModelPart& rModelPart, const Parameters& MethodParameters)
{
	const auto specific_parameters = MethodParameters["method_specific_settings"];
	const std::string& r_local_variable_name = specific_parameters["local_variable_name"].GetString();

	KRATOS_ERROR_IF_NOT(KratosComponents<ArrayVariableType>::Has(r_local_variable_name)) << "Variable " << r_local_variable_name << " not known" << std::endl;
	const ArrayVariableType& r_variable = KratosComponents<ArrayVariableType>::Get(r_local_variable_name);

	Vector3 global_vector;

	//global_vector is the global fiber direction given by user
	//read this global direction
	CheckAndReadVectors(specific_parameters, "global_fiber_direction", global_vector);

	//Normalize global_vector
	global_vector /=  norm_2(global_vector);

	const auto& r_process_info = rModelPart.GetProcessInfo();

	// Declare working variables
	Matrix local_coordinate_orientation;


	// Loop over all elements in part
	for (auto &element : rModelPart.Elements())
	{

		// get local axis in cartesian coordinates
		element.Calculate(LOCAL_ELEMENT_ORIENTATION, local_coordinate_orientation, r_process_info);

		Vector local_axis_1 = ZeroVector(3);
		Vector local_axis_2 = ZeroVector(3);
		Vector local_axis_3 = ZeroVector(3);

		for (size_t i=0;i<3;++i)
		{
			local_axis_1[i] = local_coordinate_orientation(i,0);
			local_axis_2[i] = local_coordinate_orientation(i,1);
			local_axis_3[i] = local_coordinate_orientation(i,2);
		}

		// normalise local axis vectors (global cartesian)
		local_axis_1 /= norm_2(local_axis_1);
		local_axis_2 /= norm_2(local_axis_2);
		local_axis_3 /= norm_2(local_axis_3);

		// (Abaqus default projection)
		// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
		// Shell local axis 1 is the projection of Global X vector onto the shell surface.
		// If the Global X vector is normal to the shell surface,
		// the shell local 1-direction is the projection of the
		// Global Z vector onto the shell surface

		// First, check if specified global_vector is normal to the shell surface
		if (std::abs(inner_prod(global_vector, local_axis_1)) < std::numeric_limits<double>::epsilon() && std::abs(inner_prod(global_vector, local_axis_2)) < std::numeric_limits<double>::epsilon())
		{
			KRATOS_ERROR << "Global direction is perpendicular to element " << element.GetId() << " please define a different projection plane or use another type of projection "
				<< ", available: radial,spherical" << std::endl;
		}
		else
		{
			// Second, project the global vector onto the shell surface
			// http://www.euclideanspace.com/maths/geometry/elements/plane/lineOnPlane/index.htm
			// vector to be projected = vec_a
			// Surface normal = vec_b
			const Vector& vec_a = global_vector;
			const Vector& vec_b = local_axis_3;

			Vector a_cross_b = ZeroVector(3);
			Vector projected_result = ZeroVector(3);

			MathUtils<double>::CrossProduct(a_cross_b, vec_a, vec_b);
			MathUtils<double>::CrossProduct(projected_result, vec_b, a_cross_b);
			//noramlize projected result
			projected_result /= MathUtils<double>::Norm(projected_result);

			element.SetValue(r_variable, projected_result);
		}
	}
}

void ProjectVectorOnSurfaceUtility::RadialProjection(ModelPart& rModelPart, const Parameters& MethodParameters)
{
	const auto specific_parameters = MethodParameters["method_specific_settings"];
	const std::string& r_local_variable_name = specific_parameters["local_variable_name"].GetString();

	KRATOS_ERROR_IF_NOT(KratosComponents<ArrayVariableType>::Has(r_local_variable_name)) << "Variable " << r_local_variable_name << " not known" << std::endl;
	const ArrayVariableType& r_variable = KratosComponents<ArrayVariableType>::Get(r_local_variable_name);

	Vector3 global_vector;

	//global_vector is the global fiber direction given by user
	//read this global direction
	CheckAndReadVectors(specific_parameters, "global_fiber_direction", global_vector);

	//Normalize global_vector
	global_vector /=  norm_2(global_vector);

	const auto& r_process_info = rModelPart.GetProcessInfo();

	// Declare working variables
	Matrix local_coordinate_orientation;

	// Loop over all elements in part
	for (auto &element : rModelPart.Elements())
	{

		// get local axis in cartesian coordinates
		element.Calculate(LOCAL_ELEMENT_ORIENTATION, local_coordinate_orientation, r_process_info);

		Vector local_axis_1 = ZeroVector(3);
		Vector local_axis_2 = ZeroVector(3);
		Vector local_axis_3 = ZeroVector(3);

		for (size_t i=0;i<3;++i)
		{
			local_axis_1[i] = local_coordinate_orientation(i,0);
			local_axis_2[i] = local_coordinate_orientation(i,1);
			local_axis_3[i] = local_coordinate_orientation(i,2);
		}

		// normalise local axis vectors (global cartesian)
		local_axis_1 /= norm_2(local_axis_1);
		local_axis_2 /= norm_2(local_axis_2);
		local_axis_3 /= norm_2(local_axis_3);

		// (Abaqus default projection)
		// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
		// Shell local axis 1 is the projection of Global X vector onto the shell surface.
		// If the Global X vector is normal to the shell surface,
		// the shell local 1-direction is the projection of the
		// Global Z vector onto the shell surface

		// First, check if specified global_vector is normal to the shell surface
		if (std::abs(inner_prod(global_vector, local_axis_1)) < std::numeric_limits<double>::epsilon() && std::abs(inner_prod(global_vector, local_axis_2)) < std::numeric_limits<double>::epsilon())
		{
			KRATOS_ERROR << "Global direction is perpendicular to element " << element.GetId() << " please define a different projection plane or use another type of projection "
				<< ", available: planar,spherical" << std::endl;
		}
		else
		{
			Vector projected_result = ZeroVector(3);

			MathUtils<double>::CrossProduct(projected_result, global_vector, local_axis_3);

			//noramlize projected result
			projected_result /= MathUtils<double>::Norm(projected_result);

			element.SetValue(r_variable, projected_result);
		}
	}
}

void ProjectVectorOnSurfaceUtility::SphericalProjection(ModelPart& rModelPart, const Parameters& MethodParameters)
{
	KRATOS_ERROR << "SphericalProjection not implemented" << std::endl;
}

void ProjectVectorOnSurfaceUtility::CheckAndReadVectors(Parameters ThisParameters, const std::string KeyName, Vector3 &rVector)
{
	if (ThisParameters[KeyName].size() != 3)
	{
		KRATOS_ERROR << "\" " << KeyName << "\" is not of size 3" << std::endl;
	}


	rVector[0] = ThisParameters[KeyName][0].GetDouble();
	rVector[1] = ThisParameters[KeyName][1].GetDouble();
	rVector[2] = ThisParameters[KeyName][2].GetDouble();

	if (inner_prod(rVector, rVector) < std::numeric_limits<double>::epsilon())
	{
		KRATOS_ERROR << "Vector \" " << KeyName << "\" has zero length" << std::endl;
	}
}

} // namespace Kratos.
