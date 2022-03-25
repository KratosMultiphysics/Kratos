// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//

#ifndef MASS_RESPONSE_FUNCTION_UTILITY_H
#define MASS_RESPONSE_FUNCTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "processes/find_nodal_neighbours_process.h"
#include "structural_mechanics_application_variables.h"
#include "custom_processes/total_structural_mass_process.h"

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

class MassResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of MassResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(MassResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	MassResponseFunctionUtility(ModelPart& model_part, Parameters responseSettings)
	: mrModelPart(model_part)
	{
		std::string gradient_mode = responseSettings["gradient_settings"]["gradient_mode"].GetString();
		if (gradient_mode.compare("finite_differencing") == 0)
		{
			double delta = responseSettings["gradient_settings"]["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;

		for (auto& control_type : responseSettings["control_types"]){
			if (control_type.GetString()=="thickness")
				thickness_grad = true;
			if (control_type.GetString()=="shape")
				shape_grad = true;
			if (control_type.GetString()=="topology")
				density_grad = true;
		}


	}

	/// Destructor.
	virtual ~MassResponseFunctionUtility()
	{
	}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	// ==============================================================================
	void Initialize()
	{}

	// --------------------------------------------------------------------------
	double CalculateValue()
	{
		KRATOS_TRY;

		// Variables
		double total_mass = 0.0;
		const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

		if(!density_grad)
		{
			// Incremental computation of total mass
			for (auto& elem_i : mrModelPart.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active)
					total_mass += TotalStructuralMassProcess::CalculateElementMass(elem_i, domain_size);
			}
		}
		else
		{
			for (auto& elem_i : mrModelPart.Elements())
			{
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active)
				{
					const auto& r_geom = elem_i.GetGeometry();	
					const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
					const auto& integration_points = r_geom.IntegrationPoints(integration_method);
					const auto& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);
					for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
					{
						Matrix J0,InvJ0;
						GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[i_point], J0);
						double detJ0;
						MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
						const double IntToReferenceWeight = integration_points[i_point].Weight() * detJ0;
						const auto& rN = row(Ncontainer,i_point);
						int node_index = 0;
						for (auto& node_i : r_geom){
							total_mass += node_i.FastGetSolutionStepValue(DENSITY) * rN[node_index] * IntToReferenceWeight;
							node_index++;
						}						
					}
				}				
			}
		}

		return total_mass;

		KRATOS_CATCH("");
	}
	// --------------------------------------------------------------------------
	void CalculateShapeGradient()
	{
		KRATOS_TRY;

		VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

		const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

		// Start process to identify element neighbors for every node
		FindNodalNeighboursProcess neighorFinder(mrModelPart);
		neighorFinder.Execute();

		for(auto& node_i : mrModelPart.Nodes())
		{
			// Get all neighbor elements of current node
			GlobalPointersVector<Element >& ng_elem = node_i.GetValue(NEIGHBOUR_ELEMENTS);

			// Compute total mass of all neighbor elements before finite differencing
			double mass_before_fd = 0.0;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				const bool element_is_active = ng_elem_i.IsDefined(ACTIVE) ? ng_elem_i.Is(ACTIVE) : true;

				// Compute mass according to element dimension
				if(element_is_active)
					mass_before_fd += TotalStructuralMassProcess::CalculateElementMass(ng_elem_i, domain_size);
			}
			
			// Compute sensitivities using finite differencing in the three spatial direction
			array_3d gradient(3, 0.0);

			// Apply pertubation in X-direction and recompute total mass of all neighbor elements
			double mass_after_fd = 0.0;
			node_i.X() += mDelta;
			node_i.X0() += mDelta;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				const bool element_is_active = ng_elem_i.IsDefined(ACTIVE) ? ng_elem_i.Is(ACTIVE) : true;

				// Compute mass according to element dimension
				if(element_is_active)
					mass_after_fd += TotalStructuralMassProcess::CalculateElementMass(ng_elem_i, domain_size);
			}
			gradient[0] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i.X() -= mDelta;
			node_i.X0() -= mDelta;

			// Apply pertubation in Y-direction and recompute total mass of all neighbor elements
			mass_after_fd = 0.0;
			node_i.Y() += mDelta;
			node_i.Y0() += mDelta;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				const bool element_is_active = ng_elem_i.IsDefined(ACTIVE) ? ng_elem_i.Is(ACTIVE) : true;

				// Compute mass according to element dimension
				if(element_is_active)
					mass_after_fd += TotalStructuralMassProcess::CalculateElementMass(ng_elem_i, domain_size);
			}
			gradient[1] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i.Y() -= mDelta;
			node_i.Y0() -= mDelta;

			// Apply pertubation in Z-direction and recompute total mass of all neighbor elements
			mass_after_fd = 0.0;
			node_i.Z() += mDelta;
			node_i.Z0() += mDelta;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				const bool element_is_active = ng_elem_i.IsDefined(ACTIVE) ? ng_elem_i.Is(ACTIVE) : true;

				// Compute mass according to element dimension
				if(element_is_active)
					mass_after_fd += TotalStructuralMassProcess::CalculateElementMass(ng_elem_i, domain_size);
			}
			gradient[2] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i.Z() -= mDelta;
			node_i.Z0() -= mDelta;

			// Compute sensitivity
			noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) = gradient;

		}
		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateThicknessGradient()
	{
		KRATOS_TRY;
		VariableUtils().SetHistoricalVariableToZero(THICKNESS_SENSITIVITY, mrModelPart.Nodes());
		const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

		for(auto& elem_i : mrModelPart.Elements())
		{				
			const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
			if(element_is_active)
			{
				Properties& elem_i_prop = elem_i.GetProperties();
				auto original_thickness = elem_i_prop.GetValue(THICKNESS);
				elem_i_prop.SetValue(THICKNESS,1.0);
				double elem_gradient = TotalStructuralMassProcess::CalculateElementMass(elem_i, domain_size);
				elem_i_prop.SetValue(THICKNESS,original_thickness);								
				const auto& r_geom = elem_i.GetGeometry();	
				const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
				const auto& integration_points = r_geom.IntegrationPoints(integration_method);
				const unsigned int NumGauss = integration_points.size();
				Vector GaussPtsJDet = ZeroVector(NumGauss);
				r_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);
				const auto& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);
				double elem_area = r_geom.Area(); 			
				for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
				{
					const double IntToReferenceWeight = integration_points[i_point].Weight() * GaussPtsJDet[i_point];
					const auto& rN = row(Ncontainer,i_point);
					int node_index = 0;
					for (auto& node_i : r_geom){
						node_i.FastGetSolutionStepValue(THICKNESS_SENSITIVITY) += elem_gradient * IntToReferenceWeight * rN[node_index]/elem_area;
						node_index++;
					}						
				}
			}
		}

		KRATOS_CATCH("");

	}

	// --------------------------------------------------------------------------
	void CalculateDensityGradient()
	{
		KRATOS_TRY;
		VariableUtils().SetHistoricalVariableToZero(YOUNG_MODULUS_SENSITIVITY, mrModelPart.Nodes());

		for (auto& elem_i : mrModelPart.Elements())
		{
			const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
			if(element_is_active)
			{
				const auto& r_geom = elem_i.GetGeometry();	
				const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
				const auto& integration_points = r_geom.IntegrationPoints(integration_method);
				const auto& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);
				for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
				{

					Matrix J0,InvJ0;
					GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[i_point], J0);
					double detJ0;
					MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);

					const double IntToReferenceWeight = integration_points[i_point].Weight() * detJ0;
					const auto& rN = row(Ncontainer,i_point);
					int node_index = 0;
					for (auto& node_i : r_geom){
						node_i.FastGetSolutionStepValue(YOUNG_MODULUS_SENSITIVITY) += rN[node_index] * IntToReferenceWeight;
						node_index++;
					}						
				}
			}				
		}

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateGradient()
	{
		KRATOS_TRY;

		if(shape_grad)
			CalculateShapeGradient();
			
		if(thickness_grad)
			CalculateThicknessGradient();
			
		if(density_grad)
			CalculateDensityGradient();


		KRATOS_CATCH("");
	}

	// ==============================================================================

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
	std::string Info() const
	{
		return "MassResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "MassResponseFunctionUtility";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
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

	// ==============================================================================

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

	ModelPart &mrModelPart;
	double mDelta;
	bool shape_grad = false;
	bool thickness_grad = false;
	bool density_grad = false;	

	///@}
///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

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
	//      MassResponseFunctionUtility& operator=(MassResponseFunctionUtility const& rOther);

	/// Copy constructor.
	//      MassResponseFunctionUtility(MassResponseFunctionUtility const& rOther);

	///@}

}; // Class MassResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // MASS_RESPONSE_FUNCTION_UTILITY_H
