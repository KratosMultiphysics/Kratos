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
		std::string gradient_mode = responseSettings["gradient_mode"].GetString();
		if (gradient_mode.compare("finite_differencing") == 0)
		{
			double delta = responseSettings["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;

		mConsiderDiscretization =  responseSettings["consider_discretization"].GetBool();
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

		// Incremental computation of total mass
		for (auto& elem_i : mrModelPart.Elements())
			total_mass += GetElementMass(elem_i);

		return total_mass;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateGradient()
	{
		KRATOS_TRY;

		// Formula computed in general notation:
		// \frac{dm_{total}}{dx}

		// First gradients are initialized
		VariableUtils().SetToZero_VectorVar(SHAPE_SENSITIVITY, mrModelPart.Nodes());


		// Start process to identify element neighbors for every node
		FindNodalNeighboursProcess neighorFinder = FindNodalNeighboursProcess(mrModelPart, 10, 10);
		neighorFinder.Execute();

		for(auto& node_i : mrModelPart.Nodes())
		{
			// Get all neighbor elements of current node
			WeakPointerVector<Element >& ng_elem = node_i.GetValue(NEIGHBOUR_ELEMENTS);

			// Compute total mass of all neighbor elements before finite differencing
			double mass_before_fd = 0.0;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];

				// Compute mass according to element dimension
				mass_before_fd += GetElementMass(ng_elem_i);
			}

			// Compute sensitivities using finite differencing in the three spatial direction
			array_3d gradient(3, 0.0);

			// Apply pertubation in X-direction and recompute total mass of all neighbor elements
			double mass_after_fd = 0.0;
			node_i.X() += mDelta;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];

				// Compute mass according to element dimension
				mass_after_fd += GetElementMass(ng_elem_i);
			}
			gradient[0] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i.X() -= mDelta;

			// Apply pertubation in Y-direction and recompute total mass of all neighbor elements
			mass_after_fd = 0.0;
			node_i.Y() += mDelta;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];

				// Compute mass according to element dimension
				mass_after_fd += GetElementMass(ng_elem_i);
			}
			gradient[1] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i.Y() -= mDelta;

			// Apply pertubation in Z-direction and recompute total mass of all neighbor elements
			mass_after_fd = 0.0;
			node_i.Z() += mDelta;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];

				// Compute mass according to element dimension
				mass_after_fd += GetElementMass(ng_elem_i);
			}
			gradient[2] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i.Z() -= mDelta;

			// Compute sensitivity
			noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) = gradient;

		}

		if (mConsiderDiscretization)
			this->ConsiderDiscretization();

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
  	virtual void ConsiderDiscretization()
	{

		std::cout<< "> Considering discretization size!" << std::endl;
		for(auto& node_i : mrModelPart.Nodes())
		{
			// Get all neighbor elements of current node
			WeakPointerVector<Element >& ng_elem = node_i.GetValue(NEIGHBOUR_ELEMENTS);

			// Compute total mass of all neighbor elements before finite differencing
			double scaling_factor = 0.0;
			for(std::size_t i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();

				// Compute mass according to element dimension
				scaling_factor += element_geometry.DomainSize();
			}

			// apply scaling
			node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY) /= scaling_factor;
		}
	}

	// --------------------------------------------------------------------------
	double GetElementMass(Element& element)
	{
		Element::GeometryType& geometry = element.GetGeometry();
		if (geometry.LocalSpaceDimension() == 3)
			return geometry.Volume()*element.GetProperties()[DENSITY];
		else if (geometry.LocalSpaceDimension() == 2)
			return geometry.Area()*element.GetProperties()[THICKNESS]*element.GetProperties()[DENSITY];
		else if (geometry.LocalSpaceDimension() == 1)
			return geometry.Length()*element.GetProperties()[CROSS_AREA]*element.GetProperties()[DENSITY];
		else if (geometry.LocalSpaceDimension() == 0){
			return element.GetValue(NODAL_MASS);
		}
		else
			KRATOS_ERROR << "Invalid local dimension found in element!" << std::endl;
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
	bool mConsiderDiscretization;

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
