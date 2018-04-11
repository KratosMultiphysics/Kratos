// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MASS_RESPONSE_FUNCTION_H
#define MASS_RESPONSE_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/find_nodal_neighbours_process.h"
#include "structural_mechanics_application_variables.h"
#include "response_function.h"

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

class MassResponseFunction : public ResponseFunction
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of MassResponseFunction
	KRATOS_CLASS_POINTER_DEFINITION(MassResponseFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	MassResponseFunction(ModelPart& model_part, Parameters responseSettings)
	: mr_model_part(model_part)
	{
		std::string gradientMode = responseSettings["gradient_mode"].GetString();
		if (gradientMode.compare("finite_differencing") == 0)
		{
			double delta = responseSettings["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified gradient_mode not recognized. The only option is: finite_differencing. Specified gradient_mode: ", gradientMode);

		mConsiderDiscretization =  responseSettings["consider_discretization"].GetBool();
	}

	/// Destructor.
	virtual ~MassResponseFunction()
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
		for (auto& elem_i : mr_model_part.Elements())
		{
			double elem_density = elem_i.GetProperties()[DENSITY];

			// Compute mass according to element dimension
			double elem_volume = GetElementVolume(elem_i);
			total_mass +=  elem_density*elem_volume;
		}

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
		array_3d zeros_array(3, 0.0);
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			noalias(node_i->FastGetSolutionStepValue(SHAPE_SENSITIVITY)) = zeros_array;

		// Start process to identify element neighbors for every node
		FindNodalNeighboursProcess neighorFinder = FindNodalNeighboursProcess(mr_model_part, 10, 10);
		neighorFinder.Execute();

		for(ModelPart::NodeIterator node_i=mr_model_part.NodesBegin(); node_i!=mr_model_part.NodesEnd(); node_i++)
		{
			// Get all neighbor elements of current node
			WeakPointerVector<Element >& ng_elem = node_i->GetValue(NEIGHBOUR_ELEMENTS);

			// Compute total mass of all neighbor elements before finite differencing
			double mass_before_fd = 0.0;
			for(unsigned int i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				double elem_density = ng_elem_i.GetProperties()[DENSITY];

				// Compute mass according to element dimension
				double elem_volume = GetElementVolume(ng_elem_i);
				mass_before_fd +=  elem_density*elem_volume;
			}

			// Compute sensitivities using finite differencing in the three spatial direction
			array_3d gradient(3, 0.0);

			// Apply pertubation in X-direction and recompute total mass of all neighbor elements
			double mass_after_fd = 0.0;
			node_i->X() += mDelta;
			for(unsigned int i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				double elem_density = ng_elem_i.GetProperties()[DENSITY];

				// Compute mass according to element dimension
				double elem_volume = GetElementVolume(ng_elem_i);
				mass_after_fd +=  elem_density*elem_volume;
			}
			gradient[0] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i->X() -= mDelta;

			// Apply pertubation in Y-direction and recompute total mass of all neighbor elements
			mass_after_fd = 0.0;
			node_i->Y() += mDelta;
			for(unsigned int i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				double elem_density = ng_elem_i.GetProperties()[DENSITY];

				// Compute mass according to element dimension
				double elem_volume = GetElementVolume(ng_elem_i);
				mass_after_fd +=  elem_density*elem_volume;
			}
			gradient[1] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i->Y() -= mDelta;

			// Apply pertubation in Z-direction and recompute total mass of all neighbor elements
			mass_after_fd = 0.0;
			node_i->Z() += mDelta;
			for(unsigned int i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				double elem_density = ng_elem_i.GetProperties()[DENSITY];

				// Compute mass according to element dimension
				double elem_volume = GetElementVolume(ng_elem_i);
				mass_after_fd +=  elem_density*elem_volume;
			}
			gradient[2] = (mass_after_fd - mass_before_fd) / mDelta;
			node_i->Z() -= mDelta;

			// Compute sensitivity
			noalias(node_i->FastGetSolutionStepValue(SHAPE_SENSITIVITY)) = gradient;

		}

		if (mConsiderDiscretization)
			this->ConsiderDiscretization();

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
  	virtual void ConsiderDiscretization(){

		std::cout<< "> Considering discretization size!" << std::endl;
		for(ModelPart::NodeIterator node_i=mr_model_part.NodesBegin(); node_i!=mr_model_part.NodesEnd(); node_i++)
		{
			// Get all neighbor elements of current node
			WeakPointerVector<Element >& ng_elem = node_i->GetValue(NEIGHBOUR_ELEMENTS);

			// Compute total mass of all neighbor elements before finite differencing
			double scaling_factor = 0.0;
			for(unsigned int i = 0; i < ng_elem.size(); i++)
			{
				Element& ng_elem_i = ng_elem[i];
				Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();

				// Compute mass according to element dimension
				scaling_factor += element_geometry.DomainSize();
			}

			// apply scaling
			node_i->FastGetSolutionStepValue(SHAPE_SENSITIVITY) /= scaling_factor;
		}
	}

	// --------------------------------------------------------------------------
	double GetElementVolume(Element& element)
	{
		Element::GeometryType& geometry = element.GetGeometry();
		if (geometry.LocalSpaceDimension() == 3)
			return geometry.Volume();
		else if (geometry.LocalSpaceDimension() == 2)
			return geometry.Area()*element.GetProperties()[THICKNESS];
		else if (geometry.LocalSpaceDimension() == 1)
			return geometry.Length()*element.GetProperties()[CROSS_AREA];
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
	virtual std::string Info() const
	{
		return "MassResponseFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "MassResponseFunction";
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

	ModelPart &mr_model_part;
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
	//      MassResponseFunction& operator=(MassResponseFunction const& rOther);

	/// Copy constructor.
	//      MassResponseFunction(MassResponseFunction const& rOther);

	///@}

}; // Class MassResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // MASS_RESPONSE_FUNCTION_H
