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
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "processes/find_nodal_neighbours_process.h"
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

class MassResponseFunction : ResponseFunction
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
		// Set gradient mode
		std::string gradientMode = responseSettings["gradient_mode"].GetString();

		// Strings for comparison

		// Mode 3: global finite differencing
		if (gradientMode.compare("finite_differencing") == 0)
		{
			m_gradient_mode = 3;
			double delta = responseSettings["step_size"].GetDouble();
			mDelta = delta;
		}
		// Throw error message in case of wrong specification
		else
			KRATOS_ERROR << "Specified gradient_mode not recognized. Options are: finite_differencing. Specified gradient_mode: " << gradientMode << std::endl;

		mConsiderDiscretization =  responseSettings["discretization_weighting"].GetBool();

		// Initialize member variables to NULL
		m_initial_value = 0.0;
		m_initial_value_defined = false;
		m_total_mass = 0.0;
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
	void CalculateValue()
	{
		KRATOS_TRY;

		// Variables
		m_total_mass = 0.0;

		// Incremental computation of total mass
		for (ModelPart::ElementsContainerType::iterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		{
			Element::GeometryType& element_geometry = elem_i->GetGeometry();
			double elem_density = elem_i->GetProperties()[DENSITY];

			// Compute mass according to element dimension
			double elem_volume = 0.0;
			if( IsElementOfTypeShell(element_geometry) )
				elem_volume = element_geometry.Area()*elem_i->GetProperties()[THICKNESS];
			else
				elem_volume = element_geometry.Volume();
			m_total_mass +=  elem_density*elem_volume;
		}

		// Set initial value if not done yet
		if(!m_initial_value_defined)
		{
			m_initial_value = m_total_mass;
			m_initial_value_defined = true;
		}

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
			noalias(node_i->FastGetSolutionStepValue(MASS_SHAPE_GRADIENT)) = zeros_array;

		switch (m_gradient_mode)
		{
		// Global finite differencing
    	case 3:
		{
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
					Kratos::Element& ng_elem_i = ng_elem[i];
					Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( IsElementOfTypeShell(element_geometry) )
						elem_volume = element_geometry.Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = element_geometry.Volume();
					mass_before_fd +=  elem_density*elem_volume;
				}

				// Compute sensitivities using finite differencing in the three spatial direction
				array_3d gradient(3, 0.0);

				// Apply pertubation in X-direction and recompute total mass of all neighbor elements
				double mass_after_fd = 0.0;
				node_i->X() += mDelta;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element& ng_elem_i = ng_elem[i];
					Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( IsElementOfTypeShell(element_geometry) )
						elem_volume = element_geometry.Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = element_geometry.Volume();
					mass_after_fd +=  elem_density*elem_volume;
				}
				gradient[0] = (mass_after_fd - mass_before_fd) / mDelta;
				node_i->X() -= mDelta;

				// Apply pertubation in Y-direction and recompute total mass of all neighbor elements
				mass_after_fd = 0.0;
				node_i->Y() += mDelta;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element& ng_elem_i = ng_elem[i];
					Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( IsElementOfTypeShell(element_geometry) )
						elem_volume = element_geometry.Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = element_geometry.Volume();
					mass_after_fd +=  elem_density*elem_volume;
				}
				gradient[1] = (mass_after_fd - mass_before_fd) / mDelta;
				node_i->Y() -= mDelta;

				// Apply pertubation in Z-direction and recompute total mass of all neighbor elements
				mass_after_fd = 0.0;
				node_i->Z() += mDelta;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element& ng_elem_i = ng_elem[i];
					Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( IsElementOfTypeShell(element_geometry) )
						elem_volume = element_geometry.Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = element_geometry.Volume();
					mass_after_fd +=  elem_density*elem_volume;
				}
				gradient[2] = (mass_after_fd - mass_before_fd) / mDelta;
				node_i->Z() -= mDelta;

				// Compute sensitivity
				noalias(node_i->FastGetSolutionStepValue(MASS_SHAPE_GRADIENT)) = gradient;
			}

			if (mConsiderDiscretization)
				this->ConsiderDiscretization();
		}

		}


		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double GetInitialValue()
	{
		KRATOS_TRY;

		if(!m_initial_value_defined)
			KRATOS_ERROR << "Initial value not yet defined! First compute it by calling \"calculate_value()\"!" << std::endl;

		return m_initial_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double GetValue()
	{
		KRATOS_TRY;

		return m_total_mass;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	boost::python::dict GetGradient()
	{
		KRATOS_TRY;

		// Dictionary to store all sensitivities along with Ids of corresponding nodes
		boost::python::dict dFdX;

		// Fill dictionary with gradient information
		for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
			dFdX[node_i->Id()] = node_i->FastGetSolutionStepValue(MASS_SHAPE_GRADIENT);

		return dFdX;

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
				Kratos::Element& ng_elem_i = ng_elem[i];
				Element::GeometryType& element_geometry = ng_elem_i.GetGeometry();

				// Compute mass according to element dimension
				if( IsElementOfTypeShell(element_geometry) )
					scaling_factor += element_geometry.Area();
				else
					scaling_factor += element_geometry.Volume();
			}

			// apply scaling
			node_i->FastGetSolutionStepValue(MASS_SHAPE_GRADIENT) /= scaling_factor;
		}
	}

	// --------------------------------------------------------------------------
	bool IsElementOfTypeShell( Element::GeometryType& given_element_geometry )
	{
		if(given_element_geometry.WorkingSpaceDimension() != given_element_geometry.LocalSpaceDimension())
			return true;
		else
		    return false;
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
	unsigned int m_gradient_mode;
	double m_total_mass;
	double mDelta;
	double m_initial_value;
	bool m_initial_value_defined;
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
