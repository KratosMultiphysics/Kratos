// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//                        $Author:           armin.geiser@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
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
	MassResponseFunction(ModelPart &model_part, boost::python::dict response_settings)
	: mr_model_part(model_part)
	{
		// Set gradient mode
		boost::python::extract<const char *> grad_mode(response_settings["gradient_mode"]);

		// Mode 3: global finite differencing
		if (std::strcmp(grad_mode, "finite_differencing") == 0)
		{
			m_gradient_mode = 3;
			boost::python::extract<double> delta(response_settings["step_size"]);
			m_delta = delta;
		}

		// Throw error message in case of wrong specification
		else
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified gradient_mode not recognized. Options are: finite_differencing ", grad_mode);

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
	void initialize()
	{}

	// --------------------------------------------------------------------------
	void calculate_value()
	{
		KRATOS_TRY;

		// Variables
		m_total_mass = 0.0;

		// Incremental computation of total mass
		for (ModelPart::ElementsContainerType::iterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
		{
			const unsigned int elem_dimension = elem_i->GetGeometry().WorkingSpaceDimension();
			double elem_density = elem_i->GetProperties()[DENSITY];

			// Compute mass according to element dimension
			double elem_volume = 0.0;
			if( elem_dimension == 2 )
				elem_volume = elem_i->GetGeometry().Area()*elem_i->GetProperties()[THICKNESS];
			else
				elem_volume = elem_i->GetGeometry().Volume();
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
	void calculate_gradient()
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

			for(NodesContainerType::iterator node_i=mr_model_part.NodesBegin(); node_i!=mr_model_part.NodesEnd(); node_i++)
			{
				// Get all neighbor elements of current node
				WeakPointerVector<Element >& ng_elem = node_i->GetValue(NEIGHBOUR_ELEMENTS);

				// Compute total mass of all neighbor elements before finite differencing
				double mass_before_fd = 0.0;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element ng_elem_i = ng_elem[i];
					const unsigned int elem_dimension = ng_elem_i.GetGeometry().WorkingSpaceDimension();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( elem_dimension == 2 )
						elem_volume = ng_elem_i.GetGeometry().Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = ng_elem_i.GetGeometry().Volume();
					mass_before_fd +=  elem_density*elem_volume;
				}

				// Compute sensitivities using finite differencing in the three spatial direction
				array_3d gradient(3, 0.0);

				// Apply pertubation in X-direction and recompute total mass of all neighbor elements
				double mass_after_fd = 0.0;
				node_i->X() += m_delta;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element ng_elem_i = ng_elem[i];
					const unsigned int elem_dimension = ng_elem_i.GetGeometry().WorkingSpaceDimension();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( elem_dimension == 2 )
						elem_volume = ng_elem_i.GetGeometry().Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = ng_elem_i.GetGeometry().Volume();
					mass_after_fd +=  elem_density*elem_volume;
				}
				gradient[0] = (mass_after_fd - mass_before_fd) / m_delta;
				node_i->X() -= m_delta;

				// Apply pertubation in Y-direction and recompute total mass of all neighbor elements
				mass_after_fd = 0.0;
				node_i->Y() += m_delta;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element ng_elem_i = ng_elem[i];
					const unsigned int elem_dimension = ng_elem_i.GetGeometry().WorkingSpaceDimension();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( elem_dimension == 2 )
						elem_volume = ng_elem_i.GetGeometry().Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = ng_elem_i.GetGeometry().Volume();
					mass_after_fd +=  elem_density*elem_volume;
				}
				gradient[1] = (mass_after_fd - mass_before_fd) / m_delta;
				node_i->Y() -= m_delta;

				// Apply pertubation in Z-direction and recompute total mass of all neighbor elements
				mass_after_fd = 0.0;
				node_i->Z() += m_delta;
				for(unsigned int i = 0; i < ng_elem.size(); i++)
				{
					Kratos::Element ng_elem_i = ng_elem[i];
					const unsigned int elem_dimension = ng_elem_i.GetGeometry().WorkingSpaceDimension();
					double elem_density = ng_elem_i.GetProperties()[DENSITY];

					// Compute mass according to element dimension
					double elem_volume = 0.0;
					if( elem_dimension == 2 )
						elem_volume = ng_elem_i.GetGeometry().Area()*ng_elem_i.GetProperties()[THICKNESS];
					else
						elem_volume = ng_elem_i.GetGeometry().Volume();
					mass_after_fd +=  elem_density*elem_volume;
				}
				gradient[2] = (mass_after_fd - mass_before_fd) / m_delta;
				node_i->Z() -= m_delta;

				// Compute sensitivity
				noalias(node_i->FastGetSolutionStepValue(MASS_SHAPE_GRADIENT)) = gradient;
			}

			break;
		}

		}

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double get_initial_value()
	{
		KRATOS_TRY;

		if(!m_initial_value_defined)
			KRATOS_THROW_ERROR(std::logi:error, "Initial value not yet defined! First compute it by calling \"calculate_value()\"", m_initial_value_defined);

		return m_initial_value;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	double get_value()
	{
		KRATOS_TRY;

		return m_total_mass;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	boost::python::dict get_gradient()
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
	double m_delta;
	double m_initial_value;
	bool m_initial_value_defined;

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
