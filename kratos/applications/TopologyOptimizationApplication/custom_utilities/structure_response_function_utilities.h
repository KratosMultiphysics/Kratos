// ==============================================================================
/*
 KratosTopologyOptimizationApplication
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
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
//
// ==============================================================================

#if !defined(KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED)
#define  KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"

// Application includes
#include "topology_optimization_application.h"


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

/// Solution utility to compute structural analysis responses.
/** Detail class definition.

 */

class StructureResponseFunctionUtilities
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of StructureResponseFunctionUtilities
	KRATOS_CLASS_POINTER_DEFINITION(StructureResponseFunctionUtilities);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	StructureResponseFunctionUtilities( ModelPart& model_part )
	: mr_structure_model_part(model_part)
	{
	}

	/// Destructor.
	virtual ~StructureResponseFunctionUtilities()
	{
	}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- COMPUTE STRAIN ENERGY -------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// Computes the strain energy as the objective function of the optimization problem.
	double ComputeStrainEnergy()
	{
		KRATOS_TRY;

		clock_t begin = clock();
		std::cout<<"  Start calculating strain energy."<<std::endl;

		double Out = 0.0;
		double Global_Strain_Energy = 0.0;

		// Loop over all elements to calculate their local objective function and sum it into the global objective function (Global Strain Energy)
		for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			element_i->Calculate(LOCAL_STRAIN_ENERGY, Out, mr_structure_model_part.GetProcessInfo());
			Global_Strain_Energy += element_i->GetValue(LOCAL_STRAIN_ENERGY);
		}

		clock_t end = clock();
		std::cout << "  Strain energy calculated                  [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		// Return this obtained Global Strain Energy value as the objective function of the complete system
		return Global_Strain_Energy;

		KRATOS_CATCH("");
	}

	double ComputeVolumeFraction()
	{
		KRATOS_TRY;

		clock_t begin = clock();
		std::cout<<"  Start calculating volume fraction."<<std::endl;

		int number_elements = 0;
		double Global_Volume_Fraction = 0.0;

		// Loop over all elements to obtain their X_PHYS and know how many elements the model has
		for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			Global_Volume_Fraction += element_i->GetValue(X_PHYS);
			number_elements++;
		}

		// Calculate and return the Global Volume Fraction by knowing how many elements the model has
		Global_Volume_Fraction = Global_Volume_Fraction/number_elements;

		clock_t end = clock();
		std::cout << "  Volume fraction calculated                 [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		return Global_Volume_Fraction;

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
		return "StructureResponseFunctionUtilities";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "StructureResponseFunctionUtilities";
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

	ModelPart& mr_structure_model_part;

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
	//StructureResponseFunctionUtilities& operator=(StructureResponseFunctionUtilities const& rOther);

	/// Copy constructor.
	//StructureResponseFunctionUtilities(StructureResponseFunctionUtilities const& rOther);


	///@}

}; // Class StructureResponseFunctionUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED */
