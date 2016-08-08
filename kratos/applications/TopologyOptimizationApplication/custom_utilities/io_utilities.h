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

#if !defined(KRATOS_IO_UTILITIES_H_INCLUDED)
#define  KRATOS_IO_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

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

/// Solution utility to filter results.
/** Detail class definition.

 */

class IOUtilities
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of IOUtilities
	KRATOS_CLASS_POINTER_DEFINITION(IOUtilities);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	IOUtilities( )
	{
	}

	/// Destructor.
	virtual ~IOUtilities()
	{
	}

	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- Save optimization results in a restart file (mdpa) --------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	void SaveOptimizationResults( const char* RestartInputFile, ModelPart& rModelPart, const char* RestartOutputFile )
	{
		KRATOS_TRY;

		std::cout<<"\n::[Saving optimization results as restart file]::"<<std::endl;

		// Create an empty .mdpa restart file
		std::ofstream FileToBeCreated;
		FileToBeCreated.open(RestartOutputFile, std::ios::trunc);

		// Read the original .mdpa file
		std::ifstream FileToBeRead(RestartInputFile);
		if(FileToBeRead==0)
			KRATOS_THROW_ERROR(std::invalid_argument, "Specified restart input file does not exist: ",RestartInputFile);

		// Write the given input file except the block covering X_PHYS (this is to be replaced by the current optimization results)
		bool DensityBlockActive = false;
		std::string LineString;
		while (std::getline(FileToBeRead, LineString))
		{
			if(LineString.find("Begin ElementalData X_PHYS")!=std::string::npos)
				DensityBlockActive = true;
			else if(LineString.find("End ElementalData")!=std::string::npos and DensityBlockActive)
			{
				DensityBlockActive = false;
				continue;
			}
			else if(DensityBlockActive == false)
				FileToBeCreated << LineString << "\n";
		}

		// Write the actual X_PHYS elemental data
		FileToBeCreated << "\nBegin ElementalData X_PHYS\n";
		for(ModelPart::ElementsContainerType::iterator elem_i = rModelPart.ElementsBegin(); elem_i!=rModelPart.ElementsEnd(); elem_i++)
		{
			FileToBeCreated << "    " << elem_i->Id() << " " << elem_i->GetValue(X_PHYS) << "\n";
		}
		FileToBeCreated << "End ElementalData\n";

		// Close files
		FileToBeCreated.close();
		FileToBeRead.close();

		std::cout<<"  Restart File succesfully generated under the name " << RestartOutputFile <<std::endl;

		KRATOS_CATCH("");
	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- WRITE STL FILE FROM SURFACE MESH  -------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// Generates a .STL file from a a provided surface mesh
	void WriteSurfaceAsSTLFile(const char* file_name, ModelPart& rSurfaceModelPart)
	{
		KRATOS_TRY;

		std::cout<<"\n::[Generating STL file]::"<<std::endl;

		// Write stl of surface model part
		std::ofstream myfile;
		myfile.open (file_name);
		myfile << "solid Layer0\n";
		for ( ModelPart::ConditionIterator cond_i = rSurfaceModelPart.ConditionsBegin(); cond_i != rSurfaceModelPart.ConditionsEnd(); ++cond_i )
		{
			array_1d<double,3> area_normal = cond_i->GetValue(NORMAL);
			myfile << "  facet normal " << area_normal[0] <<" " << area_normal[1] << " " << area_normal[2] <<"\n";
			myfile << "    outer loop\n";

			Element::GeometryType& rNodes = cond_i->GetGeometry();
			for(unsigned int i = 0; i<rNodes.size(); i++)
				myfile << "      vertex "<< rNodes[i].X() <<" " << rNodes[i].Y() << " " << rNodes[i].Z() <<"\n";

			myfile << "    end loop\n";
			myfile << "  end facet\n";

		}
		myfile << "endsolid Layer0\n";
		myfile.close();

		std::cout<<"  STL File succesfully generated under the name " << file_name <<std::endl;

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
		return "IOUtilities";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "IOUtilities";
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
	//IOUtilities& operator=(IOUtilities const& rOther);

	/// Copy constructor.
	//IOUtilities(IOUtilities const& rOther);


	///@}

}; // Class IOUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_IO_UTILITIES_H_INCLUDED */
