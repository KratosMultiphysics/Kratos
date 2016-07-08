//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <sstream>

// External includes


// Project includes
#include "input_output/logger_output.h"


namespace Kratos
{
	std::string LoggerOutput::Info() const
	{
		return "LoggerOutput";
	}

	/// Print information about this object.
	void LoggerOutput::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void LoggerOutput::PrintData(std::ostream& rOStream) const
	{
		rOStream << "Max Level : " << mMaxLevel;
	}

	/// char stream function
	LoggerOutput& LoggerOutput::operator << (const char * rString)
	{
		mrStream << rString;

		return *this;
	}

	LoggerOutput& LoggerOutput::operator << (std::ostream& (*pf)(std::ostream&))
	{
		std::stringstream buffer;
		pf(buffer);

		mrStream << buffer.str();

		return *this;
	}


	/// output stream function
	std::ostream& operator << (std::ostream& rOStream,
		const LoggerOutput& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}


}  // namespace Kratos.
