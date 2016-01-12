//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Pooyan Dadvand 
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "utilities/logger.h"


namespace Kratos
{

	Logger::Logger()
	{
	}

	Logger::~Logger()
	{
	}

    std::string Logger::Info() const
	{
		return "Logger";
	}
      
      /// Print information about this object.
    void Logger::PrintInfo(std::ostream& rOStream) const
	{
	}
      /// Print object's data.
    void Logger::PrintData(std::ostream& rOStream) const
	{
	}
	  
	std::string Logger::CleanFunctionName(const std::string& FunctionName, const std::string& FileName, int LineNumber)
	{
		std::stringstream buffer;
		buffer << Filter(FunctionName) + " [ " + FileName +  " , Line " << LineNumber << " ] ";
		return buffer.str();
	}
      
	std::string Logger::Filter(const std::string& ThisString)
	{
		std::string buffer(ThisString);

		ReplaceAll(buffer, "Kratos::", ""); 
		ReplaceAll(buffer, "__cdecl", "");
		ReplaceAll(buffer, "class", "");
		ReplaceAll(buffer, "Dof<double>", "Dof");
		ReplaceAll(buffer, "Node<3, Dof >", "Node");
		ReplaceAll(buffer, "Point<3,double>", "Point");
		ReplaceAll(buffer, "boost::", "");
		ReplaceAll(buffer, "numeric::", "");
		ReplaceAll(buffer, "std::allocator<double>", "");
		ReplaceAll(buffer, "std::allocator< Point >", "");
		ReplaceAll(buffer, "<double,  >", "<double>");


		return buffer;
	}
 	std::string Logger::ReplaceAll(std::string& ThisString, const std::string& FromString, const std::string& ToString)
	{
		std::size_t start_position = 0;
		while((start_position = ThisString.find(FromString,start_position)) != std::string::npos)
		{
			ThisString.replace(start_position, FromString.length(), ToString);
			start_position += ToString.length(); // ...
		}
		return ThisString;
	}
      
  
}  // namespace Kratos.


