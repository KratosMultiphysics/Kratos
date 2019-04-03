//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Pooyan Dadvand
//
//

#if !defined(KRATOS_CODE_LOCATION_H_INCLUDED )
#define  KRATOS_CODE_LOCATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

#include "includes/kratos_export_api.h"

namespace Kratos
{
	///@addtogroup KratosCore
	///@{

	/// This class keeps a code location consist of filename, function name and line number.
	/// It also provides methods to get cleaned version of filename and function name.
	class KRATOS_API(KRATOS_CORE) CodeLocation
	{
	public:

		CodeLocation();

		CodeLocation(CodeLocation const & Other);

		CodeLocation(std::string const& FileName, std::string const& FunctionName, std::size_t LineNumber);

         ///@}
		///@name Private Operators
	    ///@{

        CodeLocation& operator=(CodeLocation const& Other) {
			mFileName = Other.mFileName;
			mFunctionName = Other.mFunctionName;
            mLineNumber = Other.mLineNumber;

            return *this;
        }


		///@name Operations
		///@{

		/// This function removes the path before the Kratos root
		std::string CleanFileName() const;

		/// This method cleans many template arguments and namespaces from the function name gives by compiler
		std::string CleanFunctionName() const;


		///@}
		///@name Access
		///@{

		const std::string& GetFileName() const;

		const std::string& GetFunctionName() const;

		int GetLineNumber() const;


		///@}

	private:
		///@name Member Variables
		///@{

		std::string mFileName;
		std::string mFunctionName;
		std::size_t mLineNumber;

		///@}
		///@name Private Operations
	    ///@{

		static void RemoveNamespace(std::string& FunctionName, const std::string& Namespace);

		static void ReduceTemplateArgumentsToFirstN(std::string& FunctionName, const std::string& TemplateName, std::size_t NumberOfArgumentsToKeep);

		static std::size_t GetNextPositionSkippingWhiteSpaces(std::string const& ThisString, std::size_t Position);

		static bool IsWhiteSpace(char C);

		static void ReplaceAll(std::string& ThisString, const std::string& FromString, const std::string& ToString);

		///@}

	}; // Class CodeLocation
	   ///@}
	   ///@name Input and output
	   ///@{

	/// output stream function
	std::ostream & operator <<(std::ostream& rOStream,
		const CodeLocation& rThis);
	///@}

#if defined(KRATOS_CODE_LOCATION)
#undef KRATOS_CODE_LOCATION
#endif

#if defined(KRATOS_CURRENT_FUNCTION)
#undef KRATOS_CURRENT_FUNCTION
#endif

#if defined(__PRETTY_FUNCTION__)
#define KRATOS_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__GNUC__)
#define KRATOS_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCTION__)
#define KRATOS_CURRENT_FUNCTION __FUNCTION__
#elif defined(__func__)
#define KRATOS_CURRENT_FUNCTION __func__
#else
#define KRATOS_CURRENT_FUNCTION "unknown function"
#endif


#define KRATOS_CODE_LOCATION Kratos::CodeLocation(__FILE__, KRATOS_CURRENT_FUNCTION, __LINE__)


  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CODE_LOCATION_H_INCLUDED  defined
