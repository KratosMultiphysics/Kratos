//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Pooyan Dadvand 
//                   
//

#include <sstream>

#include "includes/exception.h"


namespace Kratos
{

	Exception::Exception()
		: std::exception(), mMessage("Unknown Error")
		, mCallStack()
	{
	}

	Exception::Exception(const std::string& rWhat ) 
		: std::exception(), mMessage(rWhat), mCallStack()
	{
	}

	Exception::Exception(const std::string& rWhat, const CodeLocation& Location) 
		: std::exception(), mMessage(rWhat), mCallStack()

	{
		mCallStack.push_back(Location);
	}

	Exception::Exception(const Exception& Other) 
		: std::exception(Other), mMessage(Other.mMessage), mCallStack(Other.mCallStack)
	{
	}


    Exception::~Exception() throw()
    {
    }

	void Exception::append_message(std::string const& rMessage)
	{
		mMessage.append(rMessage);
	}

	void Exception::add_to_call_stack(CodeLocation const& TheLocation)
	{
		mCallStack.push_back(TheLocation);
	}

    const char* Exception::what() const throw() //noexcept
	{
		return mMessage.c_str(); // TODO: << mCallStack[0];
	}

	const CodeLocation Exception::where() const
	{
		if(mCallStack.empty())
			return CodeLocation("Unknown File", "Unknown Location", 0);
		return mCallStack[0];
	}

	const std::string& Exception::message() const
	{
		return mMessage;
	}

	std::string Exception::Info() const
	{
		return "Exception";
	}
      
      /// Print information about this object.
    void Exception::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}
      /// Print object's data.
    void Exception::PrintData(std::ostream& rOStream) const
	{
		rOStream << "Error: " << mMessage << std::endl;
		rOStream << "   in: " << where();
	}


	Exception& Exception::operator << (CodeLocation const& TheLocation)
	{
		add_to_call_stack(TheLocation);

		return *this;
	}

	/// input stream function
	std::istream& operator >> (std::istream& rIStream,
		Exception& rThis)
	{
		return rIStream;
	}


	/// char stream function
    Exception& Exception::operator << (const char * rString)
	{
        append_message(rString);

        return *this;
	}

    Exception& Exception::operator << (std::ostream& (*pf)(std::ostream&))
	{
		std::stringstream buffer;
		pf(buffer);

        append_message(buffer.str());

        return *this;
	}

	std::ostream& operator << (std::ostream& rOStream, const Exception& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}


}  // namespace Kratos.


