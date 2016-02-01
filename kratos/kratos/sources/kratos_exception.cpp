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
#include "includes/kratos_exception.h"


namespace Kratos
{

	KratosException::KratosException() 
		: std::exception(), mWhat("Unknown Error"), mWhere("Unknown origin"), mMessage("Unknown Error")
	{
		mWhat.append("\n");
		mWhat.append(mWhere);
	}

	KratosException::KratosException(const std::string& rWhat ) 
		: std::exception(), mWhat(rWhat), mWhere("Unknown origin"), mMessage(rWhat)
	{
		mWhat.append("\n");
		mWhat.append(mWhere);
	}

	KratosException::KratosException(const std::string& rWhat, const std::string& rWhere) 
		: std::exception(), mWhat(rWhat), mWhere(rWhere), mMessage(rWhat)
	{
		mWhat.append("\n");
		mWhat.append(mWhere);
	}

	KratosException::KratosException(const KratosException& rOther) 
		: std::exception(rOther), mWhat(rOther.mWhat), mWhere(rOther.mWhere), mMessage(rOther.mMessage)
	{
	}


    KratosException::~KratosException() throw()
    {
    }

	void KratosException::append_message(std::string const& rMessage)
	{
		mMessage.append(rMessage);
		mWhat = mMessage;
		mWhat.append("\n");
		mWhat.append(mWhere);
	}

	void KratosException::append_where(std::string const& rWhere)
	{
		mWhere.append("\n");
		mWhere.append(rWhere);
		mWhat.append("\n");
		mWhat.append(rWhere);
	}

    const char* KratosException::what() const throw() //noexcept
	{
		return mWhat.c_str();
	}

	const std::string& KratosException::where() const
	{
		return mWhere;
	}

	const std::string& KratosException::message() const
	{
		return mMessage;
	}

	std::string KratosException::Info() const
	{
		return "KratosException";
	}
      
      /// Print information about this object.
    void KratosException::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}
      /// Print object's data.
    void KratosException::PrintData(std::ostream& rOStream) const
	{
		rOStream << "Error: " << mWhat << std::endl;
		rOStream << "   in: " << mWhere;
	}

	/// input stream function
	std::istream& operator >> (std::istream& rIStream,
		KratosException& rThis)
	{
		return rIStream;
	}

	/// output stream function
  /*	inline std::ostream& operator << (std::ostream& rOStream,
		const KratosException& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
  */
	/// char stream function
    KratosException& KratosException::operator << (const char * rString)
	{
        append_message(rString);

        return *this;
	}

    KratosException& KratosException::operator << (std::ostream& (*pf)(std::ostream&))
	{
		std::stringstream buffer;
		pf(buffer);

        append_message(buffer.str());

        return *this;
	}

}  // namespace Kratos.


