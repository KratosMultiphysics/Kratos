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
//                   Carlos Roig
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "input_output/logger.h"


namespace Kratos
{

	Logger::Logger(std::string const& TheLabel) : mCurrentMessage(TheLabel), mCallStack()
	{
	}

	Logger(std::string const& TheLabel, const CodeLocation& Location) : mCurrentMessage(TheLabel), mCallStack() {
		add_to_call_stack(Location);
	}

	Logger::~Logger()
	{
		auto outputs = GetOutputsInstance();
		#pragma omp critical
		{
			for (auto i_output = outputs.begin(); i_output != outputs.end(); i_output++)
				i_output->WriteMessage(mCurrentMessage);
		}
	}

	void Logger::AddOutput(LoggerOutput const& TheOutput)
	{
		#pragma omp critical
		{
		  GetOutputsInstance().push_back(TheOutput);
		}
	}

	const CodeLocation Logger::GetCurrentLocation() const {
		if(mCallStack.empty())
			return CodeLocation("Unknown File", "Unknown Location", 0);
		return mCallStack[0];
	}

	void Logger::AddToCallStack(CodeLocation const& TheLocation)
	{
		mCallStack.push_back(TheLocation);
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

	/// Manipulator stream function
	Logger& Logger::operator << (std::ostream& (*pf)(std::ostream&))
	{
		mCurrentMessage << pf;

		return *this;
	}

	/// char stream function
	Logger& Logger::operator << (const char * rString)
	{
		mCurrentMessage << rString;

		return *this;
	}

	// Location stream function
	Logger& Logger::operator << (CodeLocation const& TheLocation)
	{
		AddToCallStack(TheLocation);

		return *this;
	}

	/// Severity stream function
	Logger& Logger::operator << (Severity const& TheSeverity)
	{
		mCurrentMessage << TheSeverity;

		return *this;
	}

	/// Category stream function
	Logger& Logger::operator << (Category const& TheCategory) {
		mCurrentMessage << TheCategory;

		return *this;
	}



}  // namespace Kratos.
