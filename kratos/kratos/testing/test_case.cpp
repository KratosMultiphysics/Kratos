//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
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


// External includes 


// Project includes
#include "testing/test_case.h"


namespace Kratos
{
	namespace Testing
	{
		TestCase::TestCase(std::string Name)
			:mName(Name), mIsEnambled(true) {}

		TestCase::~TestCase() {}

		void TestCase::Setup() {}

		void TestCase::Run()
		{
			try {
				Setup();
				TestFunction();
				TearDown();
			}
			catch (...) {
				std::cout << "Test " << Name() << " Failed!" << std::endl;
			}
		}

		void TestCase::Profile()
		{
			try {
				Setup();
				TestFunction();
				TearDown();
			}
			catch (...) {
				std::cout << "Test " << Name() << " Failed!" << std::endl;
			}
		}

		void TestCase::TearDown() {}

		void TestCase::Enable()
		{
			mIsEnambled = true;
		}

		void TestCase::Disable() 
		{
			mIsEnambled = false;
		}

		const std::string& TestCase::Name()
		{
			return mName;
		}

		bool TestCase::IsEnabled()
		{
			return mIsEnambled;
		}

		bool TestCase::IsDisabled()
		{
			return !mIsEnambled;
		}

		std::string TestCase::Info() const
		{
			return "Test case " + mName;
		}

		void TestCase::PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
		}

		void TestCase::PrintData(std::ostream& rOStream) const
		{
		}




	} // manespace Testing.
}  // namespace Kratos.


