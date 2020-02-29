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


#if !defined(KRATOS_TEST_CASE_RESULT_H_INCLUDED )
#define  KRATOS_TEST_CASE_RESULT_H_INCLUDED

#include <string>
#include <iostream>

#include "includes/kratos_export_api.h"

namespace Kratos
{
	namespace Testing
	{
		///@addtogroup KratosCore
		///@{

		///@name Kratos Classes
		///@{

		/// The test case base class.
		/** Defines the interface for all the test cases and also fixtures
		*/
		class KRATOS_API(KRATOS_CORE) TestCaseResult
		{
		public:
			///@name Type Definitions
			///@{

			///@}
			///@name Life Cycle
			///@{

			/// TestCaseResult default constructors rest variables to before run
			TestCaseResult();

			/// Copy constructor
			TestCaseResult(TestCaseResult const& rOther);

			/// Destructor.
			virtual ~TestCaseResult();


			///@}
			///@name Operators
			///@{

			/// Assignment operator
			TestCaseResult& operator=(TestCaseResult const& rOther);

			///@}
			///@name Operations
			///@{

			/// Reset the results to before run state
			virtual void Reset();


			///@}
			///@name Access
			///@{

			void SetToSucceed();

			void SetToFailed();

			void SetToSkipped();

			void SetOutput(const std::string& TheOutput);

			const std::string& GetOutput() const;

			void SetErrorMessage(const std::string& TheMessage);

			const std::string& GetErrorMessage() const;

			void SetSetupElapsedTime(double ElapsedTime);

			double GetSetupElapsedTime() const;

			void SetRunElapsedTime(double ElapsedTime);

			double GetRunElapsedTime() const;

			void SetTearDownElapsedTime(double ElapsedTime);

			double GetTearDownElapsedTime() const;

			void SetElapsedTime(double ElapsedTime);

			double GetElapsedTime() const;



			///@}
			///@name Inquiry
			///@{

			bool IsSucceed() const;

			bool IsFailed() const;

			bool IsSkipped() const;

			bool IsRun() const;

			///@}
			///@name Input and output
			///@{

			/// Turn back information as a string.
			virtual std::string Info() const;

			/// Print information about this object.
			virtual void PrintInfo(std::ostream& rOStream) const;

			/// Print object's data.
			virtual void PrintData(std::ostream& rOStream) const;


			///@}
			///@name Friends
			///@{


			///@}


		private:

			///@name Type definitions
			///@{

			enum class Result {
				DidNotRun,
				Passed,
				Failed,
				Skipped
			};

			///@}
			///@name Static Member Variables
			///@{


			///@}
			///@name Member Variables
			///@{

			Result mStatus;
			std::string mOutput;
			std::string mErrorMessage;
			double mSetupElapsedTime;
			double mRunElapsedTime;
			double mTearDownElapsedTime;
			double mElapsedTime;


			///@}

		}; // Class TestCaseResult

	  ///@}

	  ///@name Input and output
	  ///@{

		/// output stream function
		inline std::ostream& operator << (std::ostream& rOStream,
			const TestCaseResult& rThis)
		{
			rThis.PrintInfo(rOStream);
			rOStream << std::endl;
			rThis.PrintData(rOStream);

			return rOStream;
		}
		///@}

		///@} addtogroup block
	} // manespace Testing.
}  // namespace Kratos.

#endif // KRATOS_TEST_CASE_RESULT_H_INCLUDED  defined
