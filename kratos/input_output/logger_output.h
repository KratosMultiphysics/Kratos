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


#if !defined(KRATOS_LOGGER_OUTPUT_H_INCLUDED )
#define  KRATOS_LOGGER_OUTPUT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <map>


// External includes


// Project includes
#include "includes/define.h"
#include "input_output/logger_message.h"

namespace Kratos
{
		///@addtogroup KratosCore
		///@{

		///@name Kratos Classes
		///@{

		/// LoggerOutput is the base class for all logger outputs.
		/** LoggerOutput defines the interface for all logger outputs
			and also provides the basic (and default) functionalities
			to be extended in other outputs
		*/
		class KRATOS_API(KRATOS_CORE) LoggerOutput
		{
		public:
			///@name Type Definitions
			///@{


            /// Pointer definition of LoggerOutput
            KRATOS_CLASS_POINTER_DEFINITION(LoggerOutput);

			///@}
			///@name Enums
			///@{

			///@}
			///@name Life Cycle
			///@{

			LoggerOutput(std::ostream& rOutputStream) 
				: mrStream(rOutputStream), mMaxLevel(1), mSeverity(LoggerMessage::Severity::INFO), mCategory(LoggerMessage::Category::STATUS) {}

			LoggerOutput(LoggerOutput const& Other) 
				: mrStream(Other.mrStream), mMaxLevel(Other.mMaxLevel), mSeverity(Other.mSeverity), mCategory(Other.mCategory) {}

			/// Destructor.
			virtual ~LoggerOutput() {}


			///@}
			///@name Operators
			///@{

			LoggerOutput& operator=(LoggerOutput const& Other) = delete;

			///@}
			///@name Operations
			///@{


			///@}
			///@name Access
			///@{

			virtual void WriteHeader();

			virtual void WriteMessage(LoggerMessage const& TheMessage);

			virtual void Flush();

			void SetMaxLevel(std::size_t TheLevel) {
				mMaxLevel = TheLevel;
			}

			std::size_t GetMaxLevel() const {
				return mMaxLevel;
			}

			void SetSeverity(LoggerMessage::Severity const& TheSeverity) {
				mSeverity = TheSeverity;
			}

			LoggerMessage::Severity GetSeverity() const {
				return mSeverity;
			}

			void SetCategory(LoggerMessage::Category const& TheCategory) {
				mCategory = TheCategory;
			}

			LoggerMessage::Category GetCategory() const {
				return mCategory;
			}

			///@}
			///@name Inquiry
			///@{

			///@}
			///@name Input and output
			///@{

			/// Turn back information as a string.
			virtual std::string Info() const;

			/// Print information about this object.
			virtual void PrintInfo(std::ostream& rOStream) const;

			/// Print object's data.
			virtual void PrintData(std::ostream& rOStream) const;

			
			/// string stream function
			template<class StreamValueType>
			LoggerOutput& operator << (StreamValueType const& rValue)
			{
				std::stringstream buffer;
				buffer << rValue;

				mrStream << buffer.str();

				return *this;
			}

			/// Manipulator stream function
			LoggerOutput& operator << (std::ostream& (*pf)(std::ostream&));

			/// char stream function
			LoggerOutput& operator << (const char * rString);



			///@}
        protected:

            std::ostream& GetStream() {return mrStream;}

		private:
			///@name Life Cycle
			///@{

			///@}
			///@name Member Variables
			///@{

			std::ostream& mrStream;
			std::size_t mMaxLevel;
			LoggerMessage::Severity mSeverity;
			LoggerMessage::Category mCategory;

			///@}
		}; // Class LoggerOutput

	  ///@}

	  ///@name Input and output
	  ///@{

		/// output stream function
		std::ostream& operator << (std::ostream& rOStream,
			const LoggerOutput& rThis);

		///@}
		///@name macros
		///@{


		///@}

		///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_LOGGER_OUTPUT_H_INCLUDED  defined
