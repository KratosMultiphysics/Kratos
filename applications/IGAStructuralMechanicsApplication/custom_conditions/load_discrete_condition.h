#if !defined(KRATOS_MESHLESS_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_LOAD_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "structural_mechanics_application_variables.h"

#include "iga_structural_mechanics_application_variables.h"
#include "custom_conditions/meshless_base_condition.h"


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

	/// Short class definition.
	/** Detail class definition.
	*/
	class LoadDiscreteCondition
		: public MeshlessBaseCondition
	{
	public:
		///@name Type Definitions
		///@{

		/// Counted pointer of LoadDiscreteCondition
		KRATOS_CLASS_POINTER_DEFINITION(LoadDiscreteCondition);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		LoadDiscreteCondition(
			IndexType NewId, 
			GeometryType::Pointer pGeometry);
		LoadDiscreteCondition(
			IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties);

		/// Destructor.
		virtual ~LoadDiscreteCondition();


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{

		Condition::Pointer Create(IndexType NewId, NodesArrayType const&
			ThisNodes, PropertiesType::Pointer pProperties) const override;

		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType&
			rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

		void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo&
			rCurrentProcessInfo) override;
		//virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo&
			rCurrentProcessInfo) override;

		void GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo&
			CurrentProcessInfo) override;

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
		//      virtual String Info() const;

		/// Print information about this object.
		//      virtual void PrintInfo(std::ostream& rOStream) const;

		/// Print object's data.
		//      virtual void PrintData(std::ostream& rOStream) const;


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

		friend class Serializer;

		// A private default constructor necessary for serialization
		LoadDiscreteCondition() : MeshlessBaseCondition()
		{
		}

		virtual void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
		}

		virtual void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
		}

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
		//LoadDiscreteCondition& operator=(const LoadDiscreteCondition& rOther);

		/// Copy constructor.
		//LoadDiscreteCondition(const LoadDiscreteCondition& rOther);

		///@}

	}; // Class LoadDiscreteCondition

	   ///@}

	   ///@name Type Definitions
	   ///@{

	   ///@}
	   ///@name Input and output
	   ///@{

	   /// input stream function
	   /*  inline std::istream& operator >> (std::istream& rIStream,
	   LoadDiscreteCondition& rThis);
	   */
	   /// output stream function
	   /*  inline std::ostream& operator << (std::ostream& rOStream,
	   const LoadDiscreteCondition& rThis)
	   {
	   rThis.PrintInfo(rOStream);
	   rOStream << std::endl;
	   rThis.PrintData(rOStream);

	   return rOStream;
	   }*/
	   ///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_LOAD_CONDITION_H_INCLUDED  defined 