#if !defined(KRATOS_MESHLESS_FORCE_INTERFACE_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_FORCE_INTERFACE_CONDITION_H_INCLUDED



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

	///@name Kratos Classes
	///@{

	/// Interface condition for IGA surface towards disrecte elements.
	/** The condition is designed to allow the interface between DEM-particles and
	*   the IGA-surface models.
	*/
	class MeshlessForceInterfaceCondition
		: public MeshlessBaseCondition
	{
	public:
		///@name Type Definitions
		///@{
		/// Counted pointer of MeshlessForceInterfaceCondition
		KRATOS_CLASS_POINTER_DEFINITION(MeshlessForceInterfaceCondition);

		///@}
		///@name Life Cycle
		///@{
		/// Default constructor.
		MeshlessForceInterfaceCondition(
			IndexType NewId, 
			GeometryType::Pointer pGeometry);
		MeshlessForceInterfaceCondition(
			IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties);

		/// Destructor.
		virtual ~MeshlessForceInterfaceCondition();
		///@}
		///@name Operations
		///@{

		Condition::Pointer Create(IndexType NewId, NodesArrayType const&
			ThisNodes, PropertiesType::Pointer pProperties) const override;

		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType&
			rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

		void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo&
			rCurrentProcessInfo) override;

		/**
		* @brief This function provides the place to perform checks on the completeness of the input.
		* @details It is designed to be called only once (or anyway, not often) typically at the beginning
		* of the calculations, so to verify that nothing is missing from the input
		* or that no common error is found.
		* @param rCurrentProcessInfo the current process info instance
		*/
		int Check(const ProcessInfo& rCurrentProcessInfo) override;

		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo&
			rCurrentProcessInfo) override;

		void GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo&
			CurrentProcessInfo) override;


		/** Sets on rValues the nodal velocities
		*   @param rValues: The values of velocities
		*   @param Step: The step to be computed
		*/
		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0
		) override;

		///@}
	protected:

	private:
		///@name Member Variables
		///@{
		friend class Serializer;

		// A private default constructor necessary for serialization
		MeshlessForceInterfaceCondition() : MeshlessBaseCondition()
		{}

		virtual void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
		}

		virtual void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
		}
		///@}
	}; // Class MeshlessForceInterfaceCondition
}  // namespace Kratos.

#endif // KRATOS_MESHLESS_FORCE_INTERFACE_CONDITION_H_INCLUDED  defined 