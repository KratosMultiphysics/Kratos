#if !defined(KRATOS_MESHLESS_POINT_CONDITION_H_INCLUDED )
#define  KRATOS_MESHLESS_POINT_CONDITION_H_INCLUDED



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

	/// Condition to check values at specific points.
	/** The condition is designed to allow the checking of values concerning the
	*   displacements at the specific point
	*/
	class MeshlessPointCondition
		: public MeshlessBaseCondition
	{
	public:
		///@name Type Definitions
		///@{
		/// Counted pointer of MeshlessPointCondition
		KRATOS_CLASS_POINTER_DEFINITION(MeshlessPointCondition);

		///@}
		///@name Life Cycle
		///@{
		// Constructor void
		MeshlessPointCondition() : MeshlessBaseCondition()
		{};

		// Constructor using an array of nodes
		MeshlessPointCondition(
			IndexType NewId, 
			GeometryType::Pointer pGeometry)
			: MeshlessBaseCondition(NewId, pGeometry){};

		// Constructor using an array of nodes with properties
		MeshlessPointCondition(
			IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties)
			: MeshlessBaseCondition(NewId, pGeometry, pProperties){};

		/// Destructor.
		virtual ~MeshlessPointCondition() {};
		///@}
		///@name Operations
		///@{
		Condition::Pointer Create(IndexType NewId, NodesArrayType const&
			ThisNodes, PropertiesType::Pointer pProperties) const override
		{
			return MeshlessBaseCondition::Pointer(new MeshlessPointCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
		};

		/**
		* @brief Sets on rResult the ID's of the element degrees of freedom
		* @param rResult The vector containing the equation id
		* @param rCurrentProcessInfo The current process info instance
		*/
		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo&
			rCurrentProcessInfo) override;

		/**
		* @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
		* @param rElementalDofList The vector containing the dof of the element
		* @param rCurrentProcessInfo The current process info instance
		*/
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


		/**
		* @brief Calculate a Vector Variable on the Element Constitutive Law
		* @param rVariable The variable we want to get
		* @param rOutput The values obtained int the integration points
		* @param rCurrentProcessInfo the current process info instance
		*/
		void CalculateOnIntegrationPoints(
			const Variable<array_1d<double, 3>>& rVariable,
			std::vector<array_1d<double, 3>>& rOutput,
			const ProcessInfo& rCurrentProcessInfo
		) override;


		/**
		* @brief Get on rVariable a Vector Value from the Condition
		* @param rVariable The variable we want to get
		* @param rValues The results in the integration points
		* @param rCurrentProcessInfo the current process info instance
		*/
		void GetValueOnIntegrationPoints(
			const Variable<array_1d<double, 3>>& rVariable,
			std::vector<array_1d<double, 3>>& rValues,
			const ProcessInfo& rCurrentProcessInfo
		) override;
		///@}
	protected:

	private:
		///@name Member Variables
		///@{
		friend class Serializer;

		virtual void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
		}

		virtual void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
		}
		///@}
	}; // Class MeshlessPointCondition
}  // namespace Kratos.

#endif // KRATOS_MESHLESS_POINT_CONDITION_H_INCLUDED  defined 