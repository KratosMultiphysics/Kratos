//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined( CONVDIFF_INTERFACE_ELEMENT_H_INCLUDED )
#define  CONVDIFF_INTERFACE_ELEMENT_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"

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

	/** \brief ConvDiffInterfaceElement
	*
	* Da aggiungere
	*/
	class ConvDiffInterfaceElement : public Element
	{
	public:

		///@name Type Definitions
		///@{

		KRATOS_CLASS_POINTER_DEFINITION(ConvDiffInterfaceElement);

		///@}

		///@name Classes
		///@{

		struct InterfaceIndexPermutation
		{
			bool HasPermutation_disp;
			std::vector<SizeType> Permutation_disp;
			bool HasPermutation_temp;
			std::vector<SizeType> Permutation_temp;
		};

    ///@}

		///@name Life Cycle
		///@{

		ConvDiffInterfaceElement(IndexType NewId,
			GeometryType::Pointer pGeometry);

		ConvDiffInterfaceElement(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties);

		virtual ~ConvDiffInterfaceElement();

		///@}

		///@name Operations
		///@{

		// Basic

		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

		void Initialize();

		void ResetConstitutiveLaw();

		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

		void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo);

		int Check(const ProcessInfo& rCurrentProcessInfo);

		void CleanMemory();

		void GetValuesVector(Vector& values, int Step = 0);

		void GetFirstDerivativesVector(Vector& values, int Step = 0);

		void GetSecondDerivativesVector(Vector& values, int Step = 0);

		void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

		void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

		void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

		void CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);

		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateRightHandSide(VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		// Results calculation on integration points

		void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
										  std::vector<double>& rOutput,
										  const ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable< array_1d< double, 3 > >& rVariable,
										  std::vector< array_1d<double, 3 > >& rOutput,
										  const ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable< Vector >& rVariable,
										  std::vector< Vector >& rOutput,
										  const ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
										  std::vector< Matrix >& rOutput,
										  const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, 
										 std::vector<double>& rValues, 
										 const ProcessInfo& rCurrentProcessInfo);
	
		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, 
										 std::vector<Vector>& rValues, 
										 const ProcessInfo& rCurrentProcessInfo);
	
		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, 
										 std::vector<Matrix>& rValues, 
										 const ProcessInfo& rCurrentProcessInfo);
	
		void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, 
										 std::vector<array_1d<double,3> >& rValues, 
										 const ProcessInfo& rCurrentProcessInfo);
	
		void GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable, 
										 std::vector<array_1d<double,6> >& rValues, 
										 const ProcessInfo& rCurrentProcessInfo);
	
	///@}

	protected:

		///@name Protected Lyfe Cycle
		///@{

		/**
		* Protected empty constructor
		*/
		ConvDiffInterfaceElement() : Element()
		{
		}

		///@}

	private:

		///@name Private Operations
		///@{

		void DecimalCorrection(Vector& a);

		void InitializeContactData();



		virtual void CalculatePermutation_disp(InterfaceIndexPermutation& p);

		virtual void CalculatePermutation_temp(InterfaceIndexPermutation& p);

		virtual void CalculateDeltaPosition(Matrix& rDeltaPosition);

		virtual void CalculateJacobianAndTransformationMatrix(const SizeType pointID,
			Matrix& delta_position,
			Matrix& jacobian,
			double& J,
			Matrix& iR);

		virtual double CalculateIntegrationWeight(double J,
			double iw);

		virtual void CalculateLocalDisplacementVector(const InterfaceIndexPermutation& P,
			const Matrix& R,
			const Vector& UG,
			Vector& UL);

		virtual void CalculateLocalTemperature(const InterfaceIndexPermutation& P,
			const Matrix& R,
			const Vector& TG,
			Vector& TL);

		virtual void TransformToGlobalAndAdd(const InterfaceIndexPermutation& P,
			const Matrix& R,
			const Matrix& LHS_local,
			const Vector& RHS_local,
			Matrix& LHS_global,
			Vector& RHS_global,
			const bool LHSrequired,
			const bool RHSrequired);

		virtual void CalculateBMatrix(const SizeType pointID,
			Matrix& B);

		virtual void CalculateBMatrixU(const SizeType pointID,
			Matrix& B);

		virtual void CalculateGeneralizedTemperatureStrains(const SizeType pointID,
			const Matrix& B,
			const Vector& T,
			Vector& generalizedTemperatureStrains);

		virtual void CalculateGeneralizedDeltaDisplacement(const SizeType pointID,
			const Matrix& B,
			const Vector& U,
			Vector& deltaDispl);

		virtual void GetDisplacementVector(Vector& values);

	    void CalculateAll(MatrixType& rLeftHandSideMatrix,
	                      VectorType& rRightHandSideVector,
	                      ProcessInfo& rCurrentProcessInfo,
	                      const bool LHSrequired,
	                      const bool RHSrequired);

	    ///@}

	    ///@name Static Member Variables
	    ///@{
	    ///@}
    
	    ///@name Member Variables
	    ///@{

		bool mInitialized;
		std::vector< ConstitutiveLaw::Pointer > mConstitutiveLawVector;

		std::vector<Vector> m_init_gradT;
		Vector m_delta_disp;
		Vector mgeneralizedTemperatureStresses;

		///@}

		///@name Serialization
		///@{

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const;

		virtual void load(Serializer& rSerializer);

		///@}

		///@name Private  Access
		///@{
		///@}

		///@name Private Inquiry
		///@{
		///@}

		///@name Un accessible methods
		///@{
		///@}

	};

}
#endif // CONVDIFF_INTERFACE_ELEMENT_H_INCLUDED
