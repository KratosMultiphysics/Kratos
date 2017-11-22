// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder
//

#if !defined(SHELL_THIN_ADJOINT_ELEMENT_3D4N_H_INCLUDED )
#define  SHELL_THIN_ADJOINT_ELEMENT_3D4N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_elements/shell_thin_element_3D4N.hpp"
//#include "custom_utilities/shell_cross_section.hpp"
#include "utilities/quaternion.h"
//#include "custom_utilities/shellq4_local_coordinate_system.hpp"

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


	class ShellThinAdjointElement3D4N : public ShellThinElement3D4N
	{
	public:

		///@name Type Definitions
		///@{
		KRATOS_CLASS_POINTER_DEFINITION(ShellThinAdjointElement3D4N);

     	typedef Element::PropertiesType PropertiesType;

     	typedef Element::DofsArrayType DofsArrayType;

		typedef Element::NodesArrayType NodesArrayType;

		typedef Element::IndexType IndexType;

		///@}

		///@name Life Cycle
		///@{
		ShellThinAdjointElement3D4N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			bool NLGeom = false);

		ShellThinAdjointElement3D4N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties,
			bool NLGeom = false);

		ShellThinAdjointElement3D4N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties,
			CoordinateTransformationBasePointerType pCoordinateTransformation);

		~ShellThinAdjointElement3D4N() override;

		///@}

		///@name Operations
		///@{
		// Basic

		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

		Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

		void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

		int Check(const ProcessInfo& rCurrentProcessInfo) override;

		void GetValuesVector(Vector& values, int Step = 0) override;


		double GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rVariable);

		double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable);

    	void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, 
											const ProcessInfo& rCurrentProcessInfo) override;
	
    	void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput, 
											const ProcessInfo& rCurrentProcessInfo) override;

    	void Calculate(const Variable<Vector >& rVariable, Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;                                         

		void Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;   

    	void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);    

    	void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, const Variable<Vector>& rStressVariable,
                                        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);
	
    	void CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable, 
                                            const Variable<Vector>& rStressVariable,
                                             Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);                       

		void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;    

		// Results calculation on integration points
		void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
			std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


		void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
			std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		//void Initialize() override;
		///@}



		///@name Public specialized Access - Temporary
		///@{
		///@}

	protected:

		///@name Protected Lyfe Cycle
		///@{
		/**
		* Protected empty constructor
		*/
		ShellThinAdjointElement3D4N() : ShellThinElement3D4N()
		{
		}

		///@}

	private:

		///@name Private Classes
		///@{


		///@}

		///@name Private Operations
		///@{
		

		///@}

		///@name Static Member Variables
		///@{
		///@}

		///@name Member Variables
		///@{


		///@}

		///@name Serialization
		///@{
		friend class Serializer;

		void save(Serializer& rSerializer) const override;

		void load(Serializer& rSerializer) override;

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
#endif // SHELL_THIN_ADJOINT_ELEMENT_3D4N_H_INCLUDED