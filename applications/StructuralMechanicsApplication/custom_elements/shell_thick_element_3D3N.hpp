// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//       contact:    A.Winterstein@tum.de
//					 

#if !defined(SHELL_THICK_ELEMENT_3D3N_H_INCLUDED )
#define  SHELL_THICK_ELEMENT_3D3N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "custom_utilities/shellt3_local_coordinate_system.hpp"
#include "utilities/quaternion.h"

namespace Kratos
{


	///@name Kratos Globals
	///@{
	///@}

	///@name Type Definitions
	///@{
	///@}

	class ShellT3_CoordinateTransformation;

	///@name  Enum's
	///@{
	///@}

	///@name  Functions
	///@{
	///@}

	///@name Kratos Classes
	///@{

	/** \brief ShellThickElement3D3N
	*
	* This element represents a 3-node Shell element
	* based on the Discrete Shear Gap theory (DSG) by Bletzinger.
	* This element is formulated for small strains,
	* but can be used in Geometrically nonlinear problems
	* involving large displacements and rotations
	* using a Corotational Coordinate Transformation.
	* Material nonlinearity is handled by means of the cross section object.
	*/

	/*
	Shell formulation reference:
	1.	Bletzinger, K.U., Bischoff, M. and Ramm, E., 2000. A unified approach for 
		shear-locking-free triangular and rectangular shell finite elements. 
		Computers & Structures, 75(3), pp.321-334.
	2.	Rama, G.,  Marinkovic, D.,  Zehn, M., 2016. Efficient co-rotational 
		3-node shell element.
		American Journal of Engineering and Applied Sciences, Volume 9, Issue 2, 
		Pages 420-431.
	*/


	class ShellThickElement3D3N : public Element
	{
	public:

		///@name Type Definitions
		///@{

		KRATOS_CLASS_POINTER_DEFINITION(ShellThickElement3D3N);

		typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

		typedef ShellT3_CoordinateTransformation CoordinateTransformationBaseType;

		typedef boost::shared_ptr<CoordinateTransformationBaseType> CoordinateTransformationBasePointerType;

		typedef array_1d<double, 3> Vector3Type;

		typedef Quaternion<double> QuaternionType;

		///@}

		///@name Classes
		///@{

		///@}

		///@name Life Cycle
		///@{

		ShellThickElement3D3N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			bool NLGeom = false);

		ShellThickElement3D3N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties,
			bool NLGeom = false);

		ShellThickElement3D3N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties,
			CoordinateTransformationBasePointerType pCoordinateTransformation);

		~ShellThickElement3D3N() override;

		///@}

		///@name Operations
		///@{

		// Basic

		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

		IntegrationMethod GetIntegrationMethod() const override;

		void Initialize() override;

		void ResetConstitutiveLaw() override;

		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

		void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;
		
		int Check(const ProcessInfo& rCurrentProcessInfo) override;

		void GetValuesVector(Vector& values, int Step = 0) override;

		void GetFirstDerivativesVector(Vector& values, int Step = 0) override;

		void GetSecondDerivativesVector(Vector& values, int Step = 0) override;

		void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;	//corotational formulation

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override; //corotational formulation

		void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override; //corotational formulation

		void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override; //corotational formulation

		void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

		void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;
		

		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		
		void CalculateRightHandSide(VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		// Results calculation on integration points

		
		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
		
		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6> >& rVariable, std::vector<array_1d<double, 6> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		// More results calculation on integration points to interface with python
		void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
			std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
			std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateOnIntegrationPoints(const Variable<array_1d<double,
			3> >& rVariable, std::vector<array_1d<double, 3> >& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateOnIntegrationPoints(const Variable<array_1d<double,
			6> >& rVariable, std::vector<array_1d<double, 6> >& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		// Calculate functions
		void Calculate(const Variable<Matrix >& rVariable,
			Matrix& Output,
			const ProcessInfo& rCurrentProcessInfo) override;

		void Calculate(const Variable<double>& rVariable,
			double& Output,
			const ProcessInfo& rCurrentProcessInfo) override;

		///@}

		///@name Public specialized Access - Temporary
		///@{
		
		void SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections);
		
		///@}

	protected:

		///@name Protected Lyfe Cycle
		///@{

		/**
		* Protected empty constructor
		*/
		ShellThickElement3D3N() : Element()
		{
		}

		///@}

	private:

		///@name Private Classes
		///@{

		class CalculationData
		{

		public:

			// ---------------------------------------
			// calculation-constant data
			// ----------------------------------------
			// these data are allocated and constructed
			// at the beginning of the calculation

			ShellT3_LocalCoordinateSystem LCS0; /*!< reference coordinate system */
			ShellT3_LocalCoordinateSystem LCS;  /*!< current coordinate system */

			double dA;
			double hMean;
			double TotalArea;

			std::vector< array_1d<double, 3> > gpLocations;

			MatrixType dNxy = ZeroMatrix(3, 2);  /*!< shape function cartesian derivatives */
			VectorType N = ZeroVector(3); /*!< shape function vector at the current integration point */

			VectorType globalDisplacements = ZeroVector(18); /*!< global displacement vector */
			VectorType localDisplacements = ZeroVector(18);  /*!< local displacement vector */

			bool CalculateRHS; /*!< flag for the calculation of the right-hand-side vector */
			bool CalculateLHS; /*!< flag for the calculation of the left-hand-side vector */

			const bool parabolic_composite_transverse_shear_strains = false;

			// ---------------------------------------
			// Testing flags
			// ---------------------------------------
			// These should all be FALSE unless you are testing, or 
			// investigating the effects of element enhancements!

			const bool basicTriCST = false;	// bool to use basic CST 
			// displacement-based shear formulation. This should be FALSE unless 
			// you are testing

			const bool ignore_shear_stabilization = false; // bool to 
			// ignore stabilizing the transverse shear part of the material 
			// matrix. This should be false unless you are testing

			const bool smoothedDSG = false; // bool to use smoothed DSG 
			// formulation according to [Nguyen-Thoi et al., 2013]. 
			// This should be false unless you are testing

			const bool specialDSGc3 = false; // bool to use experimental
			// DSGc3 formulation not yet complete.
			// This should be false unless you are testing

			// ---------------------------------------
			// calculation-variable data
			// ---------------------------------------
			// these data are updated during the
			// calculations

			size_t gpIndex;

			// ---------------------------------------
			// calculation-variable data
			// ---------------------------------------
			// these data are updated during the
			// calculations, but they are allocated
			// only once(the first time they are used)
			// to avoid useless re-allocations

			MatrixType B = ZeroMatrix(8, 18);   /*!< total strain-displacement matrix at the current integration point */

			double h_e;		/*!< longest edge of triangle */
			double alpha = 0.1;	// modifier of shear material matrix stabilization parameter
								// refer Lyly(1993)
			double shearStabilisation;

			Matrix D = ZeroMatrix(8, 8);		/*!< section constitutive matrix at the current integration point */

			VectorType generalizedStrains = ZeroVector(8);  /*!< generalized strain vector at the current integration point */

			VectorType generalizedStresses = ZeroVector(8); /*!< generalized stress vector at the current integration point */

			ShellCrossSection::SectionParameters SectionParameters; /*!< parameters for cross section calculations */

			std::vector<VectorType> rlaminateStrains;

			std::vector<VectorType> rlaminateStresses;

		public:

			const ProcessInfo& CurrentProcessInfo;

		public:

			CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
				const ProcessInfo& rCurrentProcessInfo);

		};

		///@}

		///@name Private Operations
		///@{
		void CalculateStressesFromForceResultants
			(VectorType& rstresses, 
				const double& rthickness);

		void CalculateLaminaStrains(CalculationData& data);

		void CalculateLaminaStresses(CalculationData& data);

		double CalculateTsaiWuPlaneStress(const CalculationData& data, const Matrix& rLamina_Strengths, const unsigned int& rCurrent_Ply);

		void CalculateVonMisesStress(const CalculationData& data, const Variable<double>& rVariable, double& rVon_Mises_Result);

		void CalculateShellElementEnergy(const CalculationData& data, const Variable<double>& rVariable, double& rEnergy_Result);

		void CheckGeneralizedStressOrStrainOutput(const Variable<Matrix>& rVariable, int& iJob, bool& bGlobal);

		void DecimalCorrection(Vector& a);

		void SetupOrientationAngles();

		void CalculateSectionResponse(CalculationData& data);

		void InitializeCalculationData(CalculationData& data);

		void CalculateDSGc3Contribution(CalculationData& data, MatrixType& rLeftHandSideMatrix);

		void CalculateSmoothedDSGBMatrix(CalculationData& data);

		void CalculateDSGShearBMatrix(Matrix& shearBMatrix, const double& a, const double& b, const double& c, const double& d, const double& A);

		void AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector);

		void CalculateAll(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo,
			const bool LHSrequired,
			const bool RHSrequired);

		bool TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double, 3> >& rVariable,
			std::vector<array_1d<double, 3> >& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		bool TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
			std::vector<Matrix>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		///@}

		///@name Static Member Variables
		///@{
		///@}

		///@name Member Variables
		///@{

		CoordinateTransformationBasePointerType mpCoordinateTransformation; /*!< The Coordinate Transformation */

		CrossSectionContainerType mSections; /*!< Container for cross section associated to each integration point */

		IntegrationMethod mThisIntegrationMethod; /*!< Currently selected integration method */

		double mOrthotropicSectionRotation = 0.0; /*!< In-plane rotation angle for orthotropic section */

												  ///@}

												  ///@name Serialization
												  ///@{
		void printMatrix(Matrix& matrixIn, std::string stringIn);


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
#endif // SHELL_THICK_ELEMENT_3D3N_H_INCLUDED