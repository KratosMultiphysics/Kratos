/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(SHELL_THIN_ELEMENT_3D3N_H_INCLUDED )
#define  SHELL_THIN_ELEMENT_3D3N_H_INCLUDED


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

/** \brief ShellThinElement3D3N
*
* This element represents a 3-node Shell element
* based on the Assumed Natural DEviatoric Strain (ANDES) by Felippa.
* This element is formulated for small strains, 
* but can be used in Geometrically nonlinear problems
* involving large displacements and rotations 
* using a Corotational Coordinate Transformation.
* Material nonlinearity is handled by means of the cross section object.
*/
class ShellThinElement3D3N : public Element
{
public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(ShellThinElement3D3N);
    
    typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

    typedef ShellT3_CoordinateTransformation CoordinateTransformationBaseType;

	typedef boost::shared_ptr<CoordinateTransformationBaseType> CoordinateTransformationBasePointerType;

    typedef array_1d<double, 3> Vector3Type;

    typedef Quaternion<double> QuaternionType;

    ///@}

	///@name Classes
    ///@{

	// TODO: Add Calulation Data
	
	///@}

    ///@name Life Cycle
    ///@{

    ShellThinElement3D3N(IndexType NewId, 
                         GeometryType::Pointer pGeometry, 
                         bool NLGeom = false);
    
    ShellThinElement3D3N(IndexType NewId, 
                         GeometryType::Pointer pGeometry, 
                         PropertiesType::Pointer pProperties, 
                         bool NLGeom = false);

    ShellThinElement3D3N(IndexType NewId, 
                         GeometryType::Pointer pGeometry, 
                         PropertiesType::Pointer pProperties, 
                         CoordinateTransformationBasePointerType pCoordinateTransformation);

    virtual ~ShellThinElement3D3N();

    ///@}
    
    ///@name Operations
    ///@{

    // Basic

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    IntegrationMethod GetIntegrationMethod() const;
    
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

    void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

    void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo);

    // Results calculation on integration points

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

	void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

	void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo);

	void GetValueOnIntegrationPoints(const Variable<array_1d<double,6>>& rVariable, std::vector<array_1d<double,6>>& rValues, const ProcessInfo& rCurrentProcessInfo);

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
    ShellThinElement3D3N() : Element()
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

		MatrixType L;

		MatrixType Q1;
		MatrixType Q2;
		MatrixType Q3;

		MatrixType Te;
		MatrixType TTu;

		double dA;
		double hMean;
		double TotalArea;
		double TotalVolume;
		std::vector< array_1d<double,3> > gpLocations;

		MatrixType dNxy; /*!< shape function cartesian derivatives */

		VectorType globalDisplacements; /*!< global displacement vector */
        VectorType localDisplacements;  /*!< local displacement vector */

		bool CalculateRHS; /*!< flag for the calculation of the right-hand-side vector */
		bool CalculateLHS; /*!< flag for the calculation of the left-hand-side vector */

		// ---------------------------------------
		// calculation-variable data
		// ---------------------------------------
		// these data are updated during the
		// calculations

		double beta0;
		size_t gpIndex;

		// ---------------------------------------
		// calculation-variable data
		// ---------------------------------------
		// these data are updated during the
		// calculations, but they are allocated
		// only once(the first time they are used)
		// to avoid useless re-allocations

		MatrixType B;   /*!< total strain-displacement matrix at the current integration point */
		MatrixType D;   /*!< section constitutive matrix at the current integration point */
		MatrixType BTD; /*!< auxiliary matrix to store the product B'*D */

		VectorType generalizedStrains;  /*!< generalized strain vector at the current integration point */
        VectorType generalizedStresses; /*!< generalized stress vector at the current integration point */

		VectorType N; /*!< shape function vector at the current integration point */

		MatrixType Q; /*!< 3x3 - stores the weighted sum of Q1, Q2 and Q3 */
		MatrixType Qh; /*!< 3x9 - the higher order B matrix */
		MatrixType TeQ; /*!< 3x3 - stores the product Te*Q */

		VectorType H1;
		VectorType H2;
		VectorType H3;
		VectorType H4;
		MatrixType Bb;

		ShellCrossSection::Parameters SectionParameters; /*!< parameters for cross section calculations */

		array_1d< Vector3Type, 3 > Sig;

	public:

		const ProcessInfo& CurrentProcessInfo;

	public:

		CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
			            const ProcessInfo& rCurrentProcessInfo);

	};

	///@}

    ///@name Private Operations
    ///@{

	void DecimalCorrection(Vector& a);

    void SetupOrientationAngles();

	void InitializeCalculationData(CalculationData& data);

	void CalculateBMatrix(CalculationData& data);

	void CalculateBeta0(CalculationData& data);

	void CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS);

	void ApplyCorrectionToRHS(CalculationData& data, VectorType& RHS);

	void AddBodyForces(ShellT3_LocalCoordinateSystem& LCS, VectorType& rRightHandSideVector);

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      ProcessInfo& rCurrentProcessInfo,
                      const bool LHSrequired,
                      const bool RHSrequired);

    bool TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double,3>>& rVariable,
		                                                    std::vector<array_1d<double,3>>& rValues, 
														    const ProcessInfo& rCurrentProcessInfo);

	bool TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<array_1d<double,3>>& rVariable,
		                                                             std::vector<array_1d<double,3>>& rValues, 
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
#endif // SHELL_THIN_ELEMENT_3D3N_H_INCLUDED
