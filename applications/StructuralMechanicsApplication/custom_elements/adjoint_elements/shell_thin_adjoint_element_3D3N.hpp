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

#if !defined(SHELL_THIN_ADJOINT_ELEMENT_3D3N_H_INCLUDED )
#define  SHELL_THIN_ADJOINT_ELEMENT_3D3N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "custom_utilities/shellt3_local_coordinate_system.hpp"
#include "utilities/quaternion.h"
#include "custom_elements/shell_thin_element_3D3N.hpp"

namespace Kratos
{


///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

//class ShellT3_CoordinateTransformation;

///@name  Enum's
///@{
///@}

///@name  Functions
///@{
///@}

///@name Kratos Classes
///@{

/** \brief ShellThinAdjointElement3D3N
 *
 * This element is inherited from ShellThinElement3D3N. It is the corresponding 
 * element to it and is used for solving the adjoint problem and for computing sensitivities.
 */
class ShellThinAdjointElement3D3N : public ShellThinElement3D3N 
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ShellThinAdjointElement3D3N);

     typedef Element::PropertiesType PropertiesType;

     typedef Element::DofsArrayType DofsArrayType;

   // typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

   // typedef ShellT3_CoordinateTransformation CoordinateTransformationBaseType;

    //typedef boost::shared_ptr<CoordinateTransformationBaseType> CoordinateTransformationBasePointerType;

   // typedef array_1d<double, 3> Vector3Type;

   // typedef Quaternion<double> QuaternionType;

    ///@}

    ///@name Classes
    ///@{

    // TODO: Add Calulation Data

    ///@}

    ///@name Life Cycle
    ///@{

    ShellThinAdjointElement3D3N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         bool NLGeom = false);

    ShellThinAdjointElement3D3N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties,
                         bool NLGeom = false);

    ShellThinAdjointElement3D3N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties,
                         CoordinateTransformationBasePointerType pCoordinateTransformation);

    ~ShellThinAdjointElement3D3N() override;

    ///@}

    ///@name Operations
    ///@{

    // Basic
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    //IntegrationMethod GetIntegrationMethod() const override;

    void InitializeAgainAfterDisturbing();

    //void ResetConstitutiveLaw() override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    //void CleanMemory() override;

    void GetValuesVector(Vector& values, int Step = 0) override;

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, 
											const ProcessInfo& rCurrentProcessInfo) override;
	
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput, 
											const ProcessInfo& rCurrentProcessInfo) override;

	void Calculate(const Variable<Matrix >& rVariable, Matrix& Output,
                           const ProcessInfo& rCurrentProcessInfo) override;    

    void Calculate(const Variable<Vector >& rVariable, Vector& Output,
                           const ProcessInfo& rCurrentProcessInfo) override;                                                                

    //void GetFirstDerivativesVector(Vector& values, int Step = 0) override;

    //void GetSecondDerivativesVector(Vector& values, int Step = 0) override;

    //void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    //void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    //void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    //void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    //void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    //void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;

    //void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                            //  VectorType& rRightHandSideVector,
                            //  ProcessInfo& rCurrentProcessInfo) override;

    //void CalculateRightHandSide(VectorType& rRightHandSideVector,
                               // ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;    
                      
    // Results calculation on integration points

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput,
					                    const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //void GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable, std::vector<array_1d<double,6> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

    ///@name Public specialized Access - Temporary
    ///@{

    //void SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections);

    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{

    /**
     * Protected empty constructor
     */
    ShellThinAdjointElement3D3N() : ShellThinElement3D3N()
    {
    }

    ///@}

private:

    ///@name Private Classes
    ///@{

    //class CalculationData
    //{

    //public:

        // ---------------------------------------
        // calculation-constant data
        // ----------------------------------------
        // these data are allocated and constructed
        // at the beginning of the calculation

       // ShellT3_LocalCoordinateSystem LCS0; /*!< reference coordinate system */
       // ShellT3_LocalCoordinateSystem LCS;  /*!< current coordinate system */

       /* MatrixType L;

        MatrixType Q1;
        MatrixType Q2;
        MatrixType Q3;

        MatrixType Te;
        MatrixType TTu;

        double dA;
        double hMean;
        double TotalArea;
        double TotalVolume;
        std::vector< array_1d<double,3> > gpLocations;*/

        //MatrixType dNxy; /*!< shape function cartesian derivatives */

        //VectorType globalDisplacements; /*!< global displacement vector */
        //VectorType localDisplacements;  /*!< local displacement vector */

        //bool CalculateRHS; /*!< flag for the calculation of the right-hand-side vector */
        //bool CalculateLHS; /*!< flag for the calculation of the left-hand-side vector */

        // ---------------------------------------
        // calculation-variable data
        // ---------------------------------------
        // these data are updated during the
        // calculations

        //double beta0;
        //size_t gpIndex;

        // ---------------------------------------
        // calculation-variable data
        // ---------------------------------------
        // these data are updated during the
        // calculations, but they are allocated
        // only once(the first time they are used)
        // to avoid useless re-allocations

        //MatrixType B;   /*!< total strain-displacement matrix at the current integration point */
        //MatrixType D;   /*!< section constitutive matrix at the current integration point */
        //MatrixType BTD; /*!< auxiliary matrix to store the product B'*D */

        //VectorType generalizedStrains;  /*!< generalized strain vector at the current integration point */
        //VectorType generalizedStresses; /*!< generalized stress vector at the current integration point */

        //VectorType N; /*!< shape function vector at the current integration point */

        //MatrixType Q; /*!< 3x3 - stores the weighted sum of Q1, Q2 and Q3 */
        //MatrixType Qh; /*!< 3x9 - the higher order B matrix */
        //MatrixType TeQ; /*!< 3x3 - stores the product Te*Q */

        //VectorType H1;
        //VectorType H2;
        //VectorType H3;
        //VectorType H4;
        //MatrixType Bb;

        //ShellCrossSection::Parameters SectionParameters; /*!< parameters for cross section calculations */

        //array_1d< Vector3Type, 3 > Sig;

    //public:

       // const ProcessInfo& CurrentProcessInfo;

    //public:

      //  CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
        //                const ProcessInfo& rCurrentProcessInfo);

    //};

    ///@}

    ///@name Private Operations
    ///@{

   /*void DecimalCorrection(Vector& a);

    void SetupOrientationAngles();

    void InitializeCalculationData(CalculationData& data);

    void CalculateBMatrix(CalculationData& data);

    void CalculateBeta0(CalculationData& data);

    void CalculateSectionResponse(CalculationData& data);

    void CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS);

    void ApplyCorrectionToRHS(CalculationData& data, VectorType& RHS);

    void AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector);*/

    /*void CalculateAllPrimal(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      ProcessInfo& rCurrentProcessInfo,
                      const bool LHSrequired,
                      const bool RHSrequired);*/

    /*bool TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double,3> >& rVariable,
            std::vector<array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo);

    bool TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo);*/

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{

    //CoordinateTransformationBasePointerType mpCoordinateTransformation; /*!< The Coordinate Transformation */

    //CrossSectionContainerType mSections; /*!< Container for cross section associated to each integration point */

    //IntegrationMethod mThisIntegrationMethod; /*!< Currently selected integration method */

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
#endif // SHELL_THIN_ADJOINT_ELEMENT_3D3N_H_INCLUDED
