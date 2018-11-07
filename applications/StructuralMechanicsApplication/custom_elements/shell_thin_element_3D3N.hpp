// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//

#if !defined(SHELL_THIN_ELEMENT_3D3N_H_INCLUDED )
#define  SHELL_THIN_ELEMENT_3D3N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_elements/base_shell_element.h"
#include "custom_utilities/shellt3_local_coordinate_system.hpp"

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ShellThinElement3D3N : public BaseShellElement
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ShellThinElement3D3N);

    typedef ShellT3_CoordinateTransformation CoordinateTransformationBaseType;

    typedef Kratos::shared_ptr<CoordinateTransformationBaseType> CoordinateTransformationBasePointerType;

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

    ~ShellThinElement3D3N() override;

    ///@}

    ///@name Operations
    ///@{

    // Basic

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    void Initialize() override;

    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

	// More results calculation on integration points to interface with python
	void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateOnIntegrationPoints(const Variable<array_1d<double,
		3> >& rVariable, std::vector<array_1d<double, 3> >& rOutput,
		const ProcessInfo& rCurrentProcessInfo) override;

    // Calculate functions
    void Calculate(const Variable<Matrix >& rVariable,
        Matrix& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

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
    ShellThinElement3D3N() : BaseShellElement()
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
        SizeType gpIndex;

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
		std::vector<VectorType> rlaminateStrains;	/*!< laminate strain vector at all surfaces at the current integration point */
		std::vector<VectorType> rlaminateStresses;	/*!< laminate stress vector at all surfaces at the current integration point */

        VectorType N; /*!< shape function vector at the current integration point */

        MatrixType Q; /*!< 3x3 - stores the weighted sum of Q1, Q2 and Q3 */
        MatrixType Qh; /*!< 3x9 - the higher order B matrix */
        MatrixType TeQ; /*!< 3x3 - stores the product Te*Q */

        VectorType H1;
        VectorType H2;
        VectorType H3;
        VectorType H4;
        MatrixType Bb;

        ShellCrossSection::SectionParameters SectionParameters; /*!< parameters for cross section calculations */

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

	void CheckGeneralizedStressOrStrainOutput(const Variable<Matrix>& rVariable, int& ijob, bool& bGlobal);

	void CalculateStressesFromForceResultants(VectorType& rstresses,
		const double& rthickness);

	void CalculateLaminaStrains(CalculationData& data);

	void CalculateLaminaStresses(CalculationData& data);

	double CalculateTsaiWuPlaneStress(const CalculationData& data, const Matrix& rLamina_Strengths, const unsigned int& rCurrent_Ply);

	void CalculateVonMisesStress(const CalculationData& data, const Variable<double>& rVariable, double& rVon_Mises_Result);

    void DecimalCorrection(Vector& a);

    void SetupOrientationAngles() override;

    void InitializeCalculationData(CalculationData& data);

    void CalculateBMatrix(CalculationData& data);

    void CalculateBeta0(CalculationData& data);

    void CalculateSectionResponse(CalculationData& data);

    void CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS);

    void ApplyCorrectionToRHS(CalculationData& data, VectorType& RHS);

    void AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector);

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) override;

    bool TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo);

    /**
    * Returns the behavior of this shell (thin/thick)
    * @return the shell behavior
    */
    ShellCrossSection::SectionBehaviorType GetSectionBehavior() override;

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{

    CoordinateTransformationBasePointerType mpCoordinateTransformation; /*!< The Coordinate Transformation */

    SizeType mStrainSize = 6;

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
#endif // SHELL_THIN_ELEMENT_3D3N_H_INCLUDED
