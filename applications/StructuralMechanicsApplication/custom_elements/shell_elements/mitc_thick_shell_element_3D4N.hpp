// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Massimo Petracca
//                   Alejandro Cornejo
//

#pragma once


// System includes
#include <type_traits>

// External includes

// Project includes
#include "base_shell_element.h"
#include "custom_utilities/shell_utilities.h"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"
#include "custom_utilities/shellq4_local_coordinate_system.hpp"

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

/** \brief MITCThickShellElement3D4N
 *
 * This element represents a 4-node bilinear Shell element
 * based on the Enhanced Assumed Strain Method (E.A.S.) for the membrane part
 * and on the Mixed Interpolation of Tensorial Components (M.I.T.C.)
 * for the transverse shear part.
 * This element is formulated for small strains,
 * but can be used in Geometrically nonlinear problems
 * involving large displacements and rotations
 * using a Corotational Coordinate Transformation.
 * Material nonlinearity is handled by means of the plane stress constitutive laws.
 */
template <ShellKinematics TKinematics>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) MITCThickShellElement3D4N :
    public BaseShellElement<typename std::conditional<TKinematics==ShellKinematics::NONLINEAR_COROTATIONAL,
        ShellQ4_CorotationalCoordinateTransformation,
        ShellQ4_CoordinateTransformation>::type>
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MITCThickShellElement3D4N);

    using BaseType = BaseShellElement<typename std::conditional<TKinematics==ShellKinematics::NONLINEAR_COROTATIONAL,
        ShellQ4_CorotationalCoordinateTransformation,
        ShellQ4_CoordinateTransformation>::type>;

    using QuaternionType = Quaternion<double>;
    using GeometryType = Element::GeometryType;
    using PropertiesType = Element::PropertiesType;
    using NodesArrayType = Element::NodesArrayType;
    using MatrixType = Element::MatrixType;
    using VectorType = Element::VectorType;
    using SizeType = Element::SizeType;
    using Element::GetGeometry;
    using Element::GetProperties;
    using Vector3Type = typename BaseType::Vector3Type;
    static constexpr double drilling_factor = 1.0e-1;

    ///@}

    ///@name Classes
    ///@{

    /** \brief MITC4Params
     *
     * This class performs some operations and stores some data to compute
     * the transverse shear contribution of the stiffness matrix using the
     * M.I.T.C. formulation.
     *
     * References:
     * -   Dvorkin,Bathe, "A continuum mechanics based four node shell
     *     element for general nonlinear analysis",
     *     Eng.Comput.,vol. 1, 77-88, 1984
     * -   Bathe, Dvorkin, "Short communication A four-node plate bending element
     *     based on Mindlin/Reissner plate theory and a Mixed Interpolation",
     *     International Journal for Numerical Methods in Eng.,
     *     vol. 21, 367-383, 1985
     */
    struct MITC4Params {

        double Ax;
        double Ay;
        double Bx;
        double By;
        double Cx;
        double Cy;
        Matrix Transformation;
        Matrix ShearStrains;

        MITC4Params(const ShellQ4_LocalCoordinateSystem& LCS);
    };

    class EASOperator; // forward declaration

    /** \brief EASOperatorStorage
     *
     * This class is meant to store persistent data for the EAS calculations.
     * This class is stored in the element and used by the EASOperator.
     */
    class EASOperatorStorage
    {

    public:

        friend class EASOperator;

        typedef Element::GeometryType GeometryType;

    public:

        EASOperatorStorage();

        inline void Initialize(const GeometryType& geom);

        inline void InitializeSolutionStep();

        inline void FinalizeSolutionStep();

        inline void FinalizeNonLinearIteration(const Vector& displacementVector);

    private:
        array_1d<double, 5> alpha;           // (trial) vector containing the 5 enhanced strain parameters
        array_1d<double, 5> alpha_converged; // (converged) vector containing the 5 enhanced strain parameters

        array_1d<double, 24> displ;           // (trial) vector containing the displacement vector
        array_1d<double, 24> displ_converged; // (converged) vector containing the displacement vector

        array_1d<double, 5> residual;     // vector containing the 5 residuals for the 5 enhanced strain parameters
        BoundedMatrix<double, 5, 5> Hinv; // 5x5 matrix that stores H^-1
        BoundedMatrix<double, 5, 24> L;   // 5x24 coupling matrix

        bool mInitialized; // Initialization flag

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const;

        virtual void load(Serializer& rSerializer);

    };

    /** \brief EASOperator
     *
     * This class performs some operations and stores some data to compute
     * the membrane contribution of the stiffness matrix using the
     * Enhanced Assumed Strain formulation.
     *
     * References:
     * -   J.C.Simo,M.S.Rifai, "A class of mixed assumed strain methods
     *     and the method of incompatible modes",
     *     International Journal for Numerical Methods in Eng.,
     *     vol. 29, 1595-1638, 1990
     */
    class EASOperator
    {

    public:

        /**
         * Constructor
         */
        EASOperator(const ShellQ4_LocalCoordinateSystem& LCS, EASOperatorStorage& storage);

    public:

        /**
         * this method should be called in the Gauss Loop
         * after the standard strain-displacement matrix has been computed, as well as the standard
         * generalized strains, but before the computation of the constitutive laws.
         */
        inline void GaussPointComputation_Step1(double xi, double eta, const ShellUtilities::JacobianOperator& jac,
                                                Vector& generalizedStrains,
                                                EASOperatorStorage& storage);

        /**
         * this method should be called in the Gauss Loop
         * after the standard computation of the constitutive laws.
         * note:
         * the input algorithmic tangent and generalized stress vector are assumed to be already multiplied
         * by the integration weight, and their size is assumed to be those of a standard shell element
         * (i.e. algorithmicTangent(8x8), generalizedStresses(8x1))
         */
        inline void GaussPointComputation_Step2(const Matrix& D,
                                                const Matrix& B,
                                                const Vector& S,
                                                EASOperatorStorage& storage);

        /**
         * this method should be called at the end of the Gauss Loop,
         * when the integration is terminated, but before transforming everything
         * to the global system: Here we are still operating in the element local
         * coordinate system.
         */
        inline void ComputeModfiedTangentAndResidual(Matrix& rLeftHandSideMatrix,
                Vector& rRightHandSideVector,
                EASOperatorStorage& storage);

    private:

        Matrix mF0inv;           // 3x3 inverse deformation matrix at the element center
        double mJ0;              // determinant of the jacobian at the element center
        Vector mEnhancedStrains; // vector of 3 enhanced strains [e.xx, e.yy, 2e.xy]
        Matrix mG;               // 3x5 interpolation matrix in cartesian coordinates
    };

    ///@}

    ///@name Life Cycle
    ///@{

    MITCThickShellElement3D4N(IndexType NewId,
                              GeometryType::Pointer pGeometry);

    MITCThickShellElement3D4N(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties);

    // ~MITCThickShellElement3D4N() override = default;

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


    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief this is called for non-linear analysis at the end of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    // More results calculation on integration points to interface with python

    using BaseType::CalculateOnIntegrationPoints;

    /**
     * @brief Calculate values of type double on integration points (stub)
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate values of type Vector on integration points (stub)
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * This method provides the place to perform checks on the completeness of the input
    * and the compatibility with the problem options as well as the contitutive laws selected
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    * this method is: MANDATORY
    */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;


    /**
     * @brief This method gets a value directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @tparam TType The type considered
     */
    template <class TType>
    void GetValueOnConstitutiveLaw(
        const Variable<TType> &rVariable,
        std::vector<TType> &rOutput)
    {
        const auto& r_integration_points = GetGeometry().IntegrationPoints();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            mConstitutiveLawVector[point_number]->GetValue(rVariable, rOutput[point_number]);
        }
    }

    /**
     * @brief This method computes directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @param rProcessInfo the current process info instance
     * @tparam TType The type considered
     */
    template<class TType>
    void CalculateOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const auto& r_props = GetProperties();
        const auto& r_geom = GetGeometry();
        const Matrix& shapeFunctions = r_geom.ShapeFunctionsValues();
        Vector iN(shapeFunctions.size2());

        // Compute the local coordinate system.
        ShellQ4_LocalCoordinateSystem localCoordinateSystem(
            this->mpCoordinateTransformation->CreateLocalCoordinateSystem());

        ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
            this->mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        // Prepare all the parameters needed for the MITC formulation.
        // This is to be done here outside the Gauss Loop.
        MITC4Params shearParameters(referenceCoordinateSystem);

        // Instantiate the Jacobian Operator.
        // This will store:
        // the jacobian matrix, its inverse, its determinant
        // and the derivatives of the shape functions in the local
        // coordinate system
        ShellUtilities::JacobianOperator jacOp;
        array_1d<double, 4> dArea;

        // Instantiate all strain-displacement matrices.

        Matrix B(8, 24, 0.0);
        Vector Bdrilling(24, 0.0);

        // Instantiate all section tangent matrices.
        Matrix D(8, 8, 0.0);

        // Instantiate strain and stress-resultant vectors
        Vector generalizedStrains(8);
        Vector generalizedStresses(8);

        // Get the current displacements in global coordinate system

        Vector globalDisplacements(24);
        this->GetValuesVector(globalDisplacements, 0);

        // Get the current displacements in local coordinate system
        Vector localDisplacements(this->mpCoordinateTransformation->CalculateLocalDisplacements(localCoordinateSystem, globalDisplacements));

        // Instantiate the EAS Operator.
        // This will apply the Enhanced Assumed Strain Method for the calculation
        // of the membrane contribution.
        EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

        // Initialize parameters for the cross section calculation
        ShellCrossSection::SectionParameters parameters(r_geom, r_props, rCurrentProcessInfo);
        parameters.SetGeneralizedStrainVector(generalizedStrains);
        parameters.SetGeneralizedStressVector(generalizedStresses);
        parameters.SetConstitutiveMatrix(D);
        Flags& r_options = parameters.GetOptions();
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Gauss Loop.
        for (SizeType integration_point = 0; integration_point < 4; integration_point++) {

            // get a reference of the current integration point and shape functions
            const GeometryType::IntegrationPointType& ip = r_geom.IntegrationPoints()[integration_point];
            noalias(iN) = row(shapeFunctions, integration_point);

            // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
            // and Shape functions derivatives in the local coordinate system
            jacOp.Calculate(referenceCoordinateSystem, r_geom.ShapeFunctionLocalGradient(integration_point));

            // Compute all strain-displacement matrices
            CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

            // Calculate strain vectors in local coordinate system
            noalias(generalizedStrains) = prod(B, localDisplacements);

            // Apply the EAS method to modify the membrane part of the strains computed above.
            EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains, mEASStorage);

            // Calculate the response of the Cross Section
            parameters.SetShapeFunctionsValues(iN);
            parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());

            ConstitutiveLaw::Parameters cl_values(r_geom, r_props, rCurrentProcessInfo);
            cl_values.SetStrainVector(generalizedStrains);
            cl_values.SetStressVector(generalizedStresses);
            cl_values.SetConstitutiveMatrix(D);

            mConstitutiveLawVector[integration_point]->CalculateValue(cl_values, rVariable, rOutput[integration_point]);
        }
    }

    /**
     * @brief This method calculates the shear modulus to be used.
     * In elasticity the shear modulus is G = E / (2*(1+nu)).
     * However, in some cases, such as in plasticity, the constitutive matrix may not have this form.
     * In nonlinear cases, we compute an equivalent isotropic D and retrieve its shear modulus.
    */
    double CalculateEquivalentShearModulus(const Matrix& rConstitutiveMatrix);


    /**
     * @brief This method computes the mass matrix
     */
    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

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
    MITCThickShellElement3D4N() : BaseType()
    {
    }

    ///@}

private:

    ///@name Private Operations
    ///@{

    /**
     * @brief This method calculates the B matrix for the bending, shear and membrane behavior of the element
     */
    void CalculateBMatrix(double xi, double eta,
                          const ShellUtilities::JacobianOperator& Jac, const MITC4Params& params,
                          const Vector& N,
                          Matrix& B, Vector& Bdrill);

    /**
     * @brief This method computes the element stiffness matrix and residual vector
     */
    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      const bool CalculateStiffnessMatrixFlag,
                      const bool CalculateResidualVectorFlag) override;

    /**
     * @brief This method computes the body forces contribution to the right hand side vector
     */
    void AddBodyForces(const array_1d<double,4>& dA, VectorType& rRightHandSideVector);


    /**
    * Returns the behavior of this shell (thin/thick)
    * @return the shell behavior
    */
    ShellCrossSection::SectionBehaviorType GetSectionBehavior() const override;


    /**
     * @brief This method gets a calculates the material response directly in the CL
     */
    void CalculateMaterialResponse(
        ShellCrossSection::SectionParameters& rSectionParameters,
        const SizeType& rPointNumber,
        const bool CalculateStress,
        const bool CalculateConstitutive,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief This method gets a Finalizes the material response directly in the CL
     */
    void FinalizeMaterialResponse(
        ShellCrossSection::SectionParameters& rSectionParameters,
        const SizeType& rPointNumber,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief This method gets a Initializes the material response directly in the CL
     */
    void InitializeMaterialResponse(
        ShellCrossSection::SectionParameters& rSectionParameters,
        const SizeType& rPointNumber,
        const ProcessInfo& rProcessInfo);

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{

    EASOperatorStorage mEASStorage; /// The storage instance for the EAS Operator
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

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
