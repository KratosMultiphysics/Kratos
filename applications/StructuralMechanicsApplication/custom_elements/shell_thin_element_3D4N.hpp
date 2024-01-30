// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Peter Wilson
//       Contact:    A.Winterstein[at]tum.de
//

#pragma once

// System includes
#include <type_traits>

// External includes

// Project includes
#include "custom_elements/base_shell_element.h"
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
/** \brief ShellThinElement3D4N
*
* This element represents a 4-node Shell element.
* The membrane part is Felippa's assumed Natural DEviatoric Strain (ANDES)
* formulation, while the bending part is the Discrete Kirchhoff
* Quadrilateral.
* This element is formulated for small strains,
* but can be used in Geometrically nonlinear problems
* involving large displacements and rotations
* using a Corotational Coordinate Transformation.
* Material nonlinearity is handled by means of the cross section object.
*/

/*
Shell formulation references:

ANDES formulation:
Bjorn Haugen. "Buckling and Stability Problems for Thin Shell Structures
Using High Performance Finite Elements". Dissertation. Colorado: University
of Colorado, 1994.

ANDES filter matrix H modification as per:
Carlos A Felippa. "Supernatural QUAD4: a template formulation".    In: Computer
methods in applied mechanics and engineering 195.41 (2006), pp. 5316-5342.

DKQ original formulation:
Jean-Louis Batoz and Mabrouk Ben Tahar. "Evaluation of a new quadrilateral
thin plate bending element". In: International Journal for Numerical Methods
in Engineering 18.11 (1982), pp. 1655-1677.

Clearly presented DKQ formulation:
Fabian Rojas Barrales. "Development of a nonlinear quadrilateral layered
membrane element with drilling degrees of freedom and a nonlinear
quadrilateral thin flat layered shell element for the modeling of reinforced
concrete walls". Dissertation. Los Angeles, California: University of
Southern California, 2012. */

template <ShellKinematics TKinematics>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ShellThinElement3D4N : public
    BaseShellElement<typename std::conditional<TKinematics==ShellKinematics::NONLINEAR_COROTATIONAL,
        ShellQ4_CorotationalCoordinateTransformation,
        ShellQ4_CoordinateTransformation>::type>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ShellThinElement3D4N);

    using BaseType = BaseShellElement<typename std::conditional<TKinematics==ShellKinematics::NONLINEAR_COROTATIONAL,
        ShellQ4_CorotationalCoordinateTransformation,
        ShellQ4_CoordinateTransformation>::type>;

    typedef Quaternion<double> QuaternionType;

    using GeometryType = Element::GeometryType;

    using PropertiesType = Element::PropertiesType;

    using NodesArrayType = Element::NodesArrayType;

    using MatrixType = Element::MatrixType;

    using VectorType = Element::VectorType;

    using SizeType = Element::SizeType;

    using Element::GetGeometry;

    using Element::GetProperties;

    using Vector3Type = typename BaseType::Vector3Type;

    ///@}

    ///@name Life Cycle
    ///@{
    ShellThinElement3D4N(IndexType NewId,
                         GeometryType::Pointer pGeometry);

    ShellThinElement3D4N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties);

    ~ShellThinElement3D4N() override = default;

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


    // More results calculation on integration points to interface with python

    using BaseType::CalculateOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

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
    ShellThinElement3D4N() : BaseType()
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

        ShellQ4_LocalCoordinateSystem LCS;  /*!< current coordinate system */
        ShellQ4_LocalCoordinateSystem LCS0; /*!< reference coordinate system */

        // Unit vectors (in cartesian coords)
        Vector s_xi = ZeroVector(3);    /*!< xi unit vector in cartesian coords */
        Vector s_eta = ZeroVector(3);    /*!< xi unit vector in cartesian coords */

        // Geometry data
        array_1d<Vector, 4> r_cartesian;    /*!< array of cartesian point positions */
        array_1d<double, 4> dA;    /*!< array of integration areas (incl. weight) */

        // Displacements
        VectorType globalDisplacements = ZeroVector(24); /*!< global displacement vector */
        VectorType localDisplacements = ZeroVector(24);  /*!< local displacement vector */

        // Element flags
        bool CalculateRHS; /*!< flag for the calculation of the right-hand-side vector */
        bool CalculateLHS; /*!< flag for the calculation of the left-hand-side vector */

        // ---------------------------------------
        // Testing flags
        // ---------------------------------------
        // These should both be FALSE unless you are testing, or
        // investigating the effects of element enhancements!

        const bool basicQuad = false;    /*!< flag for using basic membrane
                                        formulation - should be FALSE unless
                                        you are testing */

        // ---------------------------------------
        // calculation-variable data
        // ---------------------------------------
        // these data are updated during the
        // calculations

        array_1d<double, 4> N;    /*!< SF values at parametric point */
        SizeType gpIndex;    /*!< Index of current Gauss Point (zero based) */

        // ---------------------------------------
        // calculation-variable data
        // ---------------------------------------
        // these data are updated during the
        // calculations, but they are allocated
        // only once(the first time they are used)
        // to avoid useless re-allocations

        // ANDES membrane data
        const double alpha = 1.5;
        MatrixType L_mem = ZeroMatrix(3, 12); /*!< basic membrane lumping matrix */
        MatrixType H_mem_mod = ZeroMatrix(7, 12);    /*!< higher order membrane filter matrix */
        MatrixType Z = ZeroMatrix(12, 12);    /*!< transformation from Haugen to Kratos DOFs */
        MatrixType B_h_1 = ZeroMatrix(3, 7);    /*!< higher order membrane B1 matrix */
        MatrixType B_h_2 = ZeroMatrix(3, 7);    /*!< higher order membrane B2 matrix */
        MatrixType B_h_3 = ZeroMatrix(3, 7);    /*!< higher order membrane B3 matrix */
        MatrixType B_h_4 = ZeroMatrix(3, 7);    /*!< higher order membrane B4 matrix */
        MatrixType B_h_bar = ZeroMatrix(3, 7);    /*!< higher order membrane B_bar matrix */
        MatrixType T_13 = ZeroMatrix(3, 3);
        MatrixType T_24 = ZeroMatrix(3, 3);

        // DKQ bending data
        array_1d<double, 4> DKQ_a;
        array_1d<double, 4> DKQ_b;
        array_1d<double, 4> DKQ_c;
        array_1d<double, 4> DKQ_d;
        array_1d<double, 4> DKQ_e;
        MatrixType DKQ_indices = ZeroMatrix(4, 2);
        array_1d<Matrix, 4> DKQ_invJac;
        array_1d<Matrix, 4> DKQ_jac;
        array_1d<double, 4> DKQ_jac_det;


        // General total element data
        MatrixType B = ZeroMatrix(6, 24);   /*!< total strain-displacement matrix at the current integration point */
        MatrixType D = ZeroMatrix(6, 6);    /*!< section constitutive matrix at the current integration point */
        MatrixType BTD = ZeroMatrix(24, 6);  /*!< auxiliary matrix to store the product B'*D */

        VectorType generalizedStrains = ZeroVector(6);  /*!< generalized strain vector at the current integration point */
        VectorType generalizedStresses = ZeroVector(6); /*!< generalized stress vector at the current integration point */
        std::vector<VectorType> rlaminateStrains;    /*!< laminate strain vector at all surfaces at the current integration point */
        std::vector<VectorType> rlaminateStresses;    /*!< laminate stress vector at all surfaces at the current integration point */

        ShellUtilities::JacobianOperator jacOp;
        ShellCrossSection::SectionParameters SectionParameters; /*!< parameters for cross section calculations */

    public:

        const ProcessInfo& CurrentProcessInfo;

    public:

        CalculationData(const ShellQ4_LocalCoordinateSystem& localcoordsys,
                        const ShellQ4_LocalCoordinateSystem& refcoordsys,
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

    void InitializeCalculationData(CalculationData& data);

    void CalculateBMatrix(CalculationData& data);

    void CalculateSectionResponse(CalculationData& data);

    void CalculateGaussPointContribution(CalculationData& data,
                                         MatrixType& LHS, VectorType& RHS);

    void AddBodyForces(CalculationData& data,
                       VectorType& rRightHandSideVector); //not required for dyn

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      const bool CalculateStiffnessMatrixFlag,
                      const bool CalculateResidualVectorFlag) override;

    bool TryCalculateOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo);

    /**
    * Returns the behavior of this shell (thin/thick)
    * @return the shell behavior
    */
    ShellCrossSection::SectionBehaviorType GetSectionBehavior() const override;

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
