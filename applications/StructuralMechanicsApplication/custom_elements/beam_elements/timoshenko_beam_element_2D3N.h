// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "timoshenko_beam_element_2D2N.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

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

/**
 * @class LinearTimoshenkoBeamElement2D3N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the Timoshenko beam element of 3 nodes. 5th order degree Hermitic polynomials
 * for deflection and 2nd degree for axial displacements. The ordering of the local shape functions
 * assume a " 0--1--2 " ordering and then we swap the components when moving to global coordinates
 *
 *                                    ^ y, v
 *                                    |
 * Global Ordering of the nodes:      0 ------ 2 ------- 1 --> x, u      and rotation theta node-wise
 * Reference: Felippa and Oñate,
 * "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability",
 * Archives of Comp. Methods in Eng. (2021) 28:2021-2080
 * DOI: https://doi.org/10.1007/s11831-020-09515-0; Then adapted to the 3 noded straight beam element.
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearTimoshenkoBeamElement2D3N
    : public LinearTimoshenkoBeamElement2D2N
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = LinearTimoshenkoBeamElement2D2N;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LinearTimoshenkoBeamElement2D3N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    LinearTimoshenkoBeamElement2D3N() = default;

    // Constructor using an array of nodes
    LinearTimoshenkoBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
    }

    // Constructor using an array of nodes with properties
    LinearTimoshenkoBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
    }

    // Copy constructor
    LinearTimoshenkoBeamElement2D3N(LinearTimoshenkoBeamElement2D2N const& rOther)
        : BaseType(rOther)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LinearTimoshenkoBeamElement2D3N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<LinearTimoshenkoBeamElement2D3N>(NewId, pGeom, pProperties);
    }

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns a 9 component vector including the values of the DoFs
     * in LOCAL axes
     */
    void GetNodalValuesVector(VectorType& rNodalValue) const override;

    /**
     * @brief Computes the axial strain (El), shear strain (gamma_xy) and bending curvature (kappa)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
     * @param rNodalValues The vector containing the nodal values in local axes
     */
    double CalculateAxialStrain     (const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const override;
    double CalculateShearStrain     (const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const override;
    double CalculateBendingCurvature(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const override;

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector including nul values to the axial u terms
     * @param rLocalSizeVector The 6 local components of v and theta
     */
    void GlobalSizeVector(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) const override
    {
        // The ordering of the local vector nodes is different from global
        rGlobalSizeVector.clear();
        rGlobalSizeVector[1] = rLocalSizeVector[0];
        rGlobalSizeVector[2] = rLocalSizeVector[1];

        rGlobalSizeVector[7] = rLocalSizeVector[2];
        rGlobalSizeVector[8] = rLocalSizeVector[3];
        rGlobalSizeVector[4] = rLocalSizeVector[4];
        rGlobalSizeVector[5] = rLocalSizeVector[5];
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector including only the axial u terms
     * @param rLocalSizeVector The 3 local components of u
     */
    void GlobalSizeAxialVector(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) override
    {
        // The ordering of the local vector nodes is different from global
        rGlobalSizeVector.clear();
        rGlobalSizeVector[0] = rLocalSizeVector[0];
        rGlobalSizeVector[6] = rLocalSizeVector[1];
        rGlobalSizeVector[3] = rLocalSizeVector[2];
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief This function returns the 4 shape functions used for interpolating the transverse displacement v. (denoted as N)
     * Also its derivatives
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFirstDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetSecondDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetThirdDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFourthDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;

    /**
     * @brief This function returns the 4 shape functions used for interpolating the total rotation Theta (N_theta)
     * Also its derivative
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFirstDerivativesNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;

    /**
     * @brief This function returns the 2 shape functions used for interpolating the axial displacement u0
     * Also its derivatives
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFirstDerivativesNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;

    /**
     * @brief This function rotates the LHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateLHS(
        MatrixType &rLHS,
        const GeometryType &rGeometry) override;

    /**
     * @brief This function rotates the RHS from local to global coordinates
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateRHS(
        VectorType &rRHS,
        const GeometryType &rGeometry) override;

    /**
     * @brief This function rotates the LHS and RHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateAll(
        MatrixType &rLHS,
        VectorType &rRHS,
        const GeometryType &rGeometry) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Linear Timoshenko 3N straight Beam Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;

}; // class LinearTimoshenkoBeamElement2D2N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
