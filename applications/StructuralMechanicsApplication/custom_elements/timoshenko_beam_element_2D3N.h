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
#include "custom_elements/timoshenko_beam_element_2D2N.h"
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
 * @class TimoshenkoBeamElement2D3N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the Timoshenko beam element of 3 nodes. 5th order degree Hermitic polynomials
 * for deflection and 3rd degree for axial displacements
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TimoshenkoBeamElement2D3N
    : public TimoshenkoBeamElement2D2N
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = TimoshenkoBeamElement2D2N;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TimoshenkoBeamElement2D3N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    TimoshenkoBeamElement2D3N()
    {
    }

    // Constructor using an array of nodes
    TimoshenkoBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
    }

    // Constructor using an array of nodes with properties
    TimoshenkoBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
    }

    // Copy constructor
    TimoshenkoBeamElement2D3N(TimoshenkoBeamElement2D2N const& rOther)
        : BaseType(rOther)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TimoshenkoBeamElement2D3N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<TimoshenkoBeamElement2D3N>(NewId, pGeom, pProperties);
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
    void GetNodalValuesVector(VectorType& rNodalValue) override;

    /**
     * @brief This method returns the length of the FE
     */
    double GetAngle() override
    {
        return StructuralMechanicsElementUtilities::GetReferenceRotationAngle2D3NBeam(GetGeometry());
    }

    /**
     * @brief Computes the axial strain (El), shear strain (gamma_xy) and bending curvature (kappa)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
     * @param rNodalValues The vector containing the nodal values in local axes
     */
    double CalculateAxialStrain     (const double Length, const double Phi, const double xi, const VectorType& rNodalValues) override;
    double CalculateShearStrain     (const double Length, const double Phi, const double xi, const VectorType& rNodalValues) override;
    double CalculateBendingCurvature(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) override;

    /**
     * @brief Computes the length of the FE and returns it
     */
    double CalculateLength() override
    {
        return StructuralMechanicsElementUtilities::CalculateReferenceLength2D3N(*this);
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector including nul values to the axial u terms
     * @param rLocalSizeVector The 6 local components of v and theta
     */
    void GlobalSizeVector(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) override
    {
        rGlobalSizeVector.clear();
        rGlobalSizeVector[1] = rLocalSizeVector[0];
        rGlobalSizeVector[2] = rLocalSizeVector[1];
        rGlobalSizeVector[4] = rLocalSizeVector[2];
        rGlobalSizeVector[5] = rLocalSizeVector[3];
        rGlobalSizeVector[7] = rLocalSizeVector[4];
        rGlobalSizeVector[8] = rLocalSizeVector[5];
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector including only the axial u terms
     * @param rLocalSizeVector The 3 local components of u
     */
    void GlobalSizeAxialVector(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) override
    {
        rGlobalSizeVector.clear();
        rGlobalSizeVector[0] = rLocalSizeVector[0];
        rGlobalSizeVector[3] = rLocalSizeVector[1];
        rGlobalSizeVector[6] = rLocalSizeVector[2];
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
    void GetShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;
    void GetFirstDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;
    void GetSecondDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;
    void GetThirdDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;
    void GetFourthDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;

    /**
     * @brief This function returns the 4 shape functions used for interpolating the total rotation Theta (N_theta)
     * Also its derivative
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;
    void GetFirstDerivativesNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;

    /**
     * @brief This function returns the 2 shape functions used for interpolating the axial displacement u0
     * Also its derivatives
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;
    void GetFirstDerivativesNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) override;

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
        rOStream << "Timoshenko 3N Beam Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
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

}; // class TimoshenkoBeamElement2D2N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
