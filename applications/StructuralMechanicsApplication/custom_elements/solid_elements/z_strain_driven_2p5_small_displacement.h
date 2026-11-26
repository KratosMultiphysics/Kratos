// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//                   Ignasi de Pouplana
//


#pragma once

// Project includes
#include "includes/checks.h"

// Application includes
#include "small_displacement.h"
#include "structural_mechanics_application_variables.h"

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
 * @class ZStrainDriven2p5DSmallDisplacement
 * @ingroup StructuralMechanicsApplication
 * @brief Small displacement element for 2.5D cases.
 * @details Implements a small displacement definition for structural analysis.
 * This works for arbitrary geometries in 2D domains, but takes the Z Strain as input for the calculation of the stresses.
 * @author Miguel Angel Celigueta
 * @author Ignasi de Pouplana
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ZStrainDriven2p5DSmallDisplacement
    : public SmallDisplacement
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef SmallDisplacement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of ZStrainDriven2p5DSmallDisplacement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ZStrainDriven2p5DSmallDisplacement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ZStrainDriven2p5DSmallDisplacement(IndexType NewId, GeometryType::Pointer pGeometry);
    ZStrainDriven2p5DSmallDisplacement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    ZStrainDriven2p5DSmallDisplacement(ZStrainDriven2p5DSmallDisplacement const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~ZStrainDriven2p5DSmallDisplacement() override = default;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

     /**
      * @brief Set a double Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValuesOnIntegrationPoints(
        const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Small Displacement Solid Element #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()
            && BaseType::mConstitutiveLawVector[0] != nullptr) {
          buffer << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          buffer << " (no constitutive law)";
        }
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Small Displacement Solid Element #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()
            && BaseType::mConstitutiveLawVector[0] != nullptr) {
          rOStream << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          rOStream << " (no constitutive law)";
        }
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<double> mImposedZStrainVector; /// The vector containing the imposed z strains

    ///@}
    ///@name Protected Operators
    ///@{

    ZStrainDriven2p5DSmallDisplacement() : SmallDisplacement()
    {
    }

    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        ) override;

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
    ///@}

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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}

}; // Class ZStrainDriven2p5DSmallDisplacement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
