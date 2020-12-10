// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_SMALL_DISPLACEMENT_SIMP_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_SIMP_ELEMENT_H_INCLUDED

// Project includes
//#include "solid_mechanics_application.h"
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "includes/variables.h"
//#include "small_displacement_element.hpp"



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
 * @class SmallDisplacement
 * @ingroup StructuralMechanicsApplication
 * @brief Small displacement element for 2D and 3D geometries.
 * @details Implements a small displacement definition for structural analysis. This works for arbitrary geometries in 2D and 3D
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementSIMPElement
    : public BaseSolidElement
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


    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of SmallDisplacement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementSIMPElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry ):BaseSolidElement(NewId,pGeometry)
    {};
    SmallDisplacementSIMPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):BaseSolidElement(NewId,pGeometry,pProperties)
    {};

    // Copy constructor
    SmallDisplacementSIMPElement(SmallDisplacementSIMPElement const& rOther)
    {};

    /// Destructor.
    ~SmallDisplacementSIMPElement() override
    {};

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

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
 // =============================================================================================================================================
    // STARTING / ENDING METHODS
    // =============================================================================================================================================

    /// Function that gets the value on the Integration Point (For printing purposes in the output GiD)
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /// Function to calculate the sensitivities and the objective function
    void Calculate(const Variable<double> &rVariable, double &rOutput, const ProcessInfo &rCurrentProcessInfo);

    /// Function that overwrites the CalculateOnIntegrationPoints, to insert the X_PHYS into all Gauss Points of the given element
    /// That allows printing X_PHYS as elemental value in GiD
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
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

    SmallDisplacementSIMPElement() : BaseSolidElement()
    {
    }

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
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class SmallDisplacementSIMPElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_H_INCLUDED  defined

