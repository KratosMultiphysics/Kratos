//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

#if !defined(KRATOS_SMALL_DISPLACEMENT_SIMP_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_SIMP_ELEMENT_H_INCLUDED

// Project includes

#include "structural_mechanics_application.h"
#include "custom_elements/small_displacement.h"


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

/// Topology Optimization Small Displacement Element for 3D geometries.

/**
 * Implements a Topology Optimization Small Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D
 */
class KRATOS_API(TOPOLOGY_OPTIMIZATION_APPLICATION) SmallDisplacementSIMPElement
    : public SmallDisplacement
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    ///Type for element variables
    ///typedef SmallDisplacement::ElementDataType ElementDataType;


    /// hinugefügt am 10.02.
    /// The base element type
    typedef SmallDisplacement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of SmallDisplacementSIMPElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SmallDisplacementSIMPElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    SmallDisplacementSIMPElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SmallDisplacementSIMPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    SmallDisplacementSIMPElement(SmallDisplacementSIMPElement const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~SmallDisplacementSIMPElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SmallDisplacementSIMPElement& operator=(SmallDisplacementSIMPElement const& rOther);

    ///@}
    ///@name Operations
    ///@{
        /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */  
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
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
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /// Function to calculate the sensitivities and the objective function
        void Calculate(
    const Variable<double>& rVariable,
    double &rOutput,
    const ProcessInfo& rCurrentProcessInfo
    ) override;  


    // =============================================================================================================================================
    // =============================================================================================================================================


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
        buffer << "Small Displacement Solid Element #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Small Displacement Solid Element #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
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

    ///@}
    ///@name Protected Operators
    ///@{

    SmallDisplacementSIMPElement() : SmallDisplacement()
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
    virtual void save(Serializer& rSerializer) const override;
    virtual void load(Serializer& rSerializer) override;


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
#endif // KRATOS_SMALL_DISPLACEMENT_TOPLOGY_OPTIMIZATION_ELEMENT_H_INCLUDED  defined
