//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/large_displacement_element.hpp"


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

/// Updated Lagrangian Element for 3D and 2D geometries.

/**
 * Implements a spatial Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) UpdatedLagrangianElement
    : public LargeDisplacementElement
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

    /// Counted pointer of UpdatedLagrangianElement
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianElement );
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    UpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianElement(UpdatedLagrangianElement const& rOther);

    /// Destructor.
    virtual ~UpdatedLagrangianElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianElement& operator=(UpdatedLagrangianElement const& rOther);

    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    /**
     * creates a  element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

  
    // //************* GETTING METHODS

    //SET

    /**
     * Set a double  Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);


    //GET:

    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);



    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize();

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //int Check(const ProcessInfo& rCurrentProcessInfo);


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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Updated Lagrangian Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Updated Lagrangian Element #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      GetGeometry().PrintData(rOStream);
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

    /**
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    std::vector< Matrix > mDeformationGradientF0;

    /**
     * Container for the total deformation gradient determinants
     */
    Vector mDeterminantF0;

    ///@}
    ///@name Protected Operators
    ///@{
    UpdatedLagrangianElement() : LargeDisplacementElement()
    {
    }

    /**
     * Initialize Element General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables& rVariables, 
					    const ProcessInfo& rCurrentProcessInfo);



    /**
     * Finalize Element Internal Variables
     */
    virtual void FinalizeStepVariables(GeneralVariables & rVariables, 
				       const double& rPointNumber );


    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables,
                                     const double& rPointNumber);

    /**
     * Calculation of the Deformation Gradient F
     */
    void CalculateDeformationGradient(const Matrix& rDN_DX,
                                      Matrix& rF,
                                      Matrix& rDeltaPosition);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
				    Matrix& rF,
				    Matrix& rDN_DX);

    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    void GetHistoricalVariables( GeneralVariables& rVariables, 
				 const double& rPointNumber );

    /**
     * Calculate on Integration Points (Calculation of internal energy)
     */
    void CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo );
    
    /**
     * Calculation of the Volume Change of the Element
     */
    virtual double& CalculateVolumeChange(double& rVolumeChange, GeneralVariables& rVariables);

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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class UpdatedLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED  defined 
