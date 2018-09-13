//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LARGE_DISPLACEMENT_BEAM_SEMC_ELEMENT_H_INCLUDED)
#define  KRATOS_LARGE_DISPLACEMENT_BEAM_SEMC_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/large_displacement_beam_element.hpp"


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

/// Beam Element for 3D space dimension

/**
 * Implements a Large Displacement definition for structural analysis.
 * This works for line geometries in 3D :: it must be extended to 2D and large displacements
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) LargeDisplacementBeamSEMCElement
    :public LargeDisplacementBeamElement
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw                         ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer     ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure        StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod           IntegrationMethod;
    ///Type definition for beam utilities
    typedef BeamMathUtils<double>                     BeamMathUtilsType;
    ///Type definition for quaternion
    typedef Quaternion<double>                           QuaternionType;
    ///Type for size
    typedef GeometryData::SizeType                             SizeType;
    ///Type for element variables
    typedef LargeDisplacementBeamElement::ElementDataType ElementDataType;

    /// Counted pointer of LargeDisplacementBeamSEMCElement
    KRATOS_CLASS_POINTER_DEFINITION( LargeDisplacementBeamSEMCElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LargeDisplacementBeamSEMCElement(IndexType NewId, GeometryType::Pointer pGeometry);

    LargeDisplacementBeamSEMCElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    LargeDisplacementBeamSEMCElement(LargeDisplacementBeamSEMCElement const& rOther);

    /// Destructor.
    ~LargeDisplacementBeamSEMCElement() override;


    ///@}
    ///@name Operators
    ///@{


   /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;



    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSVector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable, Variable<array_1d<double,3> >& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo) override;

    void AddExplicitContribution(const MatrixType& rLHSMatrix, const Variable<MatrixType>& rLHSVariable, Variable<Matrix >& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo) override;


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "Large Displacement Beam SEMC Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Large Displacement Beam SEMC Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    ///@}
    ///@name Protected Operators
    ///@{
    LargeDisplacementBeamSEMCElement() {};

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculate Element Strain Resultants
     */
    virtual void CalculateStrainResultants(Vector& rStrainResultants, ElementDataType& rVariables, double alpha);

    /**
     * Calculate Element Strain Couples
     */
    virtual void CalculateStrainCouples(Vector& rStrainCouples, ElementDataType& rVariables, double alpha);

    /**
     * Calculate Element Stress Resultants and Couples
     */
    void CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber) override;

    /**
     * Calculate current curvature
     */
    void CalculateCurrentCurvature(ElementDataType& rVariables, const Variable<array_1d<double, 3 > >& rVariable) override;


    /**
      * Calculation of the Tangent Intertia Matrix
      */
    void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   ElementDataType& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight) override;


    /**
      * Calculation of the Inertial Forces Vector
      */
    void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   ElementDataType& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight) override;



    /**
     * Get Element Strain/Stress for energy computation
     */
    void CalculateStrainEnergy(double& rEnergy, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight) override;

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


}; // Class LargeDisplacementBeamSEMCElement

} // namespace Kratos.
#endif //  KRATOS_LARGE_DISPLACEMENT_BEAM_SEMC_ELEMENT_H_INCLUDED defined

