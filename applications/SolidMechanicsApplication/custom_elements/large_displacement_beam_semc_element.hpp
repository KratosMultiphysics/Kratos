//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
// 

#if !defined(KRATOS_LARGE_DISPLACEMENT_BEAM_SEMC_ELEMENT_H_INCLUDED )
#define  KRATOS_LARGE_DISPLACEMENT_BEAM_SEMC_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/large_displacement_beam_element.hpp"


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
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, DELTA_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class LargeDisplacementBeamSEMCElement
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
    virtual ~LargeDisplacementBeamSEMCElement();


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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
 

    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize();
  
      /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
  
    //************************************************************************************
    //************************************************************************************

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

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
    void AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable, Variable<array_1d<double,3> >& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo);

    void AddExplicitContribution(const MatrixType& rLHSMatrix, const Variable<MatrixType>& rLHSVariable, Variable<Matrix >& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo);


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);
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
        buffer << "Large Displacement Beam SEMC Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Large Displacement Beam SEMC Element #" << Id();
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
     * Elemental previous curvature vectors for each integration point
     */
    std::vector<Vector>  mPreviousCurvatureVectors;

    ///@}
    ///@name Protected Operators
    ///@{
    LargeDisplacementBeamSEMCElement() {};

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Update strain current member variables
     */ 
    virtual void UpdateStrainVariables(GeneralVariables& rVariables, 
				       const unsigned int& rPointNumber);

    /**   
     * Calculate Element Strain Resultants
     */ 
    virtual void CalculateStrainResultants(Vector& rStrainResultants, GeneralVariables& rVariables, double alpha);

    /**   
     * Calculate Element Strain Couples
     */ 
    virtual void CalculateStrainCouples(Vector& rStrainCouples, GeneralVariables& rVariables, double alpha);

    /**   
     * Calculate Element Stress Resultants and Couples
     */ 
    virtual void CalculateStressResultants(GeneralVariables& rVariables, const unsigned int& rPointNumber, double alpha);

    /**   
     * Calculate current curvature
     */
    virtual void CalculateCurrentCurvature(GeneralVariables& rVariables, const Variable<array_1d<double, 3 > >& rVariable);


    /**
      * Calculation of the Tangent Intertia Matrix
      */
    virtual void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   GeneralVariables& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight);


    /**
      * Calculation of the Inertial Forces Vector
      */
    virtual void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   GeneralVariables& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight);



    /**
     * Get Element Strain/Stress for energy computation
     */
    virtual void CalculateStrainEnergy(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight);

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


    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class LargeDisplacementBeamSEMCElement

} // namespace Kratos.
#endif //  KRATOS_LARGE_DISPLACEMENT_BEAM_SEMC_ELEMENT_H_INCLUDED defined

