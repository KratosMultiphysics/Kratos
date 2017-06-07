// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_KINEMATIC_LINEAR_H_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

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

/// Total Lagrangian element for 2D and 3D geometries.

/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 */

class KinematicLinear
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

    /// Counted pointer of KinematicLinear
    KRATOS_CLASS_POINTER_DEFINITION(KinematicLinear);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KinematicLinear(IndexType NewId, GeometryType::Pointer pGeometry);
    KinematicLinear(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~KinematicLinear();

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    //TODO: ADD THE OTHER CREATE FUNCTION
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

    /**
      * This resets the constitutive law
      */
    void ResetConstitutiveLaw();

     /**
     * This function provides a more general interface to the element. 
     * It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrices: container with the output left hand side matrices
     * @param rLHSVariables: paramter describing the expected LHSs
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: the current process info instance
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Called at the end of eahc solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    /**
     * Calculate a double Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable, 
        std::vector<double>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Matrix >& rVariable, 
        std::vector< Matrix >& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
     * Set a double Value on the Element Constitutive Law
     * @param rVariable: The variable we want to set
     * @param rValues: The values to set in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void SetValueOnIntegrationPoints(
        const Variable<double>& rVariable, 
        std::vector<double>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
     /**
     * Set a Vector Value on the Element Constitutive Law
     * @param rVariable: The variable we want to set
     * @param rValues: The values to set in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
     * Set a Matrix Value on the Element Constitutive Law
     * @param rVariable: The variable we want to set
     * @param rValues: The values to set in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable, 
        std::vector<Matrix>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable, 
        std::vector<double>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable, 
        std::vector<Matrix>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This method is not defined yet...
     * @param rCurrentProcessInfo: the current process info instance
     */
    void Calculate(
        const Variable<double>& rVariable, 
        double& Output,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

    //std::string Info() const;

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
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
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
	KinematicLinear() : BaseSolidElement()
    {
    }

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);
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

    void CalculateAndAddKm(
        MatrixType& K,
        Matrix& B,
        Matrix& D,
        double weight);

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    );

    void InitializeVariables();

    virtual void InitializeMaterial();

    double CalculateIntegrationWeight(
        double GaussPointWeight,
        double DetJ0
        );

    void CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& BodyForce,
        VectorType& mResidualVector,
        const double weight
    );

    void CalculateStrain(
        const Matrix& C,
        Vector& StrainVector
        );

    /**
     * Calculation of the Deformation Matrix B
     * @param B: The deformation matrix
     * @param DN_DX: The derivatives of the shape functions
     */
    void CalculateB(
        Matrix& B,
        const Matrix& DN_DX
        );
    
    /**
     * Calculation of the equivalent deformation gradient
     * @param StrainVector: The strain tensor (Voigt notation)
     * @return The deformation gradient F
     */
    Matrix ComputeEquivalentF(const Vector& StrainVector);

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
    /// Assignment operator.
    //KinematicLinear& operator=(const KinematicLinear& rOther);
    /// Copy constructor.
    //KinematicLinear(const KinematicLinear& rOther);
    ///@}

}; // Class KinematicLinear

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_KINEMATIC_LINEAR_H_INCLUDED  defined 
