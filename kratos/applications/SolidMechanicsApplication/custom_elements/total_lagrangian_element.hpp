//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/comparison_utils.hpp"


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

class TotalLagrangianElement
    : public Element
{
protected:

    /**
      * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
      */

    struct Standard
    {
        double  detF;
        double  detF0;
        double  detJ;
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  B;
        Matrix  F;
        Matrix  F0;
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;
    };

public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of TotalLagrangianElement
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangianElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TotalLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry);
    TotalLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    ///Copy constructor
    TotalLagrangianElement(TotalLagrangianElement const& rOther);


    /// Destructor.
    virtual ~TotalLagrangianElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TotalLagrangianElement& operator=(TotalLagrangianElement const& rOther);

    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod() const;

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    void Initialize();

    void ResetConstitutiveLaw();

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

    void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);

  void SetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo );

    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValuesVector(Vector& values, int Step = 0);
    void GetFirstDerivativesVector(Vector& values, int Step = 0);
    void GetSecondDerivativesVector(Vector& values, int Step = 0);


    void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

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
	TotalLagrangianElement() : Element()
    {
    }

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateElementalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
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
    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;
    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    double mTotalDomainInitialSize;
    std::vector< Matrix > mInvJ0;
    Vector mDetJ0;
    ///@}
    ///@name Private Operators
    ///@{

    void CalculateAndAddKm(
        MatrixType& rK,
        Matrix& rB,
        Matrix& rD,
        double& IntegrationWeight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = BT * S
     */
    void CalculateAndAddKg(
        MatrixType& rK,
        Matrix& rDN_DX,
        Vector& rStressVector,
        double& rIntegrationWeight
    );


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddExternalForces(const Vector& N,
                                       const ProcessInfo& CurrentProcessInfo,
                                       Vector& BodyForce,
                                       VectorType& mResidualVector,
                                       double& rIntegrationWeight
                                      );


    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    void CalculateAndAddInternalForces(Matrix & rB,
                                       Vector& rStressVector,
                                       VectorType& rRightHandSideVector,
                                       double& rIntegrationWeight
                                      );


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    void SetStandardParameters(Standard& rVariables,
                               ConstitutiveLaw::Parameters& rValues,
                               const int & rPointNumber);


    /**
     * Initialize Standard Variables
     */ 
    void InitializeStandardVariables(Standard & rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
				  VectorType& rRightHandSideVector,
				  bool CalculateStiffnessMatrixFlag,
				  bool CalculateResidualVectorFlag);


    /**
     * Correct Precision Errors (for rigid free movements)
     */
    void DecimalCorrection(Vector& rStrainVector);


    /**
     * Initialize Variables
     */
    void InitializeVariables ();

    /**
     * Initialize Material Properties on the Constitutive Law
     */
    void InitializeMaterial ();

    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();


    /**
     * Calculation of the Green Lagrange Strain Vector
     */
    void CalculateGreenLagrangeStrain(const Matrix& C,
				      Vector& StrainVector);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& B,
				    Matrix& F,
				    Matrix& DN_DX);

    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(Standard& rVariables,
                             const double& rPointNumber);

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
//        {
//            rSerializer.save("Name", "TotalLagrangianElement");
//            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
//        }

    virtual void load(Serializer& rSerializer);
//        {
//            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
//        }



    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //TotalLagrangianElement& operator=(const TotalLagrangianElement& rOther);
    /// Copy constructor.
    //TotalLagrangianElement(const TotalLagrangianElement& rOther);
    ///@}

}; // Class TotalLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
TotalLagrangianElement& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
const TotalLagrangianElement& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}*/
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED  defined 
