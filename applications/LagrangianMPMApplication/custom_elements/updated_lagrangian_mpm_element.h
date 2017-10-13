//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//




#if !defined(KRATOS_UPDATED_LAGRANGIAN_MPM_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_MPM_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "custom_elements/meshless_base_element.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/model_part.h"

//#include "custom_geometries/meshless_geometry.h"

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

/// Short class definition.
/** Detail class definition.
*/
class UpdatedLagrangianMPMElement
    : public MeshlessBaseElement
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

    /// Counted pointer of UpdatedLagrangianMPMElement
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianMPMElement );
    ///@}

protected:

    /**
     * Flags related to the element computation
     */

    //KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    //KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    //KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    //KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */

    struct GeneralVariables
    {
    private:

        //variables including all integration points
        //const GeometryType::ShapeFunctionsGradientsType* pDN_De;
        //const Matrix* pDN_De;
//         const Vector* pNcontainer;

    public:

        StressMeasureType StressMeasure;

        //for axisymmetric use only
        double  CurrentRadius;
        double  ReferenceRadius;

        //general variables for large displacement use
        double  detF;
        double  detF0;
        double  detFT;
        double  detJ;
        double  integration_weight;
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  B;
        Matrix  F;
        Matrix  FT;
        Matrix  F0;
        Matrix  DN_DX;
        //Matrix  DN_De;
        Matrix  ConstitutiveMatrix;

        //variables including all integration points
        //GeometryType::JacobiansType J;
        //GeometryType::JacobiansType j;
        Matrix J;
        Matrix j;
        Matrix DeltaPosition;
        Matrix CurrentDisp;
        Matrix PreviousDisp;


    };


public:
    ///@name Type Definitions
    ///@{



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UpdatedLagrangianMPMElement(IndexType NewId, GeometryType::Pointer pGeometry);
    UpdatedLagrangianMPMElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~UpdatedLagrangianMPMElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;





    void  EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo );

    void  GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

    void  GetValuesVector( Vector& values, int Step );

    void  GetFirstDerivativesVector( Vector& values, int Step );

    void  GetSecondDerivativesVector( Vector& values, int Step );


    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    virtual void Initialize();

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void  CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void  CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);


    virtual void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);

    virtual void Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

protected:
    // A private default constructor necessary for serialization
    UpdatedLagrangianMPMElement() : MeshlessBaseElement()
    {
    }
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    /**
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    Matrix mDeformationGradientF0;
    /**
     * Container for the total deformation gradient determinants
     */
    double mDeterminantF0;


    /**
     * Currently selected integration methods
     */
    //IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    ConstitutiveLaw::Pointer mConstitutiveLawVector;

    /**
     * Finalize and Initialize label
     */
    bool mFinalizedStep;



    virtual void CalculateElementalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                    GeneralVariables& rVariables);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                    GeneralVariables& rVariables,
                                    Vector& rVolumeForce);


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    virtual void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables);


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            GeneralVariables& rVariables,
            Vector& rVolumeForce);


    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables);


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    virtual void SetGeneralVariables(GeneralVariables& rVariables,
                                     ConstitutiveLaw::Parameters& rValues);


    /**
     * Initialize Material Properties on the Constitutive Law
     */
    void InitializeMaterial ();


    /**
     * Reset the Constitutive Law Parameters
     */
    void ResetConstitutiveLaw();

    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);



    /**
     * Calculation of the Current Displacement
     */
    Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);


    /**
     * Initialize Element General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
      * Finalize Element Internal Variables
      */
    virtual void FinalizeStepVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
                                            Matrix& rF,
                                            Matrix& rDN_DX);

    /**
     * Calculation of the Volume Force of the Element
     */
    virtual Vector& CalculateVolumeForce(Vector& rVolumeForce, GeneralVariables& rVariables);

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


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    //friend class Serializer;


    //virtual void save(Serializer& rSerializer) const
    //{
    //KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MeshlessBaseElement);
    //}

    //virtual void load(Serializer& rSerializer)
    //{
    //KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MeshlessBaseElement);
    //}
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);
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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //UpdatedLagrangianMPMElement& operator=(const UpdatedLagrangianMPMElement& rOther);

    /// Copy constructor.
    //UpdatedLagrangianMPMElement(const UpdatedLagrangianMPMElement& rOther);


    ///@}

}; // Class UpdatedLagrangianMPMElement

///@}

///@name Type Definitions
///@{UpdatedLagrangianMPMElement


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    UpdatedLagrangianMPMElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const UpdatedLagrangianMPMElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_UPDATED_LAGRANGIAN_MPM_ELEMENT_H_INCLUDED  defined 


