//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_wP_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_U_wP_FIC_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_wP_Stab_element.hpp"

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

/// Large Displacement Lagrangian U-P Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class UpdatedLagrangianUwPFICElement
    : public UpdatedLagrangianUwPStabElement
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

    /// Counted pointer of LargeDisplacementUPElement
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUwPFICElement );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianUwPFICElement();

    /// Default constructors
    UpdatedLagrangianUwPFICElement(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianUwPFICElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianUwPFICElement(UpdatedLagrangianUwPFICElement const& rOther);


    /// Destructor.
    virtual ~UpdatedLagrangianUwPFICElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianUwPFICElement& operator=(UpdatedLagrangianUwPFICElement const& rOther);


    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    /**
     * creates a new total lagrangian updated element pointer
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

    //************* GETTING METHODS

    //SET

    /**
     * Set a double  Value on the Element Constitutive Law
     */
    //void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);


    //GET:

    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    //void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    //void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    //void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo);

    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    //void Initialize();

    /**
    * Sets on rElementalDofList the degrees of freedom of the considered element geometry
    */
    //void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    //void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rValues the nodal displacements
     */
    //void GetValuesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal velocities
     */
    //void GetFirstDerivativesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal accelerations
     */
    //void GetSecondDerivativesVector(Vector& rValues, int Step = 0);


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

    /**
     * Calculate Element Kinematics
     */
    //virtual void CalculateKinematics(GeneralVariables& rVariables,
    //                                 const double& rPointNumber);


    /**
     * Calculation of the Deformation Gradient F
     */
    //void CalculateDeformationGradient(const Matrix& rDN_DX,
    //                                  Matrix& rF,
    //                                  Matrix& rDeltaPosition);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    //virtual void CalculateDeformationMatrix(Matrix& rB,
    //                                        Matrix& rF,
    //                                        Matrix& rDN_DX);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
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

    /***
       container for don't know what 
    ***/
    //std::vector< Matrix > mDeformationGradientF0;

    /***
       container for don't know what 
    ***/
    //Vector mDeterminantF0;


    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    // TO BE DESTROYED BECAUSE IT DOES NOT MAKE ANY SENCE
    //void  GetConstants(double& rScalingConstant, double& rWaterBulk, double& rDeltaTime, double& rPermeability);

    /**
     * Calculation and addition of the matrices of the LHS
     */

    //virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
    //                                GeneralVariables& rVariables,
    //                                double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    //virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
    //                                GeneralVariables& rVariables,
    //                                Vector& rVolumeForce,
    //                                double& rIntegrationWeight);

    /**
     * Initialize Element General Variables
     */
    //virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

   /**
     * Finalize Element Internal Variables
     */
    //virtual void FinalizeStepVariables(GeneralVariables & rVariables, const double& rPointNumber);


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    //virtual void SetGeneralVariables(GeneralVariables& rVariables,
    //                                 ConstitutiveLaw::Parameters& rValues,
    //                                 const int & rPointNumber);


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
     */
    //virtual void CalculateAndAddKuum(MatrixType& rK,
    //                                 GeneralVariables & rVariables,
    //                                 double& rIntegrationWeight
    //                                );

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    //virtual void CalculateAndAddKuug(MatrixType& rK,
    //                                 GeneralVariables & rVariables,
    //                                 double& rIntegrationWeight
    //                                );

    /**
     * Calculation of the Kup matrix
     */
    //virtual void CalculateAndAddKup (MatrixType& rK,
    //                                 GeneralVariables & rVariables,
    //                                 double& rIntegrationWeight
    //                                );

    /**
     * Calculation of the Kpu matrix
     */
    /*virtual void CalculateAndAddKpu(MatrixType& rK,
                                    GeneralVariables & rVariables,
                                    double& rIntegrationWeight
                                   );*/


    /**
     * Calculation of the Kpp matrix
     */
    /*virtual void CalculateAndAddKpp(MatrixType& rK,
                                    GeneralVariables & rVariables,
                                    double& rIntegrationWeight
                                   );*/


    /**
     * Calculation of the Kpp Stabilization Term matrix
     */
    virtual void CalculateAndAddKppStab(MatrixType& rK,
                                        GeneralVariables & rVariables,
                                        double& rIntegrationWeight
                                       );
    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    //void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
    //                                   GeneralVariables& rVariables,
    //                                   Vector& rVolumeForce,
    //                                   double& rIntegrationWeight
    //                                  );


    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    /*virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            double& rIntegrationWeight
                                              );*/


    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            double& rIntegrationWeight
                                                  );

    /**
      * Calculation of the Internal Forces due to sigma. Fi = B * sigma
      */
    //virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
    //                                   GeneralVariables & rVariables,
    //                                   double& rIntegrationWeight
    //                                  );

    /**
     * Initialize System Matrices
     */
    //void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
    //                              VectorType& rRightHandSideVector,
    //                              Flags& rCalculationFlags);

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


}; // Class UpdatedLagrangianUwPFICElement



} // namespace Kratos
#endif // KRATOS_____

