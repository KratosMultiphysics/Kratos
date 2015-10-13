//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_wP_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_wP_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_wP_element.hpp"

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

/// Axisymmetric Updated Lagrangian Large Displacement Lagrangian U-Pw Element.


class AxisymUpdatedLagrangianUwPElement
    : public UpdatedLagrangianUwPElement
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
    KRATOS_CLASS_POINTER_DEFINITION( AxisymUpdatedLagrangianUwPElement );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    AxisymUpdatedLagrangianUwPElement();

    /// Default constructors
    AxisymUpdatedLagrangianUwPElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymUpdatedLagrangianUwPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymUpdatedLagrangianUwPElement(AxisymUpdatedLagrangianUwPElement const& rOther);


    /// Destructor.
    virtual ~AxisymUpdatedLagrangianUwPElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymUpdatedLagrangianUwPElement& operator=(AxisymUpdatedLagrangianUwPElement const& rOther);


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
                                      Matrix& rDeltaPosition,
                                      double& rCurrentRadius,
                                      double& rReferenceRadius);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
                                            Matrix& rF,
                                            Vector& rN,
                                            double& rCurrentRadius);


    virtual void CalculateRadius(double & rCurrentRadius, double & rReferenceRadius, const Vector& rN);

    virtual void Initialize();

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


    /**
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    //std::vector< Matrix > mDeformationGradientF0;


    /**
     * Container for the total deformation gradient determinants
     */
    //Vector mDeterminantF0;


    /**** 
       the time step (requiered). It shall be somewhere else.
    ****/    
    //double mTimeStep;

    /*** 
        Just to check a few things
     ***/
    //bool mCompressibleWater;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight);

    /**
     * Initialize Element General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

   /**
     * Finalize Element Internal Variables
     */
    //virtual void FinalizeStepVariables(GeneralVariables & rVariables, const double& rPointNumber);



    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
     */
    virtual void CalculateAndAddKuum(MatrixType& rK,
                                     GeneralVariables & rVariables,
                                     double& rIntegrationWeight
                                    );

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rK,
                                     GeneralVariables & rVariables,
                                     double& rIntegrationWeight
                                    );

    /**
     * Calculation of the Kup matrix
     */
    virtual void CalculateAndAddKup (MatrixType& rK,
                                     GeneralVariables & rVariables,
                                     double& rIntegrationWeight
                                    );

    /**
     * Calculation of the Kpu matrix
     */
    virtual void CalculateAndAddKpu(MatrixType& rK,
                                    GeneralVariables & rVariables,
                                    double& rIntegrationWeight
                                   );


    /**
     * Calculation of the Kpp matrix
     */
    virtual void CalculateAndAddKpp(MatrixType& rK,
                                    GeneralVariables & rVariables,
                                    double& rIntegrationWeight
                                   );


    /**
     * Calculation of the Kpp Stabilization Term matrix
     */
    virtual void CalculateAndAddKppStab(MatrixType& rK,
                                        GeneralVariables & rVariables,
                                        double& rIntegrationWeight
                                       );

    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            double& rIntegrationWeight
                                              );


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
    virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
                                       GeneralVariables & rVariables,
                                       double& rIntegrationWeight
                                      );

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


}; // Class UpdatedLagrangianUwPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_wP_ELEMENT_H_INCLUDED

