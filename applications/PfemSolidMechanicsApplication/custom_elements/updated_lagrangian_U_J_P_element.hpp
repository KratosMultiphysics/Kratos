//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED )
#define      KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_element.hpp"

#include "pfem_solid_mechanics_application_variables.h"

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

   /// Updated Lagrangian Large Displacement Lagrangian U-J-P Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)
   /// This is the element from Chapter 5, Vol 2. of Z


   class UpdatedLagrangianUJPElement
      : public UpdatedLagrangianUJElement
   {

      protected:
         typedef struct
         {
            unsigned int voigtsize;

            double Alpha;
            double Beta; 

            double NodalMeanStress;
            double ElementalMeanStress;

            double NodalJacobian;

            Vector StressVector;
            Vector StressVectorEC;

         } UJPElementVariables;

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
         KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUJPElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         UpdatedLagrangianUJPElement();

         /// Default constructors
         UpdatedLagrangianUJPElement(IndexType NewId, GeometryType::Pointer pGeometry);

         UpdatedLagrangianUJPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         UpdatedLagrangianUJPElement(UpdatedLagrangianUJPElement const& rOther);


         /// Destructor.
         virtual ~UpdatedLagrangianUJPElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         UpdatedLagrangianUJPElement& operator=(UpdatedLagrangianUJPElement const& rOther);


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

         //SET:

         //GET:

         /**
          * Get on rVariable a double Value from the Element Constitutive Law
          */
         void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

         void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

         void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo);

         //************* STARTING - ENDING  METHODS


         /**
          * Sets on rElementalDofList the degrees of freedom of the considered element geometry
          */
         void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

         /**
          * Sets on rResult the ID's of the element degrees of freedom
          */
         void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

         /**
          * Sets on rValues the nodal displacements
          */
         void GetValuesVector(Vector& rValues, int Step = 0);

         /**
          * Sets on rValues the nodal velocities
          */
         void GetFirstDerivativesVector(Vector& rValues, int Step = 0);

         /**
          * Sets on rValues the nodal accelerations
          */
         void GetSecondDerivativesVector(Vector& rValues, int Step = 0);


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
          * Variables to scale some equations (maybe)
          */
         double ScalingConstant1;
         double ScalingConstant2;

         /**
          * Container for historical total elastic deformation measure F0 = dx/dX
          */
         //std::vector< Matrix > mDeformationGradientF0; // LMV: Crec que no cal

         /**
          * Container for the total deformation gradient determinants
          */
         //Vector mDeterminantF0; // LMV: Crec que no cal


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
               ElementVariables& rVariables,
               double& rIntegrationWeight);

         /**
          * Initialize Element General Variables
          */
         virtual void InitializeElementVariables(ElementVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

         /**
          * Calculation and addition of the vectors of the RHS
          */

         virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementVariables& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight);


         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         virtual void CalculateAndAddKuum(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         virtual void CalculateAndAddKuug(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KuJ matrix
          */
         virtual void CalculateAndAddKuJ (MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the Kup matrix
          */
         virtual void CalculateAndAddKup (MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJu matrix
          */
         virtual void CalculateAndAddKJu(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJJ matrix
          */
         virtual void CalculateAndAddKJJ(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJp matrix
          */
         virtual void CalculateAndAddKJp(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpu matrix
          */
         virtual void CalculateAndAddKpu(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KpJ matrix
          */
         virtual void CalculateAndAddKpJ(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpp matrix
          */
         virtual void CalculateAndAddKpp(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKppStab(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJJ Stabilization Term matrix
          */
         virtual void CalculateAndAddKJJStab(MatrixType& rK,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Jacobian-Balance
          */
         virtual void CalculateAndAddJacobianForces( VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );



         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Stabilization for the Pressure
          */
         virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Stabilization for the Jacobian
          */
         virtual void CalculateAndAddStabilizedJacobian( VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               UJPElementVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Initialize System Matrices
          */
         void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
               VectorType& rRightHandSideVector,
               Flags& rCalculationFlags);

         //on integration points:
         /**
          * Calculate a double Variable on the Element Constitutive Law
          */
         void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);

         void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

         void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo);


         /*
          * Compute Some integration point variables only once
          */ 
         virtual void CalculateThisElementVariables( UJPElementVariables& rElementVariables, const ElementVariables& rVariables);
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


}; // Class UpdatedLagrangianUJPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED

