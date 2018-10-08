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


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) UpdatedLagrangianUJPElement
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

         } UJPElementData;

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
         Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

         /**
          * clones the selected element variables, creating a new one
          * @param NewId: the ID of the new element
          * @param ThisNodes: the nodes of the new element
          * @param pProperties: the properties assigned to the new element
          * @return a Pointer to the new element
          */
         Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

         //************* GETTING METHODS

         //SET:

         //GET:

         /**
          * Get on rVariable a double Value from the Element Constitutive Law
          */
         void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

         void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

         void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo) override;

         //************* STARTING - ENDING  METHODS


         /**
          * Sets on rElementalDofList the degrees of freedom of the considered element geometry
          */
         void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Sets on rResult the ID's of the element degrees of freedom
          */
         void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Sets on rValues the nodal displacements
          */
         void GetValuesVector(Vector& rValues, int Step = 0) override;

         /**
          * Sets on rValues the nodal velocities
          */
         void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

         /**
          * Sets on rValues the nodal accelerations
          */
         void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;


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

         void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
               ElementDataType& rVariables,
               double& rIntegrationWeight) override;

         /**
          * Initialize Element General Variables
          */
         void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Calculation and addition of the vectors of the RHS
          */

         void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementDataType& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight) override;


         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         void CalculateAndAddKuumElemUJP(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         void CalculateAndAddKuugElemUJP(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KuJ matrix
          */
         void CalculateAndAddKuJElemUJP (MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the Kup matrix
          */
         void CalculateAndAddKup (MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJu matrix
          */
         void CalculateAndAddKJuElemUJP(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJJ matrix
          */
         void CalculateAndAddKJJElemUJP(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJp matrix
          */
         void CalculateAndAddKJp(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpu matrix
          */
         void CalculateAndAddKpu(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KpJ matrix
          */
         void CalculateAndAddKpJ(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpp matrix
          */
         void CalculateAndAddKpp(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         void CalculateAndAddKppStab(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJJ Stabilization Term matrix
          */
         void CalculateAndAddKJJStabElemUJP(MatrixType& rK,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Jacobian-Balance
          */
         void CalculateAndAddJacobianForcesElemUJP( VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );



         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Stabilization for the Pressure
          */
         void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Stabilization for the Jacobian
          */
         void CalculateAndAddStabilizedJacobianElemUJP( VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         void CalculateAndAddInternalForcesElemUJP(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               UJPElementData& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Initialize System Matrices
          */
         void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
               VectorType& rRightHandSideVector,
               Flags& rCalculationFlags) override;

         //on integration points:
         /**
          * Calculate a double Variable on the Element Constitutive Law
          */
         void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

         void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

         void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;


         /*
          * Compute Some integration point variables only once
          */ 
         void CalculateThisElementData( UJPElementData& rElementVariables, const ElementDataType& rVariables);
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


}; // Class UpdatedLagrangianUJPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED

