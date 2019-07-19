//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:     LMonforte, MCiantia $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_STRUCTURED_SOIL_MODEL_H_INCLUDED )
#define      KRATOS_CASM_STRUCTURED_SOIL_MODEL_H_INCLUDED


// Clay And Sand Model (CASM) for structured soils

// System includes

// External includes
#include <iostream>
#include <fstream>

// Project includes
#include "custom_models/plasticity_models/non_associative_plasticity_model.hpp"


//***** the hardening law associated to this Model has ... variables
// Variables
// 0. Plastic Multiplier
// 1. Volumetric Plastic Strain
// 2. Deviatoric Plastic Strain
// 3. Abs Value Volumetric Plastic Strain
// 4. p0 (preconsolidation pressure of debonded soil)
// 5. b  (bonding)
// 6. pc (preconsolidation pressure of bonded soil)
// 7. pt (tensile strength)
// 8. NonLocal Vol Plastic Strain
// 9. NonLocal Dev Plastic Strain
//10. NonLocal Abs Vol Plastic Strain


namespace Kratos
{
   ///@addtogroup ConstitutiveModelsApplication
   ///@{

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
   template<class TElasticityModel, class TYieldSurface>
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CasmStructuredSoilModel : public NonAssociativePlasticityModel<TElasticityModel, TYieldSurface >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         typedef TElasticityModel                               ElasticityModelType;

         //yield surface
         typedef TYieldSurface                                     YieldSurfaceType;


         // derived type
         typedef NonAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>  DerivedType;

         //base type
         typedef PlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef typename BaseType::Pointer                         BaseTypePointer;
         typedef typename BaseType::SizeType                               SizeType;
         typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
         typedef typename BaseType::MatrixType                           MatrixType;
         typedef typename BaseType::VectorType                           VectorType;
         typedef typename BaseType::ModelDataType                     ModelDataType;
         typedef typename BaseType::MaterialDataType               MaterialDataType;
         typedef typename BaseType::PlasticDataType                 PlasticDataType;
         typedef typename BaseType::InternalVariablesType     InternalVariablesType;

         
         /// Pointer definition of CasmStructuredSoilModel
         KRATOS_CLASS_POINTER_DEFINITION( CasmStructuredSoilModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         CasmStructuredSoilModel() : DerivedType() { mInitialized = false; }

         /// Copy constructor.
         CasmStructuredSoilModel(CasmStructuredSoilModel const& rOther) : DerivedType(rOther), mInitialized(rOther.mInitialized) {}

         /// Assignment operator.
         CasmStructuredSoilModel& operator=(CasmStructuredSoilModel const& rOther)
         {
            DerivedType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return ( CasmStructuredSoilModel::Pointer(new CasmStructuredSoilModel(*this)) );
         }

         /// Destructor.
         virtual ~CasmStructuredSoilModel() {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{

         
         /**
          * Initialize member data
          */    
         void InitializeModel(ModelDataType& rValues) override
         {
            KRATOS_TRY

            if (mInitialized == false) {
               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               const ModelDataType& rModelData = Variables.GetModelData();
               const Properties& rMaterialProperties = rModelData.GetProperties();


               double& rP0 = Variables.Internal.Variables[4];
               double& rB  = Variables.Internal.Variables[5];
               double& rPc = Variables.Internal.Variables[6];
               double& rPt = Variables.Internal.Variables[7];
               const double& rAlphaTensile = rMaterialProperties[ALPHA_TENSILE];

               rP0 = -rMaterialProperties[PRE_CONSOLIDATION_STRESS];
               rB  = rMaterialProperties[INITIAL_BONDING];
               rPc = rP0*(1+rB);
               rPt = rP0*(rAlphaTensile*rB);

               MatrixType Stress;
               this->UpdateInternalVariables(rValues, Variables, Stress);


               mInitialized = true;
            }
            this->mElasticityModel.InitializeModel( rValues );

            KRATOS_CATCH("")
         }

         /**
          * Check
          */    
         virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
         {
            KRATOS_TRY

            //LMV: to be implemented. but should not enter in the base one

            return 0;

            KRATOS_CATCH("")
         }

         ///@}
         ///@name Access
         ///@{

         /**
          * Has Values
          */   
         virtual bool Has(const Variable<double>& rThisVariable) override
         {
            if(rThisVariable == PLASTIC_STRAIN || rThisVariable == DELTA_PLASTIC_STRAIN )
               return true;

            return false;
         }


         /**
          * Get Values
          */
         void SetValue(const Variable<double>& rVariable,
               const double& rValue,
               const ProcessInfo& rCurrentProcessInfo) override 
         {
            KRATOS_TRY

            if ( rVariable == NONLOCAL_PLASTIC_VOL_DEF) {
               this->mInternal.Variables[8] = rValue;
            }
            else if ( rVariable == NONLOCAL_PLASTIC_DEV_DEF) {
               this->mInternal.Variables[9] = rValue;
            }
            else if ( rVariable == NONLOCAL_PLASTIC_VOL_DEF_ABS) {
               this->mInternal.Variables[10] = rValue;
            } else {
               BaseType::SetValue( rVariable, rValue, rCurrentProcessInfo);
            }

            KRATOS_CATCH("")
         }


         /**
          * Get Values
          */
         virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
         {
            KRATOS_TRY

            rValue=0;

            if (rThisVariable==PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[1];
            }
            else if (rThisVariable==DELTA_PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[1]-this->mPreviousInternal.Variables[1];
            }
            else if ( rThisVariable == P0)
            {
               rValue = this->mInternal.Variables[4];
            }
            else if ( rThisVariable == BONDING)
            {
               rValue = this->mInternal.Variables[5];
            }
            else if ( rThisVariable == PC)
            {
               rValue = this->mInternal.Variables[6];
            } 
            else if ( rThisVariable == PT)
            {
               rValue = this->mInternal.Variables[7];
            } 
            else if ( rThisVariable == PLASTIC_VOL_DEF)
            {
               rValue = this->mInternal.Variables[1];
            }
            else if ( rThisVariable == PLASTIC_DEV_DEF)
            {
               rValue = this->mInternal.Variables[2];
            }
            else if ( rThisVariable == PLASTIC_VOL_DEF_ABS)
            {
               rValue = this->mInternal.Variables[3];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_VOL_DEF)
            {
               rValue = this->mPreviousInternal.Variables[8];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_DEV_DEF)
            {
               rValue = this->mPreviousInternal.Variables[9];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_VOL_DEF_ABS)
            {
               rValue = this->mPreviousInternal.Variables[10];
            }
            else {
               rValue = NonAssociativePlasticityModel<TElasticityModel, TYieldSurface>::GetValue( rThisVariable, rValue);
            }
            return rValue;

            KRATOS_CATCH("")
         }


         ///@}
         ///@name Inquiry
         ///@{


         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         virtual std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "CasmStructuredSoilModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "CasmStructuredSoilModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "CasmStructuredSoilModel Data";
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

         bool mInitialized; 

         ///@}
         ///@name Protected Operators
         ///@{


         ///@}
         ///@name Protected Operations
         ///@{

         
         
            //***************************************************************************************
            //***************************************************************************************
            // Compute Elasto Plastic Matrix
            void ComputeElastoPlasticTangentMatrix( ModelDataType & rValues, PlasticDataType & rVariables, Matrix & rEPMatrix) override
            {

               KRATOS_TRY

               // evaluate constitutive matrix and plastic flow
               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               //calculate yield function and plastic potential derivatives
               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative = this->mYieldSurface.CalculateDeltaPlasticPotential( rVariables, PlasticPotentialDerivative);

               //calcualte plastic hardening modulus H using dG/dInv
               MatrixType PlasticPotDerInvTensor;
               VectorType PlasticPotentialInvDerivative = this->mYieldSurface.CalculateDeltaStressInvPlasticPotential( rVariables, PlasticPotentialInvDerivative);
               PlasticPotDerInvTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialInvDerivative, PlasticPotDerInvTensor);
               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerInvTensor);

               VectorType AuxF = prod( trans(DeltaStressYieldCondition), rEPMatrix);
               VectorType AuxG = prod( rEPMatrix, PlasticPotentialDerivative);

               Matrix PlasticUpdateMatrix(6,6);
               noalias(PlasticUpdateMatrix) = ZeroMatrix(6,6);
               double denom = 0;
               for (unsigned int i = 0; i < 6; i++) {
                  denom += AuxF(i)*PlasticPotentialDerivative(i);
                  for (unsigned int j = 0; j < 6; j++) {
                     PlasticUpdateMatrix(i,j) = AuxF(i) * AuxG(j);
                  }
               }

               rEPMatrix -= PlasticUpdateMatrix / ( H + denom);

               KRATOS_CATCH("")
            }
            //***********************************************************************************
            //***********************************************************************************
            // Compute one step of the elasto-plastic problem
            void ComputeOneStepElastoPlasticProblem( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType & rDeltaDeformationMatrix) override
            {
               KRATOS_TRY
            
               //get constants
               const ModelDataType & rModelData = rVariables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();
               const double& rInitialP0      = rMaterialProperties[PRE_CONSOLIDATION_STRESS];
               const double& rOmega          = rMaterialProperties[PLASTIC_DEVIATORIC_STRAIN_HARDENING];
               const double& rSwellingSlope  = rMaterialProperties[SWELLING_SLOPE];
               const double& rOtherSlope     = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];
               const double& rInitialB       = rMaterialProperties[INITIAL_BONDING];
               const double& rH0             = rMaterialProperties[DEGRADATION_THRESHOLD];
               const double& rH1             = rMaterialProperties[DEGRADATION_RATE_COMPRESSION];
               const double& rH2             = rMaterialProperties[DEGRADATION_RATE_SHEAR];
               const double& rAlphaTensile   = rMaterialProperties[ALPHA_TENSILE];

               //get internal variables
               double& rPlasticMultiplier    = rVariables.Internal.Variables[0];
               double& rPlasticVolDef        = rVariables.Internal.Variables[1]; 
               double& rPlasticDevDef        = rVariables.Internal.Variables[2];
               double& rAbsPlasticDevDef     = rVariables.Internal.Variables[3];
               double& rP0                   = rVariables.Internal.Variables[4];
               double& rB                    = rVariables.Internal.Variables[5];
               double& rPc                   = rVariables.Internal.Variables[6];
               double& rPt                   = rVariables.Internal.Variables[7];

               //calculate elasticity matrix
               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               //calculate yield function and plastic potential derivatives
               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative = this->mYieldSurface.CalculateDeltaPlasticPotential( rVariables, PlasticPotentialDerivative);

               //calcualte plastic hardening modulus H using dG/dInv
               MatrixType PlasticPotDerInvTensor;
               VectorType PlasticPotentialInvDerivative = this->mYieldSurface.CalculateDeltaStressInvPlasticPotential( rVariables, PlasticPotentialInvDerivative);
               PlasticPotDerInvTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialInvDerivative, PlasticPotDerInvTensor);
               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerInvTensor);

               //calculate incremental Hencky strain vector
               MatrixType StrainMatrix = prod( rDeltaDeformationMatrix, trans( rDeltaDeformationMatrix) );
               VectorType StrainVector; 
               this->ConvertCauchyGreenTensorToHenckyVector( StrainMatrix, StrainVector);

               //calculate plastic multiplier
               VectorType AuxVector;
               AuxVector = prod( ElasticMatrix, StrainVector);
               double DeltaGamma;
               DeltaGamma = MathUtils<double>::Dot( AuxVector, DeltaStressYieldCondition);

               double Denominador = H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) );

               DeltaGamma /= Denominador;

               if ( DeltaGamma < 0)
                  DeltaGamma = 0;

               //calculate new b_e based on exponential variation of the plastic def gradient (Simo 1998)
               MatrixType UpdateMatrix;
               this->ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);
               UpdateMatrix = prod( rDeltaDeformationMatrix, UpdateMatrix);

               rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
               rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

               //calculate new stress matrix
               MatrixType StressMatrix;
               this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

               //update plastic multiplier and plastic strains
               rPlasticMultiplier += DeltaGamma;

               double VolPlasticIncr = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  VolPlasticIncr += DeltaGamma * PlasticPotentialDerivative(i);
               rPlasticVolDef += VolPlasticIncr;
               rAbsPlasticDevDef += std::fabs(VolPlasticIncr);

               double DevPlasticIncr = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  DevPlasticIncr += pow( DeltaGamma*PlasticPotentialDerivative(i) - VolPlasticIncr/3.0, 2.0);
               for (unsigned int i = 3; i < 6; i++)
                  DevPlasticIncr += 2.0 * pow( DeltaGamma*PlasticPotentialDerivative(i)/2.0 , 2.0);
               DevPlasticIncr = sqrt(2.0/3.0*DevPlasticIncr);
               rPlasticDevDef += DevPlasticIncr;

               //update P0, b, Pc, Pt
               rP0 = -rInitialP0*std::exp( (-rPlasticVolDef+rOmega*rPlasticDevDef)/(rOtherSlope-rSwellingSlope) );
               double DamageH = rH1*std::fabs(rAbsPlasticDevDef) + rH2*std::fabs(rPlasticDevDef);
               rB = rInitialB*std::exp(rH0 - DamageH);
               rPc = rP0*(1+rB);
               rPt = rP0*(rAlphaTensile*rB);

               KRATOS_CATCH("")
            }

            //********************************************************************
            //********************************************************************
            // UpdateInternalVariables
            virtual void UpdateInternalVariables(ModelDataType& rValues, PlasticDataType& rVariables, const MatrixType& rStressMatrix) override
            {
               KRATOS_TRY

               //this->mPreviousInternal.Variables[3] = this->mInternal.Variables[3];
               //this->mInternal.Variables[3] = this->mInternal.Variables[3] + fabs( rVariables.Internal.Variables[1] - this->mInternal.Variables[1]);
               for (unsigned int i = 0; i < 11; i++) {
                  //if ( i != 3) {
                     double & rCurrentPlasticVariable = rVariables.Internal.Variables[i]; 
                     double & rPreviousPlasticVariable    = this->mInternal.Variables[i];

                     this->mPreviousInternal.Variables[i] = rPreviousPlasticVariable;
                     rPreviousPlasticVariable = rCurrentPlasticVariable;
                  //}
               }

               KRATOS_CATCH("")
            }

            //***************************************************************************************
            //***************************************************************************************
            // Correct Yield Surface Drift According to 
            virtual void ReturnStressToYieldSurface( ModelDataType & rValues, PlasticDataType & rVariables) override
            {
               KRATOS_TRY

               double Tolerance = 1e-6;

               //check yield surface violation
               MatrixType StressMatrix;
               this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);
               double YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);

               if ( fabs(YieldSurface) < Tolerance)
                  return;

               //get constants
               const ModelDataType & rModelData = rVariables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();
               const double& rInitialP0      = rMaterialProperties[PRE_CONSOLIDATION_STRESS];
               const double& rOmega          = rMaterialProperties[PLASTIC_DEVIATORIC_STRAIN_HARDENING];
               const double& rSwellingSlope  = rMaterialProperties[SWELLING_SLOPE];
               const double& rOtherSlope     = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];
               const double& rInitialB       = rMaterialProperties[INITIAL_BONDING];
               const double& rH0             = rMaterialProperties[DEGRADATION_THRESHOLD];
               const double& rH1             = rMaterialProperties[DEGRADATION_RATE_COMPRESSION];
               const double& rH2             = rMaterialProperties[DEGRADATION_RATE_SHEAR];
               const double& rAlphaTensile   = rMaterialProperties[ALPHA_TENSILE];

               //get internal variables
               double& rPlasticVolDef        = rVariables.Internal.Variables[1]; 
               double& rPlasticDevDef        = rVariables.Internal.Variables[2];
               double& rAbsPlasticDevDef     = rVariables.Internal.Variables[3];
               double& rP0                   = rVariables.Internal.Variables[4];
               double& rB                    = rVariables.Internal.Variables[5];
               double& rPc                   = rVariables.Internal.Variables[6];
               double& rPt                   = rVariables.Internal.Variables[7];

               for (unsigned int i = 0; i < 150; i++) {

                  //calculate elasticity matrix
                  Matrix ElasticMatrix(6,6);
                  noalias(ElasticMatrix) = ZeroMatrix(6,6);
                  this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

                  //calculate yield function and plastic potential derivatives
                  VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
                  VectorType PlasticPotentialDerivative = this->mYieldSurface.CalculateDeltaPlasticPotential( rVariables, PlasticPotentialDerivative);

                  //calcualte plastic hardening modulus H using dG/dInv
                  MatrixType PlasticPotDerInvTensor;
                  VectorType PlasticPotentialInvDerivative = this->mYieldSurface.CalculateDeltaStressInvPlasticPotential( rVariables, PlasticPotentialInvDerivative);
                  PlasticPotDerInvTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialInvDerivative, PlasticPotDerInvTensor);
                  double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerInvTensor);

                  //calculate plastic multiplier based on yield surface violation
                  double DeltaGamma = YieldSurface;
                  DeltaGamma /= ( H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) ) );

                  //calculate new b_e based exp variation of F_p
                  MatrixType UpdateMatrix;
                  this->ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);

                  rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
                  rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

                  //calculate new stress matrix
                  this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

                  //update plastic strains
                  double VolPlasticIncr = 0.0;
                  for (unsigned int i = 0; i < 3; i++)
                     VolPlasticIncr += DeltaGamma * PlasticPotentialDerivative(i);
                  rPlasticVolDef += VolPlasticIncr;
                  rAbsPlasticDevDef += VolPlasticIncr;//?????? no fabs()

                  double DevPlasticIncr = 0.0;
                  for (unsigned int i = 0; i < 3; i++)
                     DevPlasticIncr += pow( DeltaGamma * PlasticPotentialDerivative(i) - VolPlasticIncr/3.0, 2.0);
                  for (unsigned int i = 3; i < 6; i++)
                     DevPlasticIncr += 2.0 * pow( DeltaGamma *  DeltaStressYieldCondition(i) /2.0 , 2.0);
                  DevPlasticIncr = DeltaGamma/fabs(DeltaGamma) * sqrt(2.0/3.0*DevPlasticIncr);
                  rPlasticDevDef += DevPlasticIncr;

                  //update P0, b, Pc, Pt
                  rP0 = -rInitialP0*std::exp( (-rPlasticVolDef+rOmega*rPlasticDevDef)/(rOtherSlope-rSwellingSlope) );
                  double DamageH = rH1*std::fabs(rAbsPlasticDevDef) + rH2*std::fabs(rPlasticDevDef);
                  rB = rInitialB*std::exp(rH0 - DamageH);
                  rPc = rP0*(1+rB);
                  rPt = rP0*(rAlphaTensile*rB);

                  //check yield condition
                  YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);
                  if ( fabs( YieldSurface) < Tolerance) {
                     return;
                  }
               }
               std::cout << " TheStressPointDidNotReturnCorrectly " << YieldSurface << std::endl;

               KRATOS_CATCH("")
            }
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
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Serialization
         ///@{    
         friend class Serializer;

         virtual void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
         }

         ///@}
         ///@name Un accessible methods
         ///@{


         ///@}

   }; // Class CasmStructuredSoilModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@} 
   ///@name Input and output 
   ///@{


   ///@}

   ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_CASM_STRUCTURED_SOIL_MODEL_H_INCLUDED  defined 


