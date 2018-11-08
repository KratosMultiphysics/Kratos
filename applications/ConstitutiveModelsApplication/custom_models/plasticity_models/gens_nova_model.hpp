//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                MCiantia $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_GENS_NOVA_MODEL_H_INCLUDED )
#define  KRATOS_GENS_NOVA_MODEL_H_INCLUDED

// System includes

// External includes
#include <iostream>
#include <fstream>

// Project includes
#include "custom_models/plasticity_models/non_associative_plasticity_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/gens_nova_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/gens_nova_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"
#include "custom_models/elasticity_models/tamagnini_model.hpp"



//***** the hardening law associated to this Model has ... variables
// 0. Plastic multiplier
// 1. Plastic Volumetric deformation
// 2. Plastic Deviatoric deformation
// 3. ps (mechanical)
// 4. pt (ageing)
// 5. pcSTAR = ps + (1+k) p_t
// ... (the number now is then..., xD)

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) GensNovaModel : public NonAssociativePlasticityModel<TamagniniModel, GensNovaYieldSurface<GensNovaHardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         //typedef BorjaModel                                     ElasticityModelType;
         typedef TamagniniModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef GensNovaHardeningRule                             HardeningRuleType;
         typedef GensNovaYieldSurface<HardeningRuleType>    YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

         //base type
         typedef NonAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of GensNovaModel
         KRATOS_CLASS_POINTER_DEFINITION( GensNovaModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         GensNovaModel() : BaseType() { mInitialized = false; }

         /// Copy constructor.
         GensNovaModel(GensNovaModel const& rOther) : BaseType(rOther), mInitialized(rOther.mInitialized) {}

         /// Assignment operator.
         GensNovaModel& operator=(GensNovaModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return ( GensNovaModel::Pointer(new GensNovaModel(*this)) );
         }

         /// Destructor.
         virtual ~GensNovaModel() {}


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

               const ModelDataType & rModelData = Variables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();

               double k = rMaterialProperties[KSIM];

               double & rPS     = Variables.Internal.Variables[3];
               double & rPT     = Variables.Internal.Variables[4];
               double & rPCstar = Variables.Internal.Variables[5];

               rPS = -rMaterialProperties[PS];
               rPT = -rMaterialProperties[PT];
               rPCstar = rPS + (1.0+k)*rPT;

               MatrixType Stress;
               this->UpdateInternalVariables(rValues, Variables, Stress);


               mInitialized = true;
            }
            mElasticityModel.InitializeModel( rValues );

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
         virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
         {
            KRATOS_TRY

            rValue=0;

            if (rThisVariable==PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0];
            }
            else if (rThisVariable==DELTA_PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0]-mPreviousInternal.Variables[0];
            }
            else if ( rThisVariable == PS)
            {
               rValue = this->mInternal.Variables[3];
            }
            else if ( rThisVariable == PT)
            {
               rValue = this->mInternal.Variables[4];
            }
            else if ( rThisVariable == PM)
            {
               rValue = this->mInternal.Variables[5];
            } else {
               rValue = NonAssociativePlasticityModel::GetValue( rThisVariable, rValue);
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
            buffer << "GensNovaModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "GensNovaModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "GensNovaModel Data";
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

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV


               MatrixType PlasticPotDerTensor;
               PlasticPotDerTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialDerivative, PlasticPotDerTensor);
               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerTensor);

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


               const ModelDataType & rModelData = rVariables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();
               const double & rhos = rMaterialProperties[RHOS];
               const double & rhot = rMaterialProperties[RHOT];
               double k =  rMaterialProperties[KSIM];
      
               const double & rChis = rMaterialProperties[CHIS];
               const double & rChit = rMaterialProperties[CHIT];

               MatrixType StressMatrix;
               // evaluate constitutive matrix and plastic flow
               double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
               double & rPlasticMultiplier = rVariables.Internal.Variables[0];
               double & rPlasticDevDef = rVariables.Internal.Variables[2];
               double & rPS     = rVariables.Internal.Variables[3];
               double & rPT     = rVariables.Internal.Variables[4];
               double & rPCstar = rVariables.Internal.Variables[5];

               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

               MatrixType PlasticPotDerTensor;
               PlasticPotDerTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialDerivative, PlasticPotDerTensor);
               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerTensor);

               MatrixType StrainMatrix = prod( rDeltaDeformationMatrix, trans( rDeltaDeformationMatrix) );
               VectorType StrainVector; 
               ConvertCauchyGreenTensorToHenckyVector( StrainMatrix, StrainVector);

               VectorType AuxVector;
               AuxVector = prod( ElasticMatrix, StrainVector);
               double DeltaGamma;
               DeltaGamma = MathUtils<double>::Dot( AuxVector, DeltaStressYieldCondition);

               double Denominador = H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) );

               DeltaGamma /= Denominador;

               if ( DeltaGamma < 0)
                  DeltaGamma = 0;

               MatrixType UpdateMatrix;
               ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);
               UpdateMatrix = prod( rDeltaDeformationMatrix, UpdateMatrix);


               rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
               rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

               this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

               rPlasticMultiplier += DeltaGamma;
               double VolPlasticIncr = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  VolPlasticIncr += DeltaGamma * DeltaStressYieldCondition(i);
               rPlasticVolDef += VolPlasticIncr;

               double DevPlasticIncr = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  DevPlasticIncr += pow( DeltaGamma * DeltaStressYieldCondition(i) - VolPlasticIncr/3.0, 2.0);
               for (unsigned int i = 3; i < 6; i++)
                  DevPlasticIncr += 2.0 * pow( DeltaGamma *  DeltaStressYieldCondition(i) /2.0 , 2.0);
               DevPlasticIncr = sqrt(DevPlasticIncr);
               rPlasticDevDef += DevPlasticIncr;


               double hs = rhos * ( rPS) * (     VolPlasticIncr  + rChis*sqrt(2.0/3.0) * DevPlasticIncr );
               double ht = rhot * (-rPT) * (fabs(VolPlasticIncr) + rChit*sqrt(2.0/3.0) * DevPlasticIncr );

               rPS -= hs;
               rPT -= ht;

               rPCstar = rPS + (1.0 + k)*rPT;

               KRATOS_CATCH("")
            }

            //********************************************************************
            //********************************************************************
            // UpdateInternalVariables
            virtual void UpdateInternalVariables(ModelDataType& rValues, PlasticDataType& rVariables, const MatrixType& rStressMatrix)
            {
               KRATOS_TRY

               for (unsigned int i = 0; i < 10; i++) {
                  double & plasticVolDefNew = rVariables.Internal.Variables[i]; 
                  double & plasticVolDef    = mInternal.Variables[i];

                  mPreviousInternal.Variables[i] = plasticVolDef;
                  plasticVolDef = plasticVolDefNew;
                  //std::cout <<  i << " , "  << mPreviousInternal.Variables[i] << " , " << plasticVolDef << std::endl;
               }

               KRATOS_CATCH("")
            }

            //***************************************************************************************
            //***************************************************************************************
            // Correct Yield Surface Drift According to 
            virtual void ReturnStressToYieldSurface( ModelDataType & rValues, PlasticDataType & rVariables)
            {
               KRATOS_TRY



               double Tolerance = 1e-7;

                  MatrixType StressMatrix;
                  this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);
               double YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);

               if ( fabs(YieldSurface) < Tolerance)
                  return;

               const ModelDataType & rModelData = rVariables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();
               double rhos = rMaterialProperties[RHOS];
               double rhot = rMaterialProperties[RHOT];
               double k =  rMaterialProperties[KSIM];
               // evaluate constitutive matrix and plastic flow
               double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
               //double & rPlasticMultiplier = rVariables.Internal.Variables[0];
               double & rPlasticDevDef = rVariables.Internal.Variables[2];
               double & rPS     = rVariables.Internal.Variables[3];
               double & rPT     = rVariables.Internal.Variables[4];
               double & rPCstar = rVariables.Internal.Variables[5];

               for (unsigned int i = 0; i < 150; i++) {

                  Matrix ElasticMatrix(6,6);
                  noalias(ElasticMatrix) = ZeroMatrix(6,6);
                  this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

                  VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
                  VectorType PlasticPotentialDerivative;
                  PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

                  MatrixType PlasticPotDerTensor;
                  PlasticPotDerTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialDerivative, PlasticPotDerTensor);
                  double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerTensor);

                  double DeltaGamma = YieldSurface;
                  DeltaGamma /= ( H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) ) );

                  MatrixType UpdateMatrix;
                  ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);

                  rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
                  rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

                  this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

                  double VolPlasticIncr = 0.0;
                  for (unsigned int i = 0; i < 3; i++)
                     VolPlasticIncr += DeltaGamma * DeltaStressYieldCondition(i);
                  rPlasticVolDef += VolPlasticIncr;

                  double DevPlasticIncr = 0.0;
                  for (unsigned int i = 0; i < 3; i++)
                     DevPlasticIncr += pow( DeltaGamma * DeltaStressYieldCondition(i) - VolPlasticIncr/3.0, 2.0);
                  for (unsigned int i = 3; i < 6; i++)
                     DevPlasticIncr += 2.0 * pow( DeltaGamma *  DeltaStressYieldCondition(i) /2.0 , 2.0);
                  rPlasticDevDef += sqrt(DevPlasticIncr);


                  double hs =  rhos * rPS * (VolPlasticIncr + 0.0);
                  double ht = rhot * rPT * ( -fabs(VolPlasticIncr)  + 0.0);
                  rPS -= hs;
                  rPT -= ht;

                  rPCstar = rPS + (1.0 + k)*rPT;


                  YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);


                  if ( fabs( YieldSurface) < Tolerance) {
                     return;
                  }
               }


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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Un accessible methods
         ///@{


         ///@}

   }; // Class GensNovaModel

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

#endif // KRATOS_GENS_NOVA_MODEL_H_INCLUDED  defined 


