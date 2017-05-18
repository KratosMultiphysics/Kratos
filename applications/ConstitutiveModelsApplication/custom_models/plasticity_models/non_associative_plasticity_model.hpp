//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NON_ASSOCIATIVE_PLASTICITY_MODEL_H_INCLUDED)
#define KRATOS_NON_ASSOCIATIVE_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/plasticity_model.hpp"

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
      template<class TElasticityModel, class TYieldCriterion>
      class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonAssociativePlasticityModel : public PlasticityModel<TElasticityModel,TYieldCriterion>
      {
         public:

            ///@name Type Definitions
            ///@{

            //elasticity model
            typedef TElasticityModel                               ElasticityModelType;

            //yield criterion
            typedef TYieldCriterion                                 YieldCriterionType;

            //base type
            typedef PlasticityModel<ElasticityModelType,YieldCriterionType>   BaseType;

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

            typedef ConstitutiveModelData::StrainMeasureType         StrainMeasureType;   
            typedef ConstitutiveModelData::StressMeasureType         StressMeasureType;   


            /// Pointer definition of NonAssociativePlasticityModel
            KRATOS_CLASS_POINTER_DEFINITION( NonAssociativePlasticityModel );
            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            NonAssociativePlasticityModel() : BaseType() {}

            /// Copy constructor.
            NonAssociativePlasticityModel(NonAssociativePlasticityModel const& rOther) :BaseType(rOther), mInternal(rOther.mInternal), mPreviousInternal(rOther.mPreviousInternal) {}

            /// Assignment operator.
            NonAssociativePlasticityModel& operator=(NonAssociativePlasticityModel const& rOther)
            {
               BaseType::operator=(rOther);
               mInternal = rOther.mInternal;
               mPreviousInternal = rOther.mPreviousInternal;
               return *this;
            }

            /// Clone.
            virtual ConstitutiveModel::Pointer Clone() const override
            {
               return ( NonAssociativePlasticityModel::Pointer(new NonAssociativePlasticityModel(*this)) );
            }

            /// Destructor.
            virtual ~NonAssociativePlasticityModel() {}


            ///@}
            ///@name Operators
            ///@{


            ///@}
            ///@name Operations
            ///@{

            /**
             * Calculate Stresses
             */

            virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
            {
               KRATOS_TRY

               this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);

               rValues.StressMatrix = rStressMatrix; 

               KRATOS_CATCH(" ")

            }

            /**
             * Calculate Constitutive Tensor
             */
            virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               //Initialize ConstitutiveMatrix
               rConstitutiveMatrix.clear();

               MatrixType StressMatrix = rValues.StressMatrix;
               //1.-Elastic Stress Matrix
               this->mElasticityModel.CalculateStressTensor(rValues,rValues.StressMatrix);

               // calculate elastic constitutive tensor
               this->mElasticityModel.CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);


               KRATOS_CATCH(" ")
            }

            //*******************************************************************************
            //*******************************************************************************
            // Calcualte Stress and constitutive tensor
            virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               std::cout << " mInternal Size " << mInternal.Variables << std::endl;
               double Tolerance = 1e-5;

               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               //Calculate trial stress Matrix
               this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);
               //rValues.StressMatrix = rStressMatrix;

               Variables.TrialStateFunction = this->mYieldCriterion.CalculateYieldCondition( Variables, Variables.TrialStateFunction);

               if ( Variables.TrialStateFunction  < Tolerance) {
                  // elastic loading step
                  rConstitutiveMatrix.clear();
                  this->mElasticityModel.CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);
                  return;
               }

               rConstitutiveMatrix.clear();
               this->mElasticityModel.CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);

               // elasto-plastic step. Recover Initial be
               const MatrixType & rIncrementalDeformationGradient = rValues.GetIncrementalDeformationGradientF();
               RecoverPreviousElasticLeftCauchyGreen( rIncrementalDeformationGradient, rValues.StrainMatrix );

               double InitialYieldFunction;
               this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);
               InitialYieldFunction = this->mYieldCriterion.CalculateYieldCondition( Variables, InitialYieldFunction);
               double H = this->mYieldCriterion.GetHardeningLaw().CalculateHardening( Variables, H);
               if ( (InitialYieldFunction < -Tolerance) && (Variables.TrialStateFunction > Tolerance) )
               {
                  // compute solution with change
                  ComputeSolutionWithChange( rValues, Variables, rIncrementalDeformationGradient);
               } else {
                  // compute unloading condition
                  ComputeElastoPlasticProblem( rValues, Variables, rIncrementalDeformationGradient);
                  bool UnloadingCondition = false;
                  if (UnloadingCondition) {
                     // compute solution with change
                  } else {
                     // compute plastic problem
                  }
               }

               noalias( rStressMatrix) = rValues.StressMatrix;

               if ( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES) )
                  this->UpdateInternalVariables( rValues, Variables);

               KRATOS_CATCH(" ")
            }

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
            virtual std::string Info() const override
            {
               std::stringstream buffer;
               buffer << "NonAssociativePlasticityModel" ;
               return buffer.str();
            }

            /// Print information about this object.
            virtual void PrintInfo(std::ostream& rOStream) const override
            {
               rOStream << "NonAssociativePlasticityModel";
            }

            /// Print object's data.
            virtual void PrintData(std::ostream& rOStream) const override
            {
               rOStream << "NonAssociativePlasticityModel Data";	    
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

            // internal variables:
            InternalVariablesType  mInternal;
            InternalVariablesType  mPreviousInternal;

            ///@}
            ///@name Protected Operators
            ///@{


            ///@}
            ///@name Protected Operations
            ///@{


            //***************************************************************************************
            //***************************************************************************************
            // Advance the solution first in elastic regime and then in elastoplastic
            void ComputeSolutionWithChange( ModelDataType& rValues, PlasticDataType & rVariables, const MatrixType& rIncrementalDeformationGradient)
            {
               KRATOS_TRY
               double Tolerance = 1e-5;

               double InitialTime = 0; double EndTime = 1; double HalfTime;
               double InitialStateFunction(-1), EndStateFunction(1), HalfTimeStateFunction;

               MatrixType HalfTimeDeformationGradient;
               MatrixType StressMatrix;
               MatrixType InitialLeftCauchyGreen = rValues.StrainMatrix;
               for (unsigned int i = 0; i < 150; i++)
               {
                  HalfTime = 0.5*(InitialTime + EndTime);

                  ComputeSubstepIncrementalDeformationGradient( rIncrementalDeformationGradient, 0, HalfTime, HalfTimeDeformationGradient);
                  MatrixType AuxMatrix; 
                  AuxMatrix = prod( InitialLeftCauchyGreen, trans(HalfTimeDeformationGradient));
                  AuxMatrix = prod( HalfTimeDeformationGradient, AuxMatrix);
                  rValues.StrainMatrix = AuxMatrix;

                  this->mElasticityModel.CalculateStressTensor(rValues,StressMatrix);
                  HalfTimeStateFunction = this->mYieldCriterion.CalculateYieldCondition( rVariables, HalfTimeStateFunction);

                  if ( HalfTimeStateFunction < 0.0) {
                     InitialStateFunction = HalfTimeStateFunction;
                     InitialTime = HalfTime;
                  } else {
                     EndStateFunction = HalfTimeStateFunction;
                     EndTime = HalfTime;
                  }

                  double ErrorMeasure1 = fabs( InitialStateFunction - EndStateFunction);
                  double ErrorMeasure2 = fabs( InitialTime-EndTime);

                  if ( (ErrorMeasure1 < Tolerance) && (ErrorMeasure2 < Tolerance) )
                     break;
               }

               
               // continue with plasticity
               MatrixType RemainingDeformationGradient;
               ComputeSubstepIncrementalDeformationGradient( rIncrementalDeformationGradient, HalfTime, 1,RemainingDeformationGradient);

               ComputeElastoPlasticProblem( rValues, rVariables, RemainingDeformationGradient);



               KRATOS_CATCH("")
            }

            //***********************************************************************************
            //***********************************************************************************
            // Compute the elasto-plastic problem
            void ComputeElastoPlasticProblem( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType & rIncrementalDeformationGradient)
            {
               KRATOS_TRY

               // evaluate constitutive matrix and plastic flow
               Matrix ElasticMatrix = ZeroMatrix(6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType YieldSurfaceDerivative;
               this->mYieldCriterion.CalculateYieldSurfaceDerivative( rVariables, YieldSurfaceDerivative);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = YieldSurfaceDerivative; // LMV

               double H = this->mYieldCriterion.GetHardeningLaw().CalculateDeltaHardening( rVariables, H);

               MatrixType StrainMatrix = prod( rIncrementalDeformationGradient, trans( rIncrementalDeformationGradient) );
               VectorType StrainVector; 
               ConvertCauchyGreenTensorToHenckyVector( StrainMatrix, StrainVector);

               VectorType AuxVector;
               AuxVector = prod( ElasticMatrix, StrainVector);
               double DeltaGamma;
               DeltaGamma = MathUtils<double>::Dot( AuxVector, YieldSurfaceDerivative);

               double Denominador = H + MathUtils<double>::Dot( YieldSurfaceDerivative, prod(ElasticMatrix, PlasticPotentialDerivative) );

               DeltaGamma /= Denominador;

               if ( DeltaGamma < 0)
                  DeltaGamma = 0;
               
               MatrixType UpdateMatrix;
               ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);
               UpdateMatrix = prod( rIncrementalDeformationGradient, UpdateMatrix);


               rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
               rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

               MatrixType StressMatrix;
               this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

               double & plasticVolDef = rVariables.Internal.Variables[1]; 
               for (unsigned int i = 0; i < 3; i++)
                  plasticVolDef += DeltaGamma * YieldSurfaceDerivative(i);

               KRATOS_CATCH("")
            }
            //**************************************************************************************
            //**************************************************************************************
            // Convert epsilon (vector) to b
            void ConvertHenckyVectorToCauchyGreenTensor(const VectorType&  rHenckyVector, MatrixType & rStrainMatrix)
            {
               KRATOS_TRY
              
               MatrixType HenckyTensor;
               HenckyTensor.clear();

               ConstitutiveModelUtilities::StrainVectorToTensor( rHenckyVector, HenckyTensor);
               ConvertHenckyTensorToCauchyGreenTensor( HenckyTensor, rStrainMatrix);


               KRATOS_CATCH("")
            }
            //**************************************************************************************
            //**************************************************************************************
            // Convert epsilon (matrix) to b
            void ConvertHenckyTensorToCauchyGreenTensor(const MatrixType&  rHenckyTensor, MatrixType & rStrainMatrix)
            {
               KRATOS_TRY
               
               MatrixType EigenVectors;
               EigenVectors.clear();

               rStrainMatrix.clear();
               MathUtils<double>::EigenSystem<3> ( rHenckyTensor, EigenVectors, rStrainMatrix);

               for (unsigned int i = 0; i < 3; i++)
                  rStrainMatrix(i,i) = std::exp( 2.0* rStrainMatrix(i,i));


               rStrainMatrix = prod( trans(EigenVectors), MatrixType(prod(rStrainMatrix, EigenVectors)) );


               KRATOS_CATCH("")
            }
            //**************************************************************************************
            //**************************************************************************************
            // Convert b to epsilon (Vector) // vale, m'he equivocat i no el necessito, pfpfpf
            void ConvertCauchyGreenTensorToHenckyTensor(const MatrixType&  rStrainMatrix, MatrixType & rHenckyStrain)
            {
               KRATOS_TRY
               
               MatrixType EigenVectors;
               EigenVectors.clear();

               rHenckyStrain.clear();
               MathUtils<double>::EigenSystem<3> ( rStrainMatrix, EigenVectors, rHenckyStrain);

               for (unsigned int i = 0; i < 3; i++)
                  rHenckyStrain(i,i) = std::log( rHenckyStrain(i,i))/2.0;

               rHenckyStrain = prod( trans(EigenVectors), MatrixType(prod(rHenckyStrain, EigenVectors)) );


               KRATOS_CATCH("")
            }

            //**************************************************************************************
            //**************************************************************************************
            // Convert b to epsilon (Vector)
            void ConvertCauchyGreenTensorToHenckyVector(const MatrixType&  rStrainMatrix, VectorType & rStrainVector)
            {
               KRATOS_TRY

               MatrixType HenckyStrain;
               ConvertCauchyGreenTensorToHenckyTensor( rStrainMatrix, HenckyStrain);
               rStrainVector = ConstitutiveModelUtilities::StrainTensorToVector( HenckyStrain, rStrainVector);

               KRATOS_CATCH("")
            }

            //**************************************************************************************
            //**************************************************************************************
            // divide the deformation gradient in smaller steps
            void ComputeSubstepIncrementalDeformationGradient( const MatrixType & rIncrementalDeformationGradient, const double & rReferenceConfiguration, const double & rFinalConfiguration, MatrixType & rSubstepDeformationGradient)
            {
               KRATOS_TRY

               MatrixType DeformationGradientReference;
               MatrixType DeformationGradientFinal;
               MatrixType IdentityMatrix = identity_matrix<double>(3);

               DeformationGradientReference = rReferenceConfiguration*rIncrementalDeformationGradient  + (1.0 - rReferenceConfiguration)*IdentityMatrix;
               DeformationGradientFinal     =     rFinalConfiguration*rIncrementalDeformationGradient  + (1.0 -     rFinalConfiguration)*IdentityMatrix;

               double det;
               //rSubstepDeformationGradient.clear();
               ConstitutiveModelUtilities::InvertMatrix3( DeformationGradientReference, rSubstepDeformationGradient, det);
               rSubstepDeformationGradient = prod( DeformationGradientFinal, rSubstepDeformationGradient);
               
               KRATOS_CATCH("")
            }
            //***************************************************************************************
            //***************************************************************************************
            // recalculate the elastic left cauchy n
            void RecoverPreviousElasticLeftCauchyGreen( const MatrixType & rIncrementalDeformationGradient, MatrixType & rInitialLeftCauchyGreen)
            {
               MatrixType InverseMatrix; double detMatrix;
               ConstitutiveModelUtilities::InvertMatrix3( rIncrementalDeformationGradient, InverseMatrix, detMatrix);
               rInitialLeftCauchyGreen = prod( InverseMatrix, rInitialLeftCauchyGreen);
               rInitialLeftCauchyGreen = prod( rInitialLeftCauchyGreen, trans(InverseMatrix));
            }


            /**
             * Calculate Stresses
             */
            virtual void SetWorkingMeasures(PlasticDataType& rVariables, MatrixType& rStressMatrix)
            {
               KRATOS_TRY

               const ModelDataType&  rModelData = rVariables.GetModelData();

               //working stress is Kirchhoff by default : transform stresses is working stress is PK2
               const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
               const StrainMeasureType& rStrainMeasure = rModelData.GetStrainMeasure();

               if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){	
                  KRATOS_ERROR << "calling initialize PlasticityModel .. StrainMeasure provided is inconsistent" << std::endl;
               }
               else if(  rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){

                  if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_Left ) {
                     rVariables.StrainMatrix = identity_matrix<double>(3);
                  }
                  else{
                     KRATOS_ERROR << "calling initialize PlasticityModel .. StrainMeasure provided is inconsistent" << std::endl;
                  }

               }
               else{
                  KRATOS_ERROR << "calling initialize PlasticityModel .. StressMeasure provided is inconsistent" << std::endl;
               }


               KRATOS_CATCH(" ")    
            }

            //********************************************************************
            //********************************************************************
            // Initialize Plastic Variables ( a copy from elsewhere:'P)
            virtual void InitializeVariables(ModelDataType& rValues, PlasticDataType& rVariables)
            {
               KRATOS_TRY

               //set model data pointer
               rVariables.SetModelData(rValues);

               rValues.State.Set(ConstitutiveModelData::PLASTIC_REGION,false);

               rValues.State.Set(ConstitutiveModelData::IMPLEX_ACTIVE,false);
               if( rValues.GetProcessInfo()[IMPLEX] == 1 )
                  rValues.State.Set(ConstitutiveModelData::IMPLEX_ACTIVE,true);

               rVariables.SetState(rValues.State);

               // RateFactor
               rVariables.RateFactor = 0;

               // EquivalentPlasticStrain    
               rVariables.Internal = mInternal;

               // DeltaGamma / DeltaPlasticStrain (asociative plasticity)
               rVariables.DeltaInternal.Variables.clear();

               // Flow Rule local variables
               rVariables.TrialStateFunction = 0;
               rVariables.StressNorm = 0;

               KRATOS_CATCH(" ")
            }

            //********************************************************************
            //********************************************************************
            // UpdateInternalVariables
            virtual void UpdateInternalVariables(ModelDataType& rValues, PlasticDataType& rVariables)
            {
               KRATOS_TRY

               double & plasticVolDefNew = rVariables.Internal.Variables[1]; 
               double & plasticVolDef    =mInternal.Variables[1];

               mPreviousInternal.Variables[1] = plasticVolDef;
               plasticVolDef = plasticVolDefNew;

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
            ///@name Serialization
            ///@{
            friend class Serializer;

            virtual void save(Serializer& rSerializer) const override
            {
               KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
                  rSerializer.save("InternalVariables",mInternal);
               rSerializer.save("PreviousInternalVariables",mPreviousInternal);
            }

            virtual void load(Serializer& rSerializer) override
            {
               KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
                  rSerializer.load("InternalVariables",mInternal);
               rSerializer.load("PreviousInternalVariables",mPreviousInternal);
            }

            ///@}
            ///@name Private Inquiry
            ///@{


            ///@}
            ///@name Un accessible methods
            ///@{

            ///@}

      }; // Class NonAssociativePlasticityModel

   ///@}

   ///@name Type Definitions
   ///@{

   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block


} // namespace Kratos

#endif // KRATOS_NON_ASSOCIATIVE_PLASTICITY_MODEL_H_INCLUDED




