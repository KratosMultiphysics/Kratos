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
#include "custom_utilities/stress_invariants_utilities.hpp"

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
      class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonAssociativePlasticityModel : public PlasticityModel<TElasticityModel,TYieldSurface>
      {
         public:

            ///@name Type Definitions
            ///@{

            //elasticity model
            typedef TElasticityModel                               ElasticityModelType;

            //yield surface
            typedef TYieldSurface                                     YieldSurfaceType;

            //base type
            typedef PlasticityModel<ElasticityModelType,YieldSurfaceType>     BaseType;

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
            NonAssociativePlasticityModel(NonAssociativePlasticityModel const& rOther) :BaseType(rOther), mInternal(rOther.mInternal), mPreviousInternal(rOther.mPreviousInternal),
               mStressMatrix(rOther.mStressMatrix), mTotalB(rOther.mTotalB) {}

            /// Assignment operator.
            NonAssociativePlasticityModel& operator=(NonAssociativePlasticityModel const& rOther)
            {
               BaseType::operator=(rOther);
               mInternal = rOther.mInternal;
               mPreviousInternal = rOther.mPreviousInternal;
               mStressMatrix = rOther.mStressMatrix;
               mTotalB = rOther.mTotalB;
               return *this;
            }

            /// Clone.
            ConstitutiveModel::Pointer Clone() const override
            {
               return Kratos::make_shared<NonAssociativePlasticityModel>(*this);
            }

            /// Destructor.
            ~NonAssociativePlasticityModel() override {}


            ///@}
            ///@name Operators
            ///@{


            ///@}
            ///@name Operations
            ///@{

            
            /**
             * Get Values
             */
            virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
            {
               KRATOS_TRY

               // do somehting
               if ( rThisVariable == STRESS_INV_P)
               {
                  double J2;
                  StressInvariantsUtilities::CalculateStressInvariants(mStressMatrix , rValue, J2);
               }
               else if ( rThisVariable == STRESS_INV_J2)
               {
                  double p;
                  StressInvariantsUtilities::CalculateStressInvariants(mStressMatrix , p, rValue);
               }
               else if ( rThisVariable == STRESS_INV_THETA)
               {
                  double p, J2;
                  StressInvariantsUtilities::CalculateStressInvariants(mStressMatrix , p, J2, rValue);
               }
               return rValue;


               KRATOS_CATCH("")
            }

            /**
             * Calculate Stresses
             */

            void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
            {
               KRATOS_TRY


               Matrix ConstitutiveMatrix(6,6);
               noalias(ConstitutiveMatrix) = ZeroMatrix(6,6);
               this->CalculateStressAndConstitutiveTensors( rValues, rStressMatrix, ConstitutiveMatrix);
               rValues.StressMatrix = rStressMatrix;

               KRATOS_CATCH(" ")

            }

            /**
             * Calculate Constitutive Tensor
             */
            void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               //Initialize ConstitutiveMatrix
               rConstitutiveMatrix.clear();

               MatrixType StressMatrix;
               this->CalculateStressAndConstitutiveTensors( rValues, StressMatrix, rConstitutiveMatrix);

               KRATOS_CATCH(" ")
            }

            //*******************************************************************************
            //*******************************************************************************
            // Calculate Stress and constitutive tensor
            void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               double Tolerance = 1e-6;

               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               //Calculate trial stress Matrix
               this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);

               Variables.TrialStateFunction = this->mYieldSurface.CalculateYieldCondition( Variables, Variables.TrialStateFunction);

               Matrix ConstitutiveMatrix(6,6);
               noalias( ConstitutiveMatrix ) = ZeroMatrix(6,6);

               if ( Variables.State().Is(ConstitutiveModelData::IMPLEX_ACTIVE) )
               {
                  const MatrixType & rDeltaDeformationMatrix = rValues.GetDeltaDeformationMatrix();
                  RecoverPreviousElasticLeftCauchyGreen( rDeltaDeformationMatrix, rValues.StrainMatrix );
                  // Calculate with implex
                  this->CalculateImplexPlasticStep(rValues, Variables, rStressMatrix, rDeltaDeformationMatrix);

                  rConstitutiveMatrix.clear();
                  this->mElasticityModel.CalculateConstitutiveTensor(rValues, ConstitutiveMatrix);
                  rConstitutiveMatrix = SetConstitutiveMatrixToTheApropiateSize( rConstitutiveMatrix, ConstitutiveMatrix, rStressMatrix);

               }
               else if ( Variables.TrialStateFunction  < Tolerance) {

                  // elastic loading step
                  rConstitutiveMatrix.clear();
                  this->mElasticityModel.CalculateConstitutiveTensor(rValues, ConstitutiveMatrix);
                  rConstitutiveMatrix = SetConstitutiveMatrixToTheApropiateSize( rConstitutiveMatrix, ConstitutiveMatrix, rStressMatrix);

               } else {

                  // elasto-plastic step. Recover Initial be
                  const MatrixType & rDeltaDeformationMatrix = rValues.GetDeltaDeformationMatrix();
                  RecoverPreviousElasticLeftCauchyGreen( rDeltaDeformationMatrix, rValues.StrainMatrix );

                  double InitialYieldFunction;
                  this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);
                  InitialYieldFunction = this->mYieldSurface.CalculateYieldCondition( Variables, InitialYieldFunction);

                  if ( (InitialYieldFunction < -Tolerance) && (Variables.TrialStateFunction > Tolerance) )
                  {
                     // compute solution with change
                     ComputeSolutionWithChange( rValues, Variables, rDeltaDeformationMatrix);
                  } else {
                     bool UnloadingCondition = false;

                     UnloadingCondition = EvaluateUnloadingCondition( rValues, Variables, rDeltaDeformationMatrix);
                     if (UnloadingCondition) {
                        // compute solution with change
                        ComputeSolutionWithChange( rValues, Variables, rDeltaDeformationMatrix);
                     } else {
                        // compute plastic problem
                        // compute unloading condition
                        ComputeSubsteppingElastoPlasticProblem( rValues, Variables, rDeltaDeformationMatrix);
                     }
                  }

                  this->ReturnStressToYieldSurface( rValues, Variables);

                  noalias( rStressMatrix) = rValues.StressMatrix;

                  this->mElasticityModel.CalculateConstitutiveTensor(rValues, ConstitutiveMatrix);
                  this->ComputeElastoPlasticTangentMatrix( rValues, Variables, ConstitutiveMatrix);
                  rConstitutiveMatrix = SetConstitutiveMatrixToTheApropiateSize( rConstitutiveMatrix, ConstitutiveMatrix, rStressMatrix );
               }

               if ( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES) ) {
                  this->UpdateInternalVariables( rValues, Variables, rStressMatrix );
                  mStressMatrix = rStressMatrix / rValues.GetTotalDeformationDet();
                  const ModelDataType&  rModelData = Variables.GetModelData();
                  const MatrixType & rTotalF = rModelData.GetTotalDeformationMatrix();
                  mTotalB = prod( rTotalF, trans(rTotalF) );
               }

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
            std::string Info() const override
            {
               std::stringstream buffer;
               buffer << "NonAssociativePlasticityModel" ;
               return buffer.str();
            }

            /// Print information about this object.
            void PrintInfo(std::ostream& rOStream) const override
            {
               rOStream << "NonAssociativePlasticityModel";
            }

            /// Print object's data.
            void PrintData(std::ostream& rOStream) const override
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
            MatrixType             mStressMatrix;
            MatrixType             mTotalB;

            ///@}
            ///@name Protected Operators
            ///@{


            ///@}
            ///@name Protected Operations
            ///@{


            //***************************************************************************************
            //***************************************************************************************
            // function to compress the tensor my way
            int AuxiliarCompressTensor( const unsigned int & rI, const unsigned int & rJ , double & rVoigtNumber)
            {

               unsigned int index;
               if ( rI == rJ) {
                  index = rI;
               } else if ( rI > rJ) {
                  index = AuxiliarCompressTensor( rJ, rI, rVoigtNumber);
               } else {
                  rVoigtNumber *= 0.5;
                  if ( rI == 0) {
                     if ( rJ == 1) {
                        index =3;
                     } else {
                        index = 5;
                     }
                  } else {
                     index = 4;
                  }
               }
               return index;

            }

            //***************************************************************************************
            //***************************************************************************************
            // Correct Yield Surface Drift According to
            Matrix & SetConstitutiveMatrixToTheApropiateSize( Matrix & rConstitutiveMatrix, Matrix & rConstMatrixBig, const MatrixType & rStressMatrix)
            {

               KRATOS_TRY

               // 1. Add what I think it is a missing term
               Matrix ExtraMatrix(6,6);
               noalias(ExtraMatrix)= ZeroMatrix(6,6);
               MatrixType Identity;
               noalias(Identity) = identity_matrix<double>(3);

               unsigned int indexi, indexj;
               for (unsigned int i = 0; i < 3; i++) {
                  for (unsigned int j = 0; j < 3; j++) {
                     for (unsigned int k = 0; k < 3; k++) {
                        for (unsigned int l = 0; l < 3; l++) {
                           double voigtNumber = 1.0;
                           indexi = AuxiliarCompressTensor( i, j, voigtNumber);
                           indexj = AuxiliarCompressTensor( k, l, voigtNumber);
                           ExtraMatrix(indexi, indexj) -= voigtNumber * (Identity(i,k)*rStressMatrix(j,l) + Identity(j,k) * rStressMatrix(i,l) );
                        }
                     }
                  }
               }
               rConstMatrixBig += ExtraMatrix;




               // 2. Set the matrix to the appropiate size

               if ( rConstitutiveMatrix.size1() == 6) {
                  noalias( rConstitutiveMatrix ) = rConstMatrixBig;
               } else if ( rConstitutiveMatrix.size1() == 3 ) {
                  rConstitutiveMatrix(0,0) = rConstMatrixBig(0,0);
                  rConstitutiveMatrix(0,1) = rConstMatrixBig(0,1);
                  rConstitutiveMatrix(0,2) = rConstMatrixBig(0,3);

                  rConstitutiveMatrix(1,0) = rConstMatrixBig(1,0);
                  rConstitutiveMatrix(1,1) = rConstMatrixBig(1,1);
                  rConstitutiveMatrix(1,2) = rConstMatrixBig(1,3);

                  rConstitutiveMatrix(2,0) = rConstMatrixBig(3,0);
                  rConstitutiveMatrix(2,1) = rConstMatrixBig(3,1);
                  rConstitutiveMatrix(2,2) = rConstMatrixBig(3,3);

               } else if ( rConstitutiveMatrix.size1() == 4 ) {
                  for (unsigned int i = 0; i < 4; i++) {
                     for (unsigned int j = 0; j < 4; j++) {
                        rConstitutiveMatrix(i,j) = rConstMatrixBig(i,j);
                     }
                  }
               }

               return rConstitutiveMatrix;

               KRATOS_CATCH("")
            }

            //***************************************************************************************
            //***************************************************************************************
            // Correct Yield Surface Drift According to 
            virtual void ReturnStressToYieldSurface( ModelDataType & rValues, PlasticDataType & rVariables)
            {
               KRATOS_TRY

               double Tolerance = 1e-6;

               double YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);

               if ( fabs(YieldSurface) < Tolerance)
                  return;

               for (unsigned int i = 0; i < 150; i++) {

                  Matrix ElasticMatrix(6,6);
                  noalias(ElasticMatrix) = ZeroMatrix(6,6);
                  this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

                  VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
                  VectorType PlasticPotentialDerivative;
                  PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

                  double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H);

                  double DeltaGamma = YieldSurface;
                  DeltaGamma /= ( H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) ) );

                  MatrixType UpdateMatrix;
                  ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);

                  rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
                  rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

                  MatrixType StressMatrix;
                  this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

                  double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
                  for (unsigned int i = 0; i < 3; i++)
                     rPlasticVolDef += DeltaGamma * DeltaStressYieldCondition(i);

                  YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);

                  if ( fabs( YieldSurface) < Tolerance) {
                     return;
                  }
               }

               KRATOS_CATCH("")
            }

            //***************************************************************************************
            //***************************************************************************************
            // Compute Elasto Plastic Matrix
            virtual void ComputeElastoPlasticTangentMatrix( ModelDataType & rValues, PlasticDataType & rVariables, Matrix & rEPMatrix) 
            {
               KRATOS_TRY

      // evaluate constitutive matrix and plastic flow
      Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H);

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

            //***************************************************************************************
            //***************************************************************************************
            // Advance the solution first in elastic regime and then in elastoplastic
            void ComputeSolutionWithChange( ModelDataType& rValues, PlasticDataType & rVariables, const MatrixType& rDeltaDeformationMatrix)
            {
               KRATOS_TRY

               double Tolerance = 1e-6;

               double InitialTime = 0; double EndTime = 1; double HalfTime;
               double InitialStateFunction(-1), EndStateFunction(1), HalfTimeStateFunction;

               MatrixType HalfTimeDeformationGradient;
               MatrixType StressMatrix;
               MatrixType InitialLeftCauchyGreen = rValues.StrainMatrix;
               for (unsigned int i = 0; i < 150; i++)
               {
                  HalfTime = 0.5*(InitialTime + EndTime);

                  ComputeSubstepIncrementalDeformationGradient( rDeltaDeformationMatrix, 0, HalfTime, HalfTimeDeformationGradient);
                  MatrixType AuxMatrix;
                  AuxMatrix = prod( InitialLeftCauchyGreen, trans(HalfTimeDeformationGradient));
                  AuxMatrix = prod( HalfTimeDeformationGradient, AuxMatrix);
                  rValues.StrainMatrix = AuxMatrix;

                  this->mElasticityModel.CalculateStressTensor(rValues,StressMatrix);
                  HalfTimeStateFunction = this->mYieldSurface.CalculateYieldCondition( rVariables, HalfTimeStateFunction);

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
               ComputeSubstepIncrementalDeformationGradient( rDeltaDeformationMatrix, HalfTime, 1,RemainingDeformationGradient);

               ComputeSubsteppingElastoPlasticProblem( rValues, rVariables, RemainingDeformationGradient);

               KRATOS_CATCH("")
            }


            //***********************************************************************************
            //***********************************************************************************
            // Evaluate the elastic unloading condition (Sloan et al, 2001)
            bool  EvaluateUnloadingCondition( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType & rDeltaDeformationMatrix)
            {
               KRATOS_TRY

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;

               MatrixType DeltaStrainMatrix;
               noalias(DeltaStrainMatrix) = prod( rDeltaDeformationMatrix, trans( rDeltaDeformationMatrix) );
               VectorType DeltaStrain;
               ConvertCauchyGreenTensorToHenckyVector( DeltaStrainMatrix, DeltaStrain);

               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType DeltaStress = prod( ElasticMatrix, DeltaStrain);


               double Norm1 = MathUtils<double>::Norm(DeltaStress);
               double Norm2 = MathUtils<double>::Norm(DeltaStressYieldCondition);

               Norm1 = Norm1*Norm2;
               if ( Norm1 < 1e-5)
                  return false;

               Norm2 = MathUtils<double>::Dot( DeltaStressYieldCondition, DeltaStress);

               if (Norm2 > 0) {
                  return false;
               } else {
                  return true;
               }

               KRATOS_CATCH("")
            }

            //***********************************************************************************
            //***********************************************************************************
            // Compute the elasto-plastic problem
            void ComputeSubsteppingElastoPlasticProblem( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType & rDeltaDeformationMatrix)
            {
               KRATOS_TRY

               double Tolerance = 1.0E-6;
               double TimeStep = 0.25;
               double MinTimeStep = 1.0e-4;
               double DoneTimeStep = 0.0;
               double MaxTimeStep = 0.5;

               MatrixType SubstepDeformationGradient;

               double ErrorMeasure;

               while (DoneTimeStep < 1.0) {

                  MatrixType InitialStress = rValues.StrainMatrix;
                  InternalVariablesType InitialInternalVariables = rVariables.Internal;

                  if ( DoneTimeStep + TimeStep >= 1.0) {
                     TimeStep = 1.0 - DoneTimeStep;
                  }

                  ComputeSubstepIncrementalDeformationGradient( rDeltaDeformationMatrix, DoneTimeStep, DoneTimeStep + TimeStep, SubstepDeformationGradient);


                  ErrorMeasure = ComputeElastoPlasticProblem( rValues, rVariables, SubstepDeformationGradient);
                  if ( ErrorMeasure < Tolerance) {
                     DoneTimeStep += TimeStep;
                  } else if ( TimeStep <= MinTimeStep) {
                     std::cout << " ExplicitStressIntegrationDidNotConvege: StressError: " << ErrorMeasure << std::endl;
                     DoneTimeStep += TimeStep;
                  } else {
                     rValues.StrainMatrix = InitialStress;
                     rVariables.Internal = InitialInternalVariables;
                  }

                  TimeStep *= pow( Tolerance / ( ErrorMeasure + 1e-8), 0.5);
                  TimeStep = std::max(TimeStep, MinTimeStep);
                  TimeStep = std::min(TimeStep, MaxTimeStep);

               }


               KRATOS_CATCH("")
            }


            //***********************************************************************************
            //***********************************************************************************
            // Compute one elasto-plastic problem with two discretizations and then compute some
            // sort of error measure
            double ComputeElastoPlasticProblem( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType &  rSubstepDeformationGradient)
            {
               KRATOS_TRY

               MatrixType InitialStrain = rValues.StrainMatrix;
               MatrixType Stress1;
               MatrixType Stress2;
               InternalVariablesType InitialInternalVariables = rVariables.Internal;

               // 1. Compute with one discretization
               this->mElasticityModel.CalculateStressTensor( rValues, Stress1);
               this->ComputeOneStepElastoPlasticProblem( rValues, rVariables, rSubstepDeformationGradient);
               Stress1 = rValues.StressMatrix;

               // 2. Compute with nSteps steps
               unsigned int nSteps = 3;
               rValues.StrainMatrix = InitialStrain;
               rVariables.Internal = InitialInternalVariables;
               this->mElasticityModel.CalculateStressTensor( rValues, Stress2);


               MatrixType IncrementalDefGradient;

               for (unsigned int i = 0; i < nSteps; i++) {
                  double tBegin = double(i)/double(nSteps);
                  double tEnd = double(i+1)/double(nSteps);
                  ComputeSubstepIncrementalDeformationGradient( rSubstepDeformationGradient, tBegin, tEnd, IncrementalDefGradient);
                  this->ComputeOneStepElastoPlasticProblem( rValues, rVariables, IncrementalDefGradient);
               }

               double ErrorMeasure = 0;
               double Denom = 0;
               for (unsigned int i = 0; i < 3; i++) {
                  for (unsigned int j = 0; j < 3; j++) {
                     ErrorMeasure += pow( Stress1(i,j) - rValues.StressMatrix(i,j) , 2);
                     Denom += pow( rValues.StressMatrix(i,j), 2);
                  }
               }

               if ( fabs(Denom) > 1.0E-5)
                  ErrorMeasure /= Denom;
               ErrorMeasure = sqrt(ErrorMeasure);


               return ErrorMeasure;

               KRATOS_CATCH("")
            }


            //***********************************************************************************
            //***********************************************************************************
            // Compute one step of the elasto-plastic problem
            virtual void ComputeOneStepElastoPlasticProblem( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType & rDeltaDeformationMatrix)
            {
               KRATOS_TRY

               MatrixType StressMatrix;
               // evaluate constitutive matrix and plastic flow
               double & rPlasticVolDef = rVariables.Internal.Variables[1];
               double & rPlasticMultiplier = rVariables.Internal.Variables[0];
               double & rPlasticDevDef = rVariables.Internal.Variables[2];

               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H);

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
               for (unsigned int i = 0; i < 3; i++)
                  rPlasticVolDef += DeltaGamma * DeltaStressYieldCondition(i);

               double update = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  update += pow( DeltaGamma * ( DeltaStressYieldCondition(i) - rPlasticVolDef/3.0) , 2.0);
               for (unsigned int i = 3; i < 6; i++)
                  update += 2.0 * pow( DeltaGamma *  DeltaStressYieldCondition(i) /2.0 , 2.0);
               rPlasticDevDef += sqrt(update);


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
               ConstitutiveModelUtilities::StrainTensorToVector( HenckyStrain, rStrainVector);

               KRATOS_CATCH("")
            }

            //**************************************************************************************
            //**************************************************************************************
            // divide the deformation gradient in smaller steps
            void ComputeSubstepIncrementalDeformationGradient( const MatrixType & rDeltaDeformationMatrix, const double & rReferenceConfiguration, const double & rFinalConfiguration, MatrixType & rSubstepDeformationGradient)
            {
               KRATOS_TRY

      MatrixType DeformationGradientReference;
               MatrixType DeformationGradientFinal;
               MatrixType Identity = IdentityMatrix(3);

               DeformationGradientReference = rReferenceConfiguration * rDeltaDeformationMatrix + (1.0 - rReferenceConfiguration) * Identity;
               DeformationGradientFinal     = rFinalConfiguration * rDeltaDeformationMatrix + (1.0 - rFinalConfiguration) * Identity;

               double det;
               rSubstepDeformationGradient.clear();
               ConstitutiveModelUtilities::InvertMatrix3( DeformationGradientReference, rSubstepDeformationGradient, det);
               rSubstepDeformationGradient = prod( DeformationGradientFinal, rSubstepDeformationGradient);

               KRATOS_CATCH("")
            }

            //***************************************************************************************
            //***************************************************************************************
            // recalculate the elastic left cauchy n
            void RecoverPreviousElasticLeftCauchyGreen( const MatrixType & rDeltaDeformationMatrix, MatrixType & rInitialLeftCauchyGreen)
            {
               KRATOS_TRY

      MatrixType InverseMatrix; double detMatrix;
               InverseMatrix.clear();
               ConstitutiveModelUtilities::InvertMatrix3( rDeltaDeformationMatrix, InverseMatrix, detMatrix);
               rInitialLeftCauchyGreen = prod( InverseMatrix, rInitialLeftCauchyGreen);
               rInitialLeftCauchyGreen = prod( rInitialLeftCauchyGreen, trans(InverseMatrix));

               KRATOS_CATCH("")
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

               if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){
                  KRATOS_ERROR << "calling initialize PlasticityModel .. StrainMeasure provided is inconsistent" << std::endl;
               }
               else if(  rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){

                  if( rStrainMeasure == ConstitutiveModelData::StrainMeasureType::CauchyGreen_Left ) {
                     rVariables.StrainMatrix = IdentityMatrix(3);
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
               if( rValues.GetProcessInfo()[IMPLEX] == 1 ) {
                  rValues.State.Set(ConstitutiveModelData::IMPLEX_ACTIVE,true);
               }

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
            virtual void UpdateInternalVariables(ModelDataType& rValues, PlasticDataType& rVariables, const MatrixType& rStressMatrix)
            {
               KRATOS_TRY

      for (unsigned int i = 0; i < 2; i++) {
         double & plasticVolDefNew = rVariables.Internal.Variables[i]; 
         double & plasticVolDef    = mInternal.Variables[i];

         mPreviousInternal.Variables[i] = plasticVolDef;
         plasticVolDef = plasticVolDefNew;
      }

               KRATOS_CATCH("")
            }

            // ****************************************************************************
            //  compute the stress state by using implex
            void  CalculateImplexPlasticStep(ModelDataType& rValues, PlasticDataType&  rVariables, MatrixType&  rStressMatrix, const MatrixType & rDeltaDeformationMatrix)
            {
               KRATOS_TRY

               // evaluate constitutive matrix and plastic flow
               double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
               double & rPlasticDevDef = rVariables.Internal.Variables[2];

               const double & rPlasticMultiplierOld = mPreviousInternal.Variables[0];
               double & rPlasticMultiplier    = rVariables.Internal.Variables[0];
               double  DeltaPlasticMultiplier = (rPlasticMultiplier - rPlasticMultiplierOld);

               if ( DeltaPlasticMultiplier < 0)
                  DeltaPlasticMultiplier = 0;

               
               this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV


               MatrixType UpdateMatrix;
               ConvertHenckyVectorToCauchyGreenTensor( -DeltaPlasticMultiplier * PlasticPotentialDerivative / 2.0, UpdateMatrix);
               UpdateMatrix = prod( rDeltaDeformationMatrix, UpdateMatrix);


               rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
               rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

               this->mElasticityModel.CalculateStressTensor( rValues, rStressMatrix);

               rPlasticMultiplier += DeltaPlasticMultiplier;
               for (unsigned int i = 0; i < 3; i++)
                  rPlasticVolDef += DeltaPlasticMultiplier * DeltaStressYieldCondition(i);

               double update = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  update += pow( DeltaPlasticMultiplier * ( DeltaStressYieldCondition(i) - rPlasticVolDef/3.0) , 2.0);
               for (unsigned int i = 3; i < 6; i++)
                  update += 2.0 * pow( DeltaPlasticMultiplier *  DeltaStressYieldCondition(i) /2.0 , 2.0);
               rPlasticDevDef += sqrt(update);

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

            void save(Serializer& rSerializer) const override
            {
               KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
                  rSerializer.save("InternalVariables",mInternal);
               rSerializer.save("PreviousInternalVariables",mPreviousInternal);
            }

            void load(Serializer& rSerializer) override
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
