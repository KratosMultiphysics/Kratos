//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:          October 2017 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/constitutive_model.hpp"

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
   class KRATOS_API(UMAT_APPLICATION) SmallStrainUmatModel : public ConstitutiveModel
   {
      protected:

         struct UmatModelData
         {
            private:

               Flags*               mpState;
               const ModelDataType* mpModelData;

            public:

               MatrixType IncrementalDeformation;
               MatrixType TotalStrainMatrix;

               //Set Data Pointers
               void SetState           (Flags& rState)                    {mpState = &rState;};
               void SetModelData       (const ModelDataType&  rModelData) {mpModelData = &rModelData;};

               //Get Data Pointers
               const ModelDataType&    GetModelData                () const {return *mpModelData;};
               const MaterialDataType& GetMaterialParameters       () const {return mpModelData->GetMaterialParameters();}; 

               //Get non const Data
               Flags& State                                        () {return *mpState;};

               //Get const Data
               const Flags&  GetState                              () const {return *mpState;};

         };


      public:

         ///@name Type Definitions
         ///@{
         typedef UmatModelData              UmatDataType;


         /// Pointer definition of SmallStrainUmatModel
         KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUmatModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.    
         SmallStrainUmatModel();

         /// Copy constructor.
         SmallStrainUmatModel(SmallStrainUmatModel const& rOther);

         /// Clone.
         virtual ConstitutiveModel::Pointer Clone() const override;

         /// Assignment operator.
         SmallStrainUmatModel& operator=(SmallStrainUmatModel const& rOther);

         /// Destructor.
         virtual ~SmallStrainUmatModel();


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{


         /**
          * Initialize member data
          */    
         virtual void InitializeModel(ModelDataType& rValues) override;


         /**
          * Finalize member data
          */    
         virtual void FinalizeModel(ModelDataType& rValues) override;


         /**
          * Calculate Strain Energy Density Functions
          */
         virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override;


         /**
          * Calculate Stresses
          */    
         virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;




         /**
          * Calculate Constitutive Tensor
          */
         virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override; 


         /**
          * Calculate Stress and Constitutive Tensor
          */
         virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;


         /**
          * Check
          */    
         virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override;

         ///@}
         ///@name Access
         ///@{

         virtual void SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
               const ProcessInfo& rCurrentProcessInfo ) override
         {
            KRATOS_TRY

      // A method to compute the initial linear strain from the stress is needed
      //if(rThisVariable == INITIAL_STRESS_VECTOR)

      // A method to compute the initial linear strain from the stress is needed
      // if(rThisVariable == INITIAL_STRAIN_VECTOR){
      //   this->mHistoryVector = rValue;
      // }

      KRATOS_CATCH(" ")
         }


         virtual void SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
               const ProcessInfo& rCurrentProcessInfo ) override
         {
            KRATOS_TRY

      // A method to compute the initial linear strain from the stress is needed
      //if(rThisVariable == INITIAL_STRESS_VECTOR)

      // A method to compute the initial linear strain from the stress is needed
      // if(rThisVariable == INITIAL_STRAIN_VECTOR){
      //   this->mHistoryVector = rValue;
      // }

      KRATOS_CATCH(" ")
         }

         /**
          * method to ask the plasticity model the list of variables (dofs)  needed from the domain
          * @param rScalarVariables : list of scalar dofs
          * @param rComponentVariables :  list of vector dofs
          */
         virtual void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
               std::vector<Variable<array_1d<double,3> > >& rComponentVariables) override
         {
            KRATOS_TRY

      rComponentVariables.push_back(DISPLACEMENT);

            KRATOS_CATCH(" ")
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
            buffer << "SmallStrainUmatModel";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "SmallStrainUmatModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "SmallStrainUmatModel Data";
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

         bool mInitializedModel;

         Vector     mStateVariablesFinalized;
         VectorType mStressVectorFinalized;
         VectorType mStrainVectorFinalized;


         ///@}
         ///@name Protected Operators
         ///@{


         ///@}
         ///@name Protected Operations
         ///@{


         //************//

         void InitializeElasticData(ModelDataType& rValues, UmatDataType& rVariables);   

         /*
            Get the dimension of StateVariables 
          */

         virtual unsigned int GetNumberOfStateVariables() {
            KRATOS_ERROR << " Calling the base GetNumberOfStateVariables " << std::endl;
         };

         /*
            Create the vector with constitutive parameters value
          */
         virtual void CreateConstitutiveParametersVector(double* & rpVector, int & rNumberParameters, const Properties & rMaterialProperties) {
            KRATOS_ERROR << " Calling the base CreateConstitutiveParametersVector " << std::endl;
         };

         /*
            Create state variables vector 
            ( tensor like should be "rotated" in large strains )
          */
         virtual void CreateStateVariablesVector(double * & rpStateVariables,int & rNumberStateVariables)
         {
            rpStateVariables = new double[rNumberStateVariables];
            for (int i = 0; i < rNumberStateVariables; i++)
            {
               rpStateVariables[i] = mStateVariablesFinalized(i);
            }
         };



         /*
            Create strain_n and incremental strain
            ( for large strains should be overrided ) 
          */

         virtual void CreateStrainsVectors( UmatDataType & rVariables, double* & rpStrain, double* & rpIncrementalStrain)
         {

            VectorType StrainVector;
            ConstitutiveModelUtilities::StrainTensorToVector( rVariables.TotalStrainMatrix, StrainVector);

            rpStrain = new double[6];
            rpIncrementalStrain = new double[6];

            for (unsigned int i = 0; i < 6; i++) {
               rpStrain[i] = mStrainVectorFinalized(i);
               rpIncrementalStrain[i] = ( StrainVector(i) - mStrainVectorFinalized(i) );

            }

         }

         /* 
            Create stress_n
            ( for large strains should be overrided )
          */
         virtual void CreateStressAtInitialState( UmatDataType & rVariables, double* & rpStressVector)
         {

            rpStressVector = new double[6];
            for (unsigned int i = 0; i < 6; i++) {
               rpStressVector[i] = mStressVectorFinalized(i);
            }

         }

         /* 
            Update constitutive model variables
          */
         virtual void UpdateVariables( UmatDataType & rVariables, double* & rpStressVector, double* & rpStateVariables, double Pressure = 0.0)
         {
            VectorType StrainVector;
            ConstitutiveModelUtilities::StrainTensorToVector( rVariables.TotalStrainMatrix, StrainVector);
            mStrainVectorFinalized = StrainVector;


            if ( fabs( Pressure) < 1e-5) {
               for (unsigned int i = 0; i < 6; i++) 
                  mStressVectorFinalized(i) = rpStressVector[i];
            } else {

               double meanP = 0;
               for (unsigned int i = 0; i < 3; i++) 
                  meanP += rpStressVector[i];
               meanP /= 3.0;

               for (unsigned int i = 0; i < 3; i++) 
                  mStressVectorFinalized(i) = rpStressVector[i] + Pressure - meanP;
               for (unsigned int i = 3; i < 6; i++) 
                  mStressVectorFinalized(i) = rpStressVector[i];
            }

            int nStateVariables = this->GetNumberOfStateVariables();
            for (int i = 0; i < nStateVariables; i++)
               mStateVariablesFinalized(i) = rpStateVariables[i];

         }

         /*
            Number of the constitutive equation in the fortran wrapper
          */
         virtual int GetConstitutiveEquationNumber()
         {
            KRATOS_ERROR << " Calling the base case of GetConstitutiveEquationNumber " << std::endl;
            return 0;
         }

         virtual void SetConstitutiveMatrix( Matrix & rC, const Matrix & rCBig, const MatrixType& rStressMatrix);
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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveModel )
            rSerializer.save("InitializedModel", mInitializedModel );
            rSerializer.save("StressVectorFinalized", mStressVectorFinalized );
            rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized );
            rSerializer.save("StateVariablesFinalized", mStateVariablesFinalized );
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveModel )
            rSerializer.load("InitializedModel", mInitializedModel );
            rSerializer.load("StressVectorFinalized", mStressVectorFinalized );
            rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized );
            rSerializer.load("StateVariablesFinalized", mStateVariablesFinalized );
         }

         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}

   }; // Class SmallStrainUmatModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED  defined 


