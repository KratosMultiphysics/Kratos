//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CAM_CLAY_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_CAM_CLAY_HARDENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_laws/hardening_law.hpp"

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CamClayHardeningLaw 
      : public HardeningLaw
      {
         protected:

            constexpr static std::size_t VarSize = 4;
         public:
            ///@name Type Definitions
            ///@{
    typedef ConstitutiveModelData::MatrixType                  MatrixType;
    typedef ConstitutiveModelData::VectorType                  VectorType;
    typedef ConstitutiveModelData::ModelData                ModelDataType;
    typedef ConstitutiveModelData::MaterialData          MaterialDataType;
    
    template<std::size_t TVarSize>
    struct InternalVariables
    {
      //internal variables
      array_1d<double, TVarSize> Variables;

      //default constructor (to initialize Variables)
      InternalVariables() { Variables.clear(); };

      const array_1d<double, TVarSize>& GetVariables() {return Variables;};

    private:

      friend class Serializer;
     
      void save(Serializer& rSerializer) const
      {
	rSerializer.save("Variables",Variables);
      };

      void load(Serializer& rSerializer)
      {
	rSerializer.load("Variables",Variables);
      };
      
    };
    

    template<std::size_t TVarSize>
    struct PlasticModelData
    {
    private:

      Flags*               mpState;
      const ModelDataType* mpModelData;
      
    public:

      //flow rule internal variables     
      double TrialStateFunction;
      double StressNorm;
      
      //hardening law internal variables
      double RateFactor;

      //internal variables
      InternalVariables<TVarSize>      Internal;
      InternalVariables<TVarSize> DeltaInternal;

      //strain matrix
      MatrixType StrainMatrix; //wildcard strain (cauchy green tensors or infinitessimal tensor)
      
      //Set Data Pointers
      void SetState           (Flags& rState)                    {mpState = &rState;};
      void SetModelData       (const ModelDataType&  rModelData) {mpModelData = &rModelData;};
      
      //Get Data Pointers
      const ModelDataType&    GetModelData                () const {return *mpModelData;};
      const MaterialDataType& GetMaterialParameters       () const {return mpModelData->GetMaterialParameters();};

      //Get non const Data
      Flags& State                                        () {return *mpState;};

      //Get const Data
      const Flags&  GetState              () const {return *mpState;};
      const double& GetTrialStateFunction () const {return TrialStateFunction;};
      const double& GetStressNorm         () const {return StressNorm;};     
      const double& GetRateFactor         () const {return RateFactor;};
      
      const InternalVariables<TVarSize>&      GetInternal () const {return Internal;};
      const InternalVariables<TVarSize>& GetDeltaInternal () const {return DeltaInternal;};
      
      const MatrixType&                   GetStrainMatrix () const {return StrainMatrix;};
      
      const array_1d<double,TVarSize>& GetInternalVariables       () const {return Internal.Variables;};
      const array_1d<double,TVarSize>& GetDeltaInternalVariables  () const {return DeltaInternal.Variables;};
     
    };

    typedef InternalVariables<VarSize>   InternalVariablesType;
    typedef PlasticModelData<VarSize>          PlasticDataType;
    

            /// Pointer definition of CamClayHardeningLaw
            KRATOS_CLASS_POINTER_DEFINITION( CamClayHardeningLaw );

            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            CamClayHardeningLaw();

            /// Copy constructor.
            CamClayHardeningLaw(CamClayHardeningLaw const& rOther);

            /// Assignment operator.
            CamClayHardeningLaw& operator=(CamClayHardeningLaw const& rOther);

            /// Clone.
            virtual HardeningLaw::Pointer Clone() const override;

            /// Destructor.
            ~CamClayHardeningLaw();

            ///@}
            ///@name Operators
            ///@{


            ///@}
            ///@name Operations
            ///@{


            /**
             * Calculate Hardening functions
             */

            virtual double& CalculateHardening(const PlasticDataType& rVariables, double& rHardening);

            /**
             * Calculate Hardening function derivatives
             */

            virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening);


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
               buffer << "CamClayHardeningLaw" ;
               return buffer.str();
            }

            /// Print information about this object.
            virtual void PrintInfo(std::ostream& rOStream) const override
            {
               rOStream << "CamClayHardeningLaw";
            }

            /// Print object's data.
            virtual void PrintData(std::ostream& rOStream) const override
            {
               rOStream << "CamClayHardeningLaw Data";
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


            /**
             * Pure isotropic hardening Theta=1;  pure kinematic hardening Theta= 0; combined isotropic-kinematic 0<Theta<1
             */
            static constexpr double mTheta = 1.0;


            ///@}
            ///@name Protected Operators
            ///@{

            ///@}
            ///@name Protected Operations
            ///@{

            /**
             * Calculate Hardening functions
             */
            virtual double& CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double& rIsotropicHardening);

            virtual double& CalculateAndAddKinematicHardening(const PlasticDataType& rVariables, double& rKinematicHardening);

            /**
             * Calculate Hardening function derivatives
             */
            virtual double& CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double& rDeltaIsotropicHardening);

            virtual double& CalculateAndAddDeltaKinematicHardening(const PlasticDataType& rVariables, double& rDeltaKinematicHardening);



            virtual double& CalculateThermalReferenceEffect(const PlasticDataType& rVariables, double& rThermalFactor);

            virtual double& CalculateThermalCurrentEffect(const PlasticDataType& rVariables, double& rThermalFactor);


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
               KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
            }

            virtual void load(Serializer& rSerializer) override
            {
               KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
            }

            ///@}
            ///@name Private Inquiry
            ///@{


            ///@}
            ///@name Un accessible methods
            ///@{


            ///@}

      }; // Class CamClayHardeningLaw

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CAM_CLAY_HARDENING_LAW_H_INCLUDED  defined 


