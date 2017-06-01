//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONSTITUTIVE_MODEL_DATA_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_MODEL_DATA_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"

#include "constitutive_models_application_variables.h"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ConstitutiveModelData
  {
  public:

    ///@name Type Definitions
    ///@{
    using VoigtIndexType = const unsigned int(*)[2];
    using SizeType       = std::size_t;
    using MatrixType     = bounded_matrix<double,3,3>;
    using VectorType     = array_1d<double,6>;
    
    //state flags
    KRATOS_DEFINE_LOCAL_FLAG( IMPLEX_ACTIVE );    
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTED_STRAIN );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTED_STRESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTED_CONSTITUTIVE_MATRIX );    
    KRATOS_DEFINE_LOCAL_FLAG( PLASTIC_REGION );
    KRATOS_DEFINE_LOCAL_FLAG( PLASTIC_RATE_REGION );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTED_RETURN_MAPPING );
    KRATOS_DEFINE_LOCAL_FLAG( UPDATE_INTERNAL_VARIABLES );
    
    enum StrainMeasureType  //supplied cauchy green strain measure
    {
      CauchyGreen_None,            //no strain measure supplied
      CauchyGreen_Left,            //left cauchy-green tensor
      CauchyGreen_Right,           //right cauchy-green tensor
    };

    
    enum StressMeasureType  //required stress measure
    {
      StressMeasure_PK1,            //stress related to reference configuration non-symmetric
      StressMeasure_PK2,            //stress related to reference configuration
      StressMeasure_Kirchhoff,      //stress related to current   configuration
      StressMeasure_Cauchy          //stress related to current   configuration
    };


    struct MaterialData
    {      
    public:
      
      //general elastic material properties 
      double PoissonCoefficient;
      double YoungModulus;
      double LameMu;
      double LameMuBar;
      double LameLambda;
      double BulkModulus;

      //general hyperelastic material properties
      std::vector<double> ModelParameters;

    public:
      
      //Get const Data
      const double& GetPoissonCoefficient() const {return PoissonCoefficient;};
      const double& GetYoungModulus      () const {return YoungModulus;};
      const double& GetLameMu            () const {return LameMu;};
      const double& GetLameMuBar         () const {return LameMuBar;};   
      const double& GetLameLambda        () const {return LameLambda;};
      const double& GetBulkModulus       () const {return BulkModulus;};
      
      const std::vector<double>& GetModelParameters   () const {return ModelParameters;};

    };


    struct ConstitutiveLawData
    {
    public:

      //elemental data
      double                       Pressure;
      double                       Temperature;
      double                       CharacteristicSize;      
      double                       DeterminantF;
      double                       DeterminantF0;
      
      //model data
      StressMeasureType            StressMeasure;       //stress measure requested
      StrainMeasureType            StrainMeasure;       //strain measure provided

      //deformation
      MatrixType                   DeformationGradientF;   //incremental F
      MatrixType                   DeformationGradientF0;  //total F
    };
    
    
    struct ModelData
    {
    private:

      const Flags*                 mpOptions;
      
      const Properties*            mpMaterialProperties;
      const ProcessInfo*           mpProcessInfo;

      SizeType                     mVoigtSize;
      VoigtIndexType               mIndexVoigtTensor;
      
      ConstitutiveLawData          mConstitutiveLawData;
      
    public:
      
      Flags                        State;
      
      MatrixType                   StressMatrix;           //wildcard stress (isochoric stress tensor)     
      MatrixType                   StrainMatrix;           //wildcard strain (cauchy green tensors or infinitessimal tensor)
      MaterialData                 MaterialParameters;

       
      //Set Data Pointers
      void SetOptions                      (const Flags&  rOptions)                 {mpOptions = &rOptions;};
      void SetMaterialProperties           (const Properties&  rMaterialProperties) {mpMaterialProperties = &rMaterialProperties;};
      void SetProcessInfo                  (const ProcessInfo& rProcessInfo)        {mpProcessInfo = &rProcessInfo;};
      void SetVoigtSize                    (const SizeType& rVoigtSize)             {mVoigtSize = rVoigtSize;};
      void SetVoigtIndexTensor             (VoigtIndexType rIndexVoigtTensor)       {mIndexVoigtTensor = rIndexVoigtTensor;};

      void  SetStressMeasure               (StressMeasureType Measure)              {mConstitutiveLawData.StressMeasure = Measure;}; 
      void  SetStrainMeasure               (StrainMeasureType Measure)              {mConstitutiveLawData.StrainMeasure = Measure;};

      //Get Data Pointers
      const Flags&          GetOptions                     () const {return *mpOptions;};
      const Properties&     GetMaterialProperties          () const {return *mpMaterialProperties;};
      const ProcessInfo&    GetProcessInfo                 () const {return *mpProcessInfo;};
      const SizeType&       GetVoigtSize                   () const {return  mVoigtSize;};
      const VoigtIndexType& GetVoigtIndexTensor            () const {return  mIndexVoigtTensor;};
      
      //Acces non const Data
      ConstitutiveLawData& rConstitutiveLawData            () {return mConstitutiveLawData;};
      MatrixType&          rStrainMatrix                   () {return StrainMatrix;};
      MatrixType&          rStressMatrix                   () {return StressMatrix;}; 
      MaterialData&        rMaterialParameters             () {return MaterialParameters;}; 

      //Get const Data
      const double&        GetPressure                     () const {return mConstitutiveLawData.Pressure;}; 
      const double&        GetTemperature                  () const {return mConstitutiveLawData.Temperature;}; 
      const double&        GetDeterminantF                 () const {return mConstitutiveLawData.DeterminantF;};
      const double&        GetDeterminantF0                () const {return mConstitutiveLawData.DeterminantF0;}; 
      const double&        GetCharacteristicSize           () const {return mConstitutiveLawData.CharacteristicSize;}; 

      const StressMeasureType& GetStressMeasure            () const {return mConstitutiveLawData.StressMeasure;}; 
      const StrainMeasureType& GetStrainMeasure            () const {return mConstitutiveLawData.StrainMeasure;};

      const MatrixType&    GetDeformationGradientF         () const {return mConstitutiveLawData.DeformationGradientF;}; 
      const MatrixType&    GetDeformationGradientF0        () const {return mConstitutiveLawData.DeformationGradientF0;};

      const ConstitutiveLawData&   GetConstitutiveLawData  () const {return mConstitutiveLawData;};
      
      const MatrixType&    GetStrainMatrix                 () const {return StrainMatrix;};
      const MatrixType&    GetStressMatrix                 () const {return StressMatrix;}; 
      const MaterialData&  GetMaterialParameters           () const {return MaterialParameters;}; 
      
    };


   
    // struct ThermalParameters
    // {
    //   //general thermal properties
    //   double ThermalExpansionCoefficient;
    //   double ReferenceTemperature;      
    // }
    
    
    /// Pointer definition of ConstitutiveModelData
    KRATOS_CLASS_POINTER_DEFINITION( ConstitutiveModelData );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.    
    ConstitutiveModelData(){}

    /// Copy constructor.
    ConstitutiveModelData(ConstitutiveModelData const& rOther){}

    /// Clone.
    ConstitutiveModelData::Pointer Clone() const
    {
      return (ConstitutiveModelData::Pointer(new ConstitutiveModelData(*this)));
    }

    /// Destructor.
    virtual ~ConstitutiveModelData(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    
    ///@}
    ///@name Access
    ///@{
    
    static inline void CalculateMaterialParameters(ModelData& rValues)
    {
      KRATOS_TRY
          
      //material properties
      const Properties& rProperties = rValues.GetMaterialProperties();
      ConstitutiveLawData& rConstitutiveLawData = rValues.rConstitutiveLawData(); 
      
      //if previously computed LameMu / LameLambda / BulkModulus
      // rValues.MaterialParameters.LameMu      = rProperties[LAME_MU];
      // rValues.MaterialParameters.LameLambda  = rProperties[LAME_LAMBDA];
      // rValues.MaterialParameters.BulkModulus = rProperties[BULK_MODULUS];
      
      // compute material properties

      // const double& YoungModulus       = rProperties[YOUNG_MODULUS];
      // const double& PoissonCoefficient = rProperties[POISSON_RATIO];
    
      // temperature dependent parameters:

      if( rProperties.HasTable(TEMPERATURE,YOUNG_MODULUS) ){
	const Table<double>& YoungModulusTable = rProperties.GetTable(TEMPERATURE,YOUNG_MODULUS);
	rValues.MaterialParameters.YoungModulus = YoungModulusTable[rConstitutiveLawData.Temperature];
      }
      else{
	rValues.MaterialParameters.YoungModulus = rProperties[YOUNG_MODULUS];
      }
    
      if( rProperties.HasTable(TEMPERATURE,POISSON_RATIO) ){
	const Table<double>& PoissonCoefficientTable = rProperties.GetTable(TEMPERATURE,POISSON_RATIO);
	rValues.MaterialParameters.PoissonCoefficient = PoissonCoefficientTable[rConstitutiveLawData.Temperature];
      }
      else{
	rValues.MaterialParameters.PoissonCoefficient = rProperties[POISSON_RATIO];
      }
           
      rValues.MaterialParameters.LameMu        = rValues.MaterialParameters.YoungModulus/(2.0*(1.0+rValues.MaterialParameters.PoissonCoefficient));
      rValues.MaterialParameters.LameLambda    = (rValues.MaterialParameters.YoungModulus*rValues.MaterialParameters.PoissonCoefficient)/((1.0+rValues.MaterialParameters.PoissonCoefficient)*(1.0-2.0*rValues.MaterialParameters.PoissonCoefficient));
      rValues.MaterialParameters.BulkModulus   = rValues.MaterialParameters.LameLambda + (2.0/3.0) * rValues.MaterialParameters.LameMu;


      //hyperelastic model parameters
      if( rProperties.Has(C10) ){
	rValues.MaterialParameters.ModelParameters.push_back(rProperties[C10]);

	//make neo-hookean consistent with the parameters:
	rValues.MaterialParameters.LameMu = 2.0 * rProperties[C10];
	rValues.MaterialParameters.BulkModulus   = rValues.MaterialParameters.LameLambda + (2.0/3.0) * rValues.MaterialParameters.LameMu;
      }

      if( rProperties.Has(C20) )
	rValues.MaterialParameters.ModelParameters.push_back(rProperties[C20]);

      if( rProperties.Has(C30) )
	rValues.MaterialParameters.ModelParameters.push_back(rProperties[C30]);
	        
      KRATOS_CATCH(" ")

    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ConstitutiveModelData";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConstitutiveModelData";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "ConstitutiveModelData Data";
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


    ///@}
    ///@name Protected Operators
    ///@{
    
    
    ///@}
    ///@name Protected Operations
    ///@{

        
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

    /// Assignment operator.
    ConstitutiveModelData& operator=(ConstitutiveModelData const& rOther);

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
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class ConstitutiveModelData

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                                    ConstitutiveModelData& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                                    const ConstitutiveModelData& rThis)
  {
      rThis.PrintInfo(rOStream);
      rOStream <<" : " << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
  }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSTITUTIVE_MODEL_DATA_H_INCLUDED  defined 


