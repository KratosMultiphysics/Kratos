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
    using MatrixType     = BoundedMatrix<double,3,3>;
    using VectorType     = array_1d<double,6>;

    //state flags
    KRATOS_DEFINE_LOCAL_FLAG( IMPLEX_ACTIVE );
    KRATOS_DEFINE_LOCAL_FLAG( STRAIN_COMPUTED );
    KRATOS_DEFINE_LOCAL_FLAG( STRESS_COMPUTED );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTITUTIVE_MATRIX_COMPUTED );
    KRATOS_DEFINE_LOCAL_FLAG( PLASTIC_REGION );
    KRATOS_DEFINE_LOCAL_FLAG( PLASTIC_RATE_REGION );
    KRATOS_DEFINE_LOCAL_FLAG( RETURN_MAPPING_COMPUTED );
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


    template <class T>
    struct VariableValue
    {
    private:

      const Variable<T> *mpVariable;
      T *mpValue;

    public:

      //constructors
      VariableValue(){};
      VariableValue(const Variable<T>& rVariable, T& rValue){mpVariable = &rVariable; mpValue = &rValue;};

      //Set Data Pointers
      void SetVariableValue            (const Variable<T>& rVariable, T& rValue) {mpVariable = &rVariable; mpValue = &rValue;};
      void SetVariable                 (const Variable<T>& rVariable) {mpVariable = &rVariable;};

      //Get-Set Data
      bool HasVariable                 (const Variable<T>& rVariable) const {return (rVariable == *mpVariable);};

      void SetValue                    (T& rValue) {*mpValue = rValue;};

    };


    struct VariableValueData
    {
    private:

      enum VarType { INTEGER, DOUBLE, VECTOR, MATRIX, ARRAY3, ARRAY6, NONE };

      VarType                                        mType;
      VariableValue<int>                    *mpIntVariable;
      VariableValue<double>              *mpDoubleVariable;
      VariableValue<Vector>              *mpVectorVariable;
      VariableValue<Matrix>              *mpMatrixVariable;
      VariableValue<array_1d<double,3> > *mpArray3Variable;
      VariableValue<array_1d<double,6> > *mpArray6Variable;

   public:

      //Constructor
      VariableValueData()
      {
	mType            = NONE;
	mpIntVariable    = nullptr;
	mpDoubleVariable = nullptr;
	mpVectorVariable = nullptr;
	mpMatrixVariable = nullptr;
	mpArray3Variable = nullptr;
	mpArray6Variable = nullptr;
      }

      // Destructor
      ~VariableValueData()
      {
	switch(mType)
	  {
	  case INTEGER:
	    delete mpIntVariable;
	    break;
	  case DOUBLE:
	    delete mpDoubleVariable;
	    break;
	  case VECTOR:
	    delete mpVectorVariable;
	    break;
	  case MATRIX:
	    delete mpMatrixVariable;
	    break;
	  case ARRAY3:
	    delete mpArray3Variable;
	    break;
	  case ARRAY6:
	    delete mpArray6Variable;
	    break;
	  case NONE:
	    break;
	  default:
	    break;
	  }
      }

      //Set Data

      void SetIntVariableValue(const Variable<int>& rVariable, int& rValue)
      {
	mType = INTEGER;
	typedef VariableValue<int> VariableValueType;
	mpIntVariable = new VariableValueType(rVariable,rValue);
      }

      void SetDoubleVariableValue(const Variable<double>& rVariable, double& rValue)
      {
	mType = DOUBLE;
	typedef VariableValue<double> VariableValueType;
	mpDoubleVariable = new VariableValueType(rVariable,rValue);
      }

      void SetVectorVariableValue(const Variable<Vector>& rVariable, Vector& rValue)
      {
	mType = VECTOR;
	typedef VariableValue<Vector> VariableValueType;
	mpVectorVariable = new VariableValueType(rVariable,rValue);
      }

      void SetMatrixVariableValue(const Variable<Matrix>& rVariable, Matrix& rValue)
      {
	mType = MATRIX;
	typedef VariableValue<Matrix> VariableValueType;
	mpMatrixVariable = new VariableValueType(rVariable,rValue);
      }

      void SetArray3VariableValue(const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& rValue)
      {
	mType = ARRAY3;
	typedef VariableValue<array_1d<double,3> > VariableValueType;
	mpArray3Variable = new VariableValueType(rVariable,rValue);
      }

      void SetArray6VariableValue(const Variable<array_1d<double,6> >& rVariable, array_1d<double,6>& rValue)
      {
	mType = ARRAY6;
	typedef VariableValue<array_1d<double,6> > VariableValueType;
	mpArray6Variable = new VariableValueType(rVariable,rValue);
      }

      //Get Data
      template<class T>
      bool GetVariableValue(VariableValue<T>& rVariableValue)
      {
	if( std::is_same<T,int>::value ){
	  if(mpIntVariable != nullptr){
	    rVariableValue = &mpIntVariable;
	    return true;
	  }
	  else{
	    return false;
	  }

	}
	else if( std::is_same<T,double>::value ){
	  if(mpDoubleVariable != nullptr){
	    rVariableValue = &mpDoubleVariable;
	    return true;
	  }
	  else{
	    return false;
	  }
	}
	else if( std::is_same<T,Vector>::value ){
	  if(mpVectorVariable != nullptr){
	    rVariableValue = &mpVectorVariable;
	    return true;
	  }
	  else{
	    return false;
	  }
	}
	else if( std::is_same<T,Matrix>::value ){
	  if(mpMatrixVariable != nullptr){
	    rVariableValue = &mpMatrixVariable;
	    return true;
	  }
	  else{
	    return false;
	  }
	}
	else if( std::is_same<T,array_1d<double,3> >::value ){
	  if(mpArray3Variable != nullptr){
	    rVariableValue = &mpArray3Variable;
	    return true;
	  }
	  else{
	    return false;
	  }
	}
	else if( std::is_same<T,array_1d<double,6> >::value ){
	  if(mpArray6Variable !=nullptr){
	    rVariableValue = &mpArray6Variable;
	    return true;
	  }
	  else{
	    return false;
	  }
	}
	else{
	  return false;
	}
      }

      void SetValue(const Variable<int>& rVariable, int& rValue)
      {
	if(mpIntVariable != nullptr)
	  if( mpIntVariable->HasVariable(rVariable) )
	    mpIntVariable->SetValue(rValue);
      }

      void SetValue(const Variable<double>& rVariable, double& rValue)
      {
	if(mpDoubleVariable != nullptr)
	  if( mpDoubleVariable->HasVariable(rVariable) )
	    mpDoubleVariable->SetValue(rValue);
      }

      void SetValue(const Variable<Vector>& rVariable, Vector& rValue)
      {
	if(mpVectorVariable != nullptr)
	  if( mpVectorVariable->HasVariable(rVariable) )
	    mpVectorVariable->SetValue(rValue);
      }

      void SetValue(const Variable<Matrix>& rVariable, Matrix& rValue)
      {
	if(mpMatrixVariable != nullptr)
	  if( mpMatrixVariable->HasVariable(rVariable) )
	    mpMatrixVariable->SetValue(rValue);

      }

      void SetValue(const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& rValue)
      {
	if(mpArray3Variable != nullptr)
	  if( mpArray3Variable->HasVariable(rVariable) )
	    mpArray3Variable->SetValue(rValue);
      }

      void SetValue(const Variable<array_1d<double,6> >& rVariable, array_1d<double,6>& rValue)
      {
	if(mpArray6Variable != nullptr)
	  if( mpArray6Variable->HasVariable(rVariable) )
	    mpArray6Variable->SetValue(rValue);
      }



    };


    struct ConstitutiveLawData
    {
    public:

      //elemental data
      double                       Pressure;
      double                       Temperature;
      double                       CharacteristicSize;
      double                       DeltaDeformationDet;  //wildcard for the determinant of the deformation increment (usually detF )
      double                       TotalDeformationDet;  //wildcard for the determinant of the total deformation     (usually detF0)

      //model data
      StressMeasureType            StressMeasure;       //stress measure requested
      StrainMeasureType            StrainMeasure;       //strain measure provided

      //deformation
      MatrixType                   DeltaDeformationMatrix;  //wildcard deformation increment (usually incremental F)
      MatrixType                   TotalDeformationMatrix;  //wildcard total deformation     (usually total F := F0)
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

      VariableValueData            InternalVariable;       //internal variable to compute and return

      //Set Data Pointers
      void SetOptions                      (const Flags&  rOptions)                 {mpOptions = &rOptions;};
      void SetMaterialProperties           (const Properties&  rMaterialProperties) {mpMaterialProperties = &rMaterialProperties;};
      void SetProcessInfo                  (const ProcessInfo& rProcessInfo)        {mpProcessInfo = &rProcessInfo;};
      void SetVoigtSize                    (const SizeType& rVoigtSize)             {mVoigtSize = rVoigtSize;};
      void SetVoigtIndexTensor             (VoigtIndexType rIndexVoigtTensor)       {mIndexVoigtTensor = rIndexVoigtTensor;};

      void SetIntVariableData              (const Variable<int>& rVariable, int& rValue) {InternalVariable.SetIntVariableValue(rVariable,rValue);};
      void SetDoubleVariableData           (const Variable<double>& rVariable, double& rValue) {InternalVariable.SetDoubleVariableValue(rVariable,rValue);};
      void SetVectorVariableData           (const Variable<Vector>& rVariable, Vector& rValue) {InternalVariable.SetVectorVariableValue(rVariable,rValue);};
      void SetMatrixVariableData           (const Variable<Matrix>& rVariable, Matrix& rValue) {InternalVariable.SetMatrixVariableValue(rVariable,rValue);};
      void SetArray3VariableData           (const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& rValue) {InternalVariable.SetArray3VariableValue(rVariable,rValue);};
      void SetArray6VariableData           (const Variable<array_1d<double,6> >& rVariable, array_1d<double,6>& rValue) {InternalVariable.SetArray6VariableValue(rVariable,rValue);};

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
      const double&        GetDeltaDeformationDet          () const {return mConstitutiveLawData.DeltaDeformationDet;};
      const double&        GetTotalDeformationDet          () const {return mConstitutiveLawData.TotalDeformationDet;};
      const double&        GetCharacteristicSize           () const {return mConstitutiveLawData.CharacteristicSize;};

      const StressMeasureType& GetStressMeasure            () const {return mConstitutiveLawData.StressMeasure;};
      const StrainMeasureType& GetStrainMeasure            () const {return mConstitutiveLawData.StrainMeasure;};

      const MatrixType&    GetDeltaDeformationMatrix       () const {return mConstitutiveLawData.DeltaDeformationMatrix;};
      const MatrixType&    GetTotalDeformationMatrix       () const {return mConstitutiveLawData.TotalDeformationMatrix;};

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
      return Kratos::make_shared<ConstitutiveModelData>(*this);
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

      //rValues.MaterialParameters.BulkModulus   = rValues.MaterialParameters.LameLambda + (2.0/3.0) * rValues.MaterialParameters.LameMu;

      //hyperelastic model parameters
      if( rProperties.Has(MATERIAL_PARAMETERS) ){
	  Vector ModelParameters = rProperties[MATERIAL_PARAMETERS];
	  for(unsigned int i=0; i<ModelParameters.size(); i++)
	  {
	      rValues.MaterialParameters.ModelParameters.push_back(ModelParameters[i]);
	  }
      }
      else{

	  if( rProperties.Has(C10) ){
	      rValues.MaterialParameters.ModelParameters.push_back(rProperties[C10]);

	      //make neo-hookean consistent with the parameters:
	      rValues.MaterialParameters.LameMu = 2.0 * rProperties[C10];
	  }

	  if( rProperties.Has(C20) )
	      rValues.MaterialParameters.ModelParameters.push_back(rProperties[C20]);

	  if( rProperties.Has(C30) )
	      rValues.MaterialParameters.ModelParameters.push_back(rProperties[C30]);

      }

      if( rProperties.Has(BULK_MODULUS) ){
	  rValues.MaterialParameters.BulkModulus = rProperties[BULK_MODULUS];

	  //make neo-hookean consistent with the parameters:
	  rValues.MaterialParameters.LameLambda = rValues.MaterialParameters.BulkModulus - (2.0/3.0) * rValues.MaterialParameters.LameMu;
      }
      else{
	  rValues.MaterialParameters.BulkModulus = rValues.MaterialParameters.LameLambda + (2.0/3.0) * rValues.MaterialParameters.LameMu;
      }

      //std::cout<<" Mu "<<rValues.MaterialParameters.LameMu<<" Lambda "<<rValues.MaterialParameters.LameLambda<<" BulkModulus "<<rValues.MaterialParameters.BulkModulus<<std::endl;


      //infinitessimal strain (plasticity mu_bar := mu)
      rValues.MaterialParameters.LameMuBar     = rValues.MaterialParameters.LameMu;

      //std::cout<<" B: Mu "<<rValues.MaterialParameters.LameMu<<" Lambda "<<rValues.MaterialParameters.LameLambda<<" BulkModulus "<<rValues.MaterialParameters.BulkModulus<<std::endl;

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
