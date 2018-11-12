//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Nelson Maireni Lafontaine
//                   Josep Maria Carbonell
//


#if !defined(KRATOS_CONSTITUTIVE_LAW )
#define  KRATOS_CONSTITUTIVE_LAW

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "containers/data_value_container.h"
#include "containers/flags.h"


namespace Kratos
{

/**
 * Base class of constitutive laws.
 */
class KRATOS_API(KRATOS_CORE) ConstitutiveLaw : public Flags
{
public:


   enum StrainMeasure
   {
        StrainMeasure_Infinitesimal,   //strain measure small displacements
        StrainMeasure_GreenLagrange,   //strain measure reference configuration
        StrainMeasure_Almansi,         //strain measure current configuration

        // True strain:
        StrainMeasure_Hencky_Material, //strain measure reference configuration
        StrainMeasure_Hencky_Spatial,  //strain measure current   configuration

        // Deformation measures:
        StrainMeasure_Deformation_Gradient, //deformation gradient as a strain measure
        StrainMeasure_Right_CauchyGreen,    //right cauchy-green tensor as a strain measure
        StrainMeasure_Left_CauchyGreen      //left  cauchy-green tensor as a strain measure
    };

    enum StressMeasure
    {
        StressMeasure_PK1,            //stress related to reference configuration non-symmetric
        StressMeasure_PK2,            //stress related to reference configuration
        StressMeasure_Kirchhoff,      //stress related to current   configuration
        StressMeasure_Cauchy          //stress related to current   configuration
    };


    /**
     * Type definitions
     * NOTE: geometries are assumed to be of type Node<3> for all problems
     */
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node < 3 > > GeometryType;

    /**
     * Counted pointer of ConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(ConstitutiveLaw);

    /**
     * Flags related to the Parameters of the Contitutive Law
     */   
    KRATOS_DEFINE_LOCAL_FLAG( USE_ELEMENT_PROVIDED_STRAIN );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_STRESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_CONSTITUTIVE_TENSOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_STRAIN_ENERGY );
    
    KRATOS_DEFINE_LOCAL_FLAG( ISOCHORIC_TENSOR_ONLY );
    KRATOS_DEFINE_LOCAL_FLAG( VOLUMETRIC_TENSOR_ONLY );
    
    KRATOS_DEFINE_LOCAL_FLAG( MECHANICAL_RESPONSE_ONLY );
    KRATOS_DEFINE_LOCAL_FLAG( THERMAL_RESPONSE_ONLY );

    KRATOS_DEFINE_LOCAL_FLAG( INCREMENTAL_STRAIN_MEASURE );

    
    ///the next two flags are designed for internal use within the constitutive law.
    ///please DO NOT use them from the API
    KRATOS_DEFINE_LOCAL_FLAG( INITIALIZE_MATERIAL_RESPONSE );    
    KRATOS_DEFINE_LOCAL_FLAG( FINALIZE_MATERIAL_RESPONSE );


    /**
     * Flags related to the Features of the Contitutive Law
     */

    KRATOS_DEFINE_LOCAL_FLAG( FINITE_STRAINS );
    KRATOS_DEFINE_LOCAL_FLAG( INFINITESIMAL_STRAINS );

    KRATOS_DEFINE_LOCAL_FLAG( THREE_DIMENSIONAL_LAW );
    KRATOS_DEFINE_LOCAL_FLAG( PLANE_STRAIN_LAW );
    KRATOS_DEFINE_LOCAL_FLAG( PLANE_STRESS_LAW );
    KRATOS_DEFINE_LOCAL_FLAG( AXISYMMETRIC_LAW );

    KRATOS_DEFINE_LOCAL_FLAG( U_P_LAW );

    KRATOS_DEFINE_LOCAL_FLAG( ISOTROPIC );
    KRATOS_DEFINE_LOCAL_FLAG( ANISOTROPIC );


    struct Features
    {

      KRATOS_CLASS_POINTER_DEFINITION(Features);
      
    /**
     * Structure "Features" to be used by the element to get the the constitutive law characteristics*
     * its variables will be used to check constitutive law and element compatibility

     * @param mOptions        flags  with the current constitutive law characteristics
     * @param mStrainSize     double with the strain vector size
     * @param mStrainMeasures vector with the strain measures accepted by the constitutive law

     */

      Flags                mOptions;
      double               mStrainSize;
      double               mSpaceDimension;
      std::vector< StrainMeasure > mStrainMeasures;

      /**
       * Constructor.
       */
      Features()
      {
      }
      
      /**
       * Destructor.
       */
      ~Features()
      {
      }

      // Set variables
      void SetOptions       (const Flags&  rOptions)      {mOptions=rOptions;};
      void SetStrainSize    (const double StrainSize)     {mStrainSize=StrainSize;};
      void SetSpaceDimension(const double SpaceDimension) {mSpaceDimension=SpaceDimension;};
      void SetStrainMeasure (const StrainMeasure Measure) {mStrainMeasures.push_back(Measure);};

      void SetStrainMeasures (const std::vector<StrainMeasure> MeasuresVector) {mStrainMeasures = MeasuresVector;};

      // Get variables
      const Flags& GetOptions () {return mOptions;};

      const double& GetStrainSize() {return mStrainSize;};
      const double& GetSpaceDimension() {return mSpaceDimension;};
      std::vector<StrainMeasure>& GetStrainMeasures() {return mStrainMeasures;};
    };



    struct Parameters
    {
        KRATOS_CLASS_POINTER_DEFINITION(Parameters);

    /**
     * Structure "Parameters" to be used by the element to pass the parameters into the constitutive law *

     * @param mOptions flags for the current Constitutive Law Parameters (input data)

     * KINEMATIC PARAMETERS:

     *** NOTE: Pointers are used only to point to a certain variable, no "new" or "malloc" can be used for this Parameters ***

     * @param mDeterminantF copy of the determinant of the Current DeformationGradient (although Current F  is also included as a matrix) (input data)
     * @param mpDeformationGradientF  pointer to the current deformation gradient (can be an empty matrix if a linear strain measure is used) (input data)
     * @param mpStrainVector pointer to the current strains (total strains) (input data) (*can be also OUTPUT with USE_ELEMENT_PROVIDED_STRAIN flag)
     * @param mpStressVector pointer to the current stresses (*OUTPUT with COMPUTE_STRESS flag)
     * @param mpConstitutiveMatrix pointer to the material tangent matrix (*OUTPUT with COMPUTE_CONSTITUTIVE_TENSOR flag)

     * GEOMETRIC PARAMETERS:
     * @param mpShapeFunctionsValues pointer to the shape functions values in the current integration point (input data)
     * @param mpShapeFunctionsDerivatives pointer to the shape functions derivaties values in the current integration point (input data)
     * @param mpElementGeometry pointer to the element's geometry (input data)

     * MATERIAL PROPERTIES:
     * @param mpMaterialProperties pointer to the material's Properties object (input data)

     * PROCESS PROPERTIES:
     * @param mpCurrentProcessInfo pointer to current ProcessInfo instance (input data)

     */


    private:

      /*** NOTE: Member Pointers are used only to point to a certain variable, no "new" or "malloc" can be used for this Parameters ***/

      Flags                mOptions;
      double               mDeterminantF;

      Vector*              mpStrainVector;
      Vector*              mpStressVector;

      const Vector*        mpShapeFunctionsValues;
      const Matrix*        mpShapeFunctionsDerivatives;

      const Matrix*        mpDeformationGradientF;
      Matrix*              mpConstitutiveMatrix;

      const ProcessInfo*   mpCurrentProcessInfo;
      const Properties*    mpMaterialProperties;
      const GeometryType*  mpElementGeometry;


    public:


      /**
       * Constructor.
       */
      Parameters ()
      {
          //Initialize pointers to NULL
          mDeterminantF=0;
          mpStrainVector=NULL;
          mpStressVector=NULL;
          mpShapeFunctionsValues=NULL;
          mpShapeFunctionsDerivatives=NULL;
          mpDeformationGradientF=NULL;
          mpConstitutiveMatrix=NULL;
          mpCurrentProcessInfo=NULL;
          mpMaterialProperties=NULL;
          mpElementGeometry=NULL;
      };


      /**
       * Constructor with Properties, Geometry and ProcessInfo
       */
      Parameters (
          const GeometryType& rElementGeometry,
          const Properties& rMaterialProperties,
          const ProcessInfo& rCurrentProcessInfo)
      :mpCurrentProcessInfo(&rCurrentProcessInfo)
      ,mpMaterialProperties(&rMaterialProperties)
      ,mpElementGeometry(&rElementGeometry)
      {
          //Initialize pointers to NULL
          mDeterminantF=0;
          mpStrainVector=NULL;
          mpStressVector=NULL;
          mpShapeFunctionsValues=NULL;
          mpShapeFunctionsDerivatives=NULL;
          mpDeformationGradientF=NULL;
          mpConstitutiveMatrix=NULL;
      };

      /**
       * Copy Constructor.
       */
      Parameters (const Parameters & rNewParameters)
        :mOptions(rNewParameters.mOptions)
        ,mDeterminantF(rNewParameters.mDeterminantF)
        ,mpStrainVector(rNewParameters.mpStrainVector)
        ,mpStressVector(rNewParameters.mpStressVector)
        ,mpShapeFunctionsValues(rNewParameters.mpShapeFunctionsValues)
        ,mpShapeFunctionsDerivatives(rNewParameters.mpShapeFunctionsDerivatives)
        ,mpDeformationGradientF(rNewParameters.mpDeformationGradientF)
        ,mpConstitutiveMatrix(rNewParameters.mpConstitutiveMatrix)
        ,mpCurrentProcessInfo(rNewParameters.mpCurrentProcessInfo)
        ,mpMaterialProperties(rNewParameters.mpMaterialProperties)
        ,mpElementGeometry(rNewParameters.mpElementGeometry)
      {
      };

      /**
       * Destructor.
       */
      ~Parameters()
      {
      }

      /**
       * Verify Parameters
       */

      bool CheckAllParameters ()
      {
        if(CheckMechanicalVariables() &&  CheckShapeFunctions() && CheckInfoMaterialGeometry ())
            return 1;
        else
            return 0;
      }


      /**
       *Check currentprocessinfo, material properties and geometry
       */

      bool CheckShapeFunctions ()
      {
        if(!mpShapeFunctionsValues)
            KRATOS_ERROR << "ShapeFunctionsValues NOT SET" << std::endl;

        if(!mpShapeFunctionsDerivatives)
            KRATOS_ERROR << "ShapeFunctionsDerivatives NOT SET" << std::endl;

        return 1;
      }

      /**
       *Check currentprocessinfo, material properties and geometry
       */

      bool CheckInfoMaterialGeometry ()
      {
        if(!mpCurrentProcessInfo)
            KRATOS_ERROR << "CurrentProcessInfo NOT SET" << std::endl;

        if(!mpMaterialProperties)
            KRATOS_ERROR << "MaterialProperties NOT SET" << std::endl;

        if(!mpElementGeometry)
            KRATOS_ERROR << "ElementGeometry NOT SET" << std::endl;

        return 1;
      }


      /**
       *Check deformation gradient, strains ans stresses assigned
       */

      bool CheckMechanicalVariables ()
      {
          if(mDeterminantF<=0)
            KRATOS_ERROR << "DeterminantF NOT SET, value <= 0" << std::endl;

          if(!mpDeformationGradientF)
            KRATOS_ERROR << "DeformationGradientF NOT SET" << std::endl;

          if(!mpStrainVector)
            KRATOS_ERROR << "StrainVector NOT SET" << std::endl;

          if(!mpStressVector)
            KRATOS_ERROR << "StressVector NOT SET" << std::endl;

          if(!mpConstitutiveMatrix)
            KRATOS_ERROR << "ConstitutiveMatrix NOT SET" << std::endl;

          return 1;
      }


      /**
       * Public Methods to access variables of the struct class
       */

      /**
       * sets the variable or the pointer of a specified variable: assigns the direction of the pointer for the mpvariables, only non const values can be modified
       */

      void Set                             (Flags ThisFlag)                           {mOptions.Set(ThisFlag);};
      void Reset                           (Flags ThisFlag)                           {mOptions.Reset(ThisFlag);};

      void SetOptions                      (const Flags&  rOptions)                   {mOptions=rOptions;};
      void SetDeterminantF                 (const double DeterminantF)                {mDeterminantF=DeterminantF;};
      //void SetDeterminantF                 (const double& rDeterminantF)              {mDeterminantF=&rDeterminantF;};

      void SetShapeFunctionsValues         (const Vector& rShapeFunctionsValues)      {mpShapeFunctionsValues=&rShapeFunctionsValues;};
      void SetShapeFunctionsDerivatives    (const Matrix& rShapeFunctionsDerivatives) {mpShapeFunctionsDerivatives=&rShapeFunctionsDerivatives;};

      void SetDeformationGradientF         (const Matrix& rDeformationGradientF)      {mpDeformationGradientF=&rDeformationGradientF;};

      void SetStrainVector                 (Vector& rStrainVector)                    {mpStrainVector=&rStrainVector;};
      void SetStressVector                 (Vector& rStressVector)                    {mpStressVector=&rStressVector;};
      void SetConstitutiveMatrix           (Matrix& rConstitutiveMatrix)              {mpConstitutiveMatrix =&rConstitutiveMatrix;};

      void SetProcessInfo                  (const ProcessInfo& rProcessInfo)          {mpCurrentProcessInfo =&rProcessInfo;};
      void SetMaterialProperties           (const Properties&  rMaterialProperties)   {mpMaterialProperties =&rMaterialProperties;};
      void SetElementGeometry              (const GeometryType& rElementGeometry)     {mpElementGeometry =&rElementGeometry;};

      /**
       * Returns the reference or the value of a specified variable: returns the value of the parameter, only non const values can be modified
       */
      Flags& GetOptions () {return mOptions;};

      const double& GetDeterminantF              () {return mDeterminantF;};
      const Vector& GetShapeFunctionsValues      () {return *mpShapeFunctionsValues;};
      const Matrix& GetShapeFunctionsDerivatives () {return *mpShapeFunctionsDerivatives;};
      const Matrix& GetDeformationGradientF      () {return *mpDeformationGradientF;};

      Vector& GetStrainVector                    () {return *mpStrainVector;};
      Vector& GetStressVector                    () {return *mpStressVector;};

      Matrix& GetConstitutiveMatrix              () {return *mpConstitutiveMatrix;};

      const ProcessInfo&  GetProcessInfo         () {return *mpCurrentProcessInfo;};
      const Properties&   GetMaterialProperties  () {return *mpMaterialProperties;};
      const GeometryType& GetElementGeometry     () {return *mpElementGeometry;};

      /**
       * Returns the reference to the value of a specified variable with not constant access
       */

      double& GetDeterminantF                  (double & rDeterminantF) {rDeterminantF=mDeterminantF; return rDeterminantF;};
      Vector& GetStrainVector                  (Vector & rStrainVector) {rStrainVector=*mpStrainVector; return rStrainVector;};
      Matrix& GetDeformationGradientF          (Matrix & rDeformationGradientF)  {rDeformationGradientF=*mpDeformationGradientF;   return rDeformationGradientF;};
      Vector& GetStressVector                  (Vector & rStressVector) {rStressVector=*mpStressVector; return rStressVector;};
      Matrix& GetConstitutiveMatrix            (Matrix & rConstitutiveMatrix) {rConstitutiveMatrix=*mpConstitutiveMatrix; return rConstitutiveMatrix;};

      /**
       * Returns if the different components has been set
       */
       
      bool IsSetDeterminantF              () {return (mDeterminantF > 0.0);};
      bool IsSetShapeFunctionsValues      () {return (mpShapeFunctionsValues != NULL);};
      bool IsSetShapeFunctionsDerivatives () {return (mpShapeFunctionsDerivatives != NULL);};
      bool IsSetDeformationGradientF      () {return (mpDeformationGradientF != NULL);};

      bool IsSetStrainVector              () {return (mpStrainVector != NULL);};
      bool IsSetStressVector              () {return (mpStressVector != NULL);};

      bool IsSetConstitutiveMatrix        () {return (mpConstitutiveMatrix != NULL);};

      bool IsSetProcessInfo               () {return (mpCurrentProcessInfo != NULL);};
      bool IsSetMaterialProperties        () {return (mpMaterialProperties != NULL);};
      bool IsSetElementGeometry           () {return (mpElementGeometry != NULL);};

    };// struct Parameters end

    /**
     * Constructor.
     */
    ConstitutiveLaw();

    /**
     * Destructor.
     */
    ~ConstitutiveLaw() override{};

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
     *      return p_clone;
     */
    virtual ConstitutiveLaw::Pointer Clone() const;

    /**
     * @return the working space dimension of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension();

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType GetStrainSize();

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<int>& rThisVariable);

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<double>& rThisVariable);

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<Vector>& rThisVariable);

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<Matrix>& rThisVariable);

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable);

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 6 > >& rThisVariable);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual int& GetValue(const Variable<int>& rThisVariable, int& rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @return the value of the specified variable
     */
    virtual Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                            array_1d<double, 3 > & rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 6 > & GetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                            array_1d<double, 6 > & rValue);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<int>& rVariable,
                          const int& Value,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<double>& rVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<Vector >& rVariable,
                          const Vector& rValue,
			  const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<Matrix >& rVariable,
                          const Matrix& rValue,
			  const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                          const array_1d<double, 3 > & rValue,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                          const array_1d<double, 6 > & rValue,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */    
    virtual int& CalculateValue(Parameters& rParameterValues, const Variable<int>& rThisVariable, int& rValue);

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */    
    virtual double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue);

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */    
    virtual Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue);

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */  
    virtual Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue);

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual array_1d<double, 3 > & CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 3 > >& rVariable,
						  array_1d<double, 3 > & rValue);

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 6 > & CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 6 > >& rVariable,
						  array_1d<double, 6 > & rValue);


     /**
      * Is called to check whether the provided material parameters in the Properties
      * match the requirements of current constitutive model.
      * @param rMaterialProperties the current Properties to be validated against.
      * @return true, if parameters are correct; false, if parameters are insufficient / faulty
      * NOTE: this has to be implemented by each constitutive model. Returns false in base class since
      * no valid implementation is contained here.
      */
    virtual bool ValidateInput(const Properties& rMaterialProperties);
    

    /**
     * returns the expected strain measure of this constitutive law (by default linear strains)
     * @return the expected strain measure
     */
    virtual StrainMeasure GetStrainMeasure();

    /**
     * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    virtual StressMeasure GetStressMeasure();

    /**
     * returns whether this constitutive model is formulated in incremental strains/stresses
     * NOTE: by default, all constitutive models should be formulated in total strains
     * @return true, if formulated in incremental strains/stresses, false otherwise
     */
    virtual bool IsIncremental();

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void InitializeSolutionStep(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                        const Vector& rShapeFunctionsValues,
                                        const ProcessInfo& rCurrentProcessInfo);

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void FinalizeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo);



    /**
     * to be called at the beginning of each step iteration
     * (e.g. from Element::InitializeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void InitializeNonLinearIteration(const Properties& rMaterialProperties,
					      const GeometryType& rElementGeometry,
					      const Vector& rShapeFunctionsValues,
					      const ProcessInfo& rCurrentProcessInfo);



    /**
     * to be called at the end of each step iteration
     * (e.g. from Element::FinalizeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
					    const GeometryType& rElementGeometry,
					    const Vector& rShapeFunctionsValues,
					    const ProcessInfo& rCurrentProcessInfo);



    /**
     * Computes the material response in terms of stresses and constitutive tensor
     * @see Parameters
     * @see StressMeasure
     */
    void CalculateMaterialResponse (Parameters& rValues,const StressMeasure& rStressMeasure);



    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    virtual void CalculateMaterialResponsePK1 (Parameters& rValues);


    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    virtual void CalculateMaterialResponsePK2 (Parameters& rValues);


    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (Parameters& rValues);


    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    virtual void CalculateMaterialResponseCauchy (Parameters& rValues);



    /**
     * Initialize the material response,  called by the element in FinalizeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */
    void InitializeMaterialResponse (Parameters& rValues,const StressMeasure& rStressMeasure);


    /**
     * Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */

    virtual void InitializeMaterialResponsePK1 (Parameters& rValues);

    /**
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */

    virtual void InitializeMaterialResponsePK2 (Parameters& rValues);

    /**
     * Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */

    virtual void InitializeMaterialResponseKirchhoff (Parameters& rValues);

    /**
     * Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */

    virtual void InitializeMaterialResponseCauchy (Parameters& rValues);


    
    /**
     * Finalize the material response,  called by the element in FinalizeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */
    void FinalizeMaterialResponse (Parameters& rValues,const StressMeasure& rStressMeasure);


    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponsePK1 (Parameters& rValues);

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponsePK2 (Parameters& rValues);

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponseKirchhoff (Parameters& rValues);

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponseCauchy (Parameters& rValues);




    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues);


    /**
     * Methods to transform strain Vectors:
     * @param rStrainVector the strain tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStrainInitial the measure of stress of the given  rStrainVector
     * @param rStrainFinal the measure of stress of the returned rStrainVector
     */
    Vector& TransformStrains        (Vector& rStrainVector,
				     const Matrix &rF,
				     StrainMeasure rStrainInitial,
				     StrainMeasure rStrainFinal);

    /**
     * Methods to transform stress Matrices:
     * @param rStressMatrix the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressInitial the measure of stress of the given  rStressMatrix
     * @param rStressFinal the measure of stress of the returned rStressMatrix
     */
    virtual Matrix& TransformStresses (Matrix& rStressMatrix,
				       const Matrix &rF,
				       const double &rdetF,
				       StressMeasure rStressInitial,
				       StressMeasure rStressFinal);


    /**
     * Methods to transform stress Vectors:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressInitial the measure of stress of the given  rStressVector
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    virtual Vector& TransformStresses (Vector& rStressVector,
				       const Matrix &rF,
				       const double &rdetF,
				       StressMeasure rStressInitial,
				       StressMeasure rStressFinal);



    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure PK1:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& TransformPK1Stresses (Vector& rStressVector,
				  const Matrix &rF,
				  const double &rdetF,
				  StressMeasure rStressFinal);

    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure PK2:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& TransformPK2Stresses (Vector& rStressVector,
				  const Matrix &rF,
				  const double &rdetF,
				  StressMeasure rStressFinal);

    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure Kirchhoff:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& TransformKirchhoffStresses (Vector& rStressVector,
					const Matrix &rF,
					const double &rdetF,
					StressMeasure rStressFinal);

    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure Cauchy:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& TransformCauchyStresses (Vector& rStressVector,
				     const Matrix &rF,
				     const double &rdetF,
				     StressMeasure rStressFinal);



    /**
     * Methods to transform Constitutive Matrices:
     * @param rConstitutiveMatrix the constitutive matrix
     * @param rF the DeformationGradientF matrix between the configurations
     */

    /**
     * This method performs a pull-back of the constitutive matrix
     */
    void PullBackConstitutiveMatrix ( Matrix& rConstitutiveMatrix,
				      const Matrix & rF );


    /**
     * This method performs a push-forward of the constitutive matrix
     */
    void PushForwardConstitutiveMatrix ( Matrix& rConstitutiveMatrix,
					 const Matrix & rF );


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    virtual void GetLawFeatures(Features& rFeatures);


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    virtual int Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo);


		virtual void SetPlasticVariables( const double& rInitialPreconPressure, const double& rInitialBonding); 
		virtual void GetHardeningParameters( double& rPreconPressure, double& rBonding); 
		virtual const double GetBonding();
		virtual const double GetPreconPressure();


    //*** OUTDATED METHODS: ***//


    /**
     * Computes the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    virtual void CalculateMaterialResponse(const Vector& StrainVector,
                                           const Matrix& DeformationGradient,
                                           Vector& StressVector,
                                           Matrix& AlgorithmicTangent,
                                           const ProcessInfo& rCurrentProcessInfo,
                                           const Properties& rMaterialProperties,
                                           const GeometryType& rElementGeometry,
                                           const Vector& rShapeFunctionsValues,
                                           bool CalculateStresses = true,
                                           int CalculateTangent = true,
                                           bool SaveInternalVariables = true);


    /**
     * Computes the volumetric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    virtual void CalculateVolumetricResponse(const double VolumetricStrain,
					     const Matrix& DeformationGradient,
					     double& VolumetricStress,
					     double& AlgorithmicBulk,
					     const ProcessInfo& rCurrentProcessInfo,
					     const Properties& rMaterialProperties,
					     const GeometryType& rElementGeometry,
					     const Vector& rShapeFunctionsValues,
					     bool CalculateStresses,
					     int CalculateTangent,
					     bool SaveInternalVariables);

    /**
     * Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables


     * TODO: add proper definition for algorithmic tangent
     */
    virtual void CalculateDeviatoricResponse(const Vector& StrainVector,
					     const Matrix& DeformationGradient,
					     Vector& StressVector,
					     Matrix& AlgorithmicTangent,
					     const ProcessInfo& rCurrentProcessInfo,
					     const Properties& rMaterialProperties,
					     const GeometryType& rElementGeometry,
					     const Vector& rShapeFunctionsValues,
					     bool CalculateStresses = true,
					     int CalculateTangent = true,
					     bool SaveInternalVariables = true);


    // VM

    virtual void CalculateCauchyStresses(Vector& Cauchy_StressVector,
                                         const Matrix& F,
                                         const Vector& PK2_StressVector,
                                         const Vector& GreenLagrangeStrainVector);


    ///@}
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
        buffer << "ConstitutiveLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ConstitutiveLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "ConstitutiveLaw has no data";
    }

    
    ///@}
    ///@name Friends
    ///@{
    ///@}
    
protected:

    ///@name Protected static Member Variables
    ///@{
    static const unsigned int msIndexVoigt3D6C [6][2];
    static const unsigned int msIndexVoigt2D4C [4][2];
    static const unsigned int msIndexVoigt2D3C [3][2];
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * This method performs a contra-variant push-forward between to tensors
     * i.e. 2nd PK stress to Kirchhoff stress
     */
    void ContraVariantPushForward( Matrix& rMatrix,
				   const Matrix& rF );

    /**
     * This method performs a contra-variant pull-back between to tensors
     * i.e. Kirchhoff stress to 2nd PK stress
     */
    void ContraVariantPullBack( Matrix& rMatrix,
				const Matrix& rF );


    /**
     * This method performs a co-variant push-forward between to tensors
     * i.e. Green-Lagrange strain to Almansi strain
     */
    void CoVariantPushForward( Matrix& rMatrix,
			       const Matrix& rF );


    /**
     * This method performs a co-variant pull-back between to tensors
     * i.e. Almansi strain to Green-Lagrange strain
     */
    void CoVariantPullBack( Matrix& rMatrix,
			    const Matrix& rF );


    /**
     * This method performs a pull-back or a push-forward between two constitutive matrices
     */
    void ConstitutiveMatrixTransformation ( Matrix& rConstitutiveMatrix,
					    const Matrix& rOriginalConstitutiveMatrix,
					    const Matrix & rF );


    /**
     * This method performs a pull-back or a push-forward between two constitutive tensor components
     */
    double& TransformConstitutiveComponent(double & rCabcd,
					   const Matrix & rConstitutiveMatrix,
					   const Matrix & rF,
					   const unsigned int& a, const unsigned int& b,
					   const unsigned int& c, const unsigned int& d);

    /**
     * This method gets the constitutive tensor components
     * from a consitutive matrix supplied in voigt notation
     */
    double& GetConstitutiveComponent(double & rCabcd,
				     const Matrix& rConstitutiveMatrix,
				     const unsigned int& a, const unsigned int& b,
				     const unsigned int& c, const unsigned int& d);

 
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


    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    }

 
    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
    
}; /* Class ConstitutiveLaw */

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  ConstitutiveLaw& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const ConstitutiveLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
 
///@}
///@} addtogroup block
 
template class KRATOS_API(KRATOS_CORE) KratosComponents<ConstitutiveLaw >;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, ConstitutiveLaw const& ThisComponent);

/**
 * Definition of ConstitutiveLaw variable
 */

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT


} /* namespace Kratos.*/
#endif /* KRATOS_CONSTITUTIVE_LAW  defined */

