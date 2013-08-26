/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

/* *********************************************************
 *
 *   Last Modified by:    $Author:   JMCarbonell$
 *   Date:                $Date:     2-06-2013$
 *   Revision:            $Revision: 1.5$
 *
 * ***********************************************************/

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
class ConstitutiveLaw : public Flags
{
public:

   
   enum StrainMeasure
    {
        StrainMeasure_Infinitesimal,  //strain measure small displacements
	StrainMeasure_GreenLagrange,  //strain measure reference configuration
	StrainMeasure_Almansi,        //strain measure current configuration
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
     * Flags related to the constitutive Law
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_STRAIN );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_STRESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_CONSTITUTIVE_TENSOR );
    
    KRATOS_DEFINE_LOCAL_FLAG( ISOCHORIC_TENSOR_ONLY );
    KRATOS_DEFINE_LOCAL_FLAG( VOLUMETRIC_TENSOR_ONLY );
      
    KRATOS_DEFINE_LOCAL_FLAG( TOTAL_TENSOR );
    
    KRATOS_DEFINE_LOCAL_FLAG( INITIAL_CONFIGURATION );
    KRATOS_DEFINE_LOCAL_FLAG( LAST_KNOWN_CONFIGURATION );
    KRATOS_DEFINE_LOCAL_FLAG( FINAL_CONFIGURATION );
    
    KRATOS_DEFINE_LOCAL_FLAG( FINALIZE_MATERIAL_RESPONSE );
    
    KRATOS_DEFINE_LOCAL_FLAG( AXISYMMETRIC );

    /**
     *Structure to be used by the element to pass the parameters into the constitutive law
     * @param mOptions flags for the current Constitutive Law Parameters
     * @param mDeterminantF copy of the determinant of the Current DeformationGradient (althought Current F is also included as a matrix)
     * @param mDeterminantF0 copy of the determinant of the Total DeformationGradient (althought Toal F0 is also included as a matrix)

     * @param mpStrainVector pointer to the current strains (total strains, input)
     * @param mpStressVector pointer to the computed stresses (output)
     * @param mpShapeFunctionsValues pointer to the shape functions values in the current integration point
     * @param mpShapeFunctionsDerivatires pointer to the shape functions derivaties values in the current integration point

     * @param mpDeformationGradientF pointer to the current deformation gradient (can be an empty matrix if a linear strain measure is used)

     * @param mpDeformationGradientF0 pointer to the total deformation gradient (can be an empty matrix if a linear strain measure is used)

     * @param mpConstitutiveMatrix pointer to the material tangent matrix (output)

     * @param mpCurrentProcessInfo pointer to current ProcessInfo instance
     * @param mpMaterialProperties pointer to the material's Properties object
     * @param mpElementGeometry pointer to the element's geometry

     */


    struct Parameters
    {
      
    private:
    
      Flags                mOptions;
      double               mDeterminantF;
      double               mDeterminantF0;
      Vector*              mpStrainVector;
      Vector*              mpStressVector; 

      const Vector*        mpShapeFunctionsValues;
      const Matrix*        mpShapeFunctionsDerivatives;

      const Matrix*        mpDeformationGradientF;
      Matrix*              mpDeformationGradientF0;
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
	mDeterminantF=-1;
	mDeterminantF0=-1;
	mpStrainVector=NULL;
	mpStressVector=NULL;
	mpShapeFunctionsValues=NULL;
	mpShapeFunctionsDerivatives=NULL;
	mpDeformationGradientF=NULL;
	mpDeformationGradientF0=NULL;
	mpConstitutiveMatrix=NULL;
	mpCurrentProcessInfo=NULL;
	mpMaterialProperties=NULL;
	mpElementGeometry=NULL;
      };
    

      /**
       * Constructor with Properties, Geometry and ProcessInfo
       */
      Parameters (const GeometryType& rElementGeometry
		  ,const Properties& rMaterialProperties
		  ,const ProcessInfo& rCurrentProcessInfo)
      :mpCurrentProcessInfo(&rCurrentProcessInfo)
      ,mpMaterialProperties(&rMaterialProperties)
      ,mpElementGeometry(&rElementGeometry)
      {  
	mDeterminantF=-1;
	mDeterminantF0=-1;
	mpStrainVector=NULL;
	mpStressVector=NULL;
	mpShapeFunctionsValues=NULL;
	mpShapeFunctionsDerivatives=NULL;
	mpDeformationGradientF=NULL;
	mpDeformationGradientF0=NULL;
	mpConstitutiveMatrix=NULL;
      };

 
      /**
       * Copy Constructor.
       */
      Parameters (const Parameters & rNewParameters)
        :mOptions(rNewParameters.mOptions)
        ,mDeterminantF(rNewParameters.mDeterminantF)
        ,mDeterminantF0(rNewParameters.mDeterminantF0)
	,mpStrainVector(rNewParameters.mpStrainVector)
	,mpStressVector(rNewParameters.mpStressVector)
	,mpShapeFunctionsValues(rNewParameters.mpShapeFunctionsValues)
	,mpShapeFunctionsDerivatives(rNewParameters.mpShapeFunctionsDerivatives)
	,mpDeformationGradientF(rNewParameters.mpDeformationGradientF)
	,mpDeformationGradientF0(rNewParameters.mpDeformationGradientF0)
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
	  KRATOS_ERROR(std::invalid_argument,"ShapeFunctionsValues NOT SET","");

	if(!mpShapeFunctionsDerivatives)
	  KRATOS_ERROR(std::invalid_argument,"ShapeFunctionsDerivatives NOT SET","");

	return 1;
      }

      /**
       *Check currentprocessinfo, material properties and geometry
       */
    

      bool CheckInfoMaterialGeometry ()
      {
	if(!mpCurrentProcessInfo)
	  KRATOS_ERROR(std::invalid_argument,"CurrentProcessInfo NOT SET","");

	if(!mpMaterialProperties)
	  KRATOS_ERROR(std::invalid_argument,"MaterialProperties NOT SET","");

	if(!mpElementGeometry)
	  KRATOS_ERROR(std::invalid_argument,"ElementGeometry NOT SET","");
     
	return 1;
      }


      /**
       *Check deformation gradient, strains ans stresses assigned
       */
      bool CheckMechanicalVariables ()
      {
	if(mDeterminantF==-1)
	  KRATOS_ERROR(std::invalid_argument,"DeterminantF NOT SET","");

	if(mDeterminantF0==-1)
	  KRATOS_ERROR(std::invalid_argument,"DeterminantF0 NOT SET","");

	if(!mpDeformationGradientF)
	  KRATOS_ERROR(std::invalid_argument,"DeformationGradientF NOT SET","");

	if(!mpDeformationGradientF0)
	  KRATOS_ERROR(std::invalid_argument,"DeformationGradientF0 NOT SET","");

	if(!mpStrainVector)
	  KRATOS_ERROR(std::invalid_argument,"StrainVector NOT SET","");

	if(!mpStressVector)
	  KRATOS_ERROR(std::invalid_argument,"StressVector NOT SET","");
     
	if(!mpConstitutiveMatrix)
	  KRATOS_ERROR(std::invalid_argument,"ConstitutiveMatrix NOT SET","");

	return 1;
      }
      

      /**
       * Public Methods to access variables of the struct class
       */

      /**
       * sets the value of a specified variable
       */

      void Set                             (Flags ThisFlag)                           {mOptions.Set(ThisFlag);};
      void Reset                           (Flags ThisFlag)                           {mOptions.Reset(ThisFlag);};

      void SetOptions                      (const Flags&  rOptions)                   {mOptions=rOptions;};
      void SetDeterminantF                 (const double& rDeterminantF)              {mDeterminantF=rDeterminantF;};
      void SetDeterminantF0                (const double& rDeterminantF0)             {mDeterminantF0=rDeterminantF0;};
      void SetStrainVector                 (Vector& rStrainVector)                    {mpStrainVector=&rStrainVector;};
      void SetStressVector                 (Vector& rStressVector)                    {mpStressVector=&rStressVector;};
 
      void SetShapeFunctionsValues         (const Vector& rShapeFunctionsValues)      {mpShapeFunctionsValues=&rShapeFunctionsValues;};
      void SetShapeFunctionsDevivatives    (const Matrix& rShapeFunctionsDerivatives) {mpShapeFunctionsDerivatives=&rShapeFunctionsDerivatives;};

      void SetDeformationGradientF         (const Matrix& rDeformationGradientF)     {mpDeformationGradientF=&rDeformationGradientF;};

      void SetDeformationGradientF0        (Matrix& rDeformationGradientF0)     {mpDeformationGradientF0=&rDeformationGradientF0;};
      
      void SetConstitutiveMatrix           (Matrix& rConstitutiveMatrix)              {mpConstitutiveMatrix =&rConstitutiveMatrix;};

      void SetProcessInfo                  (const ProcessInfo& rProcessInfo)          {mpCurrentProcessInfo =&rProcessInfo;};
      void SetMaterialProperties           (const Properties&  rMaterialProperties)   {mpMaterialProperties =&rMaterialProperties;};
      void SetElementGeometry              (const GeometryType& rElementGeometry)     {mpElementGeometry =&rElementGeometry;};

 
      /**
       * returns the value of a specified variable
       */ 
      Flags& GetOptions () {return mOptions;};
   
      const double& GetDeterminantF              () {return mDeterminantF;};     
      const Vector& GetShapeFunctionsValues      () {return *mpShapeFunctionsValues;};
      const Matrix& GetShapeFunctionsDerivatives () {return *mpShapeFunctionsDerivatives;};
      const Matrix& GetDeformationGradientF      () {return *mpDeformationGradientF;};
      Matrix& GetDeformationGradientF0           () {return *mpDeformationGradientF0;};

      double& GetDeterminantF0                   () {return mDeterminantF0;};
      Vector& GetStrainVector                    () {return *mpStrainVector;};
      Vector& GetStressVector                    () {return *mpStressVector;};

      Matrix& GetConstitutiveMatrix              () {return *mpConstitutiveMatrix;};


      const ProcessInfo&  GetProcessInfo         () {return *mpCurrentProcessInfo;};
      const Properties&   GetMaterialProperties  () {return *mpMaterialProperties;};
      const GeometryType& GetElementGeometry     () {return *mpElementGeometry;};
            


      /**
       * returns the value of a specified variable not constant access
       */ 
      
      double& GetDeterminantF                  (double & rDeterminantF) {rDeterminantF=mDeterminantF; return rDeterminantF;};
      double& GetDeterminantF0                 (double & rDeterminantF0) {rDeterminantF0=mDeterminantF0; return rDeterminantF0;};
      Vector& GetStrainVector                  (Vector & rStrainVector) {rStrainVector=*mpStrainVector; return rStrainVector;};
      Matrix& GetDeformationGradientF          (Matrix & rDeformationGradientF) {rDeformationGradientF=*mpDeformationGradientF; return rDeformationGradientF;};
      Matrix& GetDeformationGradientF0         (Matrix & rDeformationGradientF0) {rDeformationGradientF0=*mpDeformationGradientF0; return rDeformationGradientF0;};

      Vector& GetStressVector              (Vector & rStressVector) {rStressVector=*mpStressVector; return rStressVector;};
      Matrix& GetConstitutiveMatrix        (Matrix & rConstitutiveMatrix) {rConstitutiveMatrix=*mpConstitutiveMatrix; return rConstitutiveMatrix;};
  
    };
  
    /**
     * Constructor.
     */
    ConstitutiveLaw();

    /**
     * Destructor.
     */
    virtual ~ConstitutiveLaw(){};

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
    virtual void SetValue(const Variable<double>& rVariable,
                          const double& Value,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<Vector >& rVariable,
                          const Vector& Value, const ProcessInfo& rCurrentProcessInfo);
 
    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<Matrix >& rVariable,
                          const Matrix& Value, const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                          const array_1d<double, 3 > & Value,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                          const array_1d<double, 6 > & Value,
                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * Is called to check whether the provided material parameters in the Properties
     * match the requirements of current constitutive model.
     * @param props the current Properties to be validated against.
     * @return true, if parameters are correct; false, if parameters are insufficient / faulty
     * NOTE: this has to implemented by each constitutive model. Returns false in base class since
     * no valid implementation is contained here.
     */
    virtual bool ValidateInput(const Properties& props);

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
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual void InitializeMaterial(const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues);

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void InitializeSolutionStep(const Properties& props,
                                        const GeometryType& geom, //this is just to give the array of nodes
                                        const Vector& ShapeFunctionsValues,
                                        const ProcessInfo& CurrentProcessInfo);

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void FinalizeSolutionStep(const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo);
 
    virtual void InitializeNonLinearIteration(const Properties& props,
					      const GeometryType& geom,
					      const Vector& ShapeFunctionsValues,
					      const ProcessInfo& CurrentProcessInfo);


    /**
     * Computes the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param CurrentProcessInfo current ProcessInfo instance
     * @param props the material's Properties object
     * @param geom the element's geometry
     * @param ShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    virtual void CalculateMaterialResponse(const Vector& StrainVector,
                                           const Matrix& DeformationGradient,
                                           Vector& StressVector,
                                           Matrix& AlgorithmicTangent,
                                           const ProcessInfo& CurrentProcessInfo,
                                           const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues,
                                           bool CalculateStresses = true,
                                           int CalculateTangent = true,
                                           bool SaveInternalVariables = true);


    /**
     * Computes the volumetric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param CurrentProcessInfo current ProcessInfo instance
     * @param props the material's Properties object
     * @param geom the element's geometry
     * @param ShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    virtual void CalculateVolumetricResponse(const double VolumetricStrain,
					     const Matrix& DeformationGradient,
					     double& VolumetricStress,
					     double& AlgorithmicBulk,
					     const ProcessInfo& CurrentProcessInfo,
					     const Properties& props,
					     const GeometryType& geom,
					     const Vector& ShapeFunctionsValues,
					     bool CalculateStresses,
					     int CalculateTangent,
					     bool SaveInternalVariables);

    /**
     * Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param CurrentProcessInfo current ProcessInfo instance
     * @param props the material's Properties object
     * TODO: add proper definition for algorithmic tangent
     */
    virtual void CalculateDeviatoricResponse(const Vector& StrainVector,
					     const Matrix& DeformationGradient,
					     Vector& StressVector,
					     Matrix& AlgorithmicTangent,
					     const ProcessInfo& CurrentProcessInfo,
					     const Properties& props,
					     const GeometryType& geom,
					     const Vector& ShapeFunctionsValues,
					     bool CalculateStresses = true,
					     int CalculateTangent = true,
					     bool SaveInternalVariables = true);



    /**
     * Computes the material response in terms of stresses and constitutive tensor
     * @see Parameters
     * @see StressMeasures
     */

    void CalculateMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure);
    


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
     * Updates the material response,  called by the element in FinalizeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */

    void FinalizeMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure);
    

    /**
     * Updates the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponsePK1 (Parameters& rValues);
    
    /**
     * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponsePK2 (Parameters& rValues);

    /**
     * Updates the material response in terms of Kirchhoff stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponseKirchhoff (Parameters& rValues);

    /**
     * Updates the material response in terms of Cauchy stresses
     * @see Parameters
     */

    virtual void FinalizeMaterialResponseCauchy (Parameters& rValues);
 
    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void ResetMaterial(const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues);
    // VM

    virtual void CalculateCauchyStresses(Vector& Cauchy_StressVector,
                                         const Matrix& F,
                                         const Vector& PK2_StressVector,
                                         const Vector& GreenLagrangeStrainVector);
    //VM


    virtual Matrix& TransformStresses (Matrix& rStressMatrix,
				       const Matrix &rF,
				       const double &rdetF,
				       StressMeasure rStressInitial,
				       StressMeasure rStressFinal);
    


    virtual Vector& TransformStresses (Vector& rStressVector,
				       const Matrix &rF,
				       const double &rdetF,
				       StressMeasure rStressInitial,
				       StressMeasure rStressFinal);
    
    //VM


   Vector& TransformPK1Stresses (Vector& rStressVector,
				 const Matrix &rF,
				 const double &rdetF,
				 StressMeasure rStressFinal);
     

   Vector& TransformPK2Stresses (Vector& rStressVector,
				 const Matrix &rF,
				 const double &rdetF,
				 StressMeasure rStressFinal);
     

    Vector& TransformKirchhoffStresses (Vector& rStressVector,
					const Matrix &rF,
					const double &rdetF,
					StressMeasure rStressFinal);
      

    Vector& TransformCauchyStresses (Vector& rStressVector,
				     const Matrix &rF,
				     const double &rdetF,
				     StressMeasure rStressFinal);
      

    virtual Vector& TransformStrains (Vector& rStrainVector,
				       const Matrix &rF,
				       StrainMeasure rStrainInitial,
				      StrainMeasure rStrainFinal);
    
    //VM

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    virtual int Check(const Properties& props,
                      const GeometryType& geom,
                      const ProcessInfo& CurrentProcessInfo);
    

protected:


    void ContraVariantPushForward( Matrix& rMatrix,
				   const Matrix& rF );  //i.e. 2nd PK stress to Kirchhoff stress
    

    void ContraVariantPullBack( Matrix& rMatrix,
				const Matrix& rF );     //i.e. Kirchhoff stress to 2nd PK stress
        
        
    void CoVariantPushForward( Matrix& rMatrix,
			       const Matrix& rF );      //i.e. Green-Lagrange strain to Almansi strain
    

    void CoVariantPullBack( Matrix& rMatrix,
			    const Matrix& rF );         //i.e. Almansi strain to Green-Lagrange strain
 

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;


    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
// 	  rSerializer.save("",);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    }


}; /* Class ConstitutiveLaw */

/**
 * Definition of ConstitutiveLaw variable
 */
KRATOS_DEFINE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW);

} /* namespace Kratos.*/
#endif /* KRATOS_CONSTITUTIVE_LAW  defined */
