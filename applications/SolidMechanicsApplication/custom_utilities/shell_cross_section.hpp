//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             January 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(SHELL_CROSS_SECTION_H_INCLUDED)
#define SHELL_CROSS_SECTION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "containers/flags.h"
#include "custom_utilities/properties_extensions.hpp"

namespace Kratos
{

/** \brief ShellCrossSection
*
* ShellCrossSection is the base class for all shell cross sections.
* This object is meant to be used by shell elements to obtain the material response
* in terms of generalized strains (membrane strains, shear strains, curvatures) and
* generalized stresses (stress resultants, stress couples) by numerical integration
* of several through-the-thickness constitutive laws.
*
* Homogeneous / Composite Section...
*
* Constitutive Law Adaptation...
*
* References...
*
*/
class KRATOS_API(SOLID_MECHANICS_APPLICATION) ShellCrossSection : public Flags
{

public:

	class Ply;

	KRATOS_CLASS_POINTER_DEFINITION(ShellCrossSection);

	typedef Geometry<Node < 3 > > GeometryType;

	typedef std::vector< Ply > PlyCollection;

	typedef std::size_t SizeType;

	///@name Enums
    ///@{

	/** SectionBehaviorType Enum
	* Defines the supported behaviors of the cross section
	*/
	enum SectionBehaviorType
	{
		Thick, /**< Thick section (Mindlin-Reissner Plate Theory) */
		Thin /**< Thin section (Kirchhoff-Love Plate Theory) */
	};

	///@}

	///@name Classes
    ///@{

	struct Features
	{
		Flags mOptions;
		double mStrainSize;
		double mSpaceDimension;
		std::vector< ConstitutiveLaw::StrainMeasure > mStrainMeasures;
	};

	class Parameters
	{

	private:

		Flags                mOptions;

		Vector*              mpGeneralizedStrainVector;
		Vector*              mpGeneralizedStressVector;
		Matrix*              mpConstitutiveMatrix;

		const Vector*        mpShapeFunctionsValues;
		const Matrix*        mpShapeFunctionsDerivatives;
		const ProcessInfo*   mpCurrentProcessInfo;
		const Properties*    mpMaterialProperties;
		const GeometryType*  mpElementGeometry;

	public:

		Parameters()
			: mpGeneralizedStrainVector(NULL)
			, mpGeneralizedStressVector(NULL)
			, mpConstitutiveMatrix(NULL)
			, mpShapeFunctionsValues(NULL)
			, mpShapeFunctionsDerivatives(NULL)
			, mpCurrentProcessInfo(NULL)
			, mpMaterialProperties(NULL)
			, mpElementGeometry(NULL)
		{}

		Parameters (const GeometryType& rElementGeometry,
					const Properties& rMaterialProperties,
					const ProcessInfo& rCurrentProcessInfo)
			: mpGeneralizedStrainVector(NULL)
			, mpGeneralizedStressVector(NULL)
			, mpConstitutiveMatrix(NULL)
			, mpShapeFunctionsValues(NULL)
			, mpShapeFunctionsDerivatives(NULL)
			, mpCurrentProcessInfo(&rCurrentProcessInfo)
			, mpMaterialProperties(&rMaterialProperties)
			, mpElementGeometry(&rElementGeometry)
		{}

		Parameters (const Parameters & rNewParameters)
			: mOptions(rNewParameters.mOptions)
			, mpGeneralizedStrainVector(rNewParameters.mpGeneralizedStrainVector)
			, mpGeneralizedStressVector(rNewParameters.mpGeneralizedStressVector)
			, mpConstitutiveMatrix(rNewParameters.mpConstitutiveMatrix)
			, mpShapeFunctionsValues(rNewParameters.mpShapeFunctionsValues)
			, mpShapeFunctionsDerivatives(rNewParameters.mpShapeFunctionsDerivatives)
			, mpCurrentProcessInfo(rNewParameters.mpCurrentProcessInfo)
			, mpMaterialProperties(rNewParameters.mpMaterialProperties)
			, mpElementGeometry(rNewParameters.mpElementGeometry)
		{}

	public:

		/**
		*Checks shape functions and shape function derivatives
		*/
		bool CheckShapeFunctions ()
		{
			if(!mpShapeFunctionsValues)
				KRATOS_THROW_ERROR(std::invalid_argument,"ShapeFunctionsValues NOT SET","");

			if(!mpShapeFunctionsDerivatives)
				KRATOS_THROW_ERROR(std::invalid_argument,"ShapeFunctionsDerivatives NOT SET","");

			return 1;
		}

		/**
		*Checks currentprocessinfo, material properties and geometry
		*/
		bool CheckInfoMaterialGeometry ()
		{
			if(!mpCurrentProcessInfo)
				KRATOS_THROW_ERROR(std::invalid_argument,"CurrentProcessInfo NOT SET","");

			if(!mpMaterialProperties)
				KRATOS_THROW_ERROR(std::invalid_argument,"MaterialProperties NOT SET","");

			if(!mpElementGeometry)
				KRATOS_THROW_ERROR(std::invalid_argument,"ElementGeometry NOT SET","");

			return 1;
		}

		/**
		*Check deformation gradient, strains ans stresses assigned
		*/
		bool CheckMechanicalVariables ()
		{
			if(!mpGeneralizedStrainVector)
				KRATOS_THROW_ERROR(std::invalid_argument,"GenralizedStrainVector NOT SET","");

			if(!mpGeneralizedStressVector)
				KRATOS_THROW_ERROR(std::invalid_argument,"GenralizedStressVector NOT SET","");

			if(!mpConstitutiveMatrix)
				KRATOS_THROW_ERROR(std::invalid_argument,"ConstitutiveMatrix NOT SET","");

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

		void SetGeneralizedStrainVector      (Vector& rGeneralizedStrainVector)         {mpGeneralizedStrainVector=&rGeneralizedStrainVector;};
		void SetGeneralizedStressVector      (Vector& rGeneralizedStressVector)         {mpGeneralizedStressVector=&rGeneralizedStressVector;};
		void SetConstitutiveMatrix           (Matrix& rConstitutiveMatrix)              {mpConstitutiveMatrix =&rConstitutiveMatrix;};

		void SetShapeFunctionsValues         (const Vector& rShapeFunctionsValues)      {mpShapeFunctionsValues=&rShapeFunctionsValues;};
		void SetShapeFunctionsDerivatives    (const Matrix& rShapeFunctionsDerivatives) {mpShapeFunctionsDerivatives=&rShapeFunctionsDerivatives;};
		void SetProcessInfo                  (const ProcessInfo& rProcessInfo)          {mpCurrentProcessInfo =&rProcessInfo;};
		void SetMaterialProperties           (const Properties&  rMaterialProperties)   {mpMaterialProperties =&rMaterialProperties;};
		void SetElementGeometry              (const GeometryType& rElementGeometry)     {mpElementGeometry =&rElementGeometry;};

		/**
		* returns the reference or the value of a specified variable: returns the value of the parameter, only non const values can be modified
		*/

		Flags& GetOptions () {return mOptions;};

		Vector& GetGeneralizedStrainVector         () {return *mpGeneralizedStrainVector;};
		Vector& GetGeneralizedStressVector         () {return *mpGeneralizedStressVector;};
		Matrix& GetConstitutiveMatrix              () {return *mpConstitutiveMatrix;};

		const Vector& GetShapeFunctionsValues      () {return *mpShapeFunctionsValues;};
		const Matrix& GetShapeFunctionsDerivatives () {return *mpShapeFunctionsDerivatives;};
		const ProcessInfo&  GetProcessInfo         () {return *mpCurrentProcessInfo;};
		const Properties&   GetMaterialProperties  () {return *mpMaterialProperties;};
		const GeometryType& GetElementGeometry     () {return *mpElementGeometry;};
	};

	class IntegrationPoint
	{

	private:

		double mWeight;
		double mLocation;
		ConstitutiveLaw::Pointer mConstitutiveLaw;

	public:

		IntegrationPoint()
			: mWeight(0.0)
			, mLocation(0.0)
			, mConstitutiveLaw(ConstitutiveLaw::Pointer())
		{}

		IntegrationPoint(double location, double weight, const ConstitutiveLaw::Pointer pMaterial)
			: mWeight(weight)
			, mLocation(location)
			, mConstitutiveLaw(pMaterial)
		{}

		IntegrationPoint(const IntegrationPoint& other)
			: mWeight(other.mWeight)
			, mLocation(other.mLocation)
			, mConstitutiveLaw(other.mConstitutiveLaw != NULL ? other.mConstitutiveLaw->Clone() : ConstitutiveLaw::Pointer())
		{}

		IntegrationPoint & operator = (const IntegrationPoint & other) {
			if(this != &other) {
				mWeight = other.mWeight;
				mLocation = other.mLocation;
				mConstitutiveLaw = other.mConstitutiveLaw != NULL ? other.mConstitutiveLaw->Clone() : ConstitutiveLaw::Pointer();
			}
			return *this;
		}

                virtual ~IntegrationPoint(){};

	public:

		inline double GetWeight()const { return mWeight; }
		inline void SetWeight(double w) { mWeight = w; }

		inline double GetLocation()const { return mLocation; }
		inline void SetLocation(double l) { mLocation = l; }

		inline const ConstitutiveLaw::Pointer& GetConstitutiveLaw()const { return mConstitutiveLaw; }
		inline void SetConstitutiveLaw(const ConstitutiveLaw::Pointer& pLaw) { mConstitutiveLaw = pLaw; }

	private:

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			rSerializer.save("W", mWeight);
			rSerializer.save("L", mLocation);
			rSerializer.save("CLaw", mConstitutiveLaw);
		}

		virtual void load(Serializer& rSerializer)
		{
			rSerializer.load("W", mWeight);
			rSerializer.load("L", mLocation);
			rSerializer.load("CLaw", mConstitutiveLaw);
		}
	};

	class Ply
	{

	public:

		typedef std::vector< IntegrationPoint > IntegrationPointCollection;

	private:

		double mThickness;
		double mLocation;
		double mOrientationAngle;
		IntegrationPointCollection mIntegrationPoints;
		Properties::Pointer mpProperties;

	public:

		Ply()
			: mThickness(0.0)
			, mLocation(0.0)
			, mOrientationAngle(0.0)
			, mIntegrationPoints()
			, mpProperties(Properties::Pointer())
		{}

		Ply(double thickness, double location, double orientationAngle, int numPoints, const Properties::Pointer & pProperties)
			: mThickness(thickness)
			, mLocation(location)
			, mIntegrationPoints()
			, mpProperties(pProperties)
		{
			this->SetOrientationAngle(orientationAngle);
			this->SetUpIntegrationPoints(numPoints);
		}

		Ply(const Ply& other)
			: mThickness(other.mThickness)
			, mLocation(other.mLocation)
			, mOrientationAngle(other.mOrientationAngle)
			, mIntegrationPoints(other.mIntegrationPoints)
			, mpProperties(other.mpProperties)
		{}

		Ply & operator = (const Ply & other) {
			if(this != &other) {
				mThickness = other.mThickness;
				mLocation = other.mLocation;
				mOrientationAngle = other.mOrientationAngle;
				mIntegrationPoints = other.mIntegrationPoints;
				mpProperties = other.mpProperties;
			}
			return *this;
		}

                virtual ~Ply(){};

	public:

		inline double GetThickness()const { return mThickness; }
		inline void SetThickness(double thickness) { mThickness = thickness; }

		inline double GetLocation()const { return mLocation; }
		inline void SetLocation(double location) {
			if(location != mLocation) {
				for(IntegrationPointCollection::iterator it = mIntegrationPoints.begin(); it != mIntegrationPoints.end(); ++it)
					(*it).SetLocation((*it).GetLocation() + location - mLocation); // remove the last location and add the new one (this avoids to re-setup the integration points.
				mLocation = location; // update the current location
			}
		}

		inline double GetOrientationAngle()const { return mOrientationAngle; }
		inline void SetOrientationAngle(double degrees) {
			mOrientationAngle = std::fmod(degrees, 360.0);
			if(mOrientationAngle < 0.0)
				mOrientationAngle += 360.0;
		}

		inline const IntegrationPointCollection& GetIntegrationPoints()const { return mIntegrationPoints; }
		inline IntegrationPointCollection& GetIntegrationPoints() { return mIntegrationPoints; }

		inline const Properties::Pointer & GetPropertiesPointer()const { return mpProperties; }

		inline const Properties & GetProperties()const { return *mpProperties; }

		inline double CalculateMassPerUnitArea()const { return mpProperties->GetValue(DENSITY) * mThickness; }

		inline IntegrationPointCollection::size_type NumberOfIntegrationPoints()const { return mIntegrationPoints.size(); }

		inline void SetConstitutiveLawAt(IntegrationPointCollection::size_type integrationPointID, const ConstitutiveLaw::Pointer& pNewConstitutiveLaw)
		{
			if(integrationPointID < mIntegrationPoints.size())
				mIntegrationPoints[integrationPointID].SetConstitutiveLaw(pNewConstitutiveLaw);
		}

	private:

		void SetUpIntegrationPoints(int n)
		{
			KRATOS_TRY

                            const ConstitutiveLaw::Pointer & pMaterial = GetProperties()[CONSTITUTIVE_LAW];
			if(pMaterial == NULL)
                          KRATOS_THROW_ERROR(std::logic_error, "A Ply needs a constitutive law to be set. Missing constitutive law in property : ", GetProperties().Id());

			// make sure the number is greater than 0 and odd
			if(n < 0) n = -n;
			if(n == 0) n = 5;
			if(n % 2 == 0) n += 1;

			// generate the weights (composite simpson rule)
			Vector ip_w(n, 1.0);
			if(n >= 3) {
                          for(int i = 1; i < n-1; i++) {
                            double iw = (i % 2 == 0) ? 2.0 : 4.0;
                            ip_w(i) = iw;
                          }
                          ip_w /= sum( ip_w );
			}

			// generate locations (direction: top(+thickness/2) to bottom(-thickness/2)
			Vector ip_loc(n, 0.0);
			if(n >= 3) {
                          double loc_start = mLocation + 0.5 * mThickness;
                          double loc_incr = mThickness / double(n-1);
                          for(int i = 0; i < n; i++) {
                            ip_loc(i) = loc_start;
                            loc_start -= loc_incr;
                          }
			}

			// generate the integration points
			mIntegrationPoints.clear();
			mIntegrationPoints.resize(n);
			for(int i = 0; i < n; i++) {
                          IntegrationPoint& intp = mIntegrationPoints[i];
                          intp.SetWeight(ip_w(i) * mThickness);
                          intp.SetLocation(ip_loc(i));
                          intp.SetConstitutiveLaw(pMaterial->Clone());
			}

			KRATOS_CATCH("")
                            }

         private:

          friend class Serializer;

          virtual void save(Serializer& rSerializer) const
          {
            rSerializer.save("T", mThickness);
            rSerializer.save("L", mLocation);
            rSerializer.save("O", mOrientationAngle);
            rSerializer.save("IntP", mIntegrationPoints);
            rSerializer.save("Prop", mpProperties);
          }

          virtual void load(Serializer& rSerializer)
          {
            rSerializer.load("T", mThickness);
            rSerializer.load("L", mLocation);
            rSerializer.load("O", mOrientationAngle);
            rSerializer.load("IntP", mIntegrationPoints);
            rSerializer.load("Prop", mpProperties);
          }

	};

 protected:

  struct ElementVariables
  {
    double DeterminantF;
    double DeterminantF0;

    Vector StrainVector_2D;
    Vector StressVector_2D;
    Matrix ConstitutiveMatrix_2D;
    Matrix DeformationGradientF_2D;
    Matrix DeformationGradientF0_2D;

    Vector StrainVector_3D;
    Vector StressVector_3D;
    Matrix ConstitutiveMatrix_3D;
    Matrix DeformationGradientF_3D;
    Matrix DeformationGradientF0_3D;

    double GYZ;
    double GXZ;

    Matrix H;
    Matrix L;
    Matrix LT;
    Vector CondensedStressVector;
  };

  ///@}

 public:

  ///@name Life Cycle
  ///@{

  /**
   * Default constructor
   */
  ShellCrossSection();

  /**
   * Copy constructor
   * @param other the other cross section
   */
  ShellCrossSection(const ShellCrossSection & other);

  /**
   * Destructor
   */
  ~ShellCrossSection() override;

  ///@}

  ///@name Operators
  ///@{

  /**
   * Assignment operator
   * @param other the other cross section
   */
  ShellCrossSection & operator = (const ShellCrossSection & other);

  ///@}

  ///@name Public Methods
  ///@{

  /**
   * Initializes the editing of the Composite Layup.
   * After a call to this method, one or more calls to AddPly(...) can be done to create the stack.
   * After the stack is properly set it is necessary to call EndStack() to finalize the stack editing.
   */
  void BeginStack();

  /**
   * Adds a new Ply below the current one.
   * After the stack is properly set it is necessary to call EndStack() to finalize the stack editing.
   * @param thickness the thickness of the new ply.
   * @param orientationAngle the angle (degrees) between the new ply and the cross section.
   * @param numPoints the number of integration points. can be 1,3,5,7,9,... and so on.
   For numPoints = 3, the Simpson rule is used.
   For numPoints = odd number > 3, the composite Simpson rule is used.
   * @param pProperties the pointer to the properties assigned to the new ply.
   */
  void AddPly(double thickness, double orientationAngle, int numPoints, const Properties::Pointer & pProperties);

  /**
   * Finalizes the editing of the Composite Layup.
   */
  void EndStack();

  /**
   * Returns the string containing a detailed description of this object.
   * @return the string with informations
   */
  virtual std::string GetInfo()const;

  /**
   * Clone function
   * @return a pointer to a new instance of this cross section
   */
  virtual ShellCrossSection::Pointer Clone()const;

  /**
   * returns whether this cross section has specified variable
   * @param rThisVariable the variable to be checked for
   * @return true if the variable is defined in the cross section
   */
  virtual bool Has(const Variable<double>& rThisVariable);

  /**
   * returns whether this cross section has specified variable
   * @param rThisVariable the variable to be checked for
   * @return true if the variable is defined in the cross section
   */
  virtual bool Has(const Variable<Vector>& rThisVariable);

  /**
   * returns whether this cross section has specified variable
   * @param rThisVariable the variable to be checked for
   * @return true if the variable is defined in the cross section
   */
  virtual bool Has(const Variable<Matrix>& rThisVariable);

  /**
   * returns whether this cross section has specified variable
   * @param rThisVariable the variable to be checked for
   * @return true if the variable is defined in the cross section
   * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
   */
  virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable);

  /**
   * returns whether this cross section has specified variable
   * @param rThisVariable the variable to be checked for
   * @return true if the variable is defined in the cross section
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
   * Is called to check whether the provided material parameters in the Properties
   * match the requirements of current constitutive model.
   * @param rMaterialProperties the current Properties to be validated against.
   * @return true, if parameters are correct; false, if parameters are insufficient / faulty
   * NOTE: this has to be implemented by each constitutive model. Returns false in base class since
   * no valid implementation is contained here.
   */
  virtual bool ValidateInput(const Properties& rMaterialProperties);

  /**
   * This is to be called at the very beginning of the calculation
   * (e.g. from InitializeElement) in order to initialize all relevant
   * attributes of the cross section
   * @param rMaterialProperties the Properties instance of the current element
   * @param rElementGeometry the geometry of the current element
   * @param rShapeFunctionsValues the shape functions values in the current integration point
   */
  virtual void InitializeCrossSection(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues);

  /**
   * to be called at the beginning of each solution step
   * (e.g. from Element::InitializeSolutionStep)
   * @param rMaterialProperties the Properties instance of the current element
   * @param rElementGeometry the geometry of the current element
   * @param rShapeFunctionsValues the shape functions values in the current integration point
   * @param rCurrentProcessInfo the current ProcessInfo instance
   */
  virtual void InitializeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo);

  /**
   * to be called at the end of each solution step
   * (e.g. from Element::FinalizeSolutionStep)
   * @param rMaterialProperties the Properties instance of the current element
   * @param rElementGeometry the geometry of the current element
   * @param rShapeFunctionsValues the shape functions values in the current integration point
   * @param rCurrentProcessInfo the current ProcessInfo instance
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
   * @param rCurrentProcessInfo he current ProcessInfo instance
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
   * @param rCurrentProcessInfo the current ProcessInfo instance
   */
  virtual void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const Vector& rShapeFunctionsValues,
                                          const ProcessInfo& rCurrentProcessInfo);

  /**
   * Computes the section response in terms of generalized stresses and constitutive tensor
   * @param rValues the parameters for the current calculation
   * @param rStressMeasure the required stress measure
   * @see Parameters
   */
  virtual void CalculateSectionResponse(Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure);

  /**
   * Updates the section response, called by the element in FinalizeSolutionStep.
   * @param rValues the parameters for the current calculation
   * @param rStressMeasure the required stress measure
   * @see Parameters
   */
  virtual void FinalizeSectionResponse(Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure);

  /**
   * This can be used in order to reset all internal variables of the
   * cross section (e.g. if a model should be reset to its reference state)
   * @param rMaterialProperties the Properties instance of the current element
   * @param rElementGeometry the geometry of the current element
   * @param rShapeFunctionsValues the shape functions values in the current integration point
   */
  virtual void ResetCrossSection(const Properties& rMaterialProperties,
                                 const GeometryType& rElementGeometry,
                                 const Vector& rShapeFunctionsValues);

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

  /**
   * Computes the transformations matrix for shell generalized strains, given an orientation angle in radians.
   * @param radians the input angle in radians
   * @param T the output transformation matrix
   * @return
   */
  inline void GetRotationMatrixForGeneralizedStrains(double radians, Matrix & T)
  {
    double c = std::cos(radians);
    double s = std::sin(radians);

    SizeType strain_size = GetStrainSize();

    if(T.size1() != strain_size || T.size2() != strain_size)
      T.resize(strain_size, strain_size, false);
    noalias( T ) = ZeroMatrix(strain_size, strain_size);

    T(0, 0) = c * c;			T(0, 1) =   s * s;				T(0, 2) = - s * c;
    T(1, 0) = s * s;			T(1, 1) =   c * c;				T(1, 2) =   s * c;
    T(2, 0) = 2.0 * s * c;		T(2, 1) = - 2.0 * s * c;		T(2, 2) = c * c - s * s;

    project( T, range(3, 6), range(3, 6) ) = project( T, range(0, 3), range(0, 3) );

    if(strain_size == 8)
    {
      T(6, 6) =   c;		T(6, 7) = s;
      T(7, 6) = - s;		T(7, 7) = c;
    }
  }

  /**
   * Computes the transformations matrix for condensed strains, given an orientation angle in radians.
   * @param radians the input angle in radians
   * @param T the output transformation matrix
   * @return
   */
  inline void GetRotationMatrixForCondensedStrains(double radians, Matrix & T)
  {
    SizeType strain_size = GetCondensedStrainSize();

    if(T.size1() != strain_size || T.size2() != strain_size)
      T.resize(strain_size, strain_size, false);
    noalias( T ) = ZeroMatrix(strain_size, strain_size);

    T(0, 0) = 1.0; // condensed strain E.zz is always at index 0

    if(strain_size == 3) // if section is thin the condensed strains are (in order): E.zz E.yz E.xz
    {
      double c = std::cos(radians);
      double s = std::sin(radians);

      T(1, 1) =   c;		T(1, 2) = s;
      T(2, 1) = - s;		T(2, 2) = c;
    }
  }

  /**
   * Computes the transformations matrix for shell generalized stresses, given an orientation angle in radians.
   * @param radians the input angle in radians
   * @param T the output transformation matrix
   * @return
   */
  inline void GetRotationMatrixForGeneralizedStresses(double radians, Matrix & T)
  {
    double c = std::cos(radians);
    double s = std::sin(radians);

    SizeType strain_size = GetStrainSize();

    if(T.size1() != strain_size || T.size2() != strain_size)
      T.resize(strain_size, strain_size, false);
    noalias( T ) = ZeroMatrix(strain_size, strain_size);

    T(0, 0) = c * c;		T(0, 1) =   s * s;		T(0, 2) = - 2.0 * s * c;
    T(1, 0) = s * s;		T(1, 1) =   c * c;		T(1, 2) =   2.0 * s * c;
    T(2, 0) = s * c;		T(2, 1) = - s * c;		T(2, 2) = c * c - s * s;

    project( T, range(3, 6), range(3, 6) ) = project( T, range(0, 3), range(0, 3) );

    if(strain_size == 8)
    {
      T(6, 6) =   c;		T(6, 7) = s;
      T(7, 6) = - s;		T(7, 7) = c;
    }
  }

  /**
   * Computes the transformations matrix for condensed stresses, given an orientation angle in radians.
   * @param radians the input angle in radians
   * @param T the output transformation matrix
   * @return
   */
  inline void GetRotationMatrixForCondensedStresses(double radians, Matrix & T)
  {
    SizeType strain_size = GetCondensedStrainSize();

    if(T.size1() != strain_size || T.size2() != strain_size)
      T.resize(strain_size, strain_size, false);
    noalias( T ) = ZeroMatrix(strain_size, strain_size);

    T(0, 0) = 1.0; // condensed stresse S.zz is always at index 0

    if(strain_size == 3) // if section is thin the condensed stresses are (in order): S.zz S.yz S.xz
    {
      double c = std::cos(radians);
      double s = std::sin(radians);

      T(1, 1) =   c;		T(1, 2) = s;
      T(2, 1) = - s;		T(2, 2) = c;
    }
  }

  ///@}

 public:

  ///@name Public Access
  ///@{

  /**
   * Returns the total thickness of this cross section
   * @return the thickness
   */
  inline const double GetThickness()const
  {
    return mThickness;
  }

  /**
   * Returns the offset of this cross section with respect to the reference mid-surface
   * of the parent element.
   * The offset can be a positive or negative value, measured along the normal of the reference surface.
   * The default value is Zero (i.e. the center of the cross section coincides with the shell mid-surface).
   * @return the offset
   */
  inline const double GetOffset()const
  {
    return mOffset;
  }

  /**
   * Sets the offset of this cross section with respect to the reference mid-surface
   * of the parent element.
   * The offset can be a positive or negative value, measured along the normal of the reference surface.
   * The default value is Zero (i.e. the center of the cross section coincides with the shell mid-surface).
   * @param offset the offset
   */
  inline void SetOffset(double offset)
  {
    if((mOffset != offset) && (!mEditingStack))
    {
      for(PlyCollection::iterator it = mStack.begin(); it != mStack.end(); ++it)
        (*it).SetLocation((*it).GetLocation() + offset - mOffset);
      mOffset = offset;
    }
  }

  /**
   * Returns the number of plies of this cross section.
   * @return the number of plies
   */
  inline PlyCollection::size_type NumberOfPlies()const
  {
    return mStack.size();
  }

  /**
   * Returns the number of integration points in the specified ply
   * @param ply_id the 0-based index of the target ply
   * @return the number of integration points
   */
  inline SizeType NumberOfIntegrationPointsAt(SizeType ply_id)const
  {
    if(ply_id < mStack.size())
      return mStack[ply_id].NumberOfIntegrationPoints();
    return 0;
  }

  /**
   * Sets a constitutive law pointer to the specified location
   * @param ply_id the 0-based index of the target ply
   * @param point_id the 0-based index of the target integration point in the target ply
   */
  inline void SetConstitutiveLawAt(SizeType ply_id, SizeType point_id, const ConstitutiveLaw::Pointer& pNewConstitutiveLaw)
  {
    if(ply_id < mStack.size())
      mStack[ply_id].SetConstitutiveLawAt(point_id, pNewConstitutiveLaw);
  }

  /**
   * Calculates the mass per unit area of this cross section.
   * @return the mass per unit area
   */
  inline double CalculateMassPerUnitArea()const
  {
    double vol(0.0);
    for(PlyCollection::const_iterator it = mStack.begin(); it != mStack.end(); ++it)
      vol += (*it).CalculateMassPerUnitArea();
    return vol;
  }

  /**
   * Calculates the avarage mass density of this cross section.
   * @return the avarage mass density
   */
  inline double CalculateAvarageDensity()const
  {
    return CalculateMassPerUnitArea() / mThickness;
  }

  /**
   * Returns the orientation angle (in radians) of this cross section
   * with respect to the parent element.
   * @return the orientation angle in radians
   */
  inline double GetOrientationAngle()const
  {
    return mOrientation;
  }

  /**
   * Sets the orientation angle (in radians) of this cross section
   * with respect to the parent element.
   * @param radians the orientation angle in radians
   */
  inline void SetOrientationAngle(double radians)
  {
    mOrientation = radians;
  }

  /**
   * Returns the behavior of this cross section (thin/thick)
   * @return the section behavior
   */
  inline SectionBehaviorType GetSectionBehavior()const
  {
    return mBehavior;
  }

  /**
   * Sets the behavior of this cross section (thin/thick)
   * @param behavior the section behavior
   */
  inline void SetSectionBehavior(SectionBehaviorType behavior)
  {
    mBehavior = behavior;
  }

  /**
   * Returns the size of the generalized strain vector of this cross section,
   * 8 for thick sections and 6 for Thin sections
   * @return the generalized strain size
   */
  inline SizeType GetStrainSize()
  {
    return (mBehavior == Thick) ? 8 : 6;
  }

  /**
   * Returns the size of the condensed strain vector of this cross section,
   * 1 for thick sections and 3 for Thin sections
   * @return the generalized strain size
   */
  inline SizeType GetCondensedStrainSize()
  {
    return (mBehavior == Thick) ? 1 : 3;
  }

  /**
   * Returns the stiffness value to be used for the drilling part of the shell formulation
   * @return the drilling stiffness
   */
  inline double GetDrillingStiffness()const
  {
    return mDrillingPenalty;
  }

  ///@}

 private:

  ///@name Private Methods
  ///@{

  void InitializeParameters(Parameters& rValues, ConstitutiveLaw::Parameters& rMaterialValues, ElementVariables& rVariables);

  void UpdateIntegrationPointParameters(IntegrationPoint& rPoint, ConstitutiveLaw::Parameters& rMaterialValues, ElementVariables& rVariables);

  void CalculateIntegrationPointResponse(IntegrationPoint& rPoint,
                                         ConstitutiveLaw::Parameters& rMaterialValues,
                                         Parameters& rValues,
                                         ElementVariables& rVariables,
                                         const ConstitutiveLaw::StressMeasure& rStressMeasure);

  /**
   * Creates a deep copy of this cross section.
   * Note: all constitutive laws are properly cloned.
   * @param other the source cross section
   */
  void PrivateCopy(const ShellCrossSection & other);

  ///@}

 public:

  ///@name Private Methods
  ///@{

  ///@}

 private:

  ///@name Member Variables
  ///@{

  double mThickness;
  double mOffset;
  PlyCollection mStack;
  bool mEditingStack;
  bool mHasDrillingPenalty;
  double mDrillingPenalty;
  double mOrientation;
  SectionBehaviorType mBehavior;
  bool mInitialized;
  bool mNeedsOOPCondensation;
  Vector mOOP_CondensedStrains;
  Vector mOOP_CondensedStrains_converged;

  ///@}

  ///@name Serialization
  ///@{

  friend class Serializer;

  void save(Serializer& rSerializer) const override
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    rSerializer.save("th", mThickness);
    rSerializer.save("offs", mOffset);
    rSerializer.save("stack", mStack);
    rSerializer.save("edit", mEditingStack);
    rSerializer.save("dr", mHasDrillingPenalty);
    rSerializer.save("bdr", mDrillingPenalty);
    rSerializer.save("or", mOrientation);

    rSerializer.save("behav", (int)mBehavior);

    rSerializer.save("init", mInitialized);
    rSerializer.save("hasOOP", mNeedsOOPCondensation);
    rSerializer.save("OOP_eps", mOOP_CondensedStrains_converged);
  }

  void load(Serializer& rSerializer) override
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    rSerializer.load("th", mThickness);
    rSerializer.load("offs", mOffset);
    rSerializer.load("stack", mStack);
    rSerializer.load("edit", mEditingStack);
    rSerializer.load("dr", mHasDrillingPenalty);
    rSerializer.load("bdr", mDrillingPenalty);
    rSerializer.load("or", mOrientation);

    int temp;
    rSerializer.load("behav", temp);
    mBehavior = (SectionBehaviorType)temp;

    rSerializer.load("init", mInitialized);
    rSerializer.load("hasOOP", mNeedsOOPCondensation);
    rSerializer.load("OOP_eps", mOOP_CondensedStrains_converged);
  }

  ///@}

 public:

  DECLARE_ADD_THIS_TYPE_TO_PROPERTIES
  DECLARE_GET_THIS_TYPE_FROM_PROPERTIES

};

///@name Input/Output funcitons
///@{

inline std::istream & operator >> (std::istream & rIStream, ShellCrossSection & rThis);

inline std::ostream & operator << (std::ostream & rOStream, const ShellCrossSection & rThis)
{
  return rOStream << rThis.GetInfo();
}

///@}

}


#endif // SHELL_CROSS_SECTION_H_INCLUDED
