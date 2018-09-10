// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//                   Philipp Bucher
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
#include "shell_utilities.h"
#include "containers/flags.h"

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ShellCrossSection : public Flags
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


	/** \brief SectionParameters
	*
	* SectionParameters is an accessibility class for shells using the
	* ShellCrossSection class. It allows one to set and get vectors and matrices
	* associated with the shell cross section, such as strains, stresses and the
	* constitutive matrix.
	*
	* An example application is taken from shell_thick_3D4N.cpp, before it's
	* stiffness matrix gauss loop is entered:
	*
	* ShellCrossSection::SectionParameters parameters(geom, props, rCurrentProcessInfo);
	* parameters.SetGeneralizedStrainVector( generalizedStrains );
	* parameters.SetGeneralizedStressVector( generalizedStresses );
	* parameters.SetConstitutiveMatrix( D );
	* Flags& options = parameters.GetOptions();
	* options.Set(ConstitutiveLaw::COMPUTE_STRESS, RHSrequired);
	* options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, LHSrequired);
	*/
    class SectionParameters
    {

    private:

        Flags                mOptions;

        Vector*              mpGeneralizedStrainVector;
        Vector*              mpGeneralizedStressVector;
        Matrix*              mpConstitutiveMatrix;

		double				 mStenbergShearStabilization = 1.0;
		// refer https://doi.org/10.1016/j.cma.2003.12.036 section 3.1

        const Vector*        mpShapeFunctionsValues;
        const Matrix*        mpShapeFunctionsDerivatives;
        const ProcessInfo*   mpCurrentProcessInfo;
        const Properties*    mpMaterialProperties;
        const GeometryType*  mpElementGeometry;

    public:

        SectionParameters()
            : mpGeneralizedStrainVector(nullptr)
            , mpGeneralizedStressVector(nullptr)
            , mpConstitutiveMatrix(nullptr)
            , mpShapeFunctionsValues(nullptr)
            , mpShapeFunctionsDerivatives(nullptr)
            , mpCurrentProcessInfo(nullptr)
            , mpMaterialProperties(nullptr)
            , mpElementGeometry(nullptr)
        {}

        SectionParameters (const GeometryType& rElementGeometry,
                    const Properties& rMaterialProperties,
                    const ProcessInfo& rCurrentProcessInfo)
            : mpGeneralizedStrainVector(nullptr)
            , mpGeneralizedStressVector(nullptr)
            , mpConstitutiveMatrix(nullptr)
            , mpShapeFunctionsValues(nullptr)
            , mpShapeFunctionsDerivatives(nullptr)
            , mpCurrentProcessInfo(&rCurrentProcessInfo)
            , mpMaterialProperties(&rMaterialProperties)
            , mpElementGeometry(&rElementGeometry)
        {}

        SectionParameters (const SectionParameters & rNewParameters)
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

        void Set                             (Flags ThisFlag)
        {
            mOptions.Set(ThisFlag);
        };
        void Reset                           (Flags ThisFlag)
        {
            mOptions.Reset(ThisFlag);
        };

        void SetOptions                      (const Flags&  rOptions)
        {
            mOptions=rOptions;
        };

        void SetGeneralizedStrainVector      (Vector& rGeneralizedStrainVector)
        {
            mpGeneralizedStrainVector=&rGeneralizedStrainVector;
        };
        void SetGeneralizedStressVector      (Vector& rGeneralizedStressVector)
        {
            mpGeneralizedStressVector=&rGeneralizedStressVector;
        };
        void SetConstitutiveMatrix           (Matrix& rConstitutiveMatrix)
        {
            mpConstitutiveMatrix =&rConstitutiveMatrix;
        };

        void SetShapeFunctionsValues         (const Vector& rShapeFunctionsValues)
        {
            mpShapeFunctionsValues=&rShapeFunctionsValues;
        };
        void SetShapeFunctionsDerivatives    (const Matrix& rShapeFunctionsDerivatives)
        {
            mpShapeFunctionsDerivatives=&rShapeFunctionsDerivatives;
        };
        void SetProcessInfo                  (const ProcessInfo& rProcessInfo)
        {
            mpCurrentProcessInfo =&rProcessInfo;
        };
        void SetMaterialProperties           (const Properties&  rMaterialProperties)
        {
            mpMaterialProperties =&rMaterialProperties;
        };
        void SetElementGeometry              (const GeometryType& rElementGeometry)
        {
            mpElementGeometry =&rElementGeometry;
        };
		void SetStenbergShearStabilization(const double& StenbergShearStabilization)
		{
			mStenbergShearStabilization = StenbergShearStabilization;
		};

        /**
        * returns the reference or the value of a specified variable: returns the value of the parameter, only non const values can be modified
        */

        Flags& GetOptions ()
        {
            return mOptions;
        };

        Vector& GetGeneralizedStrainVector         ()
        {
            return *mpGeneralizedStrainVector;
        };
        Vector& GetGeneralizedStressVector         ()
        {
            return *mpGeneralizedStressVector;
        };
        Matrix& GetConstitutiveMatrix              ()
        {
            return *mpConstitutiveMatrix;
        };

        const Vector& GetShapeFunctionsValues      ()
        {
            return *mpShapeFunctionsValues;
        };
        const Matrix& GetShapeFunctionsDerivatives ()
        {
            return *mpShapeFunctionsDerivatives;
        };
        const ProcessInfo&  GetProcessInfo         ()
        {
            return *mpCurrentProcessInfo;
        };
        const Properties&   GetMaterialProperties  ()
        {
            return *mpMaterialProperties;
        };
        const GeometryType& GetElementGeometry     ()
        {
            return *mpElementGeometry;
        };
		double GetStenbergShearStabilization()
		{
			return mStenbergShearStabilization;
		};
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

        IntegrationPoint & operator = (const IntegrationPoint & other)
        {
            if(this != &other)
            {
                mWeight = other.mWeight;
                mLocation = other.mLocation;
                mConstitutiveLaw = other.mConstitutiveLaw != NULL ? other.mConstitutiveLaw->Clone() : ConstitutiveLaw::Pointer();
            }
            return *this;
        }

    public:

        inline double GetWeight()const
        {
            return mWeight;
        }
        inline void SetWeight(double w)
        {
            mWeight = w;
        }

        inline double GetLocation()const
        {
            return mLocation;
        }
        inline void SetLocation(double l)
        {
            mLocation = l;
        }

        inline const ConstitutiveLaw::Pointer& GetConstitutiveLaw()const
        {
            return mConstitutiveLaw;
        }
        inline void SetConstitutiveLaw(const ConstitutiveLaw::Pointer& pLaw)
        {
            mConstitutiveLaw = pLaw;
        }

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

        int mPlyIndex;
        IntegrationPointCollection mIntegrationPoints;

    public:

        Ply()
            : mPlyIndex(0)
            , mIntegrationPoints()
        {}

        Ply(const int PlyIndex, int NumIntegrationPoints, const Properties& rProps)
            : mPlyIndex(PlyIndex)
            , mIntegrationPoints()
        {
            // make sure the number is greater than 0 and odd
            KRATOS_ERROR_IF(NumIntegrationPoints < 1) << "Number of Integration points must be larger than 0!" << std::endl;
            if(NumIntegrationPoints < 0) NumIntegrationPoints = -NumIntegrationPoints;
            if(NumIntegrationPoints == 0) NumIntegrationPoints = 5;
            if(NumIntegrationPoints % 2 == 0) NumIntegrationPoints += 1;
            InitializeIntegrationPoints(rProps, NumIntegrationPoints);
        }

        Ply(const Ply& other)
            : mPlyIndex(other.mPlyIndex)
            , mIntegrationPoints(other.mIntegrationPoints)
        {}

        Ply & operator = (const Ply & other)
        {
            if(this != &other)
            {
                mPlyIndex = other.mPlyIndex;
                mIntegrationPoints = other.mIntegrationPoints;
            }
            return *this;
        }

    public:

        inline double GetThickness(const Properties& rProps) const
        {
            return ShellUtilities::GetThickness(rProps, mPlyIndex);
        }

        inline double GetLocation(const Properties& rProps) const
        {
            double my_location(0.0);

            double current_location = ShellUtilities::GetThickness(rProps) * 0.5;
            const double offset = GetOffset(rProps);

            for (int i=0; i<mPlyIndex+1; ++i)
            {
                double ply_thickness = GetThickness(rProps);
                my_location = current_location - ply_thickness*0.5 - offset;
                current_location -= ply_thickness;
            }
            return my_location;
        }

        inline double GetOrientationAngle(const Properties& rProps) const
        {
            return ShellUtilities::GetOrientationAngle(rProps, mPlyIndex);
        }

        inline double GetOffset(const Properties& rProps) const
        {
            return ShellUtilities::GetOffset(rProps);
        }

		void RecoverOrthotropicProperties(const IndexType currentPly, Properties& laminaProps);

        inline IntegrationPointCollection& GetIntegrationPoints(const Properties& rProps)
        {
            UpdateIntegrationPoints(rProps);
            return mIntegrationPoints;
        }

        inline double CalculateMassPerUnitArea(const Properties& rProps) const
        {
            return ShellUtilities::GetDensity(rProps, mPlyIndex) * GetThickness(rProps);
        }

        inline IntegrationPointCollection::size_type NumberOfIntegrationPoints() const
        {
            return mIntegrationPoints.size();
        }

        inline void SetConstitutiveLawAt(IntegrationPointCollection::size_type integrationPointID, const ConstitutiveLaw::Pointer& pNewConstitutiveLaw)
        {
            if(integrationPointID < mIntegrationPoints.size())
                mIntegrationPoints[integrationPointID].SetConstitutiveLaw(pNewConstitutiveLaw);
        }

    private:

        void InitializeIntegrationPoints(const Properties& rProps, const int NumIntegrationPoints)
        {
            KRATOS_TRY

            const ConstitutiveLaw::Pointer& pMaterial = rProps[CONSTITUTIVE_LAW];
            KRATOS_ERROR_IF(pMaterial == nullptr) << "A Ply needs a constitutive law to be set. "
                << "Missing constitutive law in property: " <<  rProps.Id() << std::endl;;

            // generate the integration points
            mIntegrationPoints.clear();
            mIntegrationPoints.resize(NumIntegrationPoints);
            for(int i=0; i<NumIntegrationPoints; ++i)
                mIntegrationPoints[i].SetConstitutiveLaw(pMaterial->Clone());

            KRATOS_CATCH("")
        }
        void UpdateIntegrationPoints(const Properties& rProps)
        {
            KRATOS_TRY

            const SizeType num_int_points = mIntegrationPoints.size();

            // generate the weights (composite simpson rule)
            Vector ip_w(num_int_points, 1.0);
            if (num_int_points >= 3)
            {
                for (IndexType i=1; i<num_int_points-1; ++i)
                {
                    double iw = (i % 2 == 0) ? 2.0 : 4.0;
                    ip_w(i) = iw;
                }
                ip_w /= sum( ip_w );
            }

            // generate locations (direction: top(+thickness/2) to bottom(-thickness/2)
            const double location = GetLocation(rProps);
            const double thickness = GetThickness(rProps);

            Vector ip_loc(num_int_points, 0.0);
            if (num_int_points >= 3)
            {
                double loc_start = location + 0.5 * thickness;
                double loc_incr = thickness / double(num_int_points-1);
                for (IndexType i=0; i<num_int_points; ++i)
                {
                    ip_loc(i) = loc_start;
                    loc_start -= loc_incr;
                }
            }

            for (IndexType i=0; i<num_int_points; ++i)
            {
                IntegrationPoint& r_int_point = mIntegrationPoints[i];
                r_int_point.SetWeight(ip_w(i) * thickness);
                r_int_point.SetLocation(ip_loc(i));
            }

            KRATOS_CATCH("")
        }

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const
        {
            rSerializer.save("idx", mPlyIndex);
            rSerializer.save("IntP", mIntegrationPoints);
        }

        virtual void load(Serializer& rSerializer)
        {
            rSerializer.load("idx", mPlyIndex);
            rSerializer.load("IntP", mIntegrationPoints);
        }

    };

protected:

    struct GeneralVariables
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
    void AddPly(const IndexType PlyIndex, int numPoints, const Properties& rProps);

    /**
    * Finalizes the editing of the Composite Layup.
    */
    void EndStack();

    /**
    * Returns the string containing a detailed description of this object.
    * @return the string with informations
    */
    virtual std::string GetInfo(const Properties& rProps);

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
    virtual double& GetValue(const Variable<double>& rThisVariable, const Properties& rProps, double& rValue);

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
    virtual void CalculateSectionResponse(SectionParameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure);

    /**
    * Updates the section response, called by the element in FinalizeSolutionStep.
    * @param rValues the parameters for the current calculation
    * @param rStressMeasure the required stress measure
    * @see Parameters
    */
    virtual void FinalizeSectionResponse(SectionParameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure);

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

        T(0, 0) = c * c;
        T(0, 1) =   s * s;
        T(0, 2) = - s * c;
        T(1, 0) = s * s;
        T(1, 1) =   c * c;
        T(1, 2) =   s * c;
        T(2, 0) = 2.0 * s * c;
        T(2, 1) = - 2.0 * s * c;
        T(2, 2) = c * c - s * s;

        project( T, range(3, 6), range(3, 6) ) = project( T, range(0, 3), range(0, 3) );

        if(strain_size == 8)
        {
            T(6, 6) =   c;
            T(6, 7) = s;
            T(7, 6) = - s;
            T(7, 7) = c;
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

            T(1, 1) =   c;
            T(1, 2) = s;
            T(2, 1) = - s;
            T(2, 2) = c;
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

        T(0, 0) = c * c;
        T(0, 1) =   s * s;
        T(0, 2) = - 2.0 * s * c;
        T(1, 0) = s * s;
        T(1, 1) =   c * c;
        T(1, 2) =   2.0 * s * c;
        T(2, 0) = s * c;
        T(2, 1) = - s * c;
        T(2, 2) = c * c - s * s;

        project( T, range(3, 6), range(3, 6) ) = project( T, range(0, 3), range(0, 3) );

        if(strain_size == 8)
        {
            T(6, 6) =   c;
            T(6, 7) = s;
            T(7, 6) = - s;
            T(7, 7) = c;
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

            T(1, 1) =   c;
            T(1, 2) = s;
            T(2, 1) = - s;
            T(2, 2) = c;
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
    inline double GetThickness(const Properties& rProps) const
    {
        double thickness = 0.0;
    	for (const auto& r_ply : mStack)
    		thickness += r_ply.GetThickness(rProps);
        return thickness;
    }

    /**
    * Returns the offset of this cross section with respect to the reference mid-surface
    * of the parent element.
    * The offset can be a positive or negative value, measured along the normal of the reference surface.
    * The default value is Zero (i.e. the center of the cross section coincides with the shell mid-surface).
    * @return the offset
    */
    inline double GetOffset(const Properties& rProps) const
    {
        KRATOS_DEBUG_ERROR_IF(mStack.size() == 0) << "no plies available!" << std::endl;
        return mStack[0].GetOffset(rProps);
    }

    /**
    * Stores the thicknesses of plies of this cross section.
    */
    void GetPlyThicknesses(const Properties& rProps, Vector& rPlyThicknesses)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mStack.size() == rPlyThicknesses.size()) << "Size mismatch!" << std::endl;
        for (IndexType i_ply=0; i_ply<mStack.size(); ++i_ply)
    		rPlyThicknesses[i_ply] = mStack[i_ply].GetThickness(rProps);
    }

    /**
    * Setup to get the integrated constitutive matrices for each ply
    */
    void SetupGetPlyConstitutiveMatrices()
    {
		// This function must be called before requesting un-integrated
		// constitutive matrices for each ply!
    	mStorePlyConstitutiveMatrices = true;
    	mPlyConstitutiveMatrices = std::vector<Matrix>(this->NumberOfPlies());

    	for (IndexType ply = 0; ply < this->NumberOfPlies(); ++ply)
    	{
    		if (mBehavior == Thick)
                mPlyConstitutiveMatrices[ply].resize(8, 8, false);
    		else
                mPlyConstitutiveMatrices[ply].resize(6, 6, false);

    		mPlyConstitutiveMatrices[ply].clear();
    	}
    }

    /**
    * Get the integrated constitutive matrices for each ply
    */
    Matrix GetPlyConstitutiveMatrix(const IndexType PlyIndex)
    {
    	return mPlyConstitutiveMatrices[PlyIndex];
    }

    /**
    * Returns the number of plies of this cross section.
    * @return the number of plies
    */
    inline SizeType NumberOfPlies() const
    {
        return mStack.size();
    }

    /**
    * Returns the number of integration points in the specified ply
    * @param PlyIndex the 0-based index of the target ply
    * @return the number of integration points
    */
    inline SizeType NumberOfIntegrationPointsAt(const IndexType PlyIndex) const
    {
        if(PlyIndex < mStack.size())
            return mStack[PlyIndex].NumberOfIntegrationPoints();
        return 0;
    }

    /**
    * Sets a constitutive law pointer to the specified location
    * @param PlyIndex the 0-based index of the target ply
    * @param point_id the 0-based index of the target integration point in the target ply
    */
    inline void SetConstitutiveLawAt(const IndexType PlyIndex, SizeType point_id, const ConstitutiveLaw::Pointer& pNewConstitutiveLaw)
    {
        if(PlyIndex < mStack.size())
            mStack[PlyIndex].SetConstitutiveLawAt(point_id, pNewConstitutiveLaw);
    }

    /**
    * Calculates the mass per unit area of this cross section.
    * @return the mass per unit area
    */
    inline double CalculateMassPerUnitArea(const Properties& rProps) const
    {
        double vol(0.0);
        for (const auto& r_ply : mStack)
            vol += r_ply.CalculateMassPerUnitArea(rProps);
        return vol;
    }

    /**
    * Calculates the avarage mass density of this cross section.
    * @return the avarage mass density
    */
    inline double CalculateAvarageDensity(const Properties& rProps) const
    {
        return CalculateMassPerUnitArea(rProps) / GetThickness(rProps);
    }

    /**
    * Returns the orientation angle (in radians) of this cross section
    * with respect to the parent element.
    * @return the orientation angle in radians
    */
    inline double GetOrientationAngle() const
    {
        return mOrientation;
    }

    /**
    * Sets the orientation angle (in radians) of this cross section
    * with respect to the parent element.
    * @param radians the orientation angle in radians
    */
    inline void SetOrientationAngle(const double Radians)
    {
        mOrientation = Radians;
    }

    /**
    * Returns the behavior of this cross section (thin/thick)
    * @return the section behavior
    */
    inline SectionBehaviorType GetSectionBehavior() const
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
    inline double GetDrillingStiffness() const
    {
        return mDrillingPenalty;
    }

    /**
    * Parses the shell orthotropic material data from properties
    */
    void ParseOrthotropicPropertyMatrix(const Properties& pProps);

    /**
    * Get orientation of laminae
    */
    void GetLaminaeOrientation(const Properties& pProps, Vector& rOrientation_Vector);

    /**
    * Get strengths of laminae
    */
    void GetLaminaeStrengths(std::vector<Matrix>& rLamina_Strengths, const Properties& rProps);
    ///@}

private:

    ///@name Private Methods
    ///@{

    void InitializeParameters(SectionParameters& rValues, ConstitutiveLaw::Parameters& rMaterialValues, GeneralVariables& rVariables);

    void UpdateIntegrationPointParameters(const IntegrationPoint& rPoint, ConstitutiveLaw::Parameters& rMaterialValues, GeneralVariables& rVariables);

    void CalculateIntegrationPointResponse(const IntegrationPoint& rPoint,
    	ConstitutiveLaw::Parameters& rMaterialValues,
    	SectionParameters& rValues,
    	GeneralVariables& rVariables,
    	const ConstitutiveLaw::StressMeasure& rStressMeasure,
    	const unsigned int& plyNumber);

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
    bool mStorePlyConstitutiveMatrices = false;
    std::vector<Matrix> mPlyConstitutiveMatrices;

    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("stack", mStack);
        rSerializer.save("edit", mEditingStack);
        rSerializer.save("dr", mHasDrillingPenalty);
        rSerializer.save("bdr", mDrillingPenalty);
        rSerializer.save("or", mOrientation);

        rSerializer.save("behav", (int)mBehavior);

        rSerializer.save("init", mInitialized);
        rSerializer.save("hasOOP", mNeedsOOPCondensation);
        rSerializer.save("OOP_eps", mOOP_CondensedStrains);
        rSerializer.save("OOP_eps_conv", mOOP_CondensedStrains_converged);
        rSerializer.save("store_ply_mat", mStorePlyConstitutiveMatrices);
        rSerializer.save("ply_mat", mPlyConstitutiveMatrices);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
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
        rSerializer.load("OOP_eps", mOOP_CondensedStrains);
        rSerializer.load("OOP_eps_conv", mOOP_CondensedStrains_converged);
        rSerializer.load("store_ply_mat", mStorePlyConstitutiveMatrices);
        rSerializer.load("ply_mat", mPlyConstitutiveMatrices);
    }

    ///@}

};

///@name Input/Output funcitons
///@{

inline std::istream & operator >> (std::istream & rIStream, ShellCrossSection & rThis);

inline std::ostream & operator << (std::ostream & rOStream, ShellCrossSection & rThis)
{
    return rOStream; // << rThis.GetInfo();
}

///@}

}


#endif // SHELL_CROSS_SECTION_H_INCLUDED