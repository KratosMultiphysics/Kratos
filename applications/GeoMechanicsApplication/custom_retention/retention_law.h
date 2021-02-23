// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_RETENTION_LAW)
#define KRATOS_RETENTION_LAW

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/process_info.h"

namespace Kratos
{

    /**
 * Base class of retention laws.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) RetentionLaw
{
public:
    /**
     * Type definitions
     * NOTE: geometries are assumed to be of type Node<3> for all problems
     */
    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node<3>> GeometryType;

    /**
     * Counted pointer of RetentionLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(RetentionLaw);

    /**
     * Flags related to the Parameters of the Contitutive Law
     */
    // KRATOS_DEFINE_LOCAL_FLAG( USE_ELEMENT_PROVIDED_DATA );
    // KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_SATURATION );

    class Parameters
    {
        KRATOS_CLASS_POINTER_DEFINITION(Parameters);

        /**
         * Structure "Parameters" to be used by the element to pass the parameters into the retention law *

        * KINEMATIC PARAMETERS:

        *** NOTE: Pointers are used only to point to a certain variable, 
        *   no "new" or "malloc" can be used for this Parameters ***

        * @param mVolumetricStrain copy of the determinant of the Current DeformationGradient (although Current F  is also included as a matrix) (input data)
        * @param mMeanStress pointer to the current stresses (*OUTPUT with COMPUTE_STRESS flag)
        * @param mpConstitutiveMatrix pointer to the material tangent matrix (*OUTPUT with COMPUTE_CONSTITUTIVE_TENSOR flag)

        * GEOMETRIC PARAMETERS:
        * @param mpElementGeometry pointer to the element's geometry (input data)

        * MATERIAL PROPERTIES:
        * @param mpMaterialProperties pointer to the material's Properties object (input data)

        * PROCESS PROPERTIES:
        * @param mpCurrentProcessInfo pointer to current ProcessInfo instance (input data)

        */

    private:
        /*** NOTE: Member Pointers are used only to point to a certain variable, 
         * no "new" or "malloc" can be used for this Parameters ***/

        double mFluidPressure;
        double mMeanStress;
        double mTemperature;
        double mVolumetricStrain;

        const ProcessInfo *mpCurrentProcessInfo;
        const Properties *mpMaterialProperties;
        const GeometryType *mpElementGeometry;

    public:
        /**
         * Constructor.
         */
        Parameters()
        {
            //Initialize pointers to NULL
            mFluidPressure = 0;
            mMeanStress = 0;
            mTemperature = 0;
            mVolumetricStrain = 0;
            mpCurrentProcessInfo = NULL;
            mpMaterialProperties = NULL;
            mpElementGeometry = NULL;
        };

        /**
         * Constructor with Properties, Geometry and ProcessInfo
         */
        Parameters(const GeometryType &rElementGeometry,
                   const Properties &rMaterialProperties,
                   const ProcessInfo &rCurrentProcessInfo) : mpCurrentProcessInfo(&rCurrentProcessInfo)
                                                           , mpMaterialProperties(&rMaterialProperties)
                                                           , mpElementGeometry(&rElementGeometry)
        {
            mFluidPressure = 0;
            mMeanStress = 0;
            mTemperature = 0;
            mVolumetricStrain = 0;
        };

        /**
         * Copy Constructor.
         */
        Parameters(const Parameters &rNewParameters) : mFluidPressure(rNewParameters.mFluidPressure)
                                                     , mMeanStress(rNewParameters.mMeanStress)
                                                     , mTemperature(rNewParameters.mTemperature)
                                                     , mVolumetricStrain(rNewParameters.mVolumetricStrain)
                                                     , mpCurrentProcessInfo(rNewParameters.mpCurrentProcessInfo)
                                                     , mpMaterialProperties(rNewParameters.mpMaterialProperties)
                                                     , mpElementGeometry(rNewParameters.mpElementGeometry){};

        /**
         * Destructor.
         */
        ~Parameters()
        {
        }

        /**
         * Public Methods to access variables of the struct class
         */

        /**
         * sets the variable or the pointer of a specified variable: assigns the direction of the pointer for the mpvariables, only non const values can be modified
         */

        void SetVolumetricStrain(const double rVolumetricStrain) { mVolumetricStrain = rVolumetricStrain; };
        void SetMeanStress      (const double rMeanStress)       { mMeanStress = rMeanStress; };
        void SetFluidPressure   (const double rFluidPressure)    { mFluidPressure = rFluidPressure; };
        void SetTemperature     (const double rTemperature)      { mTemperature = rTemperature; };

        void SetProcessInfo(const ProcessInfo &rProcessInfo) { mpCurrentProcessInfo = &rProcessInfo; };
        void SetMaterialProperties(const Properties &rMaterialProperties) { mpMaterialProperties = &rMaterialProperties; };
        void SetElementGeometry(const GeometryType &rElementGeometry) { mpElementGeometry = &rElementGeometry; };

        /**
         * Returns the reference or the value of a specified variable: returns the value of the parameter, only non const values can be modified
         */

        const double &GetVolumetricStrain() { return mVolumetricStrain; }
        const double &GetMeanStress()       { return mMeanStress;       }
        const double &GetFluidPressure()    { return mFluidPressure;    }
        const double &GetTemperature()      { return mTemperature;      }

        const ProcessInfo &GetProcessInfo()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetProcessInfo()) << "ProcessInfo is not set!" << std::endl;
            return *mpCurrentProcessInfo;
        }
        const Properties &GetMaterialProperties()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetMaterialProperties()) << "MaterialProperties is not set!" << std::endl;
            return *mpMaterialProperties;
        }
        const GeometryType &GetElementGeometry()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetElementGeometry()) << "ElementGeometry is not set!" << std::endl;
            return *mpElementGeometry;
        }

        /**
         * Returns the reference to the value of a specified variable with not constant access
         */

        double &GetVolumetricStrain(double &rVolumetricStrain)
        {
            rVolumetricStrain = mVolumetricStrain;
            return rVolumetricStrain;
        }
        double &GetMeanStress(double &rMeanStress)
        {
            rMeanStress = mMeanStress;
            return rMeanStress;
        }
        double &GetFluidPressure(double &rFluidPressure)
        {
            rFluidPressure = mFluidPressure;
            return rFluidPressure;
        }
        double &GetTemperature(double &rTemperature)
        {
            rTemperature = mTemperature;
            return rTemperature;
        }

        /**
         * Returns if the different components has been set
         */
        bool IsSetProcessInfo()        { return (mpCurrentProcessInfo != NULL); };
        bool IsSetMaterialProperties() { return (mpMaterialProperties != NULL); };
        bool IsSetElementGeometry()    { return (mpElementGeometry != NULL); };

    }; // class Parameters end

    /**
     * Constructor.
     */
    RetentionLaw();

    /**
     * Destructor.
     */
    virtual ~RetentionLaw();

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this retention law
     * @note implementation scheme:
     *      RetentionLaw::Pointer p_clone(new RetentionLaw());
     *      return p_clone;
     */
    virtual RetentionLaw::Pointer Clone() const;

    /**
     * @brief Returns whether this retention Law has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     */
    virtual bool Has(const Variable<bool> &rThisVariable);

    /**
     * @brief Returns whether this retention Law has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     */
    virtual bool Has(const Variable<int> &rThisVariable);

    /**
     * @brief Returns whether this retention Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     */
    virtual bool Has(const Variable<double> &rThisVariable);

    /**
     * @brief Returns whether this retention Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     */
    virtual bool Has(const Variable<Vector> &rThisVariable);

    /**
     * @brief Returns whether this retention Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     */
    virtual bool Has(const Variable<Matrix> &rThisVariable);

    /**
     * @brief Returns whether this retention Law has specified variable (array of 3 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     * @note Fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 3>> &rThisVariable);

    /**
     * @brief Returns whether this retention Law has specified variable (array of 6 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the retention law
     * @note Fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 6>> &rThisVariable);

    /**
     * @brief Returns the value of a specified variable (boolean)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    virtual bool &GetValue(const Variable<bool> &rThisVariable, bool &rValue);

    /**
     * Returns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    virtual int &GetValue(const Variable<int> &rThisVariable, int &rValue);

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    virtual double &GetValue(const Variable<double> &rThisVariable, double &rValue);

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    virtual Vector &GetValue(const Variable<Vector> &rThisVariable, Vector &rValue);

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @return rValue output: the value of the specified variable
     */
    virtual Matrix &GetValue(const Variable<Matrix> &rThisVariable, Matrix &rValue);

    /**
     * @brief Returns the value of a specified variable (array of 3 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    virtual array_1d<double, 3> &GetValue(const Variable<array_1d<double, 3>> &rThisVariable,
                                          array_1d<double, 3> &rValue);

    /**
     * @brief Returns the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 6> &GetValue(const Variable<array_1d<double, 6>> &rThisVariable,
                                          array_1d<double, 6> &rValue);

    /**
     * @brief Sets the value of a specified variable (boolean)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<bool> &rVariable,
                          const bool &Value,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<int> &rVariable,
                          const int &Value,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<double> &rVariable,
                          const double &rValue,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<Vector> &rVariable,
                          const Vector &rValue,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Sets the value of a specified variable (Matrix)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<Matrix> &rVariable,
                          const Matrix &rValue,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Sets the value of a specified variable (array of 3 components)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<array_1d<double, 3>> &rVariable,
                          const array_1d<double, 3> &rValue,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Sets the value of a specified variable (array of 6 components)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<array_1d<double, 6>> &rVariable,
                          const array_1d<double, 6> &rValue,
                          const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief Calculates the value of a specified variable (bool)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual bool &CalculateValue(Parameters &rParameters,
                                    const Variable<bool> &rThisVariable,
                                    bool &rValue);

    /**
     * @brief Calculates the value of a specified variable (int)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual int &CalculateValue(Parameters &rParameters,
                                const Variable<int> &rThisVariable,
                                int &rValue);

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual double &CalculateValue(Parameters &rParameters,
                                   const Variable<double> &rThisVariable,
                                   double &rValue);

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual Vector &CalculateValue(Parameters &rParameters,
                                   const Variable<Vector> &rThisVariable,
                                   Vector &rValue);

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual Matrix &CalculateValue(Parameters &rParameters,
                                    const Variable<Matrix> &rThisVariable,
                                    Matrix &rValue);

    /**
     * @brief Calculates the value of a specified variable (array of 3 components)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual array_1d<double, 3> &CalculateValue(Parameters &rParameters,
                                                const Variable<array_1d<double, 3>> &rVariable,
                                                array_1d<double, 3> &rValue);

    /**
     * returns the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 6> &CalculateValue(Parameters &rParameters,
                                                const Variable<array_1d<double, 6>> &rVariable,
                                                array_1d<double, 6> &rValue);

    virtual double CalculateSaturation(Parameters &rParameters);

    virtual double CalculateEffectiveSaturation(Parameters &rParameters);

    virtual double CalculateDerivativeOfSaturation(Parameters &rParameters);

    virtual double CalculateRelativePermeability(Parameters &rParameters);

    virtual double CalculateBishopCoefficient(Parameters &rParameters);

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the retention law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rCurrentProcessInfo process info
     */
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);

    virtual void Initialize(Parameters &rParameters);

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     */
    virtual void InitializeSolutionStep(Parameters &rParameters);

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     */
    virtual void FinalizeSolutionStep(Parameters &rParameters);

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    virtual void Finalize(Parameters &rParameters);

    /**
     * This can be used in order to reset all internal variables of the
     * retention law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void ResetMaterial(const Properties &rMaterialProperties,
                               const GeometryType &rElementGeometry,
                               const Vector &rShapeFunctionsValues);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    virtual int Check(const Properties &rMaterialProperties,
                      const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief This method is used to check that two Retention Laws are the same type (references)
     * @param rLHS The first argument
     * @param rRHS The second argument
     */
    inline static bool HasSameType(const RetentionLaw &rLHS, const RetentionLaw &rRHS)
    {
        return (typeid(rLHS) == typeid(rRHS));
    }

    /**
     * @brief This method is used to check that tow Retention Laws are the same type (pointers)
     * @param rLHS The first argument
     * @param rRHS The second argument
     */
    inline static bool HasSameType(const RetentionLaw *rLHS, const RetentionLaw *rRHS)
    {
        return RetentionLaw::HasSameType(*rLHS, *rRHS);
    }

    ///@}
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
        buffer << "RetentionLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "RetentionLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream << "RetentionLaw has no data";
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
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

    virtual void save(Serializer &rSerializer) const;

    virtual void load(Serializer &rSerializer);

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class RetentionLaw */

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                RetentionLaw &rThis);

/// output stream function

inline std::ostream &operator<<(std::ostream &rOStream,
                                const RetentionLaw &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block


} /* namespace Kratos.*/
#endif /* KRATOS_RETENTION_LAW  defined */
