//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Jordi Cotela
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

#if !defined(KRATOS_PFEM_SOLID_CONSTITUTIVE_LAW)
#define KRATOS_PFEM_SOLID_CONSTITUTIVE_LAW

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// This class contains the common infrastructure for solid constitutive laws.
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) PfemSolidConstitutiveLaw : public ConstitutiveLaw {
   public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PfemSolidConstitutiveLaw);

    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PfemSolidConstitutiveLaw();

    /// Copy constructor.
    PfemSolidConstitutiveLaw(const PfemSolidConstitutiveLaw& rOther);

    /// Destructor
    ~PfemSolidConstitutiveLaw() override;

    ///@}
    ///@name Operations
    ///@{

    /// Initialize a new instance of this type of law
    ConstitutiveLaw::Pointer Clone() const override;

    /// Calculate the response of the material for the current strain rate.
    /** This is the main method for solid constitutive laws.
     *  @note Besides computing the response, derived constitutive laws are responsible for setting the
     *  effective young modulus and poisson ratio for the material, which will be used by the element.
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /// Validate the data received by the constitutive law
    /** @param rMaterialProperties Properties of the parent Element.
     *  @param rElementGeometry Geometry of the parent Element.
     *  @param rCurrentProcessInfo ProcessInfo for the problem.
     *  @return 0 if everything is fine, other values indicate problems.
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    /// This always returns this->GetEffectiveDensity() or this->GetEffectiveYoungModulus() or
    /// this->GetEffectivePoissonRatio() and copies the result to rValue.
    /** @return The effective density
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<double>& rThisVariable,
                           double& rValue) override;

    ///@}
    ///@name Inquiry
    ///@{

    /// This lets user classes know if the constitutive law is defined for 1D, 2D or 3D.
    /** @return The number of spatial dimensions (1, 2 or 3).
     */
    SizeType WorkingSpaceDimension() override;

    /// This lets the user know the size of the strain rate vector (in Voigt notation) used by the constitutive law.
    /** @return The size of the strain rate vector.
     */
    SizeType GetStrainSize() override;

    ///@}
    ///@name Input and output
    ///@{

    /// @return A short string identifying this constitutive law instance.
    std::string Info() const override;

    /// Print basic information about this constitutive law instance.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print detailed information about this constitutive law instance and its managed data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

   protected:
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * @brief Get the Effective Young Modulus object
     * @param rParameters constitutive law parameters
     * @return double obtained effective young modulus
     */
    virtual double GetEffectiveYoungModulus(ConstitutiveLaw::Parameters& rParameters) const;

    /**
     * @brief Get the Poisson Ratio Modulus object
     * @param rParameters constitutive law parameters
     * @return double obtained effective Poisson Ratio
     */
    virtual double GetEffectivePoissonRatio(ConstitutiveLaw::Parameters& rParameters) const;

    /**
     * @brief Get the Density object
     * @param rParameters constitutive law parameters
     * @return double obtained effective density
     */
    virtual double GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const;

    /**
     * @brief Get the averaged value of the Input Variable
     * @param rParameters input variable, constitutive law parameters, step in which the variable is to evaluated
     * @return double obtained averaged variable
     */
    virtual double CalculateAveragedVariable(const Variable<double>& rVariableInput,
                                             ConstitutiveLaw::Parameters& rParameters, unsigned int step) const;

    /**
     * @brief Get the Value From Table object
     * For an table independent variable, this method returns the table dependent variable
     * value. Note that the properties container must have a table relating the two variables.
     * @param rIndependentVariable independent variable
     * @param rDependentVariable dependent variable
     * @param rParameters constitutive law parameters container
     * @return double output variable value
     */
    virtual double GetValueFromTable(const Variable<double>& rIndependentVariable, const Variable<double>& rDependentVariable,
                                     ConstitutiveLaw::Parameters& rParameters) const;

    ///@}

   private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

};  // Class PfemSolidConstitutiveLaw

}  // namespace Kratos.

#endif  // KRATOS_PFEM_SOLID_CONSTITUTIVE_LAW  defined
