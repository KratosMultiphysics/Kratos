//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined (KRATOS_FLUID_CONSTITUTIVE_LAW)
#define  KRATOS_FLUID_CONSTITUTIVE_LAW

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// This class contains the common infrastructure for fluid constitutive laws.
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidConstitutiveLaw : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FluidConstitutiveLaw);

    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FluidConstitutiveLaw();

    /// Copy constructor.
    FluidConstitutiveLaw (const FluidConstitutiveLaw& rOther);

    /// Destructor
    ~FluidConstitutiveLaw() override;

    ///@}
    ///@name Operations
    ///@{

    /// Initialize a new instance of this type of law
    ConstitutiveLaw::Pointer Clone() const override;

    /// Calculate the response of the material for the current strain rates.
    /** This is the main method for fluid constitutive laws.
     *  It should calculate, at least, the tangent matrix d(Stress)/d(Strain rate) or a suitable linearization
     *  and the stresses corresponding to the current strain rates.
     *  These are returned as rValues.GetConstitutiveMatrix() and rValues.GetStressVector(), respectively.
     *  @note Besides computing the response, derived constitutive laws are responsible for setting the 
     *  effective viscosity for the law, which will be used by the element to calculate the stabilization Tau.
     *  When implementing derived laws, please make sure to call SetEffectiveViscosity() at some point within this function.
     *  @see FluidElement, FluidElementData.
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /// Validate the data received by the constitutive law
    /** @param rMaterialProperties Properties of the parent Element.
     *  @param rElementGeometry Geometry of the parent Element.
     *  @param rCurrentProcessInfo ProcessInfo for the problem.
     *  @return 0 if everything is fine, other values indicate problems.
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    /// Useless boilerplate that's just here to avoid a compilation warning. Only the double variant is meaningful.
    int& CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<int>& rThisVariable, int& rValue) override;
    
    /// We are abusing the constitutive law interface to return the effective viscosity to the calling element through this function.
    /** it always returns this->GetEffectiveViscosity() (and copies the result to rValue). */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<double>& rThisVariable, double& rValue) override;
    
    /// Useless boilerplate that's just here to avoid a compilation warning. Only the double variant is meaningful.
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<Vector>& rThisVariable, Vector& rValue) override;
    
    /// Useless boilerplate that's just here to avoid a compilation warning. Only the double variant is meaningful.
    Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;
    
    /// Useless boilerplate that's just here to avoid a compilation warning. Only the double variant is meaningful.
    array_1d<double, 3 > & CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<array_1d<double, 3 > >& rThisVariable,array_1d<double, 3 > & rValue) override;
    
    /// Useless boilerplate that's just here to avoid a compilation warning. Only the double variant is meaningful.
    array_1d<double, 6 > & CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<array_1d<double, 6 > >& rThisVariable, array_1d<double, 6 > & rValue) override;

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

    /// Helper function to write the constitutive matrix using an effective viscosity (2D version).
    /** It returns a matrix with the same structure as for a Newtonian fluid, using the given viscosity.
     *  @param[in] EffectiveViscosity Equivalent viscosity for the fluid (dynamic units -- Pa s -- assumed).
     *  @param[out] rC Resulting constitutive matrix.
     */
    void NewtonianConstitutiveMatrix2D(double EffectiveViscosity, Matrix& rC);

    /// Helper function to write the constitutive matrix using an effective viscosity (3D version).
    /** It returns a matrix with the same structure as for a Newtonian fluid, using the given viscosity.
     *  @param[in] EffectiveViscosity Equivalent viscosity for the fluid (dynamic units -- Pa s -- assumed).
     *  @param[out] rC Resulting constitutive matrix.
     */
    void NewtonianConstitutiveMatrix3D(double EffectiveViscosity, Matrix& rC);

    ///@}
    ///@name Protected  Access
    ///@{

    /// Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
    virtual double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const;

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

}; // Class FluidConstitutiveLaw

}  // namespace Kratos.

#endif // KRATOS_FLUID_CONSTITUTIVE_LAW  defined 
