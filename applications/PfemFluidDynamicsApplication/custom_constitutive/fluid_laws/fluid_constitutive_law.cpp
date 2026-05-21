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

#include "fluid_constitutive_law.h"

namespace Kratos
{

    // Life cycle /////////////////////////////////////////////////////////////////

    PfemFluidConstitutiveLaw::PfemFluidConstitutiveLaw() : ConstitutiveLaw() {}

    PfemFluidConstitutiveLaw::PfemFluidConstitutiveLaw(const PfemFluidConstitutiveLaw &rOther) : ConstitutiveLaw(rOther) {}

    PfemFluidConstitutiveLaw::~PfemFluidConstitutiveLaw() {}

    // Public operations //////////////////////////////////////////////////////////

    ConstitutiveLaw::Pointer PfemFluidConstitutiveLaw::Clone() const
    {
        KRATOS_ERROR << "Calling base PfemFluidConstitutiveLaw::Clone method. This "
                        "class should not be instantiated. Please check your "
                        "constitutive law."
                     << std::endl;
        return Kratos::make_shared<PfemFluidConstitutiveLaw>(*this);
    }

    void PfemFluidConstitutiveLaw::CalculateMaterialResponseCauchy(Parameters &rValues)
    {
        KRATOS_ERROR << "Calling base "
                        "PfemFluidConstitutiveLaw::CalculateMaterialResponseCauchy "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
    }

    int PfemFluidConstitutiveLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                        const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_ERROR << "Calling base "
                        "PfemFluidConstitutiveLaw::Check "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
        return 999;
    }

    // Access /////////////////////////////////////////////////////////////////////

    double &PfemFluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters &rParameters,
                                                     const Variable<double> &rThisVariable, double &rValue)
    {

        if (rThisVariable.Name() == "DENSITY")
        {
            rValue = this->GetEffectiveMaterialParameter(rParameters, DENSITY);
        }
        else
        {
            KRATOS_ERROR << " An unexpected property has been passed " << std::endl;
        }
        return rValue;
    }

    // Inquiry ////////////////////////////////////////////////////////////////////

    ConstitutiveLaw::SizeType PfemFluidConstitutiveLaw::WorkingSpaceDimension()
    {
        KRATOS_ERROR << "Calling base "
                        "PfemFluidConstitutiveLaw::WorkingSpaceDimension "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
        return 0;
    }

    ConstitutiveLaw::SizeType PfemFluidConstitutiveLaw::GetStrainSize() const
    {
        KRATOS_ERROR << "Calling base "
                        "PfemFluidConstitutiveLaw::GetStrainSize "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
        return 0;
    }

    // Info ///////////////////////////////////////////////////////////////////////

    std::string PfemFluidConstitutiveLaw::Info() const { return "PfemFluidConstitutiveLaw"; }

    void PfemFluidConstitutiveLaw::PrintInfo(std::ostream &rOStream) const { rOStream << this->Info(); }

    void PfemFluidConstitutiveLaw::PrintData(std::ostream &rOStream) const { rOStream << this->Info(); }

    // Protected operations ///////////////////////////////////////////////////////

    void PfemFluidConstitutiveLaw::EffectiveViscousConstitutiveMatrix2D(double EffectiveDynamicViscosity,
                                                                        Matrix &rConstitutiveMatrix)
    {
        constexpr double two_thirds = 2.0 / 3.0;
        constexpr double four_thirds = 4.0 / 3.0;

        rConstitutiveMatrix(0, 0) = +EffectiveDynamicViscosity * four_thirds;
        rConstitutiveMatrix(0, 1) = -EffectiveDynamicViscosity * two_thirds;
        rConstitutiveMatrix(0, 2) = 0.0;
        rConstitutiveMatrix(1, 0) = -EffectiveDynamicViscosity * two_thirds;
        rConstitutiveMatrix(1, 1) = +EffectiveDynamicViscosity * four_thirds;
        rConstitutiveMatrix(1, 2) = 0.0;
        rConstitutiveMatrix(2, 0) = 0.0;
        rConstitutiveMatrix(2, 1) = 0.0;
        rConstitutiveMatrix(2, 2) = +EffectiveDynamicViscosity;
    }

    void PfemFluidConstitutiveLaw::EffectiveViscousConstitutiveMatrix3D(double EffectiveDynamicViscosity,
                                                                        Matrix &rConstitutiveMatrix)
    {
        rConstitutiveMatrix.clear();

        constexpr double two_thirds = 2.0 / 3.0;
        constexpr double four_thirds = 4.0 / 3.0;

        rConstitutiveMatrix(0, 0) = +EffectiveDynamicViscosity * four_thirds;
        rConstitutiveMatrix(0, 1) = -EffectiveDynamicViscosity * two_thirds;
        rConstitutiveMatrix(0, 2) = -EffectiveDynamicViscosity * two_thirds;

        rConstitutiveMatrix(1, 0) = -EffectiveDynamicViscosity * two_thirds;
        rConstitutiveMatrix(1, 1) = +EffectiveDynamicViscosity * four_thirds;
        rConstitutiveMatrix(1, 2) = -EffectiveDynamicViscosity * two_thirds;

        rConstitutiveMatrix(2, 0) = -EffectiveDynamicViscosity * two_thirds;
        rConstitutiveMatrix(2, 1) = -EffectiveDynamicViscosity * two_thirds;
        rConstitutiveMatrix(2, 2) = +EffectiveDynamicViscosity * four_thirds;

        rConstitutiveMatrix(3, 3) = +EffectiveDynamicViscosity;
        rConstitutiveMatrix(4, 4) = +EffectiveDynamicViscosity;
        rConstitutiveMatrix(5, 5) = +EffectiveDynamicViscosity;
    }

    double PfemFluidConstitutiveLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
    {
        KRATOS_ERROR << "Accessing base class PfemFluidConstitutiveLaw::GetEffectiveMaterialParameter." << std::endl;
        return 0.0;
    }

    double PfemFluidConstitutiveLaw::CalculateAveragedVariable(const Variable<double> &rVariableInput,
                                                               ConstitutiveLaw::Parameters &rParameters,
                                                               unsigned int step) const
    {
        const GeometryType &r_geometry = rParameters.GetElementGeometry();
        const unsigned int number_of_nodes = r_geometry.size();

        double result = 0;
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            result += r_geometry[i].FastGetSolutionStepValue(rVariableInput, step);
        }
        result /= number_of_nodes;

        return result;
    }

    double PfemFluidConstitutiveLaw::CalculateInGaussPoint(const Variable<double> &rVariableInput,
                                                           ConstitutiveLaw::Parameters &rParameters,
                                                           unsigned int step) const
    {

        const GeometryType &r_geometry = rParameters.GetElementGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const auto &r_shape_function = rParameters.GetShapeFunctionsValues();
        double result = 0;

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            result += r_shape_function[i] * r_geometry[i].FastGetSolutionStepValue(rVariableInput, step);
        }

        return result;
    }

    double PfemFluidConstitutiveLaw::GetValueFromTable(const Variable<double> &rIndependentVariable,
                                                       const Variable<double> &rDependentVariable,
                                                       ConstitutiveLaw::Parameters &rParameters) const
    {
        // Get material properties from constitutive law parameters
        const Properties &r_properties = rParameters.GetMaterialProperties();

        // Get geometry and Gauss points data
        const auto &r_geometry = rParameters.GetElementGeometry();
        const auto &r_N = rParameters.GetShapeFunctionsValues();

        // Compute the independent variable at the Gauss point
        double independent_at_gauss = 0.0;
        double dependent_at_gauss = 0.0;
        for (unsigned int i = 0; i < r_N.size(); ++i)
        {
            const double &r_val = r_geometry[i].FastGetSolutionStepValue(rIndependentVariable);
            independent_at_gauss += r_val * r_N[i];
        }

        // Retrieve the dependent variable from the table
        const auto &r_table = r_properties.GetTable(rIndependentVariable, rDependentVariable);
        dependent_at_gauss = r_table.GetValue(independent_at_gauss);

        return dependent_at_gauss;
    }

    // Serialization //////////////////////////////////////////////////////////////

    void PfemFluidConstitutiveLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    void PfemFluidConstitutiveLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

} // namespace Kratos