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

#include "solid_constitutive_law.h"

namespace Kratos
{

    // Life cycle /////////////////////////////////////////////////////////////////

    PfemSolidConstitutiveLaw::PfemSolidConstitutiveLaw() : ConstitutiveLaw() {}

    PfemSolidConstitutiveLaw::PfemSolidConstitutiveLaw(const PfemSolidConstitutiveLaw &rOther) : ConstitutiveLaw(rOther) {}

    PfemSolidConstitutiveLaw::~PfemSolidConstitutiveLaw() {}

    // Public operations //////////////////////////////////////////////////////////

    ConstitutiveLaw::Pointer PfemSolidConstitutiveLaw::Clone() const
    {
        KRATOS_ERROR << "Calling base PfemSolidConstitutiveLaw::Clone method. This "
                        "class should not be instantiated. Please check your "
                        "constitutive law."
                     << std::endl;
        return Kratos::make_shared<PfemSolidConstitutiveLaw>(*this);
    }

    void PfemSolidConstitutiveLaw::CalculateMaterialResponseCauchy(Parameters &rValues)
    {
        KRATOS_ERROR << "Calling base "
                        "PfemSolidConstitutiveLaw::CalculateMaterialResponseCauchy "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
    }

    int PfemSolidConstitutiveLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                        const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_ERROR << "Calling base "
                        "PfemSolidConstitutiveLaw::Check "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
        return 999;
    }

    // Access /////////////////////////////////////////////////////////////////////

    double &PfemSolidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters &rParameters,
                                                     const Variable<double> &rThisVariable, double &rValue)
    {
        if (rThisVariable.Name() == "DENSITY")
        {
            rValue = this->GetEffectiveMaterialParameter(rParameters, DENSITY);
        }
        else if (rThisVariable.Name() == "YOUNG_MODULUS")
        {
            rValue = this->GetEffectiveMaterialParameter(rParameters, YOUNG_MODULUS);
        }
        else if (rThisVariable.Name() == "POISSON_RATIO")
        {
            rValue = this->GetEffectiveMaterialParameter(rParameters, POISSON_RATIO);
        }
        else
        {
            KRATOS_ERROR << " An unexpected property has been passed " << std::endl;
        }
        return rValue;
    }

    double PfemSolidConstitutiveLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
    {
        KRATOS_ERROR << "Accessing base class PfemSolidConstitutiveLaw::GetEffectiveMaterialParameter." << std::endl;
        return 0.0;
    }
    // Inquiry ////////////////////////////////////////////////////////////////////

    ConstitutiveLaw::SizeType PfemSolidConstitutiveLaw::WorkingSpaceDimension()
    {
        KRATOS_ERROR << "Calling base "
                        "PfemSolidConstitutiveLaw::WorkingSpaceDimension "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
        return 0;
    }

    ConstitutiveLaw::SizeType PfemSolidConstitutiveLaw::GetStrainSize() const
    {
        KRATOS_ERROR << "Calling base "
                        "PfemSolidConstitutiveLaw::GetStrainSize "
                        "method. This class should not be instantiated. Please "
                        "check your constitutive law."
                     << std::endl;
        return 0;
    }

    // Info ///////////////////////////////////////////////////////////////////////

    std::string PfemSolidConstitutiveLaw::Info() const { return "PfemSolidConstitutiveLaw"; }

    void PfemSolidConstitutiveLaw::PrintInfo(std::ostream &rOStream) const { rOStream << this->Info(); }

    void PfemSolidConstitutiveLaw::PrintData(std::ostream &rOStream) const { rOStream << this->Info(); }

    // Protected operations ///////////////////////////////////////////////////////

    // Protected access ///////////////////////////////////////////////////////////

    double PfemSolidConstitutiveLaw::CalculateAveragedVariable(const Variable<double> &rVariableInput,
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

    double PfemSolidConstitutiveLaw::GetValueFromTable(const Variable<double> &rIndependentVariable,
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

    void PfemSolidConstitutiveLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    void PfemSolidConstitutiveLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

} // namespace Kratos
