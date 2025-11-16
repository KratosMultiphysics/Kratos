//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//


// System includes
#include <tuple>
#include <iomanip>
#include <sstream>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "rans_nut_utility.h"

namespace Kratos
{
RansNutUtility::RansNutUtility(
    ModelPart& rModelPart,
    const double RelativeTolerance,
    const double AbsoluteTolerance,
    const int EchoLevel) :
    mrModelPart(rModelPart),
    mEchoLevel(EchoLevel),
    mRelativeTolerance(RelativeTolerance),
    mAbsoluteTolerance(AbsoluteTolerance)
{
    KRATOS_TRY

    block_for_each(mrModelPart.Elements(),[&](const ElementType& rElement) {
        KRATOS_ERROR_IF_NOT(rElement.GetProperties().Has(CONSTITUTIVE_LAW)) << "CONSTITUTIVE_LAW is not found in element"
            <<" data value container. [ Element.Id() = " << rElement.Id() << " ].\n";

        const auto p_constitutive = rElement.GetProperties().GetValue(CONSTITUTIVE_LAW);
        const auto rans_cl_name = p_constitutive->Info();

        KRATOS_ERROR_IF(rans_cl_name.substr(0, 4) != "Rans")
            << "Incompatible constitutive law is used. Please use constitutive "
               "laws which starts with \"Rans*\" [ Constitutive law "
               "name = "
            << rans_cl_name << ", Element.Id() = " << rElement.Id() << " ].\n";
    });

    KRATOS_CATCH("");
}

void RansNutUtility::Initialize()
{
    KRATOS_TRY

    // Initialize elemental turbulent viscosity also as the initial values
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    block_for_each(mrModelPart.Elements(), TLSType(), [&](ElementType& rElement, TLSType& rTLS){
        rElement.SetValue(TURBULENT_VISCOSITY, CalculateTurbulentViscosity(rElement, rTLS, r_process_info));
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Following variables are properly initialized on " << mrModelPart.Name() << ":\n"
        << "       " << TURBULENT_VISCOSITY.Name() << "\n";

    KRATOS_CATCH("");
}

void RansNutUtility::InitializeCalculation()
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();
    const auto& elements = mrModelPart.Elements();
    const int number_of_elements = elements.size();

    if (static_cast<int>(mElementData.size()) < number_of_elements) {
        mElementData.resize(number_of_elements);
    }

    IndexPartition<int>(number_of_elements).for_each(TLSType(), [&](const int iElement, TLSType& rTLS){
        mElementData[iElement] = CalculateTurbulentViscosity(*(elements.begin() + iElement), rTLS, r_process_info);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Initialization of  convergence check completed on " << mrModelPart.Name() << ".\n";

    KRATOS_CATCH("");
}


bool RansNutUtility::CheckConvergence() const
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();
    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& elements = mrModelPart.Elements();
    const int number_of_elements = elements.size();

    KRATOS_ERROR_IF(static_cast<int>(mElementData.size()) < number_of_elements)
        << "Data is not properly initialized for " << mrModelPart.Name()
        << ". Please use \"InitializeCalculation\" first.\n";

    double element_dx, element_solution;
    std::tie(element_dx, element_solution) = IndexPartition<int>(number_of_elements).for_each<
        CombinedReduction<
            SumReduction<double>,
            SumReduction<double>
        >>(TLSType(), [&](const int iElement, TLSType& rTLS) -> std::tuple<double, double> {
            const double value = CalculateTurbulentViscosity(*(elements.begin() + iElement), rTLS, r_process_info);

            return std::make_tuple<double, double>(std::pow(value - mElementData[iElement], 2), std::pow(value, 2));
        });

    const std::vector<double> norm_values = {element_dx, element_solution, static_cast<double>(number_of_elements)};
    const auto& total_norm_values = r_communicator.GetDataCommunicator().SumAll(norm_values);

    const double dx = std::sqrt(total_norm_values[0]);
    double solution = std::sqrt(total_norm_values[1]);
    solution = (solution == 0.0 ? 1.0 : solution);

    const double relative_error = dx / solution;
    const double absolute_error = dx / total_norm_values[2];

    if (mEchoLevel > 0 ) {
        std::stringstream msg;
        msg << std::scientific << std::setprecision(6) << "[ "
            << "Obtained ratio: " << relative_error << "; "
            << "Expected ratio: " << mRelativeTolerance << "; "
            << "Absolute norm: " << absolute_error << ";"
            << "Expected norm: " << mAbsoluteTolerance << " ] - TURBULENT_VISCOSITY\n";
        KRATOS_INFO(this->Info()) << msg.str();
    }

    if (relative_error < mRelativeTolerance || absolute_error < mAbsoluteTolerance) {
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << " CONVERGENCE: *** CONVERGENCE ACHIEVED *** [ TURBULENT VISCOSITY ] \n";
        return true;
    }

    return false;

    KRATOS_CATCH("");
}

void RansNutUtility::UpdateTurbulentViscosity()
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();

    block_for_each(mrModelPart.Elements(), TLSType(), [&](ElementType& rElement, TLSType& rTLS){
        rElement.GetValue(TURBULENT_VISCOSITY) = CalculateTurbulentViscosity(rElement, rTLS, r_process_info);
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated turbulent viscosity on " << mrModelPart.Name() << ".\n";

    KRATOS_CATCH("");
}

double RansNutUtility::CalculateTurbulentViscosity(
    ElementType& rElement,
    TLSType& rTLS,
    const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    auto& Ws = std::get<0>(rTLS);
    auto& Ns = std::get<1>(rTLS);
    auto& dNdXs = std::get<2>(rTLS);


    // computing everything based on a fixed gauss integration rather than based
    // on the element one. This is because, in RANS there can be different elements
    // with different gauss integration methods. So in order to be consistent
    // GI_GAUSS_1 is chosen
    const auto& r_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;

    RansCalculationUtilities::CalculateGeometryData(
        rElement.GetGeometry(), r_integration_method, Ws, Ns, dNdXs);

    ConstitutiveLaw::Parameters parameters(rElement.GetGeometry(), rElement.GetProperties(), rProcessInfo);
    const Vector& N = row(Ns, 0);
    const Matrix& dNdX = dNdXs[0];
    parameters.SetShapeFunctionsValues(N);
    parameters.SetShapeFunctionsDerivatives(dNdX);

    double nu_t;
    rElement.GetProperties().GetValue(CONSTITUTIVE_LAW)->CalculateValue(parameters, TURBULENT_VISCOSITY, nu_t);

    return nu_t;

    KRATOS_CATCH("");
}

} // namespace Kratos.
