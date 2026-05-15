//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//                   
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/element_size_calculator.h"
#include "custom_constitutive/bingham_two_fluid_2d_law.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

BinghamTwoFluid2DLaw::BinghamTwoFluid2DLaw()
    : Bingham2DLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

BinghamTwoFluid2DLaw::BinghamTwoFluid2DLaw(const BinghamTwoFluid2DLaw& rOther)
    : Bingham2DLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer BinghamTwoFluid2DLaw::Clone() const {
    return Kratos::make_shared<BinghamTwoFluid2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

BinghamTwoFluid2DLaw::~BinghamTwoFluid2DLaw() {}

std::string BinghamTwoFluid2DLaw::Info() const {
    return "BinghamTwoFluid2DLaw";
}

int BinghamTwoFluid2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    for (unsigned int i = 0; i < rElementGeometry.size(); i++) {
        const Node& rNode = rElementGeometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DYNAMIC_VISCOSITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rNode);
        KRATOS_ERROR_IF(rNode.GetSolutionStepValue(DYNAMIC_VISCOSITY) <= 0.0)
            << "DYNAMIC_VISCOSITY was not correctly assigned to nodes for Constitutive Law.\n";
        KRATOS_ERROR_IF(rNode.GetSolutionStepValue(DENSITY) <= 0.0)
            << "DENSITY was not correctly assigned to nodes for Constitutive Law.\n";
    }
    return 0;
}

// double BinghamTwoFluid2DLaw::(ConstitutiveLaw::Parameters& rParameters) const
// {
//     double viscosity;
//     // EvaluateInPoint(viscosity, DYNAMIC_VISCOSITY, rParameters);
//     // const Properties& prop = rParameters.GetMaterialProperties();
//     const auto& r_geom = rParameters.GetElementGeometry();

//     // if (r_geom.Has(ARTIFICIAL_DYNAMIC_VISCOSITY)){
//     //     viscosity += r_geom.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
//     // }

//     // if (prop.Has(C_SMAGORINSKY)) {
//     //     const double csmag = prop[C_SMAGORINSKY];
//     //     if (csmag > 0.0) {
//     //         double density;
//     //         EvaluateInPoint(density, DENSITY, rParameters);
//     //         const double strain_rate = EquivalentStrainRate(rParameters);
//     //         const BoundedMatrix<double,3,2>& rDN_DX = rParameters.GetShapeFunctionsDerivatives();
//     //         const double elem_size = ElementSizeCalculator<2,3>::GradientsElementSize(rDN_DX);
//     //         double length_scale = std::pow(csmag * elem_size, 2);
//     //         viscosity += 2.0*length_scale * strain_rate * density;
//     //     }
//     // }
//     return viscosity;
// }

void BinghamTwoFluid2DLaw::EvaluateInPoint(
    double& rResult,
    const Variable<double>& rVariable,
    ConstitutiveLaw::Parameters& rParameters) const
{
    const SizeType n_nodes = 3;
    const GeometryType& r_geom = rParameters.GetElementGeometry();
    const array_1d<double,n_nodes>& rN = rParameters.GetShapeFunctionsValues();

    // Compute Gauss pt. distance value
    double dist = 0.0;
    for (unsigned int i = 0; i < n_nodes; ++i){
        dist += rN[i] * r_geom[i].FastGetSolutionStepValue(DISTANCE);
    }

    // Compute the Gauss pt. value considering the previously compute distance value
    rResult = 0.0; // Initialize result variable value
    SizeType n_avg = 0; // Number of nodes on the same side as the gauss point
    for (unsigned int i = 0; i < n_nodes; ++i) {
        if ( dist * r_geom[i].FastGetSolutionStepValue(DISTANCE) > 0.0) {
            n_avg += 1;
            rResult += r_geom[i].FastGetSolutionStepValue(rVariable);
        }
    }
    rResult /= n_avg;
}

double BinghamTwoFluid2DLaw::EquivalentStrainRate(ConstitutiveLaw::Parameters& rParameters) const{

    const Vector& S = rParameters.GetStrainVector();

    // Norm of symetric gradient (cross terms don't get the 2)
    return std::sqrt(2.0*S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2]);
}

void BinghamTwoFluid2DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Bingham2DLaw )
}

void BinghamTwoFluid2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Bingham2DLaw )
}

} // Namespace Kratos
