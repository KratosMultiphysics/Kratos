//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Daniel Diez
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "custom_utilities/element_size_calculator.h"
#include "custom_constitutive/newtonian_two_fluid_2d_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewtonianTwoFluid2DLaw::NewtonianTwoFluid2DLaw()
    : Newtonian2DLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewtonianTwoFluid2DLaw::NewtonianTwoFluid2DLaw(const NewtonianTwoFluid2DLaw& rOther)
    : Newtonian2DLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewtonianTwoFluid2DLaw::Clone() const {
    return Kratos::make_shared<NewtonianTwoFluid2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewtonianTwoFluid2DLaw::~NewtonianTwoFluid2DLaw() {}

std::string NewtonianTwoFluid2DLaw::Info() const {
    return "NewtonianTwoFluid2DLaw";
}

double NewtonianTwoFluid2DLaw::ComputeEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    double viscosity;
    EvaluateInPoint(viscosity, DYNAMIC_VISCOSITY, rParameters);
    const Properties& prop = rParameters.GetMaterialProperties();

    if (prop.Has(C_SMAGORINSKY)) {
        const double csmag = prop[C_SMAGORINSKY];
        if (csmag > 0.0) {
            double density;
            EvaluateInPoint(density, DENSITY, rParameters);
            const double strain_rate = EquivalentStrainRate(rParameters);
            const BoundedMatrix<double,3,2>& rDN_DX = rParameters.GetShapeFunctionsDerivatives();
            const double elem_size = ElementSizeCalculator<2,3>::GradientsElementSize(rDN_DX);
            double length_scale = std::pow(csmag * elem_size, 2);
            viscosity += 2.0*length_scale * strain_rate * density;
        }
    }
    return viscosity;
}

void NewtonianTwoFluid2DLaw::EvaluateInPoint(
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

double NewtonianTwoFluid2DLaw::EquivalentStrainRate(ConstitutiveLaw::Parameters& rParameters) const{

    const Vector& S = rParameters.GetStrainVector();

    // Norm of symetric gradient (cross terms don't get the 2)
    return std::sqrt(2.0*S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2]);
}

void NewtonianTwoFluid2DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NewtonianTwoFluid2DLaw )
}

void NewtonianTwoFluid2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NewtonianTwoFluid2DLaw )
}

} // Namespace Kratos
