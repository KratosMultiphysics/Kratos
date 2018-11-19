//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "custom_utilities/element_size_calculator.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewtonianTwoFluid3DLaw::NewtonianTwoFluid3DLaw()
    : Newtonian3DLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewtonianTwoFluid3DLaw::NewtonianTwoFluid3DLaw(const NewtonianTwoFluid3DLaw& rOther)
    : Newtonian3DLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewtonianTwoFluid3DLaw::Clone() const {
    return Kratos::make_shared<NewtonianTwoFluid3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewtonianTwoFluid3DLaw::~NewtonianTwoFluid3DLaw() {}


std::string NewtonianTwoFluid3DLaw::Info() const {
    return "NewtonianTwoFluid3DLaw";
}


double NewtonianTwoFluid3DLaw::ComputeEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {

    double viscosity;
    EvaluateInPoint(viscosity, DYNAMIC_VISCOSITY, rParameters);
    const Properties& prop = rParameters.GetMaterialProperties();

    if (prop.Has(C_SMAGORINSKY)) {
        const double csmag = prop[C_SMAGORINSKY];
        if (csmag > 0.0) {
            double density;
            EvaluateInPoint(density, DENSITY, rParameters);
            const double strain_rate = EquivalentStrainRate(rParameters);
            const BoundedMatrix<double, 4, 3>& DN_DX = rParameters.GetShapeFunctionsDerivatives();
            const double elem_size = ElementSizeCalculator<3,4>::GradientsElementSize(DN_DX);
            double length_scale = csmag * elem_size;
            length_scale *= length_scale;
            viscosity += 2.0*length_scale * strain_rate * density;
        }
    }
    return viscosity;
}


void NewtonianTwoFluid3DLaw::EvaluateInPoint(double& rResult,
    const Variable<double>& rVariable,
    ConstitutiveLaw::Parameters& rParameters) const {

    const SizeType nnodes = 4;
    const GeometryType& geom = rParameters.GetElementGeometry();
    const array_1d<double,nnodes>& N = rParameters.GetShapeFunctionsValues();

    //compute sign of distance on gauss point
    double dist = 0.0;
    for (unsigned int i = 0; i < nnodes; i++)
        dist += N[i] * geom[i].FastGetSolutionStepValue(DISTANCE);

    SizeType navg = 0; //number of nodes on the same side as the gauss point
    double value = 0.0;
    for (unsigned int i = 0; i < nnodes; i++) {
        if ( dist * geom[i].FastGetSolutionStepValue(DISTANCE) > 0.0) {
            navg += 1;
            value += geom[i].FastGetSolutionStepValue(rVariable);
        }
    }
    rResult = value/navg;
}


double NewtonianTwoFluid3DLaw::EquivalentStrainRate(ConstitutiveLaw::Parameters& rParameters) const {
    
    // Calculate Symetric gradient (Voigt notation)
    const SizeType n_nodes = 4;
    const SizeType voigt_vector_size = 6;
    const GeometryType& r_geom = rParameters.GetElementGeometry();
    const Matrix& rDN_DX = rParameters.GetShapeFunctionsDerivatives();

    BoundedMatrix<double, voigt_vector_size, 4 * n_nodes> StrainMatrix;

    FluidElementUtilities<n_nodes>::GetStrainMatrix( rDN_DX, StrainMatrix );
    Vector dofs(16, 0.0);
    for ( unsigned int i = 0; i < n_nodes; i++){
        dofs[4*i + 0] = r_geom[i].FastGetSolutionStepValue(VELOCITY_X);
        dofs[4*i + 1] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Y);
        dofs[4*i + 2] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Z);
    }

    Vector S = prod( StrainMatrix, dofs);

    // Norm of symetric gradient (cross terms don't get the 2)
    return std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + 2.*S[2]*S[2] + S[3]*S[3] + S[4]*S[4] + S[5]*S[5]);
}



void NewtonianTwoFluid3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

void NewtonianTwoFluid3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

} // Namespace Kratos
