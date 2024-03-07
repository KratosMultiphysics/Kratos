// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//                   Gennady Markelov
//

#pragma once

// Project includes

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos::UPwSmallStrain
{
using NodeType     = Node;
using GeometryType = Geometry<NodeType>;
using MatrixType   = Matrix;

template <unsigned int TDim, unsigned int TNumNodes>
struct ElementVariables {
    /// Properties variables
    bool   IgnoreUndrained            = false;
    bool   UseHenckyStrain            = false;
    bool   ConsiderGeometricStiffness = false;
    double DynamicViscosityInverse    = 0.0;
    double FluidDensity               = 0.0;
    double SolidDensity               = 0.0;
    double Density                    = 0.0;
    double Porosity                   = 0.0;
    double PermeabilityUpdateFactor   = 0.0;

    double                            BiotCoefficient    = 0.0;
    double                            BiotModulusInverse = 0.0;
    BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;

    /// ProcessInfo variables
    double VelocityCoefficient   = 0.0;
    double DtPressureCoefficient = 0.0;

    /// Nodal variables
    array_1d<double, TNumNodes>        PressureVector;
    array_1d<double, TNumNodes>        DtPressureVector;
    array_1d<double, TNumNodes * TDim> DisplacementVector;
    array_1d<double, TNumNodes * TDim> VelocityVector;
    array_1d<double, TNumNodes * TDim> VolumeAcceleration;

    /// General elemental variables
    Vector VoigtVector;

    /// Variables computed at each GP
    Matrix                                        B;
    BoundedMatrix<double, TDim, TNumNodes * TDim> Nu;
    array_1d<double, TDim>                        BodyAcceleration;
    array_1d<double, TDim>                        SoilGamma;

    /// Constitutive Law parameters
    Vector StrainVector;
    Vector StressVector;
    Matrix ConstitutiveMatrix;
    Vector Np;
    Matrix GradNpT;
    Matrix GradNpTInitialConfiguration;

    Matrix                                    F;
    double                                    detF = 0.0;
    Vector                                    detJContainer;
    Matrix                                    NContainer;
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

    /// Retention Law parameters
    double FluidPressure          = 0.0;
    double DegreeOfSaturation     = 0.0;
    double DerivativeOfSaturation = 0.0;
    double RelativePermeability   = 0.0;
    double BishopCoefficient      = 0.0;
    double EffectiveSaturation    = 0.0;

    // needed for updated Lagrangian:
    double detJ                                       = 0.0;
    double detJInitialConfiguration                   = 0.0;
    double IntegrationCoefficient                     = 0.0;
    double IntegrationCoefficientInitialConfiguration = 0.0;

    // Auxiliary Variables
    BoundedMatrix<double, TNumNodes * TDim, TNumNodes * TDim> UMatrix;
    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>        UPMatrix;
    BoundedMatrix<double, TNumNodes, TNumNodes * TDim>        PUMatrix;
    BoundedMatrix<double, TNumNodes, TNumNodes>               PMatrix;
    Matrix                                                    UVoigtMatrix;
    BoundedMatrix<double, TNumNodes, TDim>                    PDimMatrix;
    array_1d<double, TNumNodes * TDim>                        UVector;
    array_1d<double, TNumNodes>                               PVector;
};
} // namespace Kratos::UPwSmallStrain