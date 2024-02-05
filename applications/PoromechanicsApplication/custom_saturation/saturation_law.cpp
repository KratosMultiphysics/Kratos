//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

/* Project includes */
#include "custom_saturation/saturation_law.hpp"

namespace Kratos
{

void SaturationLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues)
{
    // nothing
}

void SaturationLaw::Initialize(Parameters &rParameters)
{
    // nothing
}

void SaturationLaw::InitializeSolutionStep(Parameters &rParameters)
{
    // nothing
}

void SaturationLaw::FinalizeSolutionStep(Parameters &rParameters)
{
    // nothing
}

void SaturationLaw::Finalize(Parameters &rParameters)
{
    // nothing
}

void SaturationLaw::ResetMaterial(const Properties &rMaterialProperties,
                                 const GeometryType &rElementGeometry,
                                 const Vector &rShapeFunctionsValues)
{
    // nothing
}

void SaturationLaw::save(Serializer& rSerializer) const
{
    // there is no member variables to be saved
}

void SaturationLaw::load(Serializer& rSerializer)
{
    // there is no member variables to be loaded
}

}