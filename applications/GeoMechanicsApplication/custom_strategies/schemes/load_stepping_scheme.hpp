// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "geomechanics_static_scheme.hpp"

template <class TSparseSpace, class TDenseSpace>
class LoadSteppingScheme : public Kratos::GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>
{

};

