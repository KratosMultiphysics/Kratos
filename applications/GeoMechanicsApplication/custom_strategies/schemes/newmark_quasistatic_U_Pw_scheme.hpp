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
//

#pragma once

#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

#include "geo_mechanics_application_variables.h"
#include "custom_strategies/schemes/geo_base_scheme.hpp"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkQuasistaticUPwScheme : public GeoBaseScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    NewmarkQuasistaticUPwScheme(double beta, double gamma, double theta)
        : GeoBaseScheme<TSparseSpace,TDenseSpace>(beta, gamma, theta)

    {
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~NewmarkQuasistaticUPwScheme() override {}


}; // Class NewmarkQuasistaticUPwScheme
}  // namespace Kratos
