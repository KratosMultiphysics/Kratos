//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//
//  Contributors:   Altug Emiroglu, https://github.com/emiroglu
//



#if !defined(KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
    typedef std::unordered_map<VariableData::KeyType, int> TMapPhiType;
    inline std::ostream& operator << (std::ostream& rOStream ,
                                  const std::unordered_map<VariableData::KeyType, int>& rThis)
    {
        for (auto const& r_pair: rThis) {
            rOStream << "{" << r_pair.first << ": " << r_pair.second << "}\n";
        }
        return rOStream;
    }

    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, int, AUX_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Matrix, ROM_BASIS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, double, HROM_WEIGHT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, TMapPhiType, MAP_PHI )

    // Modal derivative variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, unsigned int, BUILD_LEVEL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Vector, EIGENVALUE_VECTOR)
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, std::size_t, BASIS_I )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, std::size_t, BASIS_J )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, std::size_t, DERIVATIVE_INDEX )

}

#endif	/* KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED */
