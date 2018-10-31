//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Thomas Oberbichler
//

// System includes

// External includes

// Project includes
#include "NurbsBrepModeler.h"


namespace Kratos
{
    void NurbsBrepModeler::ImportGeometry(BrepJSON_IO& rBrepJSON_IO)
    {
        std::vector<BrepModel> brep_model_vector = rBrepJSON_IO.ImportGeometry(mp_model_part);
        for (auto brep_model = brep_model_vector.begin(); brep_model != brep_model_vector.end(); ++brep_model)
        {
            m_brep_model_vector.push_back(brep_model::Pointer);
        }
    }
}  // namespace Kratos.


