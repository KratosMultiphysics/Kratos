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
//                   Riccardo Rossi
//

// System includes


// External includes 


// Project includes
#include "NurbsBrepModeler.h"
#include "utilities/math_utils.h"
#include "BrepModelGeometryReader.h"


namespace Kratos
{
  void NurbsBrepModeler::SetUp(boost::python::dict cad_geometry, ModelPart& model_part)
  {
    std::cout << "Test 1" << std::endl;
    BrepModelGeometryReader& m_brep_model_geometry_reader = BrepModelGeometryReader(cad_geometry);
    std::cout << "Test 2" << std::endl;
    m_brep_model_geometry_reader.ReadGeometry(m_brep_model_vector, model_part);
  }

NurbsBrepModeler::NurbsBrepModeler()
{
}

NurbsBrepModeler::~NurbsBrepModeler()
{}

}  // namespace Kratos.


