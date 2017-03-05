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

//TODO: for testing in c++ use the testing mechanism of kratos
  //void NurbsBrepModeler::SetUp(Parameters& cad_geometry, ModelPart& model_part)
  //{
  //  std::cout << "Test 1" << std::endl;
  //  BrepModelGeometryReader& m_brep_model_geometry_reader = BrepModelGeometryReader(cad_geometry);
  //  std::cout << "Test 2" << std::endl;
  //  m_brep_model_geometry_reader.ReadGeometry(m_brep_model_vector, model_part);
  //}

NurbsBrepModeler::NurbsBrepModeler(Parameters& cad_geometry, ModelPart& model_part)
  : m_cad_geometry(cad_geometry),
    m_model_part(model_part)
{
  std::cout << "Test 1" << std::endl;
  BrepModelGeometryReader m_brep_model_geometry_reader(cad_geometry);
  std::cout << "Test 2" << std::endl;
  m_brep_model_geometry_reader.ReadGeometry(m_brep_model_vector, model_part);

}

NurbsBrepModeler::~NurbsBrepModeler()
{}

}  // namespace Kratos.


