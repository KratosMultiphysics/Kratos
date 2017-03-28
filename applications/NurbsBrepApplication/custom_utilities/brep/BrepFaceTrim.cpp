//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// System includes


// External includes 


// Project includes
#include "BrepFaceTrim.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
// --------------------------------------------------------------------------

//Constructor
BrepFaceTrim::BrepFaceTrim(unsigned int face_id, unsigned int trim_index, bool relative_direction)
  : m_face_id(face_id),
    m_trim_index(trim_index),
    m_relative_direction(relative_direction)
{
}
//Destructor
BrepFaceTrim::~BrepFaceTrim()
{}

}  // namespace Kratos.

