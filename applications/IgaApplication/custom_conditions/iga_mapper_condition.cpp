//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    TT
//

// System includes

// External includes


// Include Base h
#include "iga_mapper_condition.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Constructor.
 */


///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

void IGAMapperCondition::Calculate(const Variable<Matrix >& rVariable,
    Matrix& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    //mpModeler->CalculateIntegrationPoints(GetGeometry(), Output, GetValue(INTEGRATION_METHOD));

    //INTERFACE_EQUATION_ID

}


///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a string.

std::string IGAMapperCondition::Info() const {
  std::stringstream buffer;
  buffer << "IGAMapperCondition #" << Id();
  return buffer.str();
}

/// Print information about this object.

void IGAMapperCondition::PrintInfo(std::ostream& rOStream) const {
  rOStream << "IGAMapperCondition #" << Id();
}

/// Print object's data.

void IGAMapperCondition::PrintData(std::ostream& rOStream) const {
  pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

///@}
///@name Serialization
///@{

//void IGAMapperCondition::save(Serializer& rSerializer) const {
  //KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, @{KRATOS_CLASS_BASE} );

  // List
  // To be completed with the class member list
//}

//void IGAMapperCondition::load(Serializer& rSerializer) {
  //KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, @{KRATOS_CLASS_BASE} );

  // List
  // To be completed with the class member list
//}

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >> (std::istream& rIStream, IGAMapperCondition& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const IGAMapperCondition& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.
