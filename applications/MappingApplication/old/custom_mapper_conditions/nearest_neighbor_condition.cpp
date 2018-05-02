//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Include Base h
#include "nearest_neighbor_condition.h"
#include "mapping_application_variables.h"


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
 * Constructor using an array of nodes
 */
NearestNeighborCondition::NearestNeighborCondition(IndexType NewId, const NodesArrayType& ThisNodes)
    : Condition(NewId, ThisNodes) {
}

/**
 * Constructor using Geometry
 */
NearestNeighborCondition::NearestNeighborCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry) {
}

/**
 * Constructor using Properties
 */
NearestNeighborCondition::NearestNeighborCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties) {
}

/**
 * Copy Constructor
 */
NearestNeighborCondition::NearestNeighborCondition(NearestNeighborCondition const& rOther)
    : Condition(rOther) {
}

/**
 * Destructor
 */
NearestNeighborCondition::~NearestNeighborCondition() {
}

///@}
///@name Operators
///@{

/// Assignment operator.
NearestNeighborCondition & NearestNeighborCondition::operator=(NearestNeighborCondition const& rOther) {
  BaseType::operator=(rOther);
  Flags::operator =(rOther);
  // mpProperties = rOther.mpProperties;
  return *this;
}

///@}
///@name Operations
///@{

/**
 * CONDITIONS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param ThisNodes: the nodes of the new condition
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestNeighborCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Condition::Pointer(new NearestNeighborCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestNeighborCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Condition::Pointer(new NearestNeighborCondition(NewId, pGeom, pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new condition pointer and clones the previous condition data
 * @param NewId: the ID of the new condition
 * @param ThisNodes: the nodes of the new condition
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestNeighborCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const {
  KRATOS_TRY
  return Condition::Pointer(new NearestNeighborCondition(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
  KRATOS_CATCH("");
}

/**
 * this determines the condition equation ID vector for all conditional
 * DOFs
 * @param rResult: the condition equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestNeighborCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rResult.size() != number_of_nodes)
    rResult.resize(number_of_nodes, false);

}

/**
 * CONDITIONS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
 * they can be managed internally with a private method to do the same calculations
 * only once: MANDATORY
 */

/**
 * this is called during the assembling process in order
 * to calculate all condition contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestNeighborCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    if (this->Has(INTERFACE_INFO))
    {
        auto interface_info = this->GetValue(INTERFACE_INFO); // using auto bcs the type might change => shared to unique ptr
        if (interface_info.size() > 0)
        {
            if (rLeftHandSideMatrix.size1() != 3 || rLeftHandSideMatrix.size2() != 1)
                rLeftHandSideMatrix.resize(3,1);

            int nearest_neighbor_id = interface_info[0]->GetNeighborIds()[0];
            double nearest_neighbor_distance = interface_info[0]->GetNeighborDistances()[0];

            for (SizeType i=1; i<interface_info.size(); ++i)
            {
                const int curr_dist = interface_info[i]->GetNeighborDistances()[0];
                if (curr_dist < nearest_neighbor_distance)
                    nearest_neighbor_id = interface_info[i]->GetNeighborIds()[0];
            }

            rLeftHandSideMatrix(0,0) = 1.0;
            rLeftHandSideMatrix(1,0) = nearest_neighbor_id;
            rLeftHandSideMatrix(2,0) = this->Id();

            return;
        }
    }

    if (rLeftHandSideMatrix.size1() != 0 || rLeftHandSideMatrix.size2() != 0)
        rLeftHandSideMatrix.resize(0,0);
    std::cout << "NearestneighborConditions #" << this->Id() << "has no neighbor info" << std::endl;

}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix only
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestNeighborCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector only
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestNeighborCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int NearestNeighborCondition::Check(const ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY

  if (this->Id() < 1) {
    KRATOS_THROW_ERROR(std::logic_error, "NearestNeighborCondition found with Id 0 or negative","")
  }

  return 0;

  KRATOS_CATCH("");
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

std::string NearestNeighborCondition::Info() const {
  std::stringstream buffer;
  buffer << "NearestNeighborCondition #" << Id();
  return buffer.str();
}

/// Print information about this object.

void NearestNeighborCondition::PrintInfo(std::ostream& rOStream) const {
  rOStream << "NearestNeighborCondition #" << Id();
}

/// Print object's data.

void NearestNeighborCondition::PrintData(std::ostream& rOStream) const {
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

void NearestNeighborCondition::save(Serializer& rSerializer) const {
    KRATOS_ERROR << "This Object cannot be serialized!" << std::endl;
}

void NearestNeighborCondition::load(Serializer& rSerializer) {
    KRATOS_ERROR << "This Object cannot be serialized!" << std::endl;
}

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
inline std::istream & operator >> (std::istream& rIStream, NearestNeighborCondition& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const NearestNeighborCondition& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.
