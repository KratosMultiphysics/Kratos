//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Include Base h
#include "custom_elements/compressible_potential_flow_element.h"
#include "utilities/geometry_utilities.h"

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
CompressiblePotentialFlowElement::CompressiblePotentialFlowElement(IndexType NewId)
    : Element(NewId)  {
}

/**
 * Constructor using an array of nodes
 */
CompressiblePotentialFlowElement::CompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)  {
}

/**
 * Constructor using Geometry
 */
CompressiblePotentialFlowElement::CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)  {
}

/**
 * Constructor using Properties
 */
CompressiblePotentialFlowElement::CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)  {
}

/**
 * Copy Constructor
 */
CompressiblePotentialFlowElement::CompressiblePotentialFlowElement(CompressiblePotentialFlowElement const& rOther)
    : Element(rOther)  {
}

/**
 * Destructor
 */
CompressiblePotentialFlowElement::~CompressiblePotentialFlowElement() {
}

///@}
///@name Operators
///@{

/// Assignment operator.
CompressiblePotentialFlowElement & CompressiblePotentialFlowElement::operator=(CompressiblePotentialFlowElement const& rOther) {
  BaseType::operator=(rOther);
  Flags::operator =(rOther);
  // mpProperties = rOther.mpProperties;
  return *this;
}

///@}
///@name Operations
///@{

/**
 * ELEMENTS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer CompressiblePotentialFlowElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Element::Pointer(new CompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer CompressiblePotentialFlowElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Element::Pointer(new CompressiblePotentialFlowElement(NewId, pGeom, pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer CompressiblePotentialFlowElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const {
  KRATOS_TRY
  return Element::Pointer(new CompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
  KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void CompressiblePotentialFlowElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) 
{
  const unsigned int NumNodes = 3;
  if(this->IsNot(MARKER)) //normal element
  {
    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes, false);
    
    for (unsigned int i = 0; i < NumNodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
      
  }
  else//wake element
  {
    if (rResult.size() != 2*NumNodes)
        rResult.resize(2*NumNodes, false);
    
//     array_1d<double,NumNodes> distances = this->GetValue(ELEMENTAL_DISTANCES);
    array_1d<double,NumNodes> distances; 
    for(unsigned int i = 0; i<NumNodes; i++)
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    
    //positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if(distances[i] > 0)
            rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE,0).EquationId();
    }

    //negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if(distances[i] < 0)
            rResult[NumNodes+i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        else
            rResult[NumNodes+i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE,0).EquationId();
    }
  }


}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void CompressiblePotentialFlowElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) 
{
  const unsigned int NumNodes = 3;
  if(this->IsNot(MARKER)) //normal element
  {
    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);
    
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
  }
  else//wake element
  {
    if (rElementalDofList.size() != 2*NumNodes)
        rElementalDofList.resize(2*NumNodes);
    
//     const array_1d<double,NumNodes>& distances = this->GetValue(ELEMENTAL_DISTANCES);
    array_1d<double,NumNodes> distances; 
    for(unsigned int i = 0; i<NumNodes; i++)
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    
    //positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if(distances[i] > 0)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
    }

    //negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if(distances[i] < 0)
            rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
        else
            rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
    }
  }
}

/**
 * ELEMENTS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
 * they can be managed internally with a private method to do the same calculations
 * only once: MANDATORY
 */

/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void CompressiblePotentialFlowElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) 
{
    const unsigned int NumNodes = 3;
    const unsigned int Dim = 2;
    
    
    ElementalData<NumNodes,Dim> data;
    
    //calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
    
    //gather nodal data
    for(unsigned int i=0; i<NumNodes; i++)
    {
        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        data.distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    }
    
    if(this->Is(MARKER) == false)//normal element (non-wake)
    {        
        if(rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);
        if(rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes,false);
        rLeftHandSideMatrix.clear();
        rRightHandSideVector.clear();


        //decide if the element is split or not
//         unsigned int npos=0, nneg=0;
//         for(unsigned int i=0; i<NumNodes; i++)
//         {
//             if(data.distances[i] >= 0) npos++;
//             else nneg++;
//         }
//         bool is_split = false;
//         if(npos!=0 && nneg != 0)
//             is_split = true;
        
//         if(is_split == false) //if not split it is like a normal element
//         {
            bounded_matrix<double,NumNodes, NumNodes> lhs;
            array_1d<double,NumNodes> rhs;
            lhs.clear();
            rhs.clear();
            ComputeLHSGaussPointContribution(data.vol,lhs,data);
            ComputeRHSGaussPointContribution(data.vol,rhs,data);
            rLeftHandSideMatrix = lhs;
            rRightHandSideVector = rhs;
//         }
//         else
//             KRATOS_ERROR << "not yey implemented" ;
    }
    else //it is a wake element
    {
        //note that the lhs and rhs have double the size!!
        if(rLeftHandSideMatrix.size1() != 2*NumNodes || rLeftHandSideMatrix.size2() != 2*NumNodes)
            rLeftHandSideMatrix.resize(2*NumNodes,2*NumNodes,false);
        if(rRightHandSideVector.size() != 2*NumNodes)
            rRightHandSideVector.resize(2*NumNodes,false);
        rLeftHandSideMatrix.clear();
        rRightHandSideVector.clear();
        
//         data.distances = this->GetValue(ELEMENTAL_DISTANCES);
        
        //identify the normal to the wake. 
//         array_1d<double,Dim> n = prod(data.DN_DX, data.distances);
//         n/= norm_2(n);
        
        //compute the lhs and rhs that would correspond to it not being divided
        bounded_matrix<double,NumNodes, NumNodes> lhs;
//         array_1d<double,NumNodes> rhs;
        lhs.clear();
//         rhs.clear();
        double weight = data.vol; // * 1000.0; //applying a penalty factor
        ComputeLHSGaussPointContribution(weight,lhs,data);
//         ComputeRHSGaussPointContribution(weight,rhs,data);
        
        //assemble internally the positive part to the first Nnodes rows
        //and the negative part to the others
//         for(unsigned int i=0; i<NumNodes; i++)
//         {
//             for(unsigned int j=0; j<NumNodes; j++)
//             {
//                 rLeftHandSideMatrix(i,j) = lhs(i,j);
//                 rLeftHandSideMatrix(NumNodes+i,NumNodes+j) = lhs(i,j); //negative part
//             }
//             rRightHandSideVector(i) = rhs(i);
//             rRightHandSideVector(NumNodes+i) = rhs(i); //negative part
//         }
        
        //here impose the wake constraint
        const double k=1.0e3;
        for(unsigned int i=0; i<NumNodes; ++i)
        {
            for(unsigned int j=0; j<NumNodes; ++j)
            {
                 rLeftHandSideMatrix(i,j)                   =  (1.0+k)*lhs(i,j);
                 rLeftHandSideMatrix(i,j+NumNodes)          = -k*lhs(i,j);
                 rLeftHandSideMatrix(i+NumNodes,j)          = -k*lhs(i,j);
                 rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  (1.0+k)*lhs(i,j);
            }
        }
        Vector split_element_values(NumNodes*2);
        GetValuesOnSplitElement(split_element_values);
        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,split_element_values);
        
        
        
        
    }
    
//      KRATOS_WATCH(rLeftHandSideMatrix)
}

void CompressiblePotentialFlowElement::ComputeLHSGaussPointContribution(
                                       const double weight,
                                       bounded_matrix<double,3,3>& lhs,
                                       const ElementalData<3,2>& data)
{
    noalias(lhs) += weight*prod(data.DN_DX, trans(data.DN_DX));
}

void CompressiblePotentialFlowElement::ComputeRHSGaussPointContribution(
                                       const double weight,
                                       array_1d<double,3>& rhs,
                                       const ElementalData<3,2>& data)
{
    const unsigned int Dim = 2;
    array_1d<double,Dim> grad = prod(trans(data.DN_DX), data.phis);
    noalias(rhs) -= weight*prod(data.DN_DX, grad);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void CompressiblePotentialFlowElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) 
{
    //TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
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
int CompressiblePotentialFlowElement::Check(const ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY

  if (this->Id() < 1) {
    KRATOS_THROW_ERROR(std::logic_error, "CompressiblePotentialFlowElement found with Id 0 or negative","")
  }

  if (this->GetGeometry().Area() <= 0) {
    std::cout << "error on CompressiblePotentialFlowElement -> " << this->Id() << std::endl;
    KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0","")
  }

  for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable POSITIVE_FACE_PRESSURE on node ", this->GetGeometry()[i].Id() )
    }
    

  return 0;

  KRATOS_CATCH("");
}



void CompressiblePotentialFlowElement::GetValuesOnSplitElement(Vector& split_element_values )
{
    const unsigned int NumNodes = 3;
    array_1d<double,NumNodes> distances; 
    for(unsigned int i = 0; i<NumNodes; i++)
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    
//     Vector& distances = GetValue(ELEMENTAL_DISTANCES);
    //positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if(distances[i] > 0)
            split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        else
            split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
    }

    //negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if(distances[i] < 0)
            split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        else
            split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
    }
}

void CompressiblePotentialFlowElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        const int Dim = 2;
        const int NumNodes = 3;
        if(rValues.size() != 1) rValues.resize(1);
        
        if (rVariable == PRESSURE)
            {
                double p = 0.0;
                                    
                bool active = true;
                if ((this)->IsDefined(ACTIVE))
                    active = (this)->Is(ACTIVE);
                
                if(active && !this->Is(MARKER))
                {
                    const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY];
                    const double vinfinity_norm = norm_2(vinfinity);
                    
                    ElementalData<NumNodes,Dim> data;
        
                    //calculate shape functions
                    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
                    
                    //gather nodal data
                    for(unsigned int i=0; i<NumNodes; i++)
                    {
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                    }
                    
                    const array_1d<double,Dim> v = prod(trans(data.DN_DX), data.phis);
            
                        
                    p = (vinfinity_norm - norm_2(v))/vinfinity_norm; //0.5*(norm_2(vinfinity) - norm_2(v));
                }
                
                
                rValues[0] = p;
            }
        }

        void CompressiblePotentialFlowElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector<array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        const int Dim = 2;
        const int NumNodes = 3;
        if(rValues.size() != 1) rValues.resize(1);
        
        if (rVariable == VELOCITY)
            {
                bool active = true;
                if ((this)->IsDefined(ACTIVE))
                    active = (this)->Is(ACTIVE);
                
                array_1d<double,3> v;
                if(this->Is(MARKER) || active==false)
                    v.clear();
                else
                {
                    ElementalData<NumNodes,Dim> data;
        
                    //calculate shape functions
                    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
                    
                    //gather nodal data
                    for(unsigned int i=0; i<NumNodes; i++)
                    {
                        data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                    }
                    
                    array_1d<double,Dim> vaux = -prod(trans(data.DN_DX), data.phis);
                    v[0] = vaux[0];
                    v[1] = vaux[1];
                    v[2] = 0.0;
                }
                

                
                rValues[0] = v;
            }
        }
        
        void CompressiblePotentialFlowElement::ComputeSplitVolumes(array_1d<double,3>& distances, double& positive_vol, double& negative_vol)
        {
//             positive_vol = 0.0;
//             negative_vol = 0.0;
//             
//             const unsigned int nEdges = 3; //ATTENTION: this shall be initialized of size 6 in 3D and 3 in 2D
//             MatrixType Ncontainer;
//             VectorType gauss_volumes;
//             VectorType signs(nEdges); 
//             VectorType edge_areas(nEdges);
// 
//             // Splitting to determine the new Gauss pts.
//             //~ unsigned int ndivisions = ComputeSplitting(data, Ncontainer, gauss_volumes, DNvalues, Nenriched, signs, distances, edge_areas);
//             unsigned int ndivisions = ComputeSplitting(data, Ncontainer, gauss_volumes, signs, distances, edge_areas);
//             
//             for(unsigned int g=0; g<gauss_volumes.size(); ++g)
//             {
//                 double distance_on_gauss = inner_prod(row(Ncontainer,g), distances);
//                 if(distance_on_gauss > 0)
//                     positive_vol += gauss_volumes[g];
//                 else
//                     negative_vol += gauss_volumes[g];
//                     
//             }
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

std::string CompressiblePotentialFlowElement::Info() const {
  std::stringstream buffer;
  buffer << "CompressiblePotentialFlowElement #" << Id();
  return buffer.str();
}

/// Print information about this object.

void CompressiblePotentialFlowElement::PrintInfo(std::ostream& rOStream) const {
  rOStream << "CompressiblePotentialFlowElement #" << Id();
}

/// Print object's data.

void CompressiblePotentialFlowElement::PrintData(std::ostream& rOStream) const {
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

void CompressiblePotentialFlowElement::save(Serializer& rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

  // List
  // To be completed with the class member list
}

void CompressiblePotentialFlowElement::load(Serializer& rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

  // List
  // To be completed with the class member list
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
inline std::istream & operator >> (std::istream& rIStream, CompressiblePotentialFlowElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const CompressiblePotentialFlowElement& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.
