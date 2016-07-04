//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//




#if !defined(KRATOS_UPDATED_LAGRANGIAN_MPM_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_MPM_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "custom_elements/meshless_base_element.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

//#include "custom_geometries/meshless_geometry.h"

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

/// Short class definition.
/** Detail class definition.
*/
class UpdatedLagrangianMPMElement
    : public MeshlessBaseElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of UpdatedLagrangianMPMElement
    KRATOS_CLASS_POINTER_DEFINITION(UpdatedLagrangianMPMElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
	UpdatedLagrangianMPMElement(IndexType NewId, GeometryType::Pointer pGeometry);
	UpdatedLagrangianMPMElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~UpdatedLagrangianMPMElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;
    
    void  CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void  CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void  EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo );

    void  GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

    void  GetValuesVector( Vector& values, int Step );

    void  GetFirstDerivativesVector( Vector& values, int Step );
 
    void  GetSecondDerivativesVector( Vector& values, int Step );
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and outputElement
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
	// A private default constructor necessary for serialization
	UpdatedLagrangianMPMElement() : MeshlessBaseElement()
	{
	}

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

private:
    ///@name Static Member Variables
    ///@{
    

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;


    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MeshlessBaseElement);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MeshlessBaseElement);
    }
    
    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //UpdatedLagrangianMPMElement& operator=(const UpdatedLagrangianMPMElement& rOther);

    /// Copy constructor.
    //UpdatedLagrangianMPMElement(const UpdatedLagrangianMPMElement& rOther);


    ///@}

}; // Class UpdatedLagrangianMPMElement

///@}

///@name Type Definitions
///@{UpdatedLagrangianMPMElement


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    UpdatedLagrangianMPMElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const UpdatedLagrangianMPMElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_UPDATED_LAGRANGIAN_MPM_ELEMENT_H_INCLUDED  defined 


