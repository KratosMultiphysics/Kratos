// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#if !defined(KRATOS_IO_H_INCLUDED )
#define  KRATOS_IO_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"


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
class IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IO
    KRATOS_CLASS_POINTER_DEFINITION(IO);

	/// Local Flags
	KRATOS_DEFINE_LOCAL_FLAG( READ );
	KRATOS_DEFINE_LOCAL_FLAG( WRITE );
	KRATOS_DEFINE_LOCAL_FLAG( APPEND );
	KRATOS_DEFINE_LOCAL_FLAG( IGNORE_VARIABLES_ERROR );

    typedef Node<3> NodeType;

    typedef Mesh<NodeType, Properties, Element, Condition> MeshType;

    typedef MeshType::NodesContainerType NodesContainerType;

    typedef MeshType::PropertiesContainerType PropertiesContainerType;

    typedef MeshType::ElementsContainerType ElementsContainerType;

    typedef MeshType::ConditionsContainerType ConditionsContainerType;

    typedef std::vector<std::vector<std::size_t> > ConnectivitiesContainerType;

    typedef std::vector<std::vector<std::size_t> > PartitionIndicesContainerType;

    typedef std::vector<std::size_t> PartitionIndicesType;

    typedef std::size_t SizeType;

    typedef matrix<int> GraphType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IO() {}

    /// Destructor.
    virtual ~IO() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual bool ReadNode(NodeType& rThisNode)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class.", "")
    }

    virtual bool ReadNodes(NodesContainerType& rThisNodes)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual std::size_t ReadNodesNumber()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class.", "");
    }

    virtual void WriteNodes(NodesContainerType const& rThisNodes)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadProperties(Properties& rThisProperties)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadProperties(PropertiesContainerType& rThisProperties)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadElement(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, Element::Pointer& pThisElements)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual std::size_t  ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void WriteElements(ElementsContainerType const& rThisElements)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual std::size_t  ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadInitialValues(ModelPart& rThisModelPart)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }    
    
    virtual void ReadInitialValues(NodesContainerType& rThisNodes, ElementsContainerType& rThisElements, ConditionsContainerType& rThisConditions)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

    virtual void ReadMesh(MeshType & rThisMesh)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual void ReadModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    virtual std::size_t ReadNodalGraph(int** NodeIndices,int** NodeConnectivities)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "");
    }

    virtual void DivideInputToPartitions(SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                                         PartitionIndicesType const& NodesPartitions,
                                         PartitionIndicesType const& ElementsPartitions,
                                         PartitionIndicesType const& ConditionsPartitions,
                                         PartitionIndicesContainerType const& NodesAllPartitions,
                                         PartitionIndicesContainerType const& ElementsAllPartitions,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
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

//       /// Turn back information as a string.
//       virtual std::string Info() const
//  {
//    return "IO";
//  }

//       /// Print information about this object.
//       virtual void PrintInfo(std::ostream& rOStream) const
//  {
//    rOStream << "IO";
//  }


//       /// Print object's data.
//       virtual void PrintData(std::ostream& rOStream) const
//  {
//  }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
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
    IO& operator=(IO const& rOther);

    /// Copy constructor.
    IO(IO const& rOther);


    ///@}

}; // Class IO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  IO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const IO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@}


}  // namespace Kratos.

#endif // KRATOS_IO_H_INCLUDED  defined 
