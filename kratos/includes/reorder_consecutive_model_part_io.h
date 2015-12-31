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



#if !defined(KRATOS_REORDER_CONSECUTIVE_MODEL_PART_IO_H_INCLUDED )
#define  KRATOS_REORDER_CONSECUTIVE_MODEL_PART_IO_H_INCLUDED



// System includes
#include <string>
#include <fstream>
#include <set>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part_io.h"
#include "utilities/timer.h"


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

/// An IO class for reading and writing a modelpart
/** This class writes all modelpart data including the meshes.
*/
class ReorderConsecutiveModelPartIO : public ModelPartIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReorderConsecutiveModelPartIO
    KRATOS_CLASS_POINTER_DEFINITION(ReorderConsecutiveModelPartIO);

    typedef ModelPartIO BaseType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::MeshType MeshType;

    typedef BaseType::NodesContainerType NodesContainerType;

    typedef BaseType::PropertiesContainerType PropertiesContainerType;

    typedef BaseType::ElementsContainerType ElementsContainerType;

    typedef BaseType::ConditionsContainerType ConditionsContainerType;

    typedef BaseType::ConnectivitiesContainerType ConnectivitiesContainerType;

    typedef BaseType::OutputFilesContainerType OutputFilesContainerType;

	typedef std::map<SizeType, SizeType> IdMapType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with  filenames.
	ReorderConsecutiveModelPartIO(std::string const& Filename, const Flags Options = IO::READ|IO::NOT_IGNORE_VARIABLES_ERROR);


    /// Destructor.
    virtual ~ReorderConsecutiveModelPartIO();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

 
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


	virtual SizeType ReorderedNodeId(ModelPartIO::SizeType NodeId);
	virtual SizeType ReorderedElementId(ModelPartIO::SizeType ElementId);
	virtual SizeType ReorderedConditionId(ModelPartIO::SizeType ConditionId);

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

    ///@name Member Variables
    ///@{
	
	SizeType mNumberOfNodes;
	SizeType mNumberOfElements;
	SizeType mNumberOfConditions;

	IdMapType mNodeIdMap;
	IdMapType mElementIdMap;
	IdMapType mConditionIdMap;

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
    ReorderConsecutiveModelPartIO& operator=(ReorderConsecutiveModelPartIO const& rOther);

    /// Copy constructor.
    ReorderConsecutiveModelPartIO(ReorderConsecutiveModelPartIO const& rOther);


    ///@}

}; // Class ReorderConsecutiveModelPartIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  ReorderConsecutiveModelPartIO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const ReorderConsecutiveModelPartIO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@}


}  // namespace Kratos.

#endif // KRATOS_REORDER_CONSECUTIVE_MODEL_PART_IO_H_INCLUDED  defined 
