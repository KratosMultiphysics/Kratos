/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-10-29 14:26:54 $
//   Revision:            $Revision: 1.1 $
//
//


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
