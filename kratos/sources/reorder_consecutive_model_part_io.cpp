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


// Project includes
#include "includes/reorder_consecutive_model_part_io.h"



namespace Kratos
{
    /// Constructor with  filenames.
    ReorderConsecutiveModelPartIO::ReorderConsecutiveModelPartIO(std::string const& Filename, const Flags Options )
        : ModelPartIO(Filename, Options),
		mNumberOfNodes(0),
		mNumberOfElements(0),
		mNumberOfConditions(0),
		mNodeIdMap(),
		mElementIdMap(),
		mConditionIdMap()
    {
    }


    /// Destructor.
    ReorderConsecutiveModelPartIO::~ReorderConsecutiveModelPartIO() {}

	

	ModelPartIO::SizeType ReorderConsecutiveModelPartIO::ReorderedNodeId(ModelPartIO::SizeType NodeId)
	{
		IdMapType::iterator i = mNodeIdMap.find(NodeId);
		if(i != mNodeIdMap.end())
			return i->second;

		mNodeIdMap.insert(IdMapType::value_type(NodeId, ++mNumberOfNodes));
		return mNumberOfNodes;
	}
	ModelPartIO::SizeType ReorderConsecutiveModelPartIO::ReorderedElementId(ModelPartIO::SizeType ElementId)
	{
		IdMapType::iterator i = mElementIdMap.find(ElementId);
		if(i != mElementIdMap.end())
			return i->second;

		mElementIdMap.insert(IdMapType::value_type(ElementId, ++mNumberOfElements));
		return mNumberOfElements;
	}

	ModelPartIO::SizeType ReorderConsecutiveModelPartIO::ReorderedConditionId(ModelPartIO::SizeType ConditionId)
	{
		IdMapType::iterator i = mConditionIdMap.find(ConditionId);
		if(i != mConditionIdMap.end())
			return i->second;

		mConditionIdMap.insert(IdMapType::value_type(ConditionId, ++mNumberOfConditions));
		return mNumberOfConditions;
	}



}  // namespace Kratos.
