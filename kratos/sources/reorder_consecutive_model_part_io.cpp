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
