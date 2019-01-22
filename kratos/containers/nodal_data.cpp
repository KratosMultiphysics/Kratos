//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes


// External includes


// Project includes
#include "containers/nodal_data.h"


namespace Kratos
{

    NodalData::NodalData(IndexType TheId) 
    : mId(TheId)
    , mSolutionStepsNodalData()
    {

    }

    NodalData::NodalData(IndexType TheId, VariablesList*  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize) 
    : mId(TheId)
    , mSolutionStepsNodalData(pVariablesList,ThisData,NewQueueSize)
    {

    }

    /// Assignment operator.
    NodalData& NodalData::operator=(NodalData const& rOther){
        mId = rOther.mId;
        mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;

        return *this;
    }

    /// Turn back information as a string.
    std::string NodalData::Info() const 
    {
        return "NodalData";
    }

    /// Print information about this object.
    void NodalData::PrintInfo(std::ostream& rOStream) const 
    {
        rOStream << Info();
    }

    /// Print object's data.
    void NodalData::PrintData(std::ostream& rOStream) const 
    {
    }


    void NodalData::save(Serializer& rSerializer) const
    {
        rSerializer.save("Id",mId);
         const SolutionStepsNodalDataContainerType* pSolutionStepsNodalData = &mSolutionStepsNodalData;
        // I'm saving it as pointer so the dofs pointers will point to it as stored pointer. Pooyan.
        rSerializer.save("Solution Steps Nodal Data", pSolutionStepsNodalData);
   }

    void NodalData::load(Serializer& rSerializer)
    {
        rSerializer.load("Id",mId);
        SolutionStepsNodalDataContainerType* pSolutionStepsNodalData = &mSolutionStepsNodalData;
        rSerializer.load("Solution Steps Nodal Data", pSolutionStepsNodalData);
    }
}  // namespace Kratos.


