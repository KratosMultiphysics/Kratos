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

    NodalData::NodalData(IndexType TheId) : mId(TheId){}

    /// Assignment operator.
    NodalData& NodalData::operator=(NodalData const& rOther){
        mId = rOther.mId;

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
    }

    void NodalData::load(Serializer& rSerializer)
    {
        rSerializer.load("Id",mId);
    }
}  // namespace Kratos.


