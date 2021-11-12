//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "post_process_utilities.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{

void PostProcessUtilities::DefineAuxiliaryProperties()
{
    mWetToDryPropertiesMap.clear();
    mDryToWetPropertiesMap.clear();

    IndexType max_id = 0;
    IndexVectorType prop_ids;

    // Loop the original ids
    for (IndexType i = 0; i < mrModelPart.NumberOfProperties(); ++i)
    {
        auto prop_id = (mrModelPart.PropertiesBegin() + i)->Id();
        max_id = std::max(max_id, prop_id);
        prop_ids.push_back(prop_id);
    }

    // Creating the auxiliary ids for the dry domain
    for (auto id : prop_ids)
    {
        // Get pointers to the properties and create the dry property
        auto wet_prop = mrModelPart.pGetProperties(id);
        auto dry_prop = Kratos::make_shared<Properties>(*wet_prop);
        dry_prop->SetId(++max_id);

        // Add the new properties to the model part and add them to the maps
        mrModelPart.AddProperties(dry_prop);
        mWetToDryPropertiesMap[wet_prop->Id()] = max_id;
        mDryToWetPropertiesMap[dry_prop->Id()] = id;
    }
}

void PostProcessUtilities::AssignDryWetProperties(Flags WetStateFlag)
{
    block_for_each(mrModelPart.Elements(), [&](Element& rElem){
        if (rElem.Is(WetStateFlag))
        {
            auto search = mDryToWetPropertiesMap.find(rElem.GetProperties().Id());
            if (search != mDryToWetPropertiesMap.end()) // The element was dry
            {
                IndexType wet_prop_id = search->second;
                rElem.SetProperties(mrModelPart.pGetProperties(wet_prop_id));
            }
        }
        else
        {
            auto search = mWetToDryPropertiesMap.find(rElem.GetProperties().Id());
            if (search != mWetToDryPropertiesMap.end()) // The element was wet
            {
                IndexType dry_prop_id = search->second;
                rElem.SetProperties(mrModelPart.pGetProperties(dry_prop_id));
            }
        }
    });
}

void PostProcessUtilities::RestoreDryWetProperties()
{
    block_for_each(mrModelPart.Elements(), [&](Element& rElem){
        auto search = mDryToWetPropertiesMap.find(rElem.GetProperties().Id());
        if (search != mDryToWetPropertiesMap.end()) // The element was dry
        {
            IndexType wet_prop_id = search->second;
            rElem.SetProperties(mrModelPart.pGetProperties(wet_prop_id));
        }
    });
}

}  // namespace Kratos.
