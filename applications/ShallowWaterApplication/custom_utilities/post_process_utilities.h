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

#ifndef KRATOS_POST_PROCESS_UTILITIES_H_INCLUDED
#define KRATOS_POST_PROCESS_UTILITIES_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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

/// Auxiliary utilities for post process purpose.
/** Gid is able to representate different entities when they have diffrent properties.
 *  For each existing property, an auxiliary property is created, which means dry state.
 *  Gid will print the wet domain (original property) with a property-color, and the dry
 *  domain with another color.
 */
class PostProcessUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    /// Pointer definition of PostProcessUtilities
    KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PostProcessUtilities(ModelPart& rThisModelPart) : mrModelPart(rThisModelPart) {}

    /// Destructor.
    ~PostProcessUtilities() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void DefineAuxiliaryProperties()
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
            Properties::Pointer wet_prop = mrModelPart.pGetProperties(id);
            Properties::Pointer dry_prop(new Properties(*wet_prop));
            dry_prop->SetId(++max_id);

            // Add the new properties to the model part and add them to the maps
            mrModelPart.AddProperties(dry_prop);
            mWetToDryPropertiesMap[wet_prop->Id()] = max_id;
            mDryToWetPropertiesMap[dry_prop->Id()] = id;
        }
    }

    void AssignDryWetProperties(Flags Flag)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto elem = mrModelPart.ElementsBegin() + i;

            if (elem->Is(Flag))
            {
                auto search = mDryToWetPropertiesMap.find(elem->GetProperties().Id());
                if (search != mDryToWetPropertiesMap.end()) // The element was dry
                {
                    IndexType wet_prop_id = search->second;
                    elem->SetProperties(mrModelPart.pGetProperties(wet_prop_id));
                }
            }
            else
            {
                auto search = mWetToDryPropertiesMap.find(elem->GetProperties().Id());
                if (search != mWetToDryPropertiesMap.end()) // The element was wet
                {
                    IndexType dry_prop_id = search->second;
                    elem->SetProperties(mrModelPart.pGetProperties(dry_prop_id));
                }
            }
        }
    }

    void RestoreDryWetProperties()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto elem = mrModelPart.ElementsBegin() + i;

            auto search = mDryToWetPropertiesMap.find(elem->GetProperties().Id());
            if (search != mDryToWetPropertiesMap.end()) // The element was dry
            {
                IndexType wet_prop_id = search->second;
                elem->SetProperties(mrModelPart.pGetProperties(wet_prop_id));
            }
        }
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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


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

    ModelPart& mrModelPart;

    std::map<IndexType, IndexType> mWetToDryPropertiesMap;
    std::map<IndexType, IndexType> mDryToWetPropertiesMap;

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
    PostProcessUtilities& operator=(PostProcessUtilities const& rOther);

    /// Copy constructor.
    PostProcessUtilities(PostProcessUtilities const& rOther);


    ///@}

}; // Class PostProcessUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                PostProcessUtilities& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const PostProcessUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POST_PROCESS_UTILITIES_H_INCLUDED  defined
