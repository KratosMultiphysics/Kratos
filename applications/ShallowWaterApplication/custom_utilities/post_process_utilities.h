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
/** Gid is able to representate different entities when they have different Id and diffrent properties.
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

    void ApplyIdOffsetOnNodes()
    {
        ApplyIdOffset(mrModelPart.Nodes());
    }

    void UndoIdOffsetOnNodes()
    {
        UndoIdOffset(mrModelPart.Nodes());
    }

    void ApplyIdOffsetOnElements()
    {
        ApplyIdOffset(mrModelPart.Elements());
    }

    void UndoIdOffsetOnElements()
    {
        UndoIdOffset(mrModelPart.Elements());
    }

    void ApplyIdOffsetOnConditions()
    {
        ApplyIdOffset(mrModelPart.Conditions());
    }

    void UndoIdOffsetOnConditions()
    {
        UndoIdOffset(mrModelPart.Conditions());
    }

    void DefineAuxiliaryProperties()
    {
        mFluidToSolidPropertiesMap.clear();
        mSolidToFluidPropertiesMap.clear();
        mWetToDryPropertiesMap.clear();
        mDryToWetPropertiesMap.clear();

        // Create a two copies for each property
        const int nprop = static_cast<int>(mrModelPart.NumberOfProperties());
        ModelPart::PropertiesContainerType::iterator prop_begin = mrModelPart.PropertiesBegin();

        IndexType last_id = 0;
        IndexVectorType prop_ids;

        // Loop the original ids
        for (int i = 0; i < nprop; ++i)
        {
            auto prop = prop_begin + i;

            if (prop->Id() > last_id) {
                last_id = prop->Id();
            }
            prop_ids.push_back(prop->Id());
        }

        // Creating the auxiliary ids for the solid domain
        for (auto id : prop_ids)
        {
            // Get pointers to the properties and create the dry property
            Properties::Pointer fluid_prop = mrModelPart.pGetProperties(id);
            Properties::Pointer solid_prop(new Properties(*fluid_prop));
            solid_prop->SetId(++last_id);

            // Add the new property and add them to the maps
            mrModelPart.AddProperties(solid_prop);
            mFluidToSolidPropertiesMap[fluid_prop->Id()] = last_id;
            mSolidToFluidPropertiesMap[solid_prop->Id()] = id;
        }

        // Creating the auxiliary ids for the dry domain
        for (auto id : prop_ids)
        {
            // Get pointers to the properties and create the dry property
            Properties::Pointer wet_prop = mrModelPart.pGetProperties(id);
            Properties::Pointer dry_prop(new Properties(*wet_prop));
            dry_prop->SetId(++last_id);

            // Add the new property and add them to the maps
            mrModelPart.AddProperties(dry_prop);
            mWetToDryPropertiesMap[wet_prop->Id()] = last_id;
            mDryToWetPropertiesMap[dry_prop->Id()] = id;
        }

    }

    void AssignSolidFluidProperties()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto it_elem = mrModelPart.ElementsBegin() + i;
            IndexType solid_prop_id = mFluidToSolidPropertiesMap[it_elem->GetProperties().Id()];
            it_elem->SetProperties(mrModelPart.pGetProperties(solid_prop_id));
        }
    }

    void RestoreSolidFluidProperties()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto it_elem = mrModelPart.ElementsBegin() + i;
            IndexType fluid_prop_id = mSolidToFluidPropertiesMap[it_elem->GetProperties().Id()];
            it_elem->SetProperties(mrModelPart.pGetProperties(fluid_prop_id));
        }
    }

    void AssignDryWetProperties(Flags& rFlag)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto elem = mrModelPart.ElementsBegin() + i;

            if (elem->Is(rFlag))
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

    void SetBathymetryMeshPosition()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            it_node->Z() = it_node->FastGetSolutionStepValue(BATHYMETRY);
        }
    }

    void SetFreeSurfaceMeshPosition()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            it_node->Z() = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        }
    }

    void SetToZeroMeshPosition()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            it_node->Z() = 0.0;
        }
    }

    void ResetMeshPosition()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            it_node->Z() = it_node->Z0();
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

    std::map<IndexType, IndexType> mFluidToSolidPropertiesMap;
    std::map<IndexType, IndexType> mSolidToFluidPropertiesMap;

    std::map<IndexType, IndexType> mWetToDryPropertiesMap;
    std::map<IndexType, IndexType> mDryToWetPropertiesMap;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    template<class TContainerType>
    void ApplyIdOffset(TContainerType& rContainer)
    {
        IndexType offset = rContainer.size();
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int> (rContainer.size()); ++i)
        {
            auto it_cont = rContainer.begin() + i;
            it_cont->SetId(it_cont->Id() + offset);
        }
    }

    template<class TContainerType>
    void UndoIdOffset(TContainerType& rContainer)
    {
        int offset = static_cast<int> (rContainer.size());
        #pragma omp parallel for
        for (int i = 0; i < offset; ++i)
        {
            auto it_cont = rContainer.begin() + i;
            it_cont->SetId(it_cont->Id() - offset);
        }
    }

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
