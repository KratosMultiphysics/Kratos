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

#if !defined(KRATOS_REGISTRY_H_INCLUDED )
#define  KRATOS_REGISTRY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/registry_item.h"
#include "includes/registry_value_item.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/** Detail class definition.
*/
class Registry
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Registry
    KRATOS_CLASS_POINTER_DEFINITION(Registry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Registry(){}

    /// Destructor.
    virtual ~Registry(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    template< typename TItemType, class... TArgumentsList >
    static RegistryItem& AddItem(std::string const& ItemFullName, TArgumentsList&&... Arguments){

        auto item_path = SplitFullName(ItemFullName);

        RegistryItem* p_current_item = &GetRootRegistryItem();

        for(std::size_t i = 0 ; i < item_path.size() - 1 ; i++){
            auto& item_name = item_path[i];
            if(p_current_item->HasItem(item_name)){
                p_current_item = &p_current_item->GetItem(item_name);
            }
            else{
                p_current_item = &p_current_item->AddItem<RegistryItem>(item_name);
            }
        }

        // I am doing the last one out of the loop to create it with the given type and argument
        auto& item_name = item_path.back();
        if(p_current_item->HasItem(item_name)){
            KRATOS_ERROR << "The item \"" << ItemFullName << "\" is already registered." << std::endl;
        }
        else{
            p_current_item = &p_current_item->AddItem<TItemType>(item_name, std::forward<TArgumentsList>(Arguments)...);
        }

        return *p_current_item;
    }

    ///@}
    ///@name Access
    ///@{

    static RegistryItem& GetItem(std::string const& ItemFullName);


    static void RemoveItem(std::string const& ItemName);
    
    ///@}
    ///@name Inquiry
    ///@{

    static bool HasItem(std::string const& ItemFullName);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    std::string ToJson(std::string const& Indentation) const;


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

        static RegistryItem* mspRootRegistryItem;


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

        static RegistryItem& GetRootRegistryItem();
        static std::vector<std::string> SplitFullName(std::string const& FullName);


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Registry& operator=(Registry const& rOther);

    /// Copy constructor.
    Registry(Registry const& rOther);


    ///@}

}; // Class Registry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                Registry& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const Registry& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_REGISTRY_H_INCLUDED  defined
