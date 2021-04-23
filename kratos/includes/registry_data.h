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

#if !defined(KRATOS_REGISTRY_DATA_H_INCLUDED )
#define  KRATOS_REGISTRY_DATA_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <unordered_map>


// External includes


// Project includes
#include "includes/define.h"


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

/// The registry data to be stored by Registry class. It is the base class for some more specific ones.
/** RegistryData has a tree node structure and stores its name, an optional
 *  value, and an unorder_set of its sub data. 
 *  This structure let us to have registry of the elements and then different
 *  registries for each elements inside it.
 *  Please note that RegistryData stores a pointer to the value. 
 *  To have a copy of the value you may use the derived RegistryCopyOfData 
 *  which crates a copy in construction and delete it in its destructor 
 *  to make the memory management easier.
*/
class RegistryData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RegistryData
    KRATOS_CLASS_POINTER_DEFINITION(RegistryData);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    RegistryData() = delete;

    /// Constructor with the name
    RegistryData(std::string Name) : mName(Name), mpValue(nullptr){}

    /// Destructor.
    virtual ~RegistryData(){
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

        const std::string& Name() const{
            return mName;
        }


    ///@}
    ///@name Inquiry
    ///@{

        bool HasValue() const{
            return (mpValue != nullptr);
        }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


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

    std::string mName;
    void* mpValue;
    std::unordered_map<std::string, Kratos::unique_ptr<RegistryData>> mSubRegistryData;



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


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RegistryData& operator=(RegistryData const& rOther);

    /// Copy constructor.
    RegistryData(RegistryData const& rOther);


    ///@}

}; // Class RegistryData

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                RegistryData& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const RegistryData& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_REGISTRY_DATA_H_INCLUDED  defined
