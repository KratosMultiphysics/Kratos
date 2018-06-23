//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_SEARCH_STRUCTURE_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_STRUCTURE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_search_structure_base.h"


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

/// Rank local searching
/** This class provides features for rank local searching. It is basically a wrapper for the Bin Search Function
* It computes local neighbors and selects the best neighbor out of the search results. How these are selected is
* defined in the "EvaluateResult" function of the InterfaceObject
* If no neighbors are found, the search radius is increased by a factor ("increase_factor" in "Search")
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceSearchStructure : public InterfaceSearchStructureBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceSearchStructure
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceSearchStructure);

    using BaseType = InterfaceSearchStructureBase;

    using MapperInterfaceInfoUniquePointerType = BaseType::MapperInterfaceInfoUniquePointerType;

    using MapperLocalSystemPointerVectorPointer = BaseType::MapperLocalSystemPointerVectorPointer;

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceSearchStructure(ModelPart& rModelPartOrigin,
                             MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                             Parameters SearchSettings) :
        InterfaceSearchStructureBase(rModelPartOrigin,
                                 pMapperLocalSystems,
                                 SearchSettings)
    {
        KRATOS_INFO("InterfaceSearchStructure") << "In Ctor" << std::endl;
    }


    /// Destructor.
    virtual ~InterfaceSearchStructure() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceSearchStructure" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceSearchStructure";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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

    void InitializeSearchIteration(const Kratos::Flags& rOptions,
                                   const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo) override;

    void FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo) override;



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
    // InterfaceSearchStructure& operator=(InterfaceSearchStructure const& rOther){}

    //   /// Copy constructor.
    //   InterfaceSearchStructure(InterfaceSearchStructure const& rOther){}


    ///@}

}; // Class InterfaceSearchStructure

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_STRUCTURE_H_INCLUDED  defined
