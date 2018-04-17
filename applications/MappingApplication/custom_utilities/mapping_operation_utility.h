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

#if !defined(KRATOS_MAPPING_OPERATION_UTILITY_H)
#define  KRATOS_MAPPING_OPERATION_UTILITY_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


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
template<class TSparseSpace, class TDenseSpace>
class MappingOperationUtility
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingOperationUtility
    KRATOS_CLASS_POINTER_DEFINITION(MappingOperationUtility);

    using ModelPartPointerType = ModelPart::Pointer;
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingOperationUtility(ModelPartPointerType pInterfaceModelPart)
        : mpInterfaceModelPart(pInterfaceModelPart)
    {

    }

    /// Destructor.
    virtual ~MappingOperationUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void UpdateInterface()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void ExecuteMapping(const Variable<double>& rOriginVariable,
                                const Variable<double>& rDestinationVariable,
                                Kratos::Flags MappingOptions)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void ExecuteMapping(const Variable<array_1d<double, 3>>& rOriginVariable,
                                const Variable<array_1d<double, 3>>& rDestinationVariable,
                                Kratos::Flags MappingOptions)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
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
    virtual std::string Info() const
    {
        return "MappingOperationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPartPointerType mpInterfaceModelPart;


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
    // MappingOperationUtility& operator=(MappingOperationUtility const& rOther) {}

    /// Copy constructor.
    MappingOperationUtility(MappingOperationUtility const& rOther) {}

    ///@}

    }; // Class MappingOperationUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPING_OPERATION_UTILITY_H  defined
