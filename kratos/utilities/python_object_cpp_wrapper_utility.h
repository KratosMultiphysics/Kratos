//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(PYTHON_OBJECT_CPP_WRAPPER_UTILITY )
#define  PYTHON_OBJECT_CPP_WRAPPER_UTILITY

// System includes
#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/serializer.h"

namespace Kratos
{
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

/**
 * @class PythonObjectCppWrapperUtility
 * @ingroup KratosCore
 * @brief This is a class to wrapp python objects into c++ (using pybind)
 * @details This class is used in order to interoperate between c++ and python
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) PythonObjectCppWrapperUtility
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Counted pointer of PythonObjectCppWrapperUtility
    KRATOS_CLASS_POINTER_DEFINITION( PythonObjectCppWrapperUtility );

    /// The object type in python
    typedef pybind11::object ObjectType;
    
    /// The list [] of python
    typedef pybind11::list     ListType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PythonObjectCppWrapperUtility()= default;

    /**
     * @brief Constructor using the name of the file
     * @param rNameFile The name of the file used to construct the object
     */
    PythonObjectCppWrapperUtility(const std::string& rNameFile);

    /**
     * @brief Constructor using a list of objects
     * @param rObjectList List of objects that will be used to build the vector of objects
     */
    PythonObjectCppWrapperUtility(ListType& rObjectList);

    /**
     * @brief Constructor using just one process
     * @param rObject The process that will be added  at the begining of the vector of objects
     * @note This constructor overrides the previous ones ("everything" is an object, so here I try first to work as it is a list)
     */
    PythonObjectCppWrapperUtility(ObjectType& rObject);

    /// Destructor.
    virtual ~PythonObjectCppWrapperUtility()= default;

    ///@}
    ///@name Operators
    ///@{
    
    /// Assignment operator.
    
    PythonObjectCppWrapperUtility& operator=(PythonObjectCppWrapperUtility const& rOther)
    = default;
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief It add new objects to the existing object list
     * @param rObject The object that will be appended at the vector of objects
     */

    void AddObject(ObjectType& rObject);
    
    /**
     * @brief It add new objects to the existing objects list
     * @param rObjectList List of objects that will be appended the vector of objects
     */

    void AddObjects(ListType& rObjectList);

    /**
     * @brief It executes the method considered in the input
     * @param rNameMethod The method to be executed
     */

    void Execute(const std::string& rNameMethod);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Flags
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
        std::stringstream buffer;
        buffer << "PythonObjectCppWrapperUtility. Number of processes:" << mListPythonObjects.size();
        return buffer.str();
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PythonObjectCppWrapperUtility. Number of processes:" << mListPythonObjects.size();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "PythonObjectCppWrapperUtility. Number of processes:" << mListPythonObjects.size();
    }

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

    std::vector<ObjectType> mListPythonObjects; /// A vector containing the list of objects to be executed

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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        // TODO: Fill if necessary
//         rSerializer.save("Processes", mListPythonObjects);
    }

    virtual void load(Serializer& rSerializer)
    {
        // TODO: Fill if necessary
//         rSerializer.load("Processes", mListPythonObjects);
    }
    
    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class PythonObjectCppWrapperUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  PythonObjectCppWrapperUtility& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const PythonObjectCppWrapperUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // PYTHON_OBJECT_CPP_WRAPPER_UTILITY  defined
