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

#ifndef KRATOS_MOVE_MESH_UTILITY_H_INCLUDED
#define KRATOS_MOVE_MESH_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


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

/**
 * @class MoveMeshUtility
 * @ingroup KratosShallowWaterApplication
 * @brief Tools for lagrangian computations
 * @details Move the computational mesh over a background mesh and map data between them
 * @author Miguel Maso Sotomayor
 */
class MoveMeshUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MoveMeshUtility
    KRATOS_CLASS_POINTER_DEFINITION(MoveMeshUtility);

    typedef Node<3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * @brief Constructor
    */
    MoveMeshUtility(ModelPart& rLagrangianModelPart, ModelPart& rEulerianModelPart) :
        mrLagrangianModelPart(rLagrangianModelPart),
        mrEulerianModelPart(rEulerianModelPart)
    {}

    /**
    * @brief Destructor
    */
    virtual ~MoveMeshUtility(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    int Check();

    void MoveMesh();

    template<class TDataType>
    void MapToEulerian(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void MapToLagrangian(const Variable<TDataType>& rVariable);

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
        std::stringstream buffer;
        buffer << "MoveMeshUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << Info();}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrLagrangianModelPart;
    ModelPart& mrEulerianModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void MoveNode(NodeType& rNode);

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
    MoveMeshUtility& operator=(MoveMeshUtility const& rOther);

    /// Copy constructor.
    MoveMeshUtility(MoveMeshUtility const& rOther);

    ///@}

}; // Class MoveMeshUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                MoveMeshUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MoveMeshUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOVE_MESH_UTILITY_H_INCLUDED  defined
