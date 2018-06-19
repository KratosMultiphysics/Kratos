//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on work of Pablo Becker)
//

#if !defined(KRATOS_POINT_LOCATOR_H_INCLUDED)
#define  KRATOS_POINT_LOCATOR_H_INCLUDED


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

/**
 * @class PointLocator
 * @ingroup KratosCore
 * @brief Utility class to find an entity of a mesh based on a location
 * @details Based on the location of a point, the corresponding entity
 * (node, element or condition) in a mesh is found and it's id is returned
 * @author Philipp Bucher
 */
class KRATOS_API(KRATOS_CORE) PointLocator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointLocator
    KRATOS_CLASS_POINTER_DEFINITION(PointLocator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLocator(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

    /// Destructor.
    virtual ~PointLocator() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function finds a node based on a location
     * @param rThePoint the location to search
     * @param DistanceThreshold threshold for the distance
     * @return Id of the found node. -1 if no node was found
     */
    int FindNode(const Point& rThePoint, const double DistanceThreshold=1e-12) const;

    /**
     * @brief This function finds an element based on a location
     * @param rThePoint the location to search
     * @param rShapeFunctionValues vector containing the shape-function values for the given point
     * @return Id of the found element. -1 if no element was found
     */
    int FindElement(const Point& rThePoint, Vector& rShapeFunctionValues) const;

    /**
     * @brief This function finds a condition based on a location
     * @param rThePoint the location to search
     * @param rShapeFunctionValues vector containing the shape-function values for the given point
     * @return Id of the found condition. -1 if no condition was found
     */
    int FindCondition(const Point& rThePoint, Vector& rShapeFunctionValues) const;

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
        buffer << "PointLocator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PointLocator";}

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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function finds an object based on a location
     * @param rObjects the objects to search
     * @param rObjectType type of the object => "Element"/"Condition"
     * @param rThePoint the location to search
     * @param rObjectId Id of the found condition. -1 if no object was found
     * @param rShapeFunctionValues vector containing the shape-function values for the given point
     */
    template<typename TObjectType>
    void FindObject(const TObjectType& rObjects, const std::string& rObjectType,
                    const Point& rThePoint, int& rObjectId, Vector& rShapeFunctionValues) const;

    /**
     * @brief This function performs some checks after the search
     * @param rObjectType type of the object => "Node"/"Element"/"Condition"
     * @param rThePoint the location to search
     * @param LocalObjectsFound number of found results
     */
    void CheckResults(const std::string& rObjectType,
                      const Point& rThePoint,
                      int LocalObjectsFound) const;

    /**
     * @brief This function checks whether a node is close to a point based on a threshold
     * @param rNode the node to check
     * @param rThePoint the location to search
     * @param DistanceThreshold threshold for the distance
     * @return whether the rNode is close to rThePoint
     */
    bool NodeIsCloseEnough(const Node<3>& rNode,
                           const Point& rThePoint,
                           double DistanceThreshold) const;

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class PointLocator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const PointLocator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POINT_LOCATOR_H_INCLUDED  defined
