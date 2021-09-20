// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_NODE_CONFIGURE_FOR_NODE_SEARCH_INCLUDED_H)
#define  KRATOS_NODE_CONFIGURE_FOR_NODE_SEARCH_INCLUDED_H

// System includes
#include <string>
#include <iostream>
#include <cmath>

// Kratos includes
#include "spatial_containers/spatial_search.h"

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class NodeConfigureForNodeSearch
 * @ingroup StructuralMechanicsApplication
 * @brief Configuration file for nodes.
 * @details This class provides a configuration file for nodes to perform a node search
 * depending on the euclidean distance between two nodes.
 * @author Manuel Messmer
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NodeConfigureForNodeSearch {
public:
    /// Pointer definition of NodeConfigureForNodeSearch
    KRATOS_CLASS_POINTER_DEFINITION(NodeConfigureForNodeSearch);

    /**
     * @brief Compile time definitions
     * @param Epsilon Error tolerance for comparison operations with doubles
     * @param Dimension Dimension of the problem. Fixed to 3.
     **/
    static constexpr auto Epsilon   = std::numeric_limits<double>::epsilon();
    static constexpr auto Dimension = 3;

    ///@name Type Definitions
    ///@{
    typedef SpatialSearch                                           SearchType;

    typedef SearchType::PointType                                   PointType;
    typedef SearchType::NodesContainerType::ContainerType           ContainerType;
    typedef SearchType::NodesContainerType                          NodesContainerType;

    typedef SearchType::NodeType                                    NodeType;
    typedef ContainerType::value_type                               PointerType;
    typedef ContainerType::iterator                                 IteratorType;

    typedef SearchType::NodesContainerType::ContainerType           ResultContainerType;

    typedef ResultContainerType::iterator                           ResultIteratorType;
    typedef std::vector<double>::iterator                           DistanceIteratorType;

    typedef ContactPair<PointerType>                                ContactPairType;
    typedef std::vector<ContactPairType>                            ContainerContactType;
    typedef ContainerContactType::iterator                          IteratorContactType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    NodeConfigureForNodeSearch(){}

    /// Default destructor
    virtual ~NodeConfigureForNodeSearch(){}

    ///@}
    ///@name Operations
    ///@{

    /** @brief Calculates the bounding box for the given object.
     * @details For this configuation file, the bounding box is equal to the point given in 'rObject'.
     * @param rObject Point for which the bounding box will be calculated.
     * @param rLowPoint Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     **/
    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rHighPoint = rLowPoint  = *rObject;
    }

    /** @brief Calculates the bounding box for the given object extended with a Radius.
     * @details For this configuation file, the bounding box is equal to the point given in 'rObject' + - a radius.
     * @param rObject Point for which the bounding box will be calculated.
     * @param rLowPoint Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     * @param Radius The extension radius to be applied to the boundingbox.
     **/
    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double Radius)
    {
        auto radiusExtension = PointType(Radius, Radius, Radius);

        rLowPoint  = PointType{*rObject - radiusExtension};
        rHighPoint = PointType{*rObject + radiusExtension};
    }

    /** @brief Calculates the Center of the object.
     * @param rObject Point for which the bounding box will be calculated.
     * @param rCenter The center point of the object.
     **/
    static inline void CalculateCenter(const PointerType& rObject, PointType& rCenter)
    {
        rCenter  = *rObject;
    }

    /** @brief Tests the intersection of two objects extended with a given radius.
     * @details For this configuation file, tests if euclidean distance between the two nodes
     * is smaller than the provided 'Radius' within an Epsilon tolerance range.
     * @param rObj_1 First point of the test.
     * @param rObj_2 Second point of the test.
     * @param Radius The extension radius to be applied in the intersection.
     * @return Boolean indicating the result of the intersection test.
     **/
    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double Radius)
    {
        double distance;
        Distance(rObj_1, rObj_2, distance);

        if( distance > Epsilon + Radius){
            return false;
        }

        return true;
    }

    /** @brief Tests the intersection of one object with a boundingbox described by 'rLowPoint' and 'rHighPoint'.
     * @details For this configuation file, tests if one point is inside the boundingbox
     * described by 'rLowPoint' and 'rHighPoint' within an Epsilon tolerance range.
     * @param rObject Point of the tests.
     * @param rLowPoint Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     * @return Boolean indicating the result of the intersection test.
     **/
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        for(std::size_t i = 0; i < Dimension; i++) {
            if( (*rObject)[i] < rLowPoint[i] - Epsilon || (*rObject)[i] > rHighPoint[i] + Epsilon) {
                return false;
            }
        }

        return true;
    }

    /** @brief Tests the intersection of one object with a boundingbox described by 'rLowPoint' and 'rHighPoint'.
     * @details For this configuation file, tests if one point extended by radius is inside the boundingbox
     * described by 'rLowPoint' and 'rHighPoint' within an Epsilon tolerance range.
     * @param rObject Point of the tests.
     * @param rLowPoint Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     * @param Radius The extension radius to be applied in the intersection.
     * @return Boolean indicating the result of the intersection test.
     **/
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double Radius)
    {
        for(std::size_t i = 0; i < Dimension; i++) {
            if( ((*rObject)[i] + Radius) < rLowPoint[i] - Epsilon || ((*rObject)[i] - Radius) > rHighPoint[i] + Epsilon) {
                return false;
            }
        }

        return true;
    }

    /** @brief Calculates the distance between two objects.
     * @param rObj_1 First point.
     * @param rObj_2 Second point
     * @param distance The euclidean distance between 'rObj_1' and 'rObj_2'.
     **/
    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        double pwdDistance = 0.0f;

        for(std::size_t i = 0; i < Dimension; i++) {
            pwdDistance += std::pow((*rObj_1)[i] - (*rObj_2)[i], 2);
        }

        distance = std::sqrt(pwdDistance);
    }

    /// Turn back information as a string.
    virtual std::string Info() const {return " Spatial Containers Configure for Nodes to perform a Node Search"; }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
private:

    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NodeConfigureForNodeSearch& operator=(NodeConfigureForNodeSearch const& rOther);

    /// Copy constructor.
    NodeConfigureForNodeSearch(NodeConfigureForNodeSearch const& rOther);

    ///@}

    }; // Class NodeConfigureForNodeSearch

///@} Kratos classes

    ///@name Input and output
    ///@{

    /// input stream function

    inline std::istream& operator >> (std::istream& rIStream, NodeConfigureForNodeSearch& rThis){
        return rIStream;
        }

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const NodeConfigureForNodeSearch& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }

    ///@}

///@} addtogroup block

}   // namespace Kratos.
#endif	/* KRATOS_NODE_CONFIGURE_FOR_NODE_SEARCH_INCLUDED_H */