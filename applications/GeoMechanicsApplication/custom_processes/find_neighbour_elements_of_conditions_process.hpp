// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//


#ifndef KRATOS_FIND_NEIGHBOUR_ELEMENTS_OF_CONDITIONS_PROCESS
#define KRATOS_FIND_NEIGHBOUR_ELEMENTS_OF_CONDITIONS_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "includes/kratos_flags.h"


namespace Kratos
{
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class FindNeighbourElementsOfConditionsProcess
 * @ingroup Kratos.GeoMechanicsApplication
 * @brief Finds list of elements attached to conditions.
 * @author Vahid Galavi
 */
class FindNeighbourElementsOfConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(FindNeighbourElementsOfConditionsProcess);

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the node type
    typedef Node<3> NodeType;

    // Definition of the geometry
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for FindNeighbourElementsOfConditionsProcess Process
    /**
     * @param rModelPart The model part to check.
     */
    FindNeighbourElementsOfConditionsProcess( ModelPart& rModelPart ):  Process(),
            mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~FindNeighbourElementsOfConditionsProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY

        // Next check that the conditions are oriented accordingly
        // to do so begin by putting all of the conditions in a map
        typedef std::unordered_map<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> > hashmap;
        hashmap FacesMap;

        for (auto itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); ++itCond) {
            itCond->Set(VISITED, false);
            GeometryType& rGeometry = itCond->GetGeometry();

            DenseVector<int> Ids(rGeometry.size());

            for (IndexType i=0; i<Ids.size(); ++i) {
                rGeometry[i].Set(BOUNDARY,true);
                Ids[i] = rGeometry[i].Id();
            }

            // adds to the map
            FacesMap.insert( hashmap::value_type(Ids, std::vector<Condition::Pointer>({*itCond.base()})) );
        }

        // Now loop for all the elements and for each face of the element check if it is in the "FacesMap"
        // if it happens to be there check the orientation
        for (auto itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem) {
            // only active elements are considered
            const auto &rGeometryElement = itElem->GetGeometry();
            const auto rBoundaryGeometries = rGeometryElement.GenerateBoundariesEntities();
            const auto ElementType = rGeometryElement.GetGeometryType();

            for (IndexType iFace = 0; iFace < rBoundaryGeometries.size(); ++iFace) {
                DenseVector<int> FaceIds(rBoundaryGeometries[iFace].size());

                // faces or edges for 2D and 3D elements
                for (IndexType iNode = 0; iNode < FaceIds.size(); ++iNode) {
                    FaceIds[iNode] = rBoundaryGeometries[iFace][iNode].Id();
                }

                hashmap::iterator itFace = FacesMap.find(FaceIds);
                if (itFace != FacesMap.end()) { // It was actually found!!
                    // Mark the condition as visited. This will be useful for a check at the endif
                    std::vector<Condition::Pointer>& ListConditions = itFace->second;
                    for (Condition::Pointer pCondition : ListConditions) {
                        pCondition->Set(VISITED,true);
                    }

                    GlobalPointersVector< Element > VectorOfNeighbours;
                    VectorOfNeighbours.resize(1);
                    VectorOfNeighbours(0) = Element::WeakPointer( *itElem.base() );
                    for (Condition::Pointer pCondition : ListConditions) {
                        pCondition->SetValue(NEIGHBOUR_ELEMENTS, VectorOfNeighbours);
                    }
                }
            }
        }

        //check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
        bool AllVisited = true;
        for (auto& rCond : mrModelPart.Conditions()) {
            if (rCond.IsNot(VISITED)) {
                AllVisited = false;
                break;
            }
        }

        if (!AllVisited) {
            // Now try point loads:
            for (auto itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem) {
                // only active elements are considered
                const auto &rGeometryElement = itElem->GetGeometry();
                const auto rPointGeometries = rGeometryElement.GeneratePoints();
                const auto ElementType = rGeometryElement.GetGeometryType();

                for (IndexType iPoint = 0; iPoint < rPointGeometries.size(); ++iPoint) {
                    DenseVector<int> PointIds(rPointGeometries[iPoint].size());

                    // Points
                    for (IndexType iNode = 0; iNode < PointIds.size(); ++iNode) {
                        PointIds[iNode] = rPointGeometries[iPoint][iNode].Id();
                    }

                    hashmap::iterator itFace = FacesMap.find(PointIds);
                    if (itFace != FacesMap.end()) { // It was actually found!!
                        // Mark the condition as visited. This will be useful for a check at the endif
                        std::vector<Condition::Pointer>& ListConditions = itFace->second;
                        for (Condition::Pointer pCondition : ListConditions) {
                            pCondition->Set(VISITED,true);
                        }

                        GlobalPointersVector< Element > VectorOfNeighbours;
                        VectorOfNeighbours.resize(1);
                        VectorOfNeighbours(0) = Element::WeakPointer( *itElem.base() );
                        for (Condition::Pointer pCondition : ListConditions) {
                            pCondition->SetValue(NEIGHBOUR_ELEMENTS, VectorOfNeighbours);
                        }
                    }
                }
            }
        }

        //check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
        AllVisited = true;
        for (auto& rCond : mrModelPart.Conditions()) {
            if (rCond.IsNot(VISITED)) {
                AllVisited = false;
                KRATOS_WARNING("Condition without any corresponding element, ID ") << rCond.Id() << "\n";
            }
        }

        KRATOS_ERROR_IF_NOT(AllVisited) << "Some conditions found without any corresponding element" << std::endl;

        KRATOS_CATCH("")
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
    std::string Info() const override
    {
        return "FindNeighbourElementsOfConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNeighbourElementsOfConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


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

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindNeighbourElementsOfConditionsProcess& operator=(FindNeighbourElementsOfConditionsProcess const& rOther);

    /// Copy constructor.
    FindNeighbourElementsOfConditionsProcess(FindNeighbourElementsOfConditionsProcess const& rOther);

    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNeighbourElementsOfConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNeighbourElementsOfConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos


#endif // KRATOS_FIND_NEIGHBOUR_ELEMENTS_OF_CONDITIONS_PROCESS
