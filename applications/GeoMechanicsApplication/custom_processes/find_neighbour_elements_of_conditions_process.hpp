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


#pragma once

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
///@name Type Definitions
///@{
    using hashmap = std::unordered_multimap<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>>>;

///@}
///@name Kratos Classes
///@{

/**
 * @class FindNeighbourElementsOfConditionsProcess
 * @ingroup Kratos.GeoMechanicsApplication
 * @brief Finds list of elements attached to conditions.
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) FindNeighbourElementsOfConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(FindNeighbourElementsOfConditionsProcess);

    /// The definition of the index type
    using IndexType = std::size_t;

    /// Definition of the node type
    using NodeType = Node;

    // Definition of the geometry
    using GeometryType = Geometry<NodeType>;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for FindNeighbourElementsOfConditionsProcess Process
    /**
     * @param rModelPart The model part to check.
     */
    explicit FindNeighbourElementsOfConditionsProcess( ModelPart& rModelPart ): Process(),
            mrModelPart(rModelPart)
    {
    }

    FindNeighbourElementsOfConditionsProcess& operator=(const FindNeighbourElementsOfConditionsProcess&) = delete;
    FindNeighbourElementsOfConditionsProcess(const FindNeighbourElementsOfConditionsProcess&) = delete;
    ~FindNeighbourElementsOfConditionsProcess() override = default;

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


    /// Finds neighbour elements of conditions
    void Execute() override;

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

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
private:
    ///@name Static Member Variables
    ///@{
    hashmap::iterator FindFaceReorderingTetrahedra3D10(DenseVector<int> FaceIds, hashmap &FacesMap) const;
    hashmap::iterator FindFaceReorderingTetrahedra3D4(DenseVector<int> FaceIds, hashmap &FacesMap) const;
    hashmap::iterator FindFaceReorderingHexahedra3D8(DenseVector<int> FaceIds, hashmap &FacesMap) const;
    hashmap::iterator FindFaceReorderingHexahedra3D20(DenseVector<int> FaceIds, hashmap &FacesMap) const;


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operations
    ///@{

    bool CheckIfAllConditionsAreVisited() const;

    void CheckIf1DElementIsNeighbour(hashmap& rFacesMap);

    static void CheckForMultipleConditionsOnElement(hashmap& rFacesMap, hashmap::iterator& rItFace,
        PointerVector<Element>::iterator pItElem);

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

}