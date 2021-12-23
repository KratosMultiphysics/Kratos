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
// #include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"


namespace Kratos
{
///@name Type Definitions
///@{
    typedef std::unordered_multimap<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> > hashmap;

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
    FindNeighbourElementsOfConditionsProcess( ModelPart& rModelPart ): Process(),
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
