#if !defined(KRATOS_GENERIC_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_GENERIC_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"



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
 * @class FindNodalHProcess
 * @ingroup KratosCore
 * @brief Finds Neighbour elements of elements (by face/edge)
 * @details Finds the NEIGHBOUR_ELEMENTS of elements. For 3D , they are the elements connected by faces. For 2D, elements connected by edges
 * @details The neighbour vector has size = nfaces(volumetric geom) or nedges(planar geom). The order of the faces is the one given by the geom.GenerateFaces(). If no neigh elem is in a face, then this slot in the vector contains a nullptr
 * @author Pablo Becker
 */
class  KRATOS_API(KRATOS_CORE) GenericFindElementalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef  ModelPart::NodesContainerType NodesContainerType;
    typedef  ModelPart::ElementsContainerType ElementsContainerType;
    typedef Geometry<Node > GeometryType;


    /// Pointer definition of GenericFindElementalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(GenericFindElementalNeighboursProcess);

    typedef GlobalPointersVector<Element> ElementPointerVector;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenericFindElementalNeighboursProcess(ModelPart& model_part)
        : mrModelPart(model_part)
    {
    }

    /// Destructor.
    ~GenericFindElementalNeighboursProcess() override
    {
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;


    std::vector<bool> HasNeighboursInFaces(const Element&) ;

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
        return "FindElementalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindElementalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Private static Member Variables
    ///@{


    ///@}
    ///@name Private member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    GlobalPointer<Element> CheckForNeighbourElems (const Geometry<Node >& rBoundaryGeom,
                                                   Element & rElement,
                                                   const int Rank);


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private LifeCycle
    ///@{

    /// Assignment operator.
    GenericFindElementalNeighboursProcess& operator=(GenericFindElementalNeighboursProcess const& rOther);

    ///@}

}; // Class GenericFindElementalNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GenericFindElementalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GenericFindElementalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GENERIC_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined
