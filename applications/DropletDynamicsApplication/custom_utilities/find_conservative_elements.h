//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   
//
//

#pragma once

#include <unordered_map>

#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"
#include "modified_shape_functions/modified_shape_functions.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/rbf_shape_functions_utility.h"
namespace Kratos
{

struct FindConservativeElementsSettings
{
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};


template<bool THistorical = true>
class KRATOS_API(KRATOS_DROPLETDYNAMICSAPPLICATION) FindConservativeElementsProcess
    : public Process
{
public:

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using NodeType = Node;

    using NodeIterator = ModelPart::NodeIterator;
    
    using GeometryType = Geometry<NodeType>;

    using GeometryPointerType = Geometry<NodeType>::Pointer;

    using NodesPointerSetType = ModelPart::NodesContainerType;

    using ElementFacesMapType = std::unordered_map<
    std::vector<IndexType>,
    std::pair<bool, GeometryPointerType>,
    KeyHasherRange<std::vector<IndexType>>,
    KeyComparorRange<std::vector<IndexType>>>;

    void FINDCUTELEMENTNODES();

    /// Pointer definition of FindConservativeElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindConservativeElementsProcess);


    /// Default constructor.
    explicit FindConservativeElementsProcess(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~FindConservativeElementsProcess() override = default;


    void operator()()
    {
        Execute();
    }


    void Execute() override;


    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindConservativeElementsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindConservativeElementsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


private:

    ModelPart& mrModelPart;  /// The model part were to compute Nodal_Nighbers

    double& GetDISValue(NodeType& rNode);


    /// Assignment operator.
    FindConservativeElementsProcess& operator=(FindConservativeElementsProcess const& rOther);


};

/// input stream function
template<bool THistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  FindConservativeElementsProcess<THistorical>& rThis);
              

/// output stream function
template<bool THistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindConservativeElementsProcess<THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
