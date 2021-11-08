//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Manuel Messmer

#if !defined(KRATOS_MAP_NURBS_VOLUME_RESULTS_TO_EMBEDDED_GEOMETRY_PROCESS_H_INCLUDED )
#define  KRATOS_MAP_NURBS_VOLUME_RESULTS_TO_EMBEDDED_GEOMETRY_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/pointer_vector.h"
#include "containers/model.h"
#include "geometries/geometry.h"
#include "geometries/nurbs_volume_geometry.h"

#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class MapNurbsVolumeResultsToEmbeddedGeometryProcess
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) MapNurbsVolumeResultsToEmbeddedGeometryProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapNurbsVolumeResultsToEmbeddedGeometryProcess
    KRATOS_CLASS_POINTER_DEFINITION(MapNurbsVolumeResultsToEmbeddedGeometryProcess);

    typedef Node<3>                                             NodeType;
    typedef Geometry<NodeType>                                  GeometryType;
    typedef GeometryType::Pointer                               GeometryPointerType;
    typedef typename GeometryType::GeometriesArrayType          GeometriesArrayType;
    typedef typename GeometryType::CoordinatesArrayType         CoordinatesArrayType;
    typedef NurbsVolumeGeometry<PointerVector<NodeType>>        NurbsVolumeGeometryType;
    typedef NurbsVolumeGeometryType::Pointer                    NurbsVolumeGeometryPointerType;
    typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    MapNurbsVolumeResultsToEmbeddedGeometryProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~MapNurbsVolumeResultsToEmbeddedGeometryProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void MapNodalValues(const Variable<array_1d<double,3>>& rVariable);

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "main_model_part_name"                    : "main_model_part",
            "nurbs_volume_name"                       : "nurbs_volume",
            "embedded_model_part_name"                : "embedded_model_part",
            "nodal_results": []
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{


    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MapNurbsVolumeResultsToEmbeddedGeometryProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapNurbsVolumeResultsToEmbeddedGeometryProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;
    Parameters mThisParameters;

    ///@}

}; // Class MapNurbsVolumeResultsToEmbeddedGeometryProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MapNurbsVolumeResultsToEmbeddedGeometryProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MapNurbsVolumeResultsToEmbeddedGeometryProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_MAP_NURBS_VOLUME_RESULTS_TO_EMBEDDED_GEOMETRY_PROCESS_H_INCLUDED  defined
