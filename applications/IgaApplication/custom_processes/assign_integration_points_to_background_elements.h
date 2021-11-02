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

#if !defined(KRATOS_ASSIGN_INTEGRATION_POINTS_TO_BACKGROUND_ELEMENTS_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_INTEGRATION_POINTS_TO_BACKGROUND_ELEMENTS_PROCESS_H_INCLUDED

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

/* @class AssignIntegrationPointsToBackgroundElementsProcess
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) AssignIntegrationPointsToBackgroundElementsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignIntegrationPointsToBackgroundElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignIntegrationPointsToBackgroundElementsProcess);

    // typedef Node<3>                                             NodeType;
    // typedef Geometry<NodeType>                                  GeometryType;
    // typedef GeometryType::Pointer                               GeometryPointerType;
    // typedef typename GeometryType::GeometriesArrayType          GeometriesArrayType;
    // typedef typename GeometryType::CoordinatesArrayType         CoordinatesArrayType;
    // typedef NurbsVolumeGeometry<PointerVector<NodeType>>        NurbsVolumeGeometryType;
    // typedef NurbsVolumeGeometryType::Pointer                    NurbsVolumeGeometryPointerType;
    // typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;
    // typedef std::size_t IndexType;
    // typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    AssignIntegrationPointsToBackgroundElementsProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~AssignIntegrationPointsToBackgroundElementsProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void AssignIntegrationPoints();

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "main_model_part_name" : "ModelPart",
            "nurbs_volume_name" : "NurbsVolume",
            "embedded_model_part_name" : "IgaModelPart",
            "integration_point_results" : []
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{


    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AssignIntegrationPointsToBackgroundElementsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignIntegrationPointsToBackgroundElementsProcess";
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

}; // Class AssignIntegrationPointsToBackgroundElementsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignIntegrationPointsToBackgroundElementsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignIntegrationPointsToBackgroundElementsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_ASSIGN_INTEGRATION_POINTS_TO_BACKGROUND_ELEMENTS_PROCESS_H_INCLUDED  defined
