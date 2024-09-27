//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

#if !defined(KRATOS_IGA_CONTACT_PROCESS_H_INCLUDED )
#define  KRATOS_IGA_CONTACT_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/pointer_vector.h"
#include "containers/model.h"
#include "geometries/geometry.h"

#include "processes/process.h"

#include "geometries/nurbs_coupling_geometry_2d.h"

#include "custom_conditions/support_contact_2D_condition.h"

#include "utilities/entities_utilities.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class IgaContactProcess
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) IgaContactProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IgaContactProcess
    KRATOS_CLASS_POINTER_DEFINITION(IgaContactProcess);

    typedef Node                                             NodeType;
    typedef Geometry<NodeType>                                  GeometryType;
    typedef GeometryType::Pointer                               GeometryPointerType;
    typedef typename GeometryType::GeometriesArrayType          GeometriesArrayType;
    typedef typename GeometryType::CoordinatesArrayType         CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    using PointType = Node;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    IgaContactProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~IgaContactProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    // void ExecuteInitialize() override {
    //     Execute();
    // };

    void ExecuteInitializeSolutionStep() override {
        Execute();
    };

    // void ExecuteBeforeSolutionLoop() override {
    //     Execute();
    // };

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "main_model_part_name" : "ModelPart",
            "nurbs_volume_name" : "NurbsVolume",
            "embedded_model_part_name" : "IgaModelPart"
        })" );

        return default_parameters;
    }


        /// Creates conditions from geometries
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pPropMaster,
        PropertiesPointerType pPropSlave) const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaContactProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaContactProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;
    SizeType mEchoLevel;


    bool GetProjection(CoordinatesArrayType& slavePoint, GeometryType &slave_geometry, GeometryType &master_geometry, double& localProjection);

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@name CAD functionalities
    ///@{

    /// Gets list of geometries from rModelPart
    void GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    ///@}

    ///@}
    ///@name Utility
    ///@{

    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;

    ///@}
    ///@name Input and output
    ///@{


   

    ///@}

}; // Class IgaContactProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  IgaContactProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IgaContactProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_CONTACT_PROCESS_H_INCLUDED  defined
