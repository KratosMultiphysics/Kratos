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

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/geometry.h"
#include "processes/process.h"
#include "utilities/entities_utilities.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/* @class IgaContactProcessGapSbm
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) IgaContactProcessGapSbm
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IgaContactProcessGapSbm
    KRATOS_CLASS_POINTER_DEFINITION(IgaContactProcessGapSbm);

    typedef Node                                             NodeType;
    typedef Geometry<NodeType>                               GeometryType;
    typedef GeometryType::Pointer                            GeometryPointerType;
    typedef typename GeometryType::GeometriesArrayType       GeometriesArrayType;
    typedef typename GeometryType::CoordinatesArrayType      CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef typename Properties::Pointer PropertiesPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    IgaContactProcessGapSbm(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~IgaContactProcessGapSbm() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitializeSolutionStep() override {
        Execute();
    };

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "analysis_model_part_name" : "ModelPart",
            "contact_sub_model_part_name" : "contact",
            "echo_level" : 0,
            "numbered_considered_neighbours" : 1,
            "contact_parameters" : {
                "slave_model_part" : {
                    "sub_model_part_name" : "",
                    "layer_name" : "",
                    "property_id" : 1
                },
                "master_model_part" : {
                    "sub_model_part_name" : "",
                    "layer_name" : "",
                    "property_id" : 1
                }
            }
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaContactProcessGapSbm";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaContactProcessGapSbm";
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
    SizeType mEchoLevel = 0;

    ModelPart* mrSlaveModelPart = nullptr; 
    ModelPart* mrMasterModelPart = nullptr; 
    ModelPart* mpContactModelPart = nullptr; 

    Properties::Pointer mpPropMaster;
    Properties::Pointer mpPropSlave;

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}

}; // Class IgaContactProcessGapSbm

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  IgaContactProcessGapSbm& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IgaContactProcessGapSbm& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
