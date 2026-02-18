//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//

# pragma once

// System includes

// External includes

// Project includes
#include "nurbs_geometry_modeler.h"

namespace Kratos {

class KRATOS_API(IGA_APPLICATION) NurbsGeometryModelerGapSbm
    : public NurbsGeometryModeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NurbsGeometryModelerGapSbm );

    using IndexType = NurbsGeometryModeler::IndexType;
    using NodeType = NurbsGeometryModeler::NodeType;

    using GeometryType = NurbsGeometryModeler::GeometryType;
    using GeometryPointerType = NurbsGeometryModeler::GeometryPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NurbsGeometryModelerGapSbm()
        : NurbsGeometryModeler() {}

    /// Constructor.
    NurbsGeometryModelerGapSbm(
        Model & rModel,
        const Parameters ModelerParameters = Parameters())
        : NurbsGeometryModeler(rModel, ModelerParameters)
    {
        mParameters.ValidateDefaults(this->GetValidParameters());
        mParameters.AddMissingParameters(this->GetDefaultParameters());

        KRATOS_ERROR_IF_NOT(mParameters.Has("gap_element_name"))
            << "NurbsGeometryModelerGapSbm: Missing \"gap_element_name\" section." << std::endl;
        
        KRATOS_ERROR_IF_NOT(mParameters["gap_interface_condition_name"].IsString())
            << "NurbsGeometryModelerGapSbm: Missing \"gap_interface_condition_name\" section." << std::endl;
        
    }

    /// Destructor.
    ~NurbsGeometryModelerGapSbm() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, 
        const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<NurbsGeometryModelerGapSbm>(rModel, ModelParameters);
    }

    // Get the default parameters
    const Parameters GetDefaultParameters() const override;
    const Parameters GetValidParameters() const;

    ///@}
    ///@name Stages
    ///@{

    ///@}

protected:

    /**
     * @brief Creates a regular grid composed out of bivariant B-splines.
     * @param PointA Lower point of bounding box.
     * @param PointB Upper point of bounding box.
     * @param Order  Polynomial degree in each direction u,v.
     * @param NumKnotSpans Number of equidistant elements/knot spans in each direction u,v.
     * @note The CP'S are defined as nodes and added to the rModelPart.
     **/
    void CreateAndAddRegularGrid2D(
        ModelPart& rModelPart, 
        const Point& rPointAXyz, 
        const Point& rPointBXyz, 
        const Point& rPointAUvw, 
        const Point& rPointBUvw,
        const std::size_t OrderU, 
        const std::size_t OrderV, 
        const std::size_t NumKnotSpansU, 
        const std::size_t NumKnotSpansV, 
        const bool AddSurfaceToModelPart) override;

    /**
     * @brief Creates a regular grid composed out of bivariant B-splines.
     * @param PointA Lower point of bounding 3D box.
     * @param PointB Upper point of bounding 3D box.
     * @param Order  Polynomial degree in each direction u,v,w.
     * @param NumKnotSpans Number of equidistant elements/knot spans in each direction u,v,w.
     * @note The CP'S are defined as nodes and added to the rModelPart.
     **/
    void CreateAndAddRegularGrid3D(
        ModelPart& rModelPart, 
        const Point& rPointAXyz, 
        const Point& rPointBXyz, 
        const Point& rPointAUvw, 
        const Point& rPointBUvw,
        const std::size_t OrderU, 
        const std::size_t OrderV, 
        const std::size_t OrderW, 
        const std::size_t NumKnotSpansU, 
        const std::size_t NumKnotSpansV, 
        const std::size_t NumKnotSpansW, 
        const bool AddVolumeToModelPart) override;

private:

    ///@name Private Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

};

} // End namesapce Kratos
