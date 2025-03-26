//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_COUPLING_GEOMETRY_COUPLING_H_INCLUDED )
#define  KRATOS_COUPLING_GEOMETRY_COUPLING_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "modeler/modeler.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.h"

namespace Kratos {

class KRATOS_API(IGA_APPLICATION) CouplingGeometryModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( CouplingGeometryModeler );

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
    typedef typename NurbsVolumeGeometryType::Pointer NurbsVolumeGeometryPointerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CouplingGeometryModeler()
        : Modeler() {}

    /// Constructor.
    CouplingGeometryModeler(
        Model & rModel,
        const Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~CouplingGeometryModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<CouplingGeometryModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    void PrepareGeometryModel() override;

    void SetupModelPart() override;

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    Model* mpModel;

    ///@}
    ///@name Private Operations
    ///@{


};

} // End namesapce Kratos
#endif // KRATOS_NURBS_GEOMETRY_MODELER_H_INCLUDED