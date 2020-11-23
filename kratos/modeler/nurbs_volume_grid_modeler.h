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

#if !defined(KRATOS_NURBS_VOLUME_GRID_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_VOLUME_GRID_MODELER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_utilities.h"
#include "includes/node.h"
#include "includes/model_part.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) NurbsVolumeGridModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NurbsVolumeGridModeler );

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;

    typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
    typedef typename NurbsVolumeGeometryType::Pointer NurbsVolumeGeometryPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NurbsVolumeGridModeler()
        : Modeler() {}

    /// Constructor.
    NurbsVolumeGridModeler(
        Model & rModel,
        const Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~NurbsVolumeGridModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<NurbsVolumeGridModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    void SetupModelPart() override;

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    Model* mpModel;
    NurbsVolumeGeometryPointerType mpGeometry;

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * @brief Creates a cartesian grid composed out of trivariant B-spline cubes.
     * @param PointA Lower point of bounding box.
     * @param PointB Upper point of bounding box.
     * @param Order  Polynomial degree in each direction u,v,w.
     * @param NumElements Number of equidistant elements/knot spans in each direction u,v,w.
     * @return Pointer to NurbsVolumeGeometry
     * @note The CP'S are defined as nodes and added to the rModelPart.
     * @todo How to deal with node Id's..
     **/
    void CreateGrid( const Point A, const Point B, SizeType OrderU, SizeType OrderV, SizeType OrderW,
        SizeType NumElementsU, SizeType NumElementsV, SizeType NumElementsW );

};

} // End namesapce Kratos
#endif // KRATOS_NURBS_VOLUME_GRID_MODELER_H_INCLUDED