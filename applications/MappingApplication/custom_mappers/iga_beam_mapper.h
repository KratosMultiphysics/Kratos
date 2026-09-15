//    Kratos Multiphysics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt

#pragma once

// Project includes
#include "mappers/mapper.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos
{

/**
 * @brief Nonlinear mapping between an IGA Bernoulli beam and a surface.
 * @details The origin is the beam and the destination is the surface. The
 * intended kinematics use three translations and one scalar twist per control
 * point. Surface loads will be transferred with the transpose of the kinematic
 * tangent. The origin elements must reference one parent NurbsCurveGeometry3D.
 * Its control points must belong to the origin and store historical DISPLACEMENT
 * and ROTATION (ROTATION_X is the scalar twist, not a global bending rotation).
 * Destination nodes must store the requested historical vector output at Map. Reference attachments
 * are initialized at construction. Map writes total historical DISPLACEMENT;
 * InverseMap writes POINT_LOAD and POINT_MOMENT_X using the kinematic transpose. The parent curve must outlive the mapper. The beam element must
 * provide reference frame rows T, N, V through LOCAL_AXES_MATRIX.
 */
template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) IgaBeamMapper
    : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(IgaBeamMapper);

    using BaseType = Mapper<TSparseSpace, TDenseSpace>;
    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;
    using NurbsCurveType = NurbsCurveGeometry<3, PointerVector<Node>>;

    struct ReferenceAttachment
    {
        Node::Pointer pSurfaceNode;
        double Parameter = 0.0;
        std::vector<Node::Pointer> ControlPoints;
        Vector ShapeFunctions;
        Matrix ShapeFunctionDerivatives;
        array_1d<double, 3> ReferencePosition = ZeroVector(3);
        array_1d<double, 3> CenterlinePosition = ZeroVector(3);
        Matrix ReferenceFrame;
        double TangentialOffset = 0.0;
        double NormalOffset = 0.0;
        double BinormalOffset = 0.0;
    };

    const std::vector<ReferenceAttachment>& GetReferenceAttachments() const
    {
        return mReferenceAttachments;
    }

    /// Constructor for the factory prototype.
    IgaBeamMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination);

    IgaBeamMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters);

    MapperUniquePointerType Clone(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters) const override;

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override;

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override;

    void Map(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions) override;

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override;

    void InverseMap(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions) override;

    /// Beam data-transfer operator interface; validates the secondary variable.
    void Map(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rRotationVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions);

    /// Beam data-transfer operator interface; validates the secondary variable.
    void InverseMap(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rMomentVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions);

    ModelPart& GetInterfaceModelPartOrigin() override;
    ModelPart& GetInterfaceModelPartDestination() override;

    std::string Info() const override;
    void PrintInfo(std::ostream& rOStream) const override;

private:
    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    const NurbsCurveType* mpReferenceCurve = nullptr;
    std::vector<ReferenceAttachment> mReferenceAttachments;
    double mProjectionTolerance = 1e-8;
    int mProjectionMaxIterations = 50;

    array_1d<double, 3> EvaluateAttachment(
        const ReferenceAttachment& rAttachment, Matrix* pTangent) const;

    void InitializeInput();
    void InitializeReferenceAttachments();
};

} // namespace Kratos
