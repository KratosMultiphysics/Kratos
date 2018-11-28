//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "nearest_element_mapper.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"
#include "utilities/geometrical_projection_utilities.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/interface_communicator_mpi.h"
#endif

namespace Kratos {

typedef std::size_t IndexType;
typedef std::size_t SizeType;

typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;
typedef Kratos::shared_ptr<GeometryType> GeometryPointerType;
typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;


bool ProjectTo1D2D(const GeometryPointerType pGeometry,
                   const Point& rPointToProject,
                   CoordinatesArrayType& rLocalCoords,
                   double& rDistance)
{
    CoordinatesArrayType local_coords_init;

    Point projected_point;

    // using the center as trial for the projection
    pGeometry->PointLocalCoordinates(local_coords_init, pGeometry->Center());

    // trying to project to the geometry
    rDistance = std::abs(GeometricalProjectionUtilities::FastProjectDirection(
        *pGeometry,
        rPointToProject,
        projected_point,
        pGeometry->UnitNormal(local_coords_init),
        pGeometry->UnitNormal(local_coords_init)));

    bool is_inside = pGeometry->IsInside(projected_point, rLocalCoords);

    return is_inside;
}

bool ProjectIntoVolume(const GeometryPointerType pGeometry,
                       const Point& rPointToProject,
                       CoordinatesArrayType& rLocalCoords,
                       double& rDistance)
{
    bool is_inside = pGeometry->IsInside(rPointToProject, rLocalCoords);

    if (is_inside) {
        // Calculate Distance
        rDistance = MapperUtilities::ComputeDistance(rPointToProject, pGeometry->Center());
        rDistance /= pGeometry->Volume(); // Normalize Distance by Volume
    }

    return is_inside;
}

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void NearestElementInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();
    const SizeType num_nodes = p_geom->PointsNumber();

    const auto geom_family = p_geom->GetGeometryFamily();

    bool is_inside;
    double proj_dist;
    CoordinatesArrayType local_coords;

    const Point point_to_proj(this->Coordinates());

    // select projection depending on type of geometry
    if ((geom_family == GeometryData::Kratos_Linear        && num_nodes == 2) || // linear line
        (geom_family == GeometryData::Kratos_Triangle      && num_nodes == 3) || // linear triangle
        (geom_family == GeometryData::Kratos_Quadrilateral && num_nodes == 4)) { // linear quad
        is_inside = ProjectTo1D2D(p_geom, point_to_proj, local_coords, proj_dist);
    }
    else if (geom_family == GeometryData::Kratos_Tetrahedra ||
             geom_family == GeometryData::Kratos_Prism ||
             geom_family == GeometryData::Kratos_Hexahedra) { // Volume projection
        is_inside = ProjectIntoVolume(p_geom, point_to_proj, local_coords, proj_dist);
    }
    else {
        is_inside = false;
        KRATOS_WARNING_ONCE("NearestElementMapper") << "Unsupported type of geometry,"
            << "trying to use an approximation (Nearest Neighbor)" << std::endl;
    }

    Vector shape_function_values;

    // if it is closer, then we update the members to make this geometry the closest projection
    if (is_inside && proj_dist < mClosestProjectionDistance) {
        SetLocalSearchWasSuccessful();
        mClosestProjectionDistance = proj_dist;
        mShapeFunctionValues.clear();
        mNodeIds.clear();

        p_geom->ShapeFunctionsValues(shape_function_values, local_coords);
        KRATOS_DEBUG_ERROR_IF_NOT(shape_function_values.size() == num_nodes)
            << "Number of SFs is different from number of nodes!" << std::endl;

        if (mShapeFunctionValues.size() != num_nodes) mShapeFunctionValues.resize(num_nodes);
        if (mNodeIds.size() != num_nodes)             mNodeIds.resize(num_nodes);
        for (IndexType i=0; i<num_nodes; ++i) {
            mShapeFunctionValues[i] = shape_function_values[i];
            KRATOS_DEBUG_ERROR_IF_NOT((*p_geom)[i].Has(INTERFACE_EQUATION_ID))
                << "Node #" << (*p_geom)[i].Id() << " does not have an Interface Id!\n" << (*p_geom)<< "\n"
                <<  (*p_geom)[i] << std::endl;
            mNodeIds[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
        }
    }
}

void NearestElementInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject::Pointer& rpInterfaceObject,
                                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();

    // looping the points of the geometry and finding the nearest neighbor
    for (const auto& r_point : p_geom->Points()) {
        const double dist = MapperUtilities::ComputeDistance(this->Coordinates(), r_point.Coordinates());

        // in case of an approximation this is the actual distance,
        // not the projected one bcs no valid projection could be found!
        if (dist < mClosestProjectionDistance) {
            mClosestProjectionDistance = dist;
            if (mNodeIds.size() != 1) mNodeIds.resize(1);
            if (mShapeFunctionValues.size() != 1) mShapeFunctionValues.resize(1);
            KRATOS_DEBUG_ERROR_IF_NOT(r_point.Has(INTERFACE_EQUATION_ID));
            mNodeIds[0] = r_point.GetValue(INTERFACE_EQUATION_ID);
            mShapeFunctionValues[0] = 1.0; // Approximation is nearest node
        }
    }

    SetIsApproximation();
}

void NearestElementLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() > 0) {
        double distance;
        double min_distance = std::numeric_limits<double>::max();
        int found_idx = -1;
        for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
            // the approximations will be processed in the next step if necessary
            if (!mInterfaceInfos[i]->GetIsApproximation()) {
                mInterfaceInfos[i]->GetValue(distance);
                if (distance < min_distance) {
                    min_distance = distance;
                    found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
                }
            }
        }

        if (found_idx == -1) { // no valid projection exists => using an approximation
            for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
                // now the approximations are being checked
                if (mInterfaceInfos[i]->GetIsApproximation()) {
                    mInterfaceInfos[i]->GetValue(distance);
                    if (distance < min_distance) {
                        min_distance = distance;
                        found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                        rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
                    }
                }
            }
        }

        KRATOS_ERROR_IF(found_idx == -1) << "Not even an approximation is found, this should not happen!"
            << std::endl; // TODO should thi sbe an error?

        std::vector<double> sf_values;

        mInterfaceInfos[found_idx]->GetValue(sf_values);

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
            rLocalMappingMatrix.resize(1, sf_values.size(), false);
        }
        for (IndexType i=0; i<sf_values.size(); ++i) {
            rLocalMappingMatrix(0,i) = sf_values[i];
        }

        mInterfaceInfos[found_idx]->GetValue(rOriginIds);

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string NearestElementLocalSystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestElementLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    buffer << " in rank " << CommRank;
    return buffer.str();
}/* Performs operations that are needed for Initialization and when the interface is updated (=> Remeshed)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    mpMapperLocalSystems->clear();

    const MapperLocalSystemPointer p_ref_local_system = Kratos::make_unique<NearestElementLocalSystem>();;

    KRATOS_ERROR_IF_NOT(mpInterfacePreprocessor) << "mpInterfacePreprocessor is a nullptr!" << std::endl;
    mpInterfacePreprocessor->CreateMapperLocalSystems(p_ref_local_system);

    BuildMappingMatrix(MappingOptions);
}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

    KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpIntefaceCommunicator is a nullptr!" << std::endl;

    const MapperInterfaceInfoUniquePointerType p_ref_interface_info = Kratos::make_unique<NearestElementInterfaceInfo>();
    const auto interface_object_construction_type_origin = InterfaceObject::Geometry_Center;

    mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                  MappingOptions,
                                                  p_ref_interface_info,
                                                  interface_object_construction_type_origin);

    KRATOS_ERROR_IF_NOT(mpMappingMatrixBuilder) << "mpMappingMatrixBuilder is a nullptr!" << std::endl;

    mpMappingMatrixBuilder->BuildMappingMatrix(mpInterfaceVectorContainerOrigin,
                                                  mpInterfaceVectorContainerDestination,
                                                  *mpMapperLocalSystems);

    // if (mEchoLevel > 0) PrintPairingInfo();
}

template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::InitializeMappingMatrixBuilder()
{
    mpMappingMatrixBuilder = Kratos::make_unique<MappingMatrixBuilder<TSparseSpace, TDenseSpace>>();
}

template<>
void NearestElementMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    Parameters search_settings(R"({})"); // TODO fill this
    search_settings.ValidateAndAssignDefaults(mMapperSettings);
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartOrigin,
                                                                      mpMapperLocalSystems,
                                                                      search_settings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void NearestElementMapper<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    Parameters search_settings(R"({})"); // TODO fill this
    search_settings.ValidateAndAssignDefaults(mMapperSettings);
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartOrigin,
                                                                         mpMapperLocalSystems,
                                                                         search_settings);
}
#endif

template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerDestination->GetVector(),
        mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
    const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
    const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

    const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
    const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
    const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

    // X-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_x_origin, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_x_destination, MappingOptions);

    // Y-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_y_origin, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_y_destination, MappingOptions);

    // Z-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_z_origin, MappingOptions);

    TSparseSpace::Mult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_z_destination, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
    const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
    const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

    const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
    const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
    const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

    // X-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_x_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_x_origin, MappingOptions);

    // Y-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_y_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_y_origin, MappingOptions);

    // Z-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_z_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        mpMappingMatrixBuilder->GetMappingMatrix(),
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_z_origin, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void NearestElementMapper<TSparseSpace, TDenseSpace>::ValidateInput(Parameters MapperSettings)
{
    MapperUtilities::CheckInterfaceModelParts(0);
    ValidateParameters(MapperSettings);

    // mEchoLevel = MapperSettings["echo_level"].GetInt();

    if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
        const double search_radius = MapperUtilities::ComputeSearchRadius(mrModelPartOrigin,
                                        mrModelPartDestination,
                                        0);
        mMapperSettings["search_radius"].SetDouble(search_radius);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class NearestElementMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class NearestElementMapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
