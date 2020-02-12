//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "move_shallow_particles_process.h"
#include "includes/checks.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

template<size_t TDim>
MoveShallowParticlesProcess<TDim>::MoveShallowParticlesProcess(
    ModelPart& rModelPart,
    ModelPart& rParticles,
    Variable<array_1d<double,3>>& rVectorVariable,
    Variable<double>& rScalarVariable,
    Parameters Settings)
 : mrModelPart(rModelPart)
 , mrParticles(rParticles)
 , mrMomentumVariable(rVectorVariable)
 , mrMassVariable(rScalarVariable)
 , mConvectionOperator(rModelPart)
{
    Check();

    Parameters default_settings = Parameters(R"(
    {
        "quadrature_order"    : 1,
        "projection_type"     : "incompressible"
    })");
    Settings.ValidateAndAssignDefaults(default_settings);

    mrParticles.AddNodalSolutionStepVariable(mrMassVariable);
    mrParticles.AddNodalSolutionStepVariable(mrMomentumVariable);
    mrParticles.AddNodalSolutionStepVariable(PARTICLE_WEIGHT);
    mrParticles.AddNodalSolutionStepVariable(PARTICLE_AREA);
    mLastParticleId = 0;
    mDryTreshold = mrModelPart.GetProcessInfo()[DRY_HEIGHT];

    const int order = Settings["quadrature_order"].GetInt();
    if (order == 1) {
        mQuadratureOrder = GeometryData::GI_GAUSS_1;
    } else if (order == 2) {
        mQuadratureOrder = GeometryData::GI_GAUSS_2;
    } else if (order == 3) {
        mQuadratureOrder = GeometryData::GI_GAUSS_3;
    } else if (order == 4) {
        mQuadratureOrder = GeometryData::GI_GAUSS_4;
    } else if (order == 5) {
        mQuadratureOrder = GeometryData::GI_GAUSS_5;
    } else {
        KRATOS_ERROR << "Unknown quadrature order" << std::endl;
    }

    const std::string& projection_type = Settings["projection_type"].GetString();
    if (projection_type == "incompressible") {
        mProjectionType = Incompressible;
    } else if (projection_type == "compressible") {
        mProjectionType = Compressible;
    } else {
        KRATOS_ERROR << "Unknown projection type" << std::endl;
    }

    InitializeVariables();
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ExecuteInitialize()
{
    SeedInitialParticles();
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ExecuteBeforeSolutionLoop()
{
    ComputeMeanSize();
    mConvectionOperator.SetMaxResults(10000);
    mConvectionOperator.UpdateSearchDatabase();
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ExecuteInitializeSolutionStep()
{
    mrParticles.GetProcessInfo()[STEP] = mrModelPart.GetProcessInfo()[STEP];
    mrParticles.GetProcessInfo()[TIME] = mrModelPart.GetProcessInfo()[TIME];
    KRATOS_WATCH("A")
    ComputeMeanVelocity();
    KRATOS_WATCH("B")
    MoveParticles();
    KRATOS_WATCH("C")
    TransferLagrangianToEulerian();
    KRATOS_WATCH("D")
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ExecuteFinalizeSolutionStep()
{
    KRATOS_WATCH("E")
    TransferEulerianToLagrangian();
    KRATOS_WATCH("F")
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::InitializeVariables()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        const auto it_elem = mrModelPart.ElementsBegin() + i;
        it_elem->SetValue(MEAN_SIZE, 0.0);
        it_elem->SetValue(MEAN_VEL_OVER_ELEM_SIZE, 0.0);
        it_elem->SetValue(NUMBER_OF_PARTICLES, 0);
    }

    if (mProjectionType == Compressible)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            it_node->SetValue(SUM_AREAS, 0.0);
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ComputeMeanSize()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        const auto it_elem = mrModelPart.ElementsBegin() + i;
        const double mean_size = std::sqrt(it_elem->GetGeometry().Area());
        it_elem->GetValue(MEAN_SIZE) = mean_size;
    }

    if (mProjectionType == Compressible)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            const auto it_elem = mrModelPart.ElementsBegin() + i;
            const double area = it_elem->GetGeometry().Area();
            for (auto& node : it_elem->GetGeometry())
            {
                #pragma omp critical
                node.GetValue(SUM_AREAS) += area;
            }
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ComputeMeanVelocity()
{
    const size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    const double nodal_weight = 1.0 / static_cast<double>(num_nodes);

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        const auto it_elem = mrModelPart.ElementsBegin() + i;
        const GeometryType& geom = it_elem->GetGeometry();

        array_1d<double, 3 > mean_velocity = ZeroVector(3);

        for (size_t j = 0; j < num_nodes ; ++j)
            mean_velocity += geom[j].FastGetSolutionStepValue(VELOCITY);
        mean_velocity *= nodal_weight;

        const double mean_size = it_elem->GetValue(MEAN_SIZE);
        const double mean_velocity_norm = norm_2(mean_velocity);
        it_elem->SetValue(MEAN_VEL_OVER_ELEM_SIZE, mean_velocity_norm / mean_size);
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::SeedInitialParticles()
{
    const size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();

    ComputeNumberOfParticlesPerElement();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        const auto it_elem = mrModelPart.ElementsBegin() + i;

        GeometryType& geom = it_elem->GetGeometry();

        int& num_of_particles = it_elem->GetValue(NUMBER_OF_PARTICLES);
        if (num_of_particles < static_cast<int>(num_nodes))
        {
            const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mQuadratureOrder);
            size_t num_gauss = geom.IntegrationPointsNumber(mQuadratureOrder);
            for (size_t g = 0; g < num_gauss; ++g)
            {
                array_1d<double,3> position;
                geom.GlobalCoordinates(position, integration_points[g]);
                const double x = position[0];
                const double y = position[1];
                const double z = position[2];
                NodeType::Pointer p_part;
                #pragma omp critical
                p_part = mrParticles.CreateNewNode(++mLastParticleId, x, y, z);
                InitializeParticle(*p_part, *it_elem);
                ++num_of_particles;

                if (p_part->FastGetSolutionStepValue(HEIGHT) < 0.0)
                {
                    p_part->Set(TO_ERASE);
                }
            }
        }
    }
    mrParticles.RemoveNodesFromAllLevels();
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ComputeNumberOfParticlesPerElement()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        (mrModelPart.ElementsBegin() + i)->GetValue(NUMBER_OF_PARTICLES) = 0;
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;
        auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        if (p_elem != nullptr) {
            int& num_of_particles = p_elem->GetValue(NUMBER_OF_PARTICLES);
            #pragma omp critical
            ++num_of_particles;
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::MoveParticles()
{
    double time_step = mrModelPart.GetProcessInfo()[DELTA_TIME];
    size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    mConvectionOperator.InitializeSearch();

    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;
        auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        Vector N(num_nodes);
        bool is_found = mConvectionOperator.Convect(time_step, *it_part, p_elem, N);
        if (is_found) {
            it_part->GetValue(CURRENT_ELEMENT) = p_elem;
        } else {
            it_part->GetValue(CURRENT_ELEMENT) = nullptr;
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::TransferLagrangianToEulerian()
{
    if (mProjectionType == Incompressible)
    {
        StandardProjection();
    }
    else if (mProjectionType == Compressible)
    {
        ProjectionWithMassConservation();
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::TransferEulerianToLagrangian()
{
    if (mProjectionType == Incompressible)
    {
        StandardParticlesUpdate();
    }
    else if (mProjectionType == Compressible)
    {
        ParticlesUpdateWithMassConservation();
    }
}

    
template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::StandardProjection()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(MASS_WEIGHT) = 0.0;
        it_node->FastGetSolutionStepValue(PROJECTED_SCALAR1) = 0.0;
        it_node->FastGetSolutionStepValue(PROJECTED_VECTOR1) = ZeroVector(3);
    }

    size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;

        // Particle contributions
        double scalar_contribution = 0;
        array_1d<double,3> vector_contribution = ZeroVector(3);
        scalar_contribution =  it_part->FastGetSolutionStepValue(mrMassVariable);
        vector_contribution =  it_part->FastGetSolutionStepValue(mrMomentumVariable);

        // Transfer the values
        auto p_current_elem = it_part->GetValue(CURRENT_ELEMENT);
        if (p_current_elem != nullptr)
        {
            auto& geom = p_current_elem->GetGeometry();

            Vector N;
            GetShapeFunctionsValues(N, geom, *it_part);

            for (size_t j = 0; j < num_nodes; ++j)
            {
                double weight = N[j]*N[j];
                geom[j].SetLock();
                geom[j].FastGetSolutionStepValue(MASS_WEIGHT) += weight;
                geom[j].FastGetSolutionStepValue(PROJECTED_SCALAR1) += weight * scalar_contribution;
                geom[j].FastGetSolutionStepValue(PROJECTED_VECTOR1) += weight * vector_contribution;
                geom[j].UnSetLock();
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        double weight = it_node->FastGetSolutionStepValue(MASS_WEIGHT);
        if (weight > 0.00001)
        {
            it_node->FastGetSolutionStepValue(PROJECTED_SCALAR1) /= weight;
            it_node->FastGetSolutionStepValue(PROJECTED_VECTOR1) /= weight;
        }
    }
}

    
template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::StandardParticlesUpdate()
{
    size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;
        auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        if (p_elem != nullptr)
        {
            const auto& geom = p_elem->GetGeometry();
            Vector N;
            GetShapeFunctionsValues(N, geom, *it_part);

            // Variables to update
            double& scalar_value = it_part->FastGetSolutionStepValue(mrMassVariable);
            double delta_scalar = 0;
            array_1d<double,3>& vector_value = it_part->FastGetSolutionStepValue(mrMomentumVariable);
            array_1d<double,3> delta_vector = ZeroVector(3);

            for (size_t j = 0; j < num_nodes; ++j)
            {
                delta_scalar += N[j] * (geom[j].FastGetSolutionStepValue(mrMassVariable) - geom[j].FastGetSolutionStepValue(mrMassVariable,1));
                delta_vector += N[j] * (geom[j].FastGetSolutionStepValue(mrMomentumVariable) - geom[j].FastGetSolutionStepValue(mrMomentumVariable,1));
            }
            scalar_value += delta_scalar;
            vector_value += delta_vector;
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ProjectionWithMassConservation()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(PROJECTED_SCALAR1) = 0.0;
        it_node->FastGetSolutionStepValue(PROJECTED_VECTOR1) = ZeroVector(3);
    }

    const size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;
        auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        if (p_elem != nullptr)
        {
            auto& geom = p_elem->GetGeometry();
            Vector N;
            GetShapeFunctionsValues(N, geom, *it_part);

            double sum_part_weights = 0.0;
            for (size_t j = 0; j < num_nodes; ++j)
            {
                sum_part_weights += N[j]*N[j];
            }

            for (size_t j = 0; j < num_nodes; ++j)
            {
                const double part_weight = N[j]*N[j] / sum_part_weights;
                const double part_area = it_part->FastGetSolutionStepValue(PARTICLE_AREA);
                const double sum_areas = geom[j].GetValue(SUM_AREAS);
                const auto scalar_value = it_part->FastGetSolutionStepValue(mrMassVariable);
                const auto vector_value = it_part->FastGetSolutionStepValue(mrMomentumVariable);
                geom[j].FastGetSolutionStepValue(PROJECTED_SCALAR1) += num_nodes * part_weight * part_area / sum_areas * scalar_value;
                geom[j].FastGetSolutionStepValue(PROJECTED_VECTOR1) += num_nodes * part_weight * part_area / sum_areas * vector_value;
            }
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ParticlesUpdateWithMassConservation()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        (mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(SUM_PARTICLES_WEIGHTS) = 0.0;
    }

    const size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        const auto it_part = mrParticles.NodesBegin() + i;
        const auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        if (p_elem != nullptr)
        {
            auto& geom = p_elem->GetGeometry();
            Vector N;
            GetShapeFunctionsValues(N, geom, *it_part);

            for (size_t j = 0; j < num_nodes; ++j)
            {
                #pragma omp critical
                geom[j].FastGetSolutionStepValue(SUM_PARTICLES_WEIGHTS) += N[j] * N[j];
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        const auto it_part = mrParticles.NodesBegin() + i;
        const auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        if (p_elem != nullptr)
        {
            const auto& geom = p_elem->GetGeometry();
            Vector N;
            GetShapeFunctionsValues(N, geom, *it_part);

            auto& part_scalar = it_part->FastGetSolutionStepValue(mrMassVariable);
            auto& part_vector = it_part->FastGetSolutionStepValue(mrMomentumVariable);

            for (size_t j = 0; j < num_nodes; ++j)
            {
                const double part_weight = N[j] * N[j] / geom[j].FastGetSolutionStepValue(SUM_PARTICLES_WEIGHTS);
                const double part_area = it_part->FastGetSolutionStepValue(PARTICLE_AREA);
                const double sum_areas = geom[j].GetValue(SUM_AREAS);
                const auto delta_scalar = geom[j].FastGetSolutionStepValue(mrMassVariable,0) - geom[j].FastGetSolutionStepValue(mrMassVariable,1);
                const auto delta_vector = geom[j].FastGetSolutionStepValue(mrMomentumVariable,0) - geom[j].FastGetSolutionStepValue(mrMomentumVariable,1);
                part_scalar += part_weight * sum_areas / part_area / num_nodes * delta_scalar;
                part_vector += part_weight * sum_areas / part_area / num_nodes * delta_vector;
            }
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::InitializeParticle(NodeType& rParticle, Element& rElement)
{
    const auto& geom = rElement.GetGeometry();
    const size_t num_nodes = geom.size();
    const size_t num_particles = geom.IntegrationPointsNumber(mQuadratureOrder);
    const double particle_weight = 1 / static_cast<double>(num_particles);
    const double elem_area = geom.Area();

    Vector N;
    GetShapeFunctionsValues(N, geom, rParticle);

    double scalar_value = 0.0;
    array_1d<double,3> vector_value = ZeroVector(3);
    for (size_t i = 0; i < num_nodes; ++i)
    {
        scalar_value += N[i] * geom[i].FastGetSolutionStepValue(mrMassVariable);
        vector_value += N[i] * geom[i].FastGetSolutionStepValue(mrMomentumVariable);
    }

    rParticle.SetValue(CURRENT_ELEMENT, &rElement);
    rParticle.FastGetSolutionStepValue(PARTICLE_WEIGHT) = particle_weight;
    rParticle.FastGetSolutionStepValue(PARTICLE_AREA) = elem_area * particle_weight;
    rParticle.FastGetSolutionStepValue(mrMassVariable) = scalar_value;
    rParticle.FastGetSolutionStepValue(mrMomentumVariable) = vector_value;
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::GetShapeFunctionsValues(
    Vector& rN,
    const GeometryType& rGeom,
    const NodeType& rParticle)
{
    array_1d<double,3> local_coordinates;
    rGeom.PointLocalCoordinates(local_coordinates, rParticle);
    rGeom.ShapeFunctionsValues(rN, local_coordinates);    
}

template<size_t TDim>
int MoveShallowParticlesProcess<TDim>::Check()
{
    const NodeType& r_node = *mrModelPart.NodesBegin();

    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrMomentumVariable, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrMassVariable, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DELTA_VECTOR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DELTA_SCALAR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MEAN_SIZE, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MASS_WEIGHT, r_node);

    KRATOS_CHECK_EQUAL(mrParticles.NumberOfNodes(), 0);

    return 0;
}

template class MoveShallowParticlesProcess<1>;
template class MoveShallowParticlesProcess<2>;

} // namespace Kratos