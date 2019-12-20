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
    Variable<double>& rScalarVariable)
 : mrModelPart(rModelPart)
 , mrParticles(rParticles)
 , mrMomentumVariable(rVectorVariable)
 , mrMassVariable(rScalarVariable)
 , mConvectionOperator(rModelPart)
{
    Check();

    mrParticles.AddNodalSolutionStepVariable(mrMassVariable);
    mrParticles.AddNodalSolutionStepVariable(mrMomentumVariable);
    mrParticles.AddNodalSolutionStepVariable(INTEGRATION_WEIGHT);
    mLastParticleId = 0;
    mDryTreshold = mrModelPart.GetProcessInfo()[DRY_HEIGHT];

    InitializeVariables();
    ComputeMeanSize();
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
    ComputeMeanVelocity();
    MoveParticles();
    TransferLagrangianToEulerian();
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::ExecuteFinalizeSolutionStep()
{
    TransferEulerianToLagrangian();
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
        if (it_elem->Is(ACTIVE))
        {
            int& num_of_particles = it_elem->GetValue(NUMBER_OF_PARTICLES);
            if (num_of_particles < num_nodes)
            {
                const GeometryType::IntegrationPointsArrayType& integration_points = it_elem->GetGeometry().IntegrationPoints();
                for (size_t g = 0; g < num_nodes; ++g)
                {
                    const double x = integration_points[g].X();
                    const double y = integration_points[g].Y();
                    const double z = integration_points[g].Z();
                    Node<3>::Pointer p_part;
                    #pragma omp critical
                    p_part = mrParticles.CreateNewNode(++mLastParticleId, x, y, z);
                    InitializeParticle(*p_part, *it_elem);
                    ++num_of_particles;
                }
            }
        }
    }
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
        int& num_of_particles = it_part->GetValue(CURRENT_ELEMENT)->GetValue(NUMBER_OF_PARTICLES);
        #pragma omp critical
        ++num_of_particles;
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::MoveParticles()
{
    double time_step = mrModelPart.GetProcessInfo()[DELTA_TIME];
    size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    mConvectionOperator.InitializeSearch();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;
        auto p_elem = it_part->GetValue(CURRENT_ELEMENT)->shared_from_this();
        Vector N(num_nodes);
        mConvectionOperator.Convect(time_step, *it_part, p_elem, N);
        it_part->GetValue(CURRENT_ELEMENT) = p_elem;
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::TransferLagrangianToEulerian()
{
    size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;

        // Particle contributions
        double particle_weight = it_part->FastGetSolutionStepValue(INTEGRATION_WEIGHT);
        double scalar_contribution = 0;
        array_1d<double,3> vector_contribution = ZeroVector(3);
        scalar_contribution = particle_weight * it_part->FastGetSolutionStepValue(mrMassVariable);
        vector_contribution = particle_weight * it_part->FastGetSolutionStepValue(mrMomentumVariable);

        // Transfer the values
        auto& geom = it_part->GetValue(CURRENT_ELEMENT)->GetGeometry();
        for (size_t j = 0; j < num_nodes; ++j)
        {
            geom[j].SetLock();
            geom[j].FastGetSolutionStepValue(mrMassVariable) = scalar_contribution;
            geom[j].FastGetSolutionStepValue(mrMomentumVariable) = vector_contribution;
            geom[j].UnSetLock();
        }
    }
}

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::TransferEulerianToLagrangian()
{
    size_t num_nodes = mrModelPart.ElementsBegin()->GetGeometry().size();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrParticles.NumberOfNodes()); ++i)
    {
        auto it_part = mrParticles.NodesBegin() + i;
        auto p_elem = it_part->GetValue(CURRENT_ELEMENT);
        GeometryType geom = p_elem->GetGeometry();
        Vector N;
        geom.ShapeFunctionsValues(N, *it_part);

        // Variables to update
        double& scalar_value = it_part->FastGetSolutionStepValue(mrMassVariable);
        double delta_scalar = 0;
        array_1d<double,3> vector_value = it_part->FastGetSolutionStepValue(mrMomentumVariable);
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

template<size_t TDim>
void MoveShallowParticlesProcess<TDim>::InitializeParticle(Node<3>& rParticle, Element& rElement)
{
    GeometryType geom = rElement.GetGeometry();
    const size_t num_nodes = geom.size();
    const double particle_weight = 1 / static_cast<double>(num_nodes);

    Vector N;
    geom.ShapeFunctionsValues(N, rParticle);

    double scalar_value = 0.0;
    array_1d<double,3> vector_value = ZeroVector(3);
    for (size_t i = 0; i < num_nodes; ++i)
    {
        scalar_value += geom[i].FastGetSolutionStepValue(mrMassVariable);
        vector_value += geom[i].FastGetSolutionStepValue(mrMomentumVariable);
    }

    rParticle.SetValue(CURRENT_ELEMENT, &rElement);
    rParticle.FastGetSolutionStepValue(INTEGRATION_WEIGHT) = particle_weight;
    rParticle.FastGetSolutionStepValue(mrMassVariable) = scalar_value;
    rParticle.FastGetSolutionStepValue(mrMomentumVariable) = vector_value;
}

template<size_t TDim>
int MoveShallowParticlesProcess<TDim>::Check()
{
    const Node<3>& r_node = *mrModelPart.NodesBegin();

    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrMomentumVariable, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mrMassVariable, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DELTA_VECTOR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DELTA_SCALAR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MEAN_SIZE, r_node);

    KRATOS_CHECK_EQUAL(mrParticles.NumberOfNodes(), 0);

    return 0;
}

template class MoveShallowParticlesProcess<1>;
template class MoveShallowParticlesProcess<2>;

} // namespace Kratos