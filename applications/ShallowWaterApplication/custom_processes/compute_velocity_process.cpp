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
#include "includes/checks.h"
#include "compute_velocity_process.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

ComputeVelocityProcess::ComputeVelocityProcess(ModelPart& rModelPart, double Epsilon)
 : mrModelPart(rModelPart)
{
    mEpsilon = Epsilon;
    this->Initialize();
}

void ComputeVelocityProcess::Initialize()
{
    if (!mInitializeWasPerformed)
    {
        if (mrModelPart.NumberOfNodes() != 0) {
            mNumNodes = mrModelPart.ElementsBegin()->GetGeometry().size();
            this->ComputeMassMatrix();
            mInitializeWasPerformed = true;
        }

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            (mrModelPart.NodesBegin() + i)->SetValue(INTEGRATION_WEIGHT, 0.0);
        }
    }
}

void ComputeVelocityProcess::Execute()
{
    this->Initialize();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        const auto it_node = mrModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
        it_node->GetValue(INTEGRATION_WEIGHT) = 0.0;
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        const auto it_elem = mrModelPart.ElementsBegin() + i;
        double height = 0.0;
        Vector nodal_discharge_x(mNumNodes);
        Vector nodal_discharge_y(mNumNodes);
        Geometry<Node<3>>& r_geom = it_elem->GetGeometry();
        for (unsigned int j = 0; j < mNumNodes; ++j)
        {
            height += r_geom[j].FastGetSolutionStepValue(HEIGHT);
            nodal_discharge_x[j] = r_geom[j].FastGetSolutionStepValue(MOMENTUM_X);
            nodal_discharge_y[j] = r_geom[j].FastGetSolutionStepValue(MOMENTUM_Y);
        }
        height /= mNumNodes;
        Matrix mass_matrix = mMassMatrix / (std::abs(height) + mEpsilon);

        Vector nodal_velocity_x(mNumNodes);
        Vector nodal_velocity_y(mNumNodes);
        nodal_velocity_x = mNumNodes * prod(mass_matrix, nodal_discharge_x);
        nodal_velocity_y = mNumNodes * prod(mass_matrix, nodal_discharge_y);
        for (unsigned int j = 0; j < mNumNodes; ++j)
        {
            #pragma omp atomic
            r_geom[j].FastGetSolutionStepValue(VELOCITY_X) += nodal_velocity_x[j];
            r_geom[j].FastGetSolutionStepValue(VELOCITY_Y) += nodal_velocity_y[j];
            r_geom[j].GetValue(INTEGRATION_WEIGHT) += 1.0;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        const auto it_node = mrModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(VELOCITY) /= it_node->GetValue(INTEGRATION_WEIGHT);
    }
}

void ComputeVelocityProcess::Clear()
{
    mMassMatrix.clear();
    mInitializeWasPerformed = false;
}

void ComputeVelocityProcess::ComputeMassMatrix()
{
    mMassMatrix.resize(mNumNodes, mNumNodes, false);
    if (mNumNodes == 2)
    {
        double one_sixth = 1.0 / 6.0;
        mMassMatrix(0,0) = 2.0;
        mMassMatrix(0,1) = 1.0;
        mMassMatrix(1,0) = 1.0;
        mMassMatrix(1,1) = 2.0;
        mMassMatrix *= one_sixth;
    }
    else if (mNumNodes == 3)
    {
        double one_twelve = 1.0 / 12.0;
        mMassMatrix(0,0) = 2.0;
        mMassMatrix(0,1) = 1.0;
        mMassMatrix(0,2) = 1.0;
        mMassMatrix(1,0) = 1.0;
        mMassMatrix(1,1) = 1.0;
        mMassMatrix(1,2) = 2.0;
        mMassMatrix(2,0) = 1.0;
        mMassMatrix(2,1) = 1.0;
        mMassMatrix(2,2) = 2.0;
        mMassMatrix *= one_twelve;
    }
    else if (mNumNodes == 4)
    {
        double one_thirty_sixth = 1.0 / 36.0;
        mMassMatrix(0,0) = 4;
        mMassMatrix(0,1) = 2;
        mMassMatrix(0,2) = 1;
        mMassMatrix(0,3) = 2;
        mMassMatrix(1,0) = 2;
        mMassMatrix(1,1) = 4;
        mMassMatrix(1,2) = 2;
        mMassMatrix(1,3) = 1;
        mMassMatrix(2,0) = 1;
        mMassMatrix(2,1) = 2;
        mMassMatrix(2,2) = 4;
        mMassMatrix(2,3) = 2;
        mMassMatrix(3,0) = 2;
        mMassMatrix(3,1) = 1;
        mMassMatrix(3,2) = 2;
        mMassMatrix(3,3) = 4;
        mMassMatrix *= one_thirty_sixth;
    }
    else
    {
        KRATOS_ERROR << this->Info() << "Method implemented for lines, triangles and quadrilaterals" << std::endl;
    }
}

}