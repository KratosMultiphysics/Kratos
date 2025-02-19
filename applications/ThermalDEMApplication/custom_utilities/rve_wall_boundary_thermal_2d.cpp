//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
// 
#include "rve_wall_boundary_thermal_2d.h"

namespace Kratos
{
    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::PreProcessGlobalResults(void) {
        RVEWallBoundary2D::PreProcessGlobalResults();
        mConductivityTensor      = ZeroMatrix(mDim,mDim);
        mConductivityTensorInner = ZeroMatrix(mDim,mDim);
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::ProcessGlobalResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        for (int i = 0; i < mNumParticles; i++) {
            ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&>(*it);

            // Particle properties
            const int id1 = particle.GetId();
            const double r1 = particle.GetRadius();
            const double x1 = particle.GetGeometry()[0][0];
            const double y1 = particle.GetGeometry()[0][1];
            std::vector<double> coords1 = {x1, y1};
            bool is_inner_particle = true;

            // Accumulate particle properties
            mAvgRadius     += r1;
            mVolSolid      += ComputeVolumeParticle(particle);
            mVolSolidInner += ComputeVolumeParticleInner(particle);

            // Loop over contacts with walls
            for (unsigned int j = 0; j < particle.mNeighbourRigidFaces.size(); j++) {
                // Set neighbor element
                if (particle.mNeighbourRigidFaces[j] == NULL || j > 1) continue;
                particle.mNeighbor_w    = dynamic_cast<DEMWall*>(particle.mNeighbourRigidFaces[j]);
                particle.mNeighborType  = WALL_NEIGHBOR_CONTACT;
                particle.mNeighborIndex = j;

                // Neighbor properties
                const int id2 = particle.mNeighbourRigidFaces[j]->GetId();
                const double x2 = particle.mNeighbourRigidFaces[j]->GetGeometry()[0][0];
                const double y2 = particle.mNeighbourRigidFaces[j]->GetGeometry()[0][1];
                const double x3 = particle.mNeighbourRigidFaces[j]->GetGeometry()[1][0];
                const double y3 = particle.mNeighbourRigidFaces[j]->GetGeometry()[1][1];

                // Compute interaction properties
                particle.ComputeInteractionProps(r_process_info);

                // Check for valid existing contact
                if (particle.mBallToRigidFaceStoredInfo.find(id2) == particle.mBallToRigidFaceStoredInfo.end() || !particle.mNeighborInContact)
                    continue;
                const double indent = particle.mBallToRigidFaceStoredInfo[id2].indentation;
                if (indent <= 0.0 || particle.CheckAdiabaticNeighbor())
                    continue;

                // Increment number of contacts
                mAvgCoordNum++;
                mNumContacts++;
                is_inner_particle = false;

                // Normal vector
                const double d  = r1 - indent;
                const double nx = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);

                // Applied wall force (normal component)
                const double fx = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[0];
                const double fy = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[1];
                std::vector<double> force = {fx, fy};
                mWallForces += std::abs(fx * nx + fy * ny);

                // Contact results
                std::vector<double> chain{x1, y1, 0.0, x1+branch[0], y1+branch[1], 0.0, fx, fy, 0.0};
                mContactChain.insert(mContactChain.end(), chain.begin(), chain.end());

                // Effective conductivity
                const double keff = particle.GetDirectConductionModel().ComputeEffectiveThermalConductivity(r_process_info, &particle);

                // Tensors
                for (unsigned int k = 0; k < mDim; k++) {
                    for (unsigned int l = 0; l < mDim; l++) {
                        mFabricTensor(k,l)       += normal[k] * normal[l];
                        mStressTensor(k,l)       += branch[k] * force[l];
                        mConductivityTensor(k,l) += normal[k] * normal[l] * keff;
                    }
                }
            }
            if (is_inner_particle) {
                mNumParticlesInner++;
            }

            // Loop over contacts with particles
            for (unsigned int j = 0; j < particle.mNeighbourElements.size(); j++) {
                // Set neighbor element  
                particle.mNeighbor_p    = dynamic_cast<ThermalSphericParticle*>(particle.mNeighbourElements[j]);
                particle.mNeighborType  = PARTICLE_NEIGHBOR;
                particle.mNeighborIndex = j;

                // Neighbor properties
                const int id2 = particle.mNeighbourElements[j]->GetId();
                const double r2 = particle.mNeighbourElements[j]->GetRadius();
                const double x2 = particle.mNeighbourElements[j]->GetGeometry()[0][0];
                const double y2 = particle.mNeighbourElements[j]->GetGeometry()[0][1];
                std::vector<double> coords2 = {x2, y2};

                // Compute interaction properties
                particle.ComputeInteractionProps(r_process_info);

                // Check for valid existing contact
                if (particle.mBallToBallStoredInfo.find(id2) == particle.mBallToBallStoredInfo.end() || !particle.mNeighborInContact)
                    continue;
                const double indent = particle.mBallToBallStoredInfo[id2].indentation;
                if (indent <= 0.0 || particle.CheckAdiabaticNeighbor())
                    continue;

                // Increment number of contacts
                mAvgCoordNum++;
                if (is_inner_particle) mAvgCoordNumInner++;

                // Normal vector
                const double d  = r1 + r2 - indent;
                const double nx = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);
                if (is_inner_particle) AddContactToRoseDiagram(mRoseDiagramInner, normal);

                // Unique contacts (each binary contact evaluated only once)
                if (id1 < id2) {
                    // Check for inner contact TODO: COMPLETE THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    const double inner_contact_len = ComputeBranchLengthInner(coords1, coords2);
                    const double inner_pipe_area = ComputePipeAreaInner();

                    const double inner_contact_ratio = inner_contact_len / d;
                    const double inner_area_ratio = inner_pipe_area / d;

                    bool is_inner_contact = (inner_contact_len != 0.0);

                    // Increment number of unique contacts
                    mNumContacts++;
                    if (is_inner_contact) mNumContactsInner++;

                    // Contact results
                    const double fx = particle.mBallToBallStoredInfo[id2].global_contact_force[0];
                    const double fy = particle.mBallToBallStoredInfo[id2].global_contact_force[1];
                    std::vector<double> force = {fx, fy};
                    std::vector<double> chain{x1, y1, 0.0, x2, y2, 0.0, fx, fy, 0.0};
                    mContactChain.insert(mContactChain.end(), chain.begin(), chain.end());

                    // Effective conductivity
                    const double keff = particle.GetDirectConductionModel().ComputeEffectiveThermalConductivity(r_process_info, &particle);

                    // Tensors
                    for (unsigned int k = 0; k < mDim; k++) {
                        for (unsigned int l = 0; l < mDim; l++) {
                            mFabricTensor(k,l)       += normal[k] * normal[l];
                            mStressTensor(k,l)       += branch[k] * force[l];
                            mConductivityTensor(k,l) += normal[k] * normal[l] * keff;
                            if (is_inner_contact) {
                                mFabricTensorInner(k,l)       += normal[k] * normal[l];
                                mStressTensorInner(k,l)       += branch[k] * force[l] * inner_contact_ratio;
                                mConductivityTensorInner(k,l) += normal[k] * normal[l] * keff * inner_contact_ratio;
                            }
                        }
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    double RVEWallBoundaryThermal2D::ComputePipeAreaInner(void) {
        return 0.0;
    }

    //------------------------------------------------------------------------------------------------------------
    bool RVEWallBoundaryThermal2D::Homogenize(void) {
        if (!RVEWallBoundary2D::Homogenize()) return false;
        HomogenizeConductivity();
        return true;
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::HomogenizeConductivity(void) {
        for (unsigned int i = 0; i < mDim; i++) {
            for (unsigned int j = 0; j < mDim; j++) {
                mConductivityTensor(i,j)      /= mVolTotal;
                mConductivityTensorInner(i,j) /= mVolInner;
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::WriteFileHeadersTensorResults(void) {
        mFileTensorResults << "R = ROW, C = COLUMN, F = FABRIC TENSOR, S = STRESS TENSOR, K = CONDUCTIVITY TENSOR" << std::endl;
        mFileTensorResults << "Ri -> C1: STEP | C2: TIME" << std::endl;
        mFileTensorResults << "Rj -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Rk -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY (INN CONTACTS)" << std::endl;
        mFileTensorResults << "Rl -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS  | C7: WALL STRESS (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Rm -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS (INN CONTACTS)" << std::endl;
        mFileTensorResults << "Rn -> C1-C4: [K11 K12 K21 K22] (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Ro -> C1-C4: [K11 K12 K21 K22] (INN CONTACTS)" << std::endl;
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::WriteResultFilesTensorResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        mFileTensorResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mFidx              << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAnisotropy        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mFidxInner              << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAnisotropyInner        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mEffStress         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mDevStress         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mWallStress        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mEffStressInner         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mDevStressInner         << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(1,1) << " "
                           << "] "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(1,1) << " "
                           << "] "
                           << std::endl;

        mFileTensorResults << std::endl;
    }
}
