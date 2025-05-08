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
                if (particle.mBallToRigidFaceStoredInfo.find(id2) == particle.mBallToRigidFaceStoredInfo.end() ||
                    particle.mBallToRigidFaceStoredInfo[id2].indentation <= 0.0 ||
                    particle.mNeighborInContact == false)
                    continue;

                // Increment number of contacts
                is_inner_particle = false;
                mAvgCoordNum++;
                mNumContacts++;

                // Normal vector
                const double d  = r1 - particle.mBallToRigidFaceStoredInfo[id2].indentation;
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
                if (particle.mBallToBallStoredInfo.find(id2) == particle.mBallToBallStoredInfo.end() ||
                    particle.mBallToBallStoredInfo[id2].indentation <= 0.0 ||
                    particle.mNeighborInContact == false)
                    continue;

                // Normal vector
                const double d  =  particle.mNeighborDistance;
                const double nx = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Check for inner contact (ATTENTION: THIS IS SPECIFIC FOR THERMAL PIPE MODEL OF HEAT TRANSFER)
                const double pipe_length = d;
                const double pipe_width = 2.0 * particle.mContactRadius;
                const double total_pipe_area = pipe_length * pipe_width;
                const double inner_pipe_area = ComputePipeAreaInner(coords1, normal, pipe_length, pipe_width, 8);
                const double inner_contact_ratio = inner_pipe_area / total_pipe_area;
                bool is_inner_contact = (inner_contact_ratio != 0.0);

                // Increment number of contacts
                mAvgCoordNum++;
                if (is_inner_particle) mAvgCoordNumInner++;

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);
                if (is_inner_contact) AddContactToRoseDiagram(mRoseDiagramInner, normal);

                // Unique contacts (each binary contact evaluated only once)
                if (id1 < id2) {
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
    // Compute overlap area between given thermal pipe and inner RVE area
    // by approximating the pipe into a set of longitudinal axes (as a simplification for calculations).
    double RVEWallBoundaryThermal2D::ComputePipeAreaInner(std::vector<double>& coords_ref, std::vector<double>& normal, double length, double width, int width_disc) {
        // Compute pipe corner coordinates based on the reference coordinates (at pipe center) and a tangential direction
        const double nx =  normal[0];
        const double ny =  normal[1];
        const double tx = -ny;
        const double ty =  nx;
        const double x  =  coords_ref[0] + tx * width/2;
        const double y  =  coords_ref[1] + ty * width/2;

        // Average pipe length inside inner RVE area
        double avg_inner_len = 0.0;
        double width_bin = width / (width_disc-1);

        for (unsigned int i = 0; i < width_disc; i++) {
            // End coordinates of current branch
            const double x1 = x  - tx * width_bin * i;
            const double y1 = y  - ty * width_bin * i;
            const double x2 = x1 + nx * length;
            const double y2 = y1 + ny * length;
            std::vector<double> coords1 = {x1,y1};
            std::vector<double> coords2 = {x2,y2};

            // Accumulate branch lengths inside inner RVE area
            avg_inner_len += RVEWallBoundary2D::ComputeBranchLengthInner(coords1, coords2);
        }
        
        // Approximated area inside inner RVE area
        return width * avg_inner_len / width_disc;
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
        RVEWallBoundary2D::WriteFileHeadersTensorResults();
        mFileTensorResults << "ROW 5: [K11 K12 K21 K22] - ALL CONTACTS" << std::endl;
        mFileTensorResults << "ROW 6: [K11 K12 K21 K22] - INN CONTACTS" << std::endl;
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::WriteResultFilesTensorResults(void) {
        RVEWallBoundary2D::WriteResultFilesTensorResults();
        mFileTensorResults << std::endl;
        mFileTensorResults << "[" 
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensor(0,0)
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensor(0,1)
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensor(1,0)
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensor(1,1)
                           << "]  "
                           << std::endl;

        mFileTensorResults << "[" 
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensorInner(0,0)
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensorInner(0,1)
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensorInner(1,0)
                           << std::setw(18) << std::left << std::fixed << std::setprecision(10) << mConductivityTensorInner(1,1)
                           << "]";
    }
}
