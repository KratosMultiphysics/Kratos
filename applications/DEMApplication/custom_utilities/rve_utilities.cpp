//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
#include "rve_utilities.h"

namespace Kratos
{
    //------------------------------------------------------------------------------------------------------------
    void RVEUtilities::Initialize(ModelPart& dem_model_part, ModelPart& fem_model_part) {
        if (!dem_model_part.GetProcessInfo()[RVE_ANALYSIS])
            return;
        InitializeVariables(dem_model_part, fem_model_part);
        AssembleWallElementVectors();
        ReadOldForces();
        OpenResultFiles();
        WriteFileHeaders();
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEUtilities::FinalizeSolutionStep(void) {
        if (!IsTimeToEvaluateRVE(mDemModelPart->GetProcessInfo()[TIME_STEPS])) return;
        const double prev_eff_stress = mEffStress;
        SetVertexCoordinates();
        SetVertexCoordinatesInner();
        PreProcessGlobalResults();
        ProcessGlobalResults();
        PostProcessGlobalResults();
        CheckEquilibrium(mEffStress, prev_eff_stress, 1.0E-10, 10);
        StopCompress();
        WriteResultFiles();
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEUtilities::Finalize(void) {
        CloseResultFiles();
    }

    //------------------------------------------------------------------------------------------------------------
    // Initialize variables that remain constant during RVE analysis, or that need initial value/memory allocation.
    void RVEUtilities::InitializeVariables(ModelPart& dem_model_part, ModelPart& fem_model_part) {
        mDemModelPart      = &dem_model_part;
        mFemModelPart      = &fem_model_part;
        mNumParticles      = mDemModelPart->GetCommunicator().LocalMesh().Elements().size();
        mNumWallElems      = mFemModelPart->GetCommunicator().LocalMesh().Conditions().size();
        mDim               = mDemModelPart->GetProcessInfo()[DOMAIN_SIZE];
        mAvgRadius         = ComputeAverageRadius();
        mVertexCoords      = ZeroMatrix(mDim,4*(mDim-1));
        mVertexCoordsInner = ZeroMatrix(mDim,4*(mDim-1));
        mIsMoving          = true;
        mIsEquilibrium     = false;
        mEquilibriumSteps  = 0;
    }

    //------------------------------------------------------------------------------------------------------------
    // Full computation of average radius of all particles (should be called only once).
    double RVEUtilities::ComputeAverageRadius(void) {
        double avgRadius = 0.0; 
        for (unsigned int i = 0; i < mNumParticles; i++) {
            ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
            avgRadius += particle.GetRadius();
        }
        return avgRadius /= mNumParticles;
    }

    //------------------------------------------------------------------------------------------------------------
    // Initialize global RVE results accumulated from individual particle and interaction contributions.
    void RVEUtilities::PreProcessGlobalResults(void) {
        mNumParticlesInner = 0;
        mNumContacts       = 0;
        mNumContactsInner  = 0;
        mAvgCoordNum       = 0.0;
        mAvgCoordNumInner  = 0.0;
        mAvgRadius         = 0.0;
        mVolSolid          = 0.0;
        mVolSolidInner     = 0.0;
        mWallForces        = 0.0;
        mFabricTensor      = ZeroMatrix(mDim,mDim);
        mFabricTensorInner = ZeroMatrix(mDim,mDim);
        mStressTensor      = ZeroMatrix(mDim,mDim);
        mStressTensorInner = ZeroMatrix(mDim,mDim);
        mRoseDiagram.assign(40,0);      // Discretized into 40 bins (angle ranges)
        mRoseDiagramInner.assign(40,0); // Discretized into 40 bins (angle ranges)
        mContactChain.clear();
    }

    //------------------------------------------------------------------------------------------------------------
    // Post-process, sum up global RVE results from individual particle/interaction contributions,
    // and homogenize discrete behavior into tensors.
    void RVEUtilities::PostProcessGlobalResults(void) {
        mAvgRadius        /= mNumParticles;
        mAvgCoordNum      /= mNumParticles;
        mAvgCoordNumInner /= mNumParticlesInner;
        mVolTotal          = ComputeVolumeRVE();
        mVolInner          = ComputeVolumeRVEInner();
        mPorosity          = ComputePorosity();
        mPorosityInner     = ComputePorosityInner();
        mWallStress        = mWallForces / ComputeSurfaceArea();
        Homogenize();
        EvaluateRoseUniformity();
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute RVE porosity considering its total volume.
    double RVEUtilities::ComputePorosity(void) {
        return 1.0 - mVolSolid / mVolTotal;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute RVE porosity considering its inner volume.
    double RVEUtilities::ComputePorosityInner(void) {
        return 1.0 - mVolSolidInner / mVolInner;
    }

    //------------------------------------------------------------------------------------------------------------
    // Perform homogenization procedures on the particle assembly to obtain tensorial variables for upscaling the discrete solution.
    bool RVEUtilities::Homogenize(void) {
        if (mNumContacts == 0 || mNumContactsInner == 0) return false;
        HomogenizeFabric();
        HomogenizeStress();
        return true;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute fabric tensor and its related properties.
    // Assumes that fabric tensor is already filled with element-to-element interaction dependent components.
    void RVEUtilities::HomogenizeFabric(void) {
        double deviatoric_fabric;
        double deviatoric_fabric_inner;
        double double_dot_product = 0.0;
        double double_dot_product_inner = 0.0;

        for (unsigned int i = 0; i < mDim; i++) {
            for (unsigned int j = 0; j < mDim; j++) {
                mFabricTensor(i,j)      /= mNumContacts;
                mFabricTensorInner(i,j) /= mNumContactsInner;
                if (i == j) {
                    deviatoric_fabric       = 4.0 * (mFabricTensor(i,j)      - (1.0/mDim));
                    deviatoric_fabric_inner = 4.0 * (mFabricTensorInner(i,j) - (1.0/mDim));
                }
                else {
                    deviatoric_fabric       = 4.0 * mFabricTensor(i,j);
                    deviatoric_fabric_inner = 4.0 * mFabricTensor(i,j);
                }
                double_dot_product       += 0.5 * deviatoric_fabric       * deviatoric_fabric;
                double_dot_product_inner += 0.5 * deviatoric_fabric_inner * deviatoric_fabric_inner;
            }
        }
        mFidx            = ComputeFabricIndex(mFabricTensor);
        mFidxInner       = ComputeFabricIndex(mFabricTensorInner);
        mAnisotropy      = sqrt(double_dot_product);
        mAnisotropyInner = sqrt(double_dot_product_inner);
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute effective stress tensor and its related properties.
    // Assumes that stress tensor is already filled with element-to-element interaction dependent components.
    void RVEUtilities::HomogenizeStress(void) {
        double deviatoric_stress;
        double deviatoric_stress_inner;
        double stress_trace = 0.0;
        double stress_trace_inner = 0.0;
        double double_dot_product = 0.0;
        double double_dot_product_inner = 0.0;

        for (unsigned int i = 0; i < mDim; i++) {
            for (unsigned int j = 0; j < mDim; j++) {
                mStressTensor(i,j)      /= mVolTotal;
                mStressTensorInner(i,j) /= mVolInner;
                if (i == j) {
                    stress_trace       += mStressTensor(i,j);
                    stress_trace_inner += mStressTensorInner(i,j);
                }
            }
        }
        mEffStress      = stress_trace       / mDim;
        mEffStressInner = stress_trace_inner / mDim;

        for (unsigned int i = 0; i < mDim; i++) {
            for (unsigned int j = 0; j < mDim; j++) {
                deviatoric_stress       = (i==j) ? mStressTensor(i,j)      - mEffStress      : mStressTensor(i,j);
                deviatoric_stress_inner = (i==j) ? mStressTensorInner(i,j) - mEffStressInner : mStressTensorInner(i,j);
                double_dot_product       += 0.5 * deviatoric_stress       * deviatoric_stress;
                double_dot_product_inner += 0.5 * deviatoric_stress_inner * deviatoric_stress_inner;
            }
        }
        mDevStress      = sqrt(double_dot_product);
        mDevStressInner = sqrt(double_dot_product_inner);
    }

    //------------------------------------------------------------------------------------------------------------
    // Check equilibrium of particles in RVE based on custom criteria and a given tolerance.
    void RVEUtilities::CheckEquilibrium(double curr_val, double prev_val, double tol, int max_eq_steps) {
        double ratio = std::abs((curr_val-prev_val) / curr_val);

        if (!mIsMoving && ratio < tol) {
            mEquilibriumSteps++;
        }
        else {
            mEquilibriumSteps = 0;
        }

        if (mEquilibriumSteps >= max_eq_steps) {
            mIsEquilibrium = true;
            mEquilibriumSteps++;
        }
        else {
            mIsEquilibrium = false;
            mEquilibriumSteps = 0;
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Check criteria to stop motion of RVE boundaries during compression (packing consolidation) stage.
    void RVEUtilities::StopCompress(void) {
        if (!mIsMoving)
            return;

        // Check if selected criteria for stopping boundary motion is satisfied
        bool check = true;
        if      (mConsolidationCriterion.compare("time") == 0)     check = (mConsolidationLimit < mDemModelPart->GetProcessInfo()[TIME]);
        else if (mConsolidationCriterion.compare("stress") == 0)   check = (mConsolidationLimit < std::abs(mEffStressInner));
        else if (mConsolidationCriterion.compare("porosity") == 0) check = (mConsolidationLimit > mPorosityInner);

        // Assign zero velocity to boundaries
        if (check)
            StopBoundaryMotion();
    }

    //------------------------------------------------------------------------------------------------------------
    // Assign zero velocity to boundaries.
    void RVEUtilities::StopBoundaryMotion(void) {
        mIsMoving = false;

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = mFemModelPart->SubModelPartsBegin(); sub_model_part != mFemModelPart->SubModelPartsEnd(); ++sub_model_part) {
            ModelPart &submp = *sub_model_part;
            array_1d<double, 3> &linear_velocity = submp[LINEAR_VELOCITY];
            linear_velocity[0] = 0.0;
            linear_velocity[1] = 0.0;
            linear_velocity[2] = 0.0;
        }

        ModelPart::ConditionsContainerType &r_conditions = mFemModelPart->GetCommunicator().LocalMesh().Conditions();
        for (unsigned int i = 0; i < mNumWallElems; i++) {
            ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
            DEMWall *p_wall = dynamic_cast<DEMWall*>(&(*it));

            for (unsigned int j = 0; j < p_wall->GetGeometry().size(); j++) {
                array_1d<double, 3> &wall_velocity = p_wall->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                noalias(wall_velocity) = ZeroVector(3);
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Read forces between neighbors saved from previous analyses.
    // File format for line i: particle id #neighbors neighbor-1_fx neighbor-1_fy neighbor-1_fz ... neighbor-j_fx neighbor-j_fy neighbor-j_fz
    void RVEUtilities::ReadOldForces(void) {
        // Open old force files with pre-defined names
        std::fstream f_old_elastic_forces_pp;
        std::fstream f_old_elastic_forces_pw;
        f_old_elastic_forces_pp.open("OLD_FORCES.txt", std::ios::in);
        f_old_elastic_forces_pw.open("OLD_FORCES_WALLS.txt", std::ios::in);

        // Read and store particle-particle (pp) forces
        if (f_old_elastic_forces_pp) {
            for (unsigned int i = 0; i < mNumParticles; i++) {
                ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
                SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
                int particle_id, n_neighbors;
                f_old_elastic_forces_pp >> particle_id >> n_neighbors;

                for (unsigned int j = 0; j < n_neighbors; j++) {
                    double fx, fy, fz;
                    f_old_elastic_forces_pp >> fx >> fy >> fz;
                    particle.mNeighbourElasticContactForces[j][0] = fx;
                    particle.mNeighbourElasticContactForces[j][1] = fy;
                    particle.mNeighbourElasticContactForces[j][2] = fz;
                }
            }
            f_old_elastic_forces_pp.close();
        }

        // Read and store particle-wall (pw) forces
        if (f_old_elastic_forces_pw) {
            for (unsigned int i = 0; i < mNumParticles; i++) {
                ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
                SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
                int particle_id, n_neighbors;
                f_old_elastic_forces_pw >> particle_id >> n_neighbors;

                for (unsigned int j = 0; j < n_neighbors; j++) {
                    double fx, fy, fz;
                    f_old_elastic_forces_pw >> fx >> fy >> fz;
                    particle.mNeighbourRigidFacesElasticContactForce[j][0] = fx;
                    particle.mNeighbourRigidFacesElasticContactForce[j][1] = fy;
                    particle.mNeighbourRigidFacesElasticContactForce[j][2] = fz;
                }
            }
            f_old_elastic_forces_pw.close();
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Open files to write selected results.
    void RVEUtilities::OpenResultFiles(void) {
        if (mWriteFreq != 0) {
            mFileGlobalResults.open("rve_results_global.txt", std::ios::out);
            KRATOS_ERROR_IF_NOT(mFileGlobalResults) << "Could not open file rve_results_global.txt!" << std::endl;
            
            mFileParticleResults.open("rve_results_particles.txt", std::ios::out);
            KRATOS_ERROR_IF_NOT(mFileParticleResults) << "Could not open file rve_results_particles.txt!" << std::endl;

            mFileContactResults.open("rve_results_contacts.txt", std::ios::out);
            KRATOS_ERROR_IF_NOT(mFileContactResults) << "Could not open file rve_results_contacts.txt!" << std::endl;

            mFileTensorResults.open("rve_results_tensors.txt", std::ios::out);
            KRATOS_ERROR_IF_NOT(mFileTensorResults) << "Could not open file rve_results_tensors.txt!" << std::endl;

            mFileRoseDiagram.open("rve_results_rose.txt", std::ios::out);
            KRATOS_ERROR_IF_NOT(mFileRoseDiagram) << "Could not open file rve_results_rose.txt!" << std::endl;
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Write file headers for selected results.
    void RVEUtilities::WriteFileHeaders(void) {
        if (mFileGlobalResults.is_open())   WriteFileHeadersGlobalResults();
        if (mFileParticleResults.is_open()) WriteFileHeadersParticleResults();
        if (mFileContactResults.is_open())  WriteFileHeadersContactResults();
        if (mFileTensorResults.is_open())   WriteFileHeadersTensorResults();
        if (mFileRoseDiagram.is_open())     WriteFileHeadersRoseDiagram();
    }

    //------------------------------------------------------------------------------------------------------------
    // Write selected results to opened files.
    void RVEUtilities::WriteResultFiles(void) {
        if (!IsTimeToPrintResults(mDemModelPart->GetProcessInfo()[TIME_STEPS])) return;
        if (mFileGlobalResults.is_open())   WriteResultFilesGlobalResults();
        if (mFileParticleResults.is_open()) WriteResultFilesParticleResults();
        if (mFileContactResults.is_open())  WriteResultFilesContactResults();
        if (mFileTensorResults.is_open())   WriteResultFilesTensorResults();
        if (mFileRoseDiagram.is_open())     WriteResultFilesRoseDiagram();
    }

    //------------------------------------------------------------------------------------------------------------
    // Close all result files
    void RVEUtilities::CloseResultFiles(void) {
        if (mFileGlobalResults.is_open())   mFileGlobalResults.close();
        if (mFileParticleResults.is_open()) mFileParticleResults.close();
        if (mFileContactResults.is_open())  mFileContactResults.close();
        if (mFileTensorResults.is_open())   mFileTensorResults.close();
        if (mFileRoseDiagram.is_open())     mFileRoseDiagram.close();
    }
}
