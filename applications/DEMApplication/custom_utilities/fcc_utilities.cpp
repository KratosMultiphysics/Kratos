//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (r.rangel@utwente.nl)
//
//  Utilities for the analysis of FCC failure mechanism under uniaxial compression:
//  * Stage 1: Isotropic compression until reaching target stress.
//  * Stage 2: Relaxation until reaching static equilibrium.
//  * Stage 3: Servo-controlled axial compression (Z direction) until failure is detected.
//
#include "fcc_utilities.h"
#define CONSOLE_PRINT_FREQ 1000
#define EVAL_FREQ_STAGE2 10000
#define FAIL_COUNT_MAX 20
namespace Kratos
{
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::Initialize(ModelPart& dem_model_part, ModelPart& fem_model_part) {
        mDemModelPart = &dem_model_part;
        mFemModelPart = &fem_model_part;
        const ProcessInfo& processInfo = mDemModelPart->GetProcessInfo();
        auto pElement = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin();
        SphericParticle& particle = static_cast<SphericParticle&>(**pElement);
        mRadius = particle.GetRadius();
        mDensity = particle.GetProperties()[PARTICLE_DENSITY];
        mNumParticles = mDemModelPart->GetCommunicator().LocalMesh().Elements().size();
        mNumWallElems = mFemModelPart->GetCommunicator().LocalMesh().Conditions().size();
        mDim = 3;
        mStage = 1;
        mRefTime = processInfo[TIME];
        mInertialNum = 0.0;
        mEnergyRatio = 0.0;
        mFailCount = 0;
        mNumParticlesInn = 0;
        mInnIds.clear();
        mWallXMin.clear();
        mWallXMax.clear();
        mWallYMin.clear();
        mWallYMax.clear();
        mWallZMin.clear();
        mWallZMax.clear();
        mFileOut.open("output.csv", std::ios::out);
        mFileOut << "Time(s),StressXX(Pa),StressYY(Pa),StressZZ(Pa),StrainZZ,MCN" << std::endl;
        AssembleWallVectors();
        IdentifyInnerParticles();
        UpdateSystem();
        UpdateWallVelocities();
        mLenZRef = mLenZ;
        KRATOS_INFO("FCCUtilities") << "FCCUtilities initialized!" << std::endl;
        KRATOS_INFO("FCCUtilities") << "Number of particles (inner): " << mNumParticles << " (" << mNumParticlesInn << ")" << std::endl;
        KRATOS_INFO("FCCUtilities") << "Number of walls: " << mNumWallElems << std::endl;
        KRATOS_INFO("FCCUtilities") << "Box length X = " << mLenX << " (" << mWallXMin[0]->GetGeometry()[0][0] << ", " << mWallXMax[0]->GetGeometry()[0][0] << ")" << std::endl;
        KRATOS_INFO("FCCUtilities") << "Box length Y = " << mLenY << " (" << mWallYMin[0]->GetGeometry()[0][1] << ", " << mWallYMax[0]->GetGeometry()[0][1] << ")" << std::endl;
        KRATOS_INFO("FCCUtilities") << "Box length Z = " << mLenZ << " (" << mWallZMin[0]->GetGeometry()[0][2] << ", " << mWallZMax[0]->GetGeometry()[0][2] << ")" << std::endl;
        KRATOS_INFO("FCCUtilities") << "Box volume = " << mLenX * mLenY * mLenZ << std::endl;
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::FinalizeSolutionStep(void) {
        const ProcessInfo& processInfo = mDemModelPart->GetProcessInfo();
        int step = processInfo[TIME_STEPS];
        if (!IsTimeToEvaluate(step)) return;
        UpdateSystem();
        if (mStage == 1) {
            if (std::abs(mStressMean) >= mTargetStress) {
                KRATOS_INFO("FCCUtilities") << "Stage 1 completed: Mean stress reached target " << std::abs(mStressMean) << std::endl << std::endl;
                mStage = 2;
                mRunFreq = mRunFreq * EVAL_FREQ_STAGE2;
                mRefTime = processInfo[TIME];
                mLenZRef = mLenZ;
            }
        }
        else if (mStage == 2) {
            if (mEnergyRatio <= mEnergyRatioMax) {
                KRATOS_INFO("FCCUtilities") << "Stage 2 completed: Energy ratio reached target " << mEnergyRatio << std::endl << std::endl;
                mStage = 3;
                mRunFreq = mRunFreq / EVAL_FREQ_STAGE2;
                mRefTime = processInfo[TIME];
                mLenZRef = mLenZ;
            }
        }
        else if (mStage == 3) {
            if (IsTimeToPrint(step))
                WriteOutputs();
            if (mInertialNum > mInertialNumMax)
                KRATOS_ERROR << "Inertial number exceeded maximum value " << mInertialNum << std::endl;
            if (CheckFail()) {
                if (mFileOut.is_open()) mFileOut.close();
                KRATOS_ERROR << "Stage 3 completed! MCN below threshold after " << mFailCount << " counts" << std::endl;
            }
        }
        UpdateWallVelocities();
        if (step % CONSOLE_PRINT_FREQ == 0) PrintStepInfo();
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::AssembleWallVectors(void) {
        ModelPart::ConditionsContainerType &rConditions = mFemModelPart->GetCommunicator().LocalMesh().Conditions();
        // First loop: determine bounding box of wall elements
        double xMin = std::numeric_limits<double>::max();
        double xMax = std::numeric_limits<double>::lowest();
        double yMin = std::numeric_limits<double>::max();
        double yMax = std::numeric_limits<double>::lowest();
        double zMin = std::numeric_limits<double>::max();
        double zMax = std::numeric_limits<double>::lowest();
        for (int i = 0; i < mNumWallElems; i++) {
            ModelPart::ConditionsContainerType::iterator it = rConditions.ptr_begin() + i;
            DEMWall *pWall = dynamic_cast<DEMWall *>(&(*it));
            Condition::GeometryType &geom = pWall->GetGeometry();
            const double x1 = geom[0][0]; const double y1 = geom[0][1]; const double z1 = geom[0][2];
            const double x2 = geom[1][0]; const double y2 = geom[1][1]; const double z2 = geom[1][2];
            const double x3 = geom[2][0]; const double y3 = geom[2][1]; const double z3 = geom[2][2];
            xMin = std::min({xMin, x1, x2, x3});
            xMax = std::max({xMax, x1, x2, x3});
            yMin = std::min({yMin, y1, y2, y3});
            yMax = std::max({yMax, y1, y2, y3});
            zMin = std::min({zMin, z1, z2, z3});
            zMax = std::max({zMax, z1, z2, z3});
        }
        // Second loop: Classify wall elements
        for (int i = 0; i < mNumWallElems; i++) {
            ModelPart::ConditionsContainerType::iterator it = rConditions.ptr_begin() + i;
            DEMWall *pWall = dynamic_cast<DEMWall *>(&(*it));
            Condition::GeometryType &geom = pWall->GetGeometry();
            const double x1 = geom[0][0]; const double y1 = geom[0][1]; const double z1 = geom[0][2];
            const double x2 = geom[1][0]; const double y2 = geom[1][1]; const double z2 = geom[1][2];
            const double x3 = geom[2][0]; const double y3 = geom[2][1]; const double z3 = geom[2][2];
            if (EqualValues(x1, x2) && EqualValues(x2, x3)) { // X plane
                if      (EqualValues(x1, xMin)) mWallXMin.push_back(pWall);
                else if (EqualValues(x1, xMax)) mWallXMax.push_back(pWall);
            }
            else if (EqualValues(y1, y2) && EqualValues(y2, y3)) { // Y plane
                if      (EqualValues(y1, yMin)) mWallYMin.push_back(pWall);
                else if (EqualValues(y1, yMax)) mWallYMax.push_back(pWall);
            }
            else if (EqualValues(z1, z2) && EqualValues(z2, z3)) { // Z plane
                if      (EqualValues(z1, zMin)) mWallZMin.push_back(pWall);
                else if (EqualValues(z1, zMax)) mWallZMax.push_back(pWall);
            }
        }
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::IdentifyInnerParticles(void) {
        // First loop: determine bounding box of particles
        double xMin = std::numeric_limits<double>::max();
        double xMax = std::numeric_limits<double>::lowest();
        double yMin = std::numeric_limits<double>::max();
        double yMax = std::numeric_limits<double>::lowest();
        double zMin = std::numeric_limits<double>::max();
        double zMax = std::numeric_limits<double>::lowest();
        for (int i = 0; i < mNumParticles; i++) {
            auto it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            SphericParticle& particle = static_cast<SphericParticle&>(**it);
            const double x = particle.GetGeometry()[0][0];
            const double y = particle.GetGeometry()[0][1];
            const double z = particle.GetGeometry()[0][2];
            xMin = std::min({xMin, x});
            xMax = std::max({xMax, x});
            yMin = std::min({yMin, y});
            yMax = std::max({yMax, y});
            zMin = std::min({zMin, z});
            zMax = std::max({zMax, z});
        }
        // Second loop: Identify inner particles
        for (int i = 0; i < mNumParticles; i++) {
            auto it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            SphericParticle& particle = static_cast<SphericParticle&>(**it);
            const double x = particle.GetGeometry()[0][0];
            const double y = particle.GetGeometry()[0][1];
            const double z = particle.GetGeometry()[0][2];
            bool on_boundary = 
            EqualValues(x, xMin) || EqualValues(x, xMax) ||
            EqualValues(y, yMin) || EqualValues(y, yMax) ||
            EqualValues(z, zMin) || EqualValues(z, zMax);
            if (!on_boundary) {
                mInnIds.push_back(particle.GetId());
                mNumParticlesInn += 1;
            }
        }
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::UpdateSystem() {
        double xMin = mWallXMin[0]->GetGeometry()[0][0];
        double xMax = mWallXMax[0]->GetGeometry()[0][0];
        double yMin = mWallYMin[0]->GetGeometry()[0][1];
        double yMax = mWallYMax[0]->GetGeometry()[0][1];
        double zMin = mWallZMin[0]->GetGeometry()[0][2];
        double zMax = mWallZMax[0]->GetGeometry()[0][2];
        mLenX = xMax - xMin;
        mLenY = yMax - yMin;
        mLenZ = zMax - zMin;
        //mStrainZ = (mLenZ - mLenZRef) / mLenZ;
        mStrainZ = std::log(mLenZRef/mLenZ);
        mVolume  = mLenX * mLenY * mLenZ;
        mMcnAll = 0.0;
        mMcnInn = 0.0;
        mStressMean = 0.0;
        mStressTensor = ZeroMatrix(mDim, mDim);
        double kineticEnergy = 0.0;
        double elasticEnergy = 0.0;
        // Loop over all particles
        for (int i = 0; i < mNumParticles; i++) {
            auto it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            SphericParticle& particle = static_cast<SphericParticle&>(**it);
            // Check if particle is inner
            const int id1 = particle.GetId();
            const double r1 = particle.GetRadius();
            bool isInner = std::find(mInnIds.begin(), mInnIds.end(), id1) != mInnIds.end();
            // Compute kinetic and elastic energies
            const double mass = particle.GetMass();
            const array_1d<double, 3>& vel = particle.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            kineticEnergy += 0.5 * mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            elasticEnergy += particle.GetElasticEnergy();
            // Loop over particle neighbors
            for (int j = 0; j < particle.mNeighbourElements.size(); j++) {
                // Check contact
                const int id2 = particle.mNeighbourElements[j]->GetId();
                const double indent = particle.mBallToBallStoredInfo[id2].indentation;
                if (particle.mBallToBallStoredInfo.find(id2) == particle.mBallToBallStoredInfo.end() || indent <= 0.0) continue;
                // Update coordination numbers
                mMcnAll += 1.0 / mNumParticles;
                if (isInner) mMcnInn += 1.0 / mNumParticlesInn;
                // Evaluate only one contact per pair
                if (id1 >= id2) continue;
                // Force vector
                const double fx = particle.mBallToBallStoredInfo[id2].global_contact_force[0];
                const double fy = particle.mBallToBallStoredInfo[id2].global_contact_force[1];
                const double fz = particle.mBallToBallStoredInfo[id2].global_contact_force[2];
                std::vector<double> force = {fx, fy, fz};
                // Branch vector
                const double r2 = particle.mNeighbourElements[j]->GetRadius();
                const double d  = r1 + r2 - indent;
                const double nx = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][1];
                const double nz = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][2];
                std::vector<double> branch = {d * nx, d * ny, d * nz};
                // Assemble stress tensor
                for (int k = 0; k < mDim; k++) {
                    for (int l = 0; l < mDim; l++) {
                        double Skl = branch[k] * force[l] / mVolume;
                        mStressTensor(k, l) += Skl;
                        if (k == l) mStressMean += Skl / mDim;
                    }
                }
            }
            // Loop over wall neighbors
            for (int j = 0; j < particle.mNeighbourRigidFaces.size(); j++) {
                // Check contact
                const int id2 = particle.mNeighbourRigidFaces[j]->GetId();
                const double indent = particle.mBallToRigidFaceStoredInfo[id2].indentation;
                if (particle.mBallToRigidFaceStoredInfo.find(id2) == particle.mBallToRigidFaceStoredInfo.end() || indent <= 0.0) continue;
                // Force vector
                const double fx = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[0];
                const double fy = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[1];
                const double fz = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[2];
                std::vector<double> force = {fx, fy, fz};
                // Branch vector
                const double d  = r1 - indent;
                const double nx = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][1];
                const double nz = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][2];
                std::vector<double> branch = {d * nx, d * ny, d * nz};
                // Assemble stress tensor
                for (int k = 0; k < mDim; k++) {
                    for (int l = 0; l < mDim; l++) {
                        double Skl = branch[k] * force[l] / mVolume;
                        mStressTensor(k, l) += Skl;
                        if (k == l) mStressMean += Skl / mDim;
                    }
                }
            }
        }
        // Update energy ratio and inertial number
        mEnergyRatio = kineticEnergy / elasticEnergy;
        mInertialNum = mStrainRate * (2.0 * mRadius) / std::sqrt(std::abs(mStressMean / mDensity));
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::UpdateWallVelocities(void) {
        double vX = 0.0;
        double vY = 0.0;
        double vZ = 0.0;
        if (mStage == 1) {
            vX = mStrainRate * mLenX / 2.0;
            vY = mStrainRate * mLenY / 2.0;
            vZ = mStrainRate * mLenZ / 2.0;
        }
        else if (mStage == 3) {
            vZ = mStrainRate * mLenZ / 2.0;
            // Servo-control in X direction
            double Sxx = std::abs(mStressTensor(0,0));
            double stressErrorX = (mTargetStress - Sxx) / mTargetStress;
            double strainRateX = mServoGain * stressErrorX;
            strainRateX = std::clamp(strainRateX, -mServoStrainRateMax, mServoStrainRateMax);
            vX = strainRateX * mLenX / 2.0;
            // Servo-control in Y direction:
            double Syy = std::abs(mStressTensor(1,1)) + (1.0);
            double stressErrorY = (mTargetStress - Syy) / mTargetStress;
            double strainRateY = mServoGain * stressErrorY;
            strainRateY = std::clamp(strainRateY, -mServoStrainRateMax, mServoStrainRateMax);
            vY = strainRateY * mLenY / 2.0;
        }
        SetWallVelocities(vX, vY, vZ);
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::SetWallVelocities(double vX, double vY, double vZ) {
        for (ModelPart::SubModelPartsContainerType::iterator subModelPart = mFemModelPart->SubModelPartsBegin(); subModelPart != mFemModelPart->SubModelPartsEnd(); ++subModelPart) {
            ModelPart &smp = *subModelPart;
            array_1d<double, 3> &linVel = smp[LINEAR_VELOCITY];
            ModelPart::ConditionsContainerType &rConditions = smp.GetCommunicator().LocalMesh().Conditions();
            for (auto iCond = rConditions.begin(); iCond != rConditions.end(); ++iCond) {
                if (&(*iCond) == mWallXMin[0]) {
                    linVel[0] = vX; linVel[1] = 0.0; linVel[2] = 0.0;
                    break;
                }
                if (&(*iCond) == mWallXMax[0]) {
                    linVel[0] = -vX; linVel[1] = 0.0; linVel[2] = 0.0;
                    break;
                }
                if (&(*iCond) == mWallYMin[0]) {
                    linVel[0] = 0.0; linVel[1] = vY; linVel[2] = 0.0;
                    break;
                }
                if (&(*iCond) == mWallYMax[0]) {
                    linVel[0] = 0.0; linVel[1] = -vY; linVel[2] = 0.0;
                    break;
                }
                if (&(*iCond) == mWallZMin[0]) {
                    linVel[0] = 0.0; linVel[1] = 0.0; linVel[2] = vZ;
                    break;
                }
                if (&(*iCond) == mWallZMax[0]) {
                    linVel[0] = 0.0; linVel[1] = 0.0; linVel[2] = -vZ;
                    break;
                }
            }
        }
    }
    //------------------------------------------------------------------------------------------------------------
    bool FCCUtilities::CheckFail(void) {
        if (mMcnInn <= mMcnFail) {
            mFailCount++;
            KRATOS_INFO("FCCUtilities") << "MCN <= " << mMcnFail << " (x" << mFailCount << ")" << std::endl;
            if (mFailCount >= FAIL_COUNT_MAX)
                return true;
        }
        else {
            mFailCount = 0;
        }
        return false;
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::PrintStepInfo(void) {
        const ProcessInfo& processInfo = mDemModelPart->GetProcessInfo();
        std::cout << "-----------------" << std::endl;
        std::cout << std::fixed << std::setprecision(0);
        std::cout << "Stage...........: " << mStage << std::endl;
        std::cout << "Step............: " << processInfo[TIME_STEPS] << std::endl;
        std::cout << std::fixed << std::setprecision(5);
        std::cout << "Time............: " << processInfo[TIME] << std::endl;
        std::cout << std::fixed << std::setprecision(0);
        std::cout << "Stresses (X,Y,Z): "
        << std::abs(mStressMean)
        << " ("
        << std::abs(mStressTensor(0,0)) << ", "
        << std::abs(mStressTensor(1,1)) << ", "
        << std::abs(mStressTensor(2,2))
        << ")" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Coord. number...: " << mMcnInn << std::endl;
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "Inert. number...: " << mInertialNum << std::endl;
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "Energy ratio....: " << mEnergyRatio << std::endl;
        if (mStage == 3) {
            std::cout << std::scientific << std::setprecision(3);
            std::cout << "Strain Z........: " << std::abs(mStrainZ) << std::endl;
        }
        std::cout << "-----------------" << std::endl << std::endl;
    }
    //------------------------------------------------------------------------------------------------------------
    void FCCUtilities::WriteOutputs(void) {
        if (mFileOut.is_open()) {
            double relTime = mDemModelPart->GetProcessInfo()[TIME] - mRefTime;
            mFileOut
            << std::scientific << std::setprecision(3)
            << relTime << ","
            << std::fixed << std::setprecision(0)
            << mStressTensor(0,0) << "," << mStressTensor(1,1) << "," << mStressTensor(2,2) << ","
            << std::scientific << std::setprecision(3)
            << std::abs(mStrainZ) << ","
            << std::fixed << std::setprecision(5)
            << mMcnInn
            << std::endl;
        }
    }
}