//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (r.rangel@utwente.nl)
//
//  Utilities for the analysis of FCC failure mechanism under uniaxial compression:
//  * Stage 1: Frictionless isotropic compression until reaching target stress.
//  * Stage 2: Relaxation until reaching static equilibrium.
//  * Stage 3: Servo-controlled axial compression (Z direction) until failure is detected.
//
#pragma once
#include <iomanip>
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_elements/rigid_body_element.h"
namespace Kratos
{
    class KRATOS_API(DEM_APPLICATION) FCCUtilities
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(FCCUtilities);
            // INPUTS
            int mRunFreq;               // Evaluation frequency in time steps
            int mWriteFreq;             // Results printing frequency in time steps
            double mFrictionPP;         // Particle-Particle friction coefficient
            double mFrictionPW;         // Particle-Wall friction coefficient
            bool mFrictionStage1;       // Flag to indicate if friction is considered in stage 1 (isotropic compression)
            double mStrainRate;         // Compression strain rate for stages 1 and 3 (isotropic compression and servo-controlled compression)
            double mInertialNumMax;     // Maximum inertial number for quasi-static compression in stage 3 (servo-controlled compression)
            double mEnergyRatioMax;     // Limit of kinetic/elastic energy ratio for reaching static equilibrium in stage 2 (relaxation)
            double mTargetStress;       // Target stress for stage 1 (isotropic compression)
            double mServoGain;          // Servo-control gain for stage 1 (isotropic compression)
            double mServoStrainRateMax; // Servo-controlled maximum strain rate
            double mMcnFail;            // Mean coordination number threshold for failure detection
            int mFailCountMax;          // Maximum number of consecutive counts with failure condition
            FCCUtilities(
                int runFreq,
                int writeFreq,
                double frictionPP,
                double frictionPW,
                bool frictionStage1,
                double strainRate,
                double inertialNumMax,
                double energyRatioMax,
                double targetStress,
                double servoGain,
                double servoStrainRateMax,
                double mcnFail,
                int failCountMax):
            mRunFreq(runFreq),
            mWriteFreq(writeFreq),
            mFrictionPP(frictionPP),
            mFrictionPW(frictionPW),
            mFrictionStage1(frictionStage1),
            mStrainRate(strainRate),
            mInertialNumMax(inertialNumMax),
            mEnergyRatioMax(energyRatioMax),
            mTargetStress(targetStress),
            mServoGain(servoGain),
            mServoStrainRateMax(servoStrainRateMax),
            mMcnFail(mcnFail),
            mFailCountMax(failCountMax) {}
            // ATTRIBUTES
            int mDim;                        // Problem dimension
            int mStage;                      // 1: isotropic compression, 2: relaxation, 3: servo-controlled compression
            int mNumParticles;               // Number of particles
            int mNumParticlesInn;            // Number of inner particles
            int mNumWallElems;               // Number of walls
            int mFailCount;                  // Number of consecutive counts with failure condition
            double mRadius;                  // Particle radius
            double mDensity;                 // Particle density
            double mMcnAll;                  // Mean coordination number (all particles)
            double mMcnInn;                  // Mean coordination number (inner particles)
            double mInertialNum;             // Inertial number
            double mRefTime;                 // Reference initial time
            double mEnergyRatio;             // Kinetic/Elastic energy ratio
            double mLenX;                    // Length of RVE in X direction
            double mLenY;                    // Length of RVE in Y direction
            double mLenZ;                    // Length of RVE in Z direction
            double mLenZRef;                 // Reference initial length of RVE in Z direction
            double mStrainZ;                 // Axial strain in the Z direction
            double mVolume;                  // Volume of the RVE
            double mStressMean;              // Mean effective stress (all contacts)
            Matrix mStressTensor;            // Cauchy stress tensor (all contacts)
            ModelPart* mDemModelPart;        // Pointer to particle model part
            ModelPart* mFemModelPart;        // Pointer to wall model part
            std::vector<int> mInnIds;        // Vector of internal particle Ids
            std::vector<DEMWall*> mWallXMin; // Vector of wall elements in negative X direction
            std::vector<DEMWall*> mWallXMax; // Vector of wall elements in positive X direction
            std::vector<DEMWall*> mWallYMin; // Vector of wall elements in negative Y direction
            std::vector<DEMWall*> mWallYMax; // Vector of wall elements in positive Y direction
            std::vector<DEMWall*> mWallZMin; // Vector of wall elements in negative Z direction
            std::vector<DEMWall*> mWallZMax; // Vector of wall elements in positive Z direction
            std::ofstream mFileOut;          // File to print results
            // METHODS
            FCCUtilities() {}
            ~FCCUtilities() {}
            void Initialize(ModelPart& dem_model_part, ModelPart& fem_model_part);
            void FinalizeSolutionStep(void);
            void AssembleWallVectors(void);
            void IdentifyInnerParticles(void);
            void UpdateFriction(void);
            void UpdateSystem(void);
            void UpdateWallVelocities(void);
            void SetWallVelocities(double vX, double vY, double vZ);
            bool CheckFail(void);
            void PrintStepInfo(void);
            void WriteOutputs(void);
            bool EqualValues(const double a, const double b) {
                return std::abs(a-b) < 1e-15;
            }
            bool IsTimeToEvaluate(int step) {
                return (step != 0 && mRunFreq != 0 && step % mRunFreq == 0);
            }
            bool IsTimeToPrint(int step) {
                return (step != 0 && mWriteFreq != 0 && step % mWriteFreq == 0);
            }
  };
}