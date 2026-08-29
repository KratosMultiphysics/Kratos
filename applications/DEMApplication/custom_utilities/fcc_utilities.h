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
            int mFailCountMax;          // Maximum number of consecutive counts with failure condition
            double mFriction;           // Friction coefficient for compression stage
            double mStrainRate;         // Compression strain rate for consolidation and compression stages
            double mInertialNumMax;     // Maximum inertial number for quasi-static compression
            double mTargetStress;       // Target stress for consolidation stage
            double mEnergyRatioMax;     // Limit kinetic/elastics energy ratio for reaching static equilibrium in relaxation stage
            double mServoGain;          // Servo-control gain for compression stage
            double mServoFilterAlpha;   // Servo-control parameter for stress filtering
            double mServoStrainRateMax; // Servo-controlled maximum strain rate
            double mMcnFail;            // Mean coordination number threshold for failure detection
            FCCUtilities(
                double friction,
                int runFreq,
                int writeFreq,
                double strainRate,
                double inertialNumMax,
                double energyRatioMax,
                double TargetStress,
                double servoGain,
                double servoFilterAlpha,
                double ServoStrainRateMax,
                double mcnFail,
                int failCountMax):
            mFriction(friction),
            mRunFreq(runFreq),
            mWriteFreq(writeFreq),
            mStrainRate(strainRate),
            mInertialNumMax(inertialNumMax),
            mEnergyRatioMax(energyRatioMax),
            mTargetStress(TargetStress),
            mServoGain(servoGain),
            mServoFilterAlpha(servoFilterAlpha),
            mServoStrainRateMax(ServoStrainRateMax),
            mMcnFail(mcnFail),
            mFailCountMax(failCountMax) {}
            // ATTRIBUTES
            int mDim;                        // Problem dimension
            int mStage;                      // 1: consolidation, 2: relaxation, 3: compression + servo-control
            int mNumParticles;               // Number of particles
            int mNumParticlesInn;            // Number of inner particles
            int mNumWallElems;               // Number of RVE wall elements
            int mFailCount;                  // Number of consecutive counts with failure condition
            double mRadius;                  // Particle radius
            double mDensity;                 // Particle density
            double mMcnAll;                  // Mean coordination number (all particles)
            double mMcnInn;                  // Mean coordination number (inner particles)
            double mInertialNum;             // Inertial number
            double mRefTime;                 // Reference time
            double mEnergyRatio;             // Kinetic/Elastic energy ratio
            double mStressXXFilt;            // Filtered stress in X direction for servo-control
            double mStressYYFilt;            // Filtered stress in Y direction for servo-control
            double mLenX;                    // Length of RVE in X direction
            double mLenY;                    // Length of RVE in Y direction
            double mLenZ;                    // Length of RVE in Z direction
            double mLenZRef;                 // Reference initial Length of RVE in Z direction
            double mStrainZ;                 // Axial strain in the Z direction
            double mVolume;                  // Volume of the RVE
            double mStressMean;              // Mean effective stress (all contacts)
            Matrix mStressTensor;            // Cauchy stress tensor (all contacts)
            ModelPart* mDemModelPart;        // Pointer to model part of DEM particles in RVE
            ModelPart* mFemModelPart;        // Pointer to model part of FEM walls in RVE boundaries
            std::vector<int> mInnIds;        // Vector of internal particle Ids
            std::vector<DEMWall*> mWallXMin; // Vector of RVE wall elements in negative X direction
            std::vector<DEMWall*> mWallXMax; // Vector of RVE wall elements in positive X direction
            std::vector<DEMWall*> mWallYMin; // Vector of RVE wall elements in negative Y direction
            std::vector<DEMWall*> mWallYMax; // Vector of RVE wall elements in positive Y direction
            std::vector<DEMWall*> mWallZMin; // Vector of RVE wall elements in negative Z direction
            std::vector<DEMWall*> mWallZMax; // Vector of RVE wall elements in positive Z direction
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
                return (step != 0 && mRunFreq != 0 && step % mRunFreq == 0.0);
            }
            bool IsTimeToPrint(int step) {
                return (step != 0 && mWriteFreq != 0 && step % mWriteFreq == 0.0);
            }
  };
}