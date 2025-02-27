//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
// Description:
// This class and its subclasses control the process of
// RVE generation, including the consolidation packing (compression) and thermal expansion stages, and
// RVE evaluation, including homogenization of discrete properties into continuum tensors.
//
// References:
// RL Rangel, Continuous-discrete numerical modeling of the thermomechanical behavior of granular media supported by experimental campaign. PhD Thesis, UPC, 2024.
// RL Rangel et al, Multiscale data-driven modeling of the thermomechanical behavior of granular media with thermal expansion effects. Computers and Geotechnics, 176:106789, 2024.
// RL Rangel et al, A continuum-discrete multiscale methodology using machine learning for thermal analysis of granular media. Computers and Geotechnics, 168:106118, 2024.
//
#pragma once
#include <iomanip>
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_elements/rigid_body_element.h"

namespace Kratos
{
    class KRATOS_API(DEM_APPLICATION) RVEUtilities
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(RVEUtilities);

            // Public attributes
            int         mEvalFreq;               // RVE evaluation frequency in time steps (input parameter)
            int         mWriteFreq;              // Results printing frequency in time steps (input parameter)
            std::string mConsolidationCriterion; // Creterion to stop consolidation phase (input parameter)
            double      mConsolidationLimit;     // Limit value used by the consolidation creterion (input parameter)
            double      mInnerVolOffset;         // Offset wrt to boundaries and relative to the average particle radius of the inner RVE volume (input parameter)

            bool   mIsMoving;          // Flag for moving boundaries in current step
            bool   mIsEquilibrium;     // Flag for static equilibrium of particles
            int    mDim;               // Problem dimension (2 or 3)
            int    mEquilibriumSteps;  // Number of consecutive steps in static equilibrium
            int    mNumWallElems;      // Number of RVE wall elements
            int    mNumParticles;      // Number of particles in RVE
            int    mNumParticlesInner; // Number of inner particles in RVE (particles not in contact with walls)
            int    mNumContacts;       // Number of unique contacts in RVE (all contacts)
            int    mNumContactsInner;  // Number of unique inner contacts in RVE (only contacts within inner volume)
            double mAvgCoordNum;       // Average coordination number (number of contacts) of all particles
            double mAvgCoordNumInner;  // Average coordination number (number of contacts) of inner particles
            double mAvgRadius;         // Average radius of particles
            double mVolSolid;          // Volume of solid (particles, including overlaps) in total volume of RVE 
            double mVolSolidInner;     // Volume of solid (particles, including overlaps) in inner volume of RVE
            double mVolTotal;          // Total volume of RVE
            double mVolInner;          // Inner volume of RVE (effective volume without boundary wall effects)
            double mPorosity;          // Porosity of RVE (total volume)
            double mPorosityInner;     // Porosity of RVE (inner volume)
            double mFidx;              // Fabric index (all contacts)
            double mFidxInner;         // Fabric index (inner contacts)
            double mAnisotropy;        // Fabric anisotropy (all particles)
            double mAnisotropyInner;   // Fabric anisotropy (inner contacts)
            double mEffStress;         // Mean effective stress (all volume)
            double mEffStressInner;    // Mean effective stress (inner volume)
            double mDevStress;         // Deviatoric stress (all volume)
            double mDevStressInner;    // Deviatoric stress (inner volume)
            double mWallForces;        // Total force applied by wall boundaries on particles (only normal component)
            double mWallStress;        // Total stress applied by wall boundaries (only normal component)
            double mRoseUnif;          // Std dev of rose diagram values for all particles
            double mRoseUnifInner;     // Std dev of rose diagram values for inner contacts

            ModelPart* mDemModelPart; // Pointer to model part of DEM particles in RVE
            ModelPart* mFemModelPart; // Pointer to model part of FEM walls in RVE boundaries

            std::vector<DEMWall*> mWallXMin; // Vector of RVE wall elements in negative X direction
            std::vector<DEMWall*> mWallXMax; // Vector of RVE wall elements in positive X direction
            std::vector<DEMWall*> mWallYMin; // Vector of RVE wall elements in negative Y direction
            std::vector<DEMWall*> mWallYMax; // Vector of RVE wall elements in positive Y direction
            std::vector<DEMWall*> mWallZMin; // Vector of RVE wall elements in negative Z direction
            std::vector<DEMWall*> mWallZMax; // Vector of RVE wall elements in positive Z direction

            Matrix mVertexCoords;      // Coordinates of RVE vertices: mVertexCoords(i,j) = coordinate i (0=x, 1=y, 2=z) of vertex j (0=Xmin-Ymin-Zmin, 1=Xmax-Ymin-Zmin, 2=Xmax-Ymax-Zmin, 3=Xmin-Ymax-Zmin, 4=Xmin-Ymin-Zmax, 5=Xmax-Ymin-Zmax, 6=Xmax-Ymax-Zmax, 7=Xmin-Ymax-Zmax)
            Matrix mVertexCoordsInner; // Coordinates of RVE inner vertices with an offset wrt original wall vertices

            std::vector<double> mContactChain;  // Vector of contact results: Coordinates (X1,Y1,Z1,X2,Y2,Z2), Forces (Fx,Fy,Fz)... of each contact
            std::vector<int> mRoseDiagram;      // Rose diagram - Polar histogram - of contacts (all particles)
            std::vector<int> mRoseDiagramInner; // Rose diagram - Polar histogram - of contacts (inner particles)
            Matrix mFabricTensor;               // Fabric tensor (all contacts)
            Matrix mFabricTensorInner;          // Fabric tensor (inner contacts)
            Matrix mStressTensor;               // Cauchy stress tensor (all volume)
            Matrix mStressTensorInner;          // Cauchy stress tensor (inner volume)

            std::ofstream mFileGlobalResults;   // File to print global RVE results (RVE dimensions, totals, averages, porosity...)
            std::ofstream mFileParticleResults; // File to print results for each particle (radius, coordinates, #contacts...)
            std::ofstream mFileContactResults;  // File to print results for each contact (forces...)
            std::ofstream mFileTensorResults;   // File to print tensorial related results (fabric, stress...)
            std::ofstream mFileRoseDiagram;     // File to print rose diagram related results

            // Public methods
            RVEUtilities() {}
            RVEUtilities(int eval_freq, int write_freq, const std::string& consolidation_criterion, double consolidation_limit, double inner_vol_offset):
            mEvalFreq(eval_freq), mWriteFreq(write_freq), mConsolidationCriterion(consolidation_criterion), mConsolidationLimit(consolidation_limit), mInnerVolOffset(inner_vol_offset) {}
            ~RVEUtilities() {}
            
            void Initialize           (ModelPart& dem_model_part, ModelPart& fem_model_part);
            void FinalizeSolutionStep (void);
            void Finalize             (void);

        protected:
            // Protected methods
            bool EqualValues(const double a, const double b) {return std::abs(a-b) < std::numeric_limits<double>::epsilon();}
            bool IsTimeToEvaluateRVE  (int time_step) {return (time_step != 0 && mEvalFreq  != 0 && time_step % mEvalFreq  == 0.0);}
            bool IsTimeToPrintResults (int time_step) {return (time_step != 0 && mWriteFreq != 0 && time_step % mWriteFreq == 0.0);}

            virtual void   InitializeVariables             (ModelPart& dem_model_part, ModelPart& fem_model_part);
            virtual double ComputeAverageRadius            (void);
            virtual void   AssembleWallElementVectors      (void) {}
            virtual void   SetVertexCoordinates            (void) {}
            virtual void   SetVertexCoordinatesInner       (void) {}
            virtual void   PreProcessGlobalResults         (void);
            virtual void   ProcessGlobalResults            (void) {}
            virtual void   PostProcessGlobalResults        (void);
            virtual double ComputeBranchLengthInner        (std::vector<double>& coords1, std::vector<double>& coords2) {return 0.0;}
            virtual double ComputeVolumeParticle           (SphericParticle& particle) {return 0.0;}
            virtual double ComputeVolumeParticleInner      (SphericParticle& particle) {return 0.0;}
            virtual double ComputeVolumeRVE                (void) {return 0.0;}
            virtual double ComputeVolumeRVEInner           (void) {return 0.0;}
            virtual double ComputePorosity                 (void);
            virtual double ComputePorosityInner            (void);
            virtual double ComputeSurfaceArea              (void) {return 0.0;}
            virtual double ComputeFabricIndex              (Matrix fabric_tensor) {return 0.0;}
            virtual bool   Homogenize                      (void);
            virtual void   HomogenizeFabric                (void);
            virtual void   HomogenizeStress                (void);
            virtual void   AddContactToRoseDiagram         (std::vector<int>& rose_diagram, std::vector<double>& normal) {}
            virtual void   EvaluateRoseUniformity          (void) {}
            virtual void   CheckEquilibrium                (double curr_val, double prev_val, double tol, int max_eq_steps);
            virtual void   StopCompress                    (void);
            virtual void   StopBoundaryMotion              (void);
            virtual void   ReadOldForces                   (void);
            virtual void   OpenResultFiles                 (void);
            virtual void   WriteFileHeaders                (void);
            virtual void   WriteFileHeadersGlobalResults   (void) {}
            virtual void   WriteFileHeadersParticleResults (void) {}
            virtual void   WriteFileHeadersContactResults  (void) {}
            virtual void   WriteFileHeadersTensorResults   (void) {}
            virtual void   WriteFileHeadersRoseDiagram     (void) {}
            virtual void   WriteResultFiles                (void);
            virtual void   WriteResultFilesGlobalResults   (void) {}
            virtual void   WriteResultFilesParticleResults (void) {}
            virtual void   WriteResultFilesContactResults  (void) {}
            virtual void   WriteResultFilesTensorResults   (void) {}
            virtual void   WriteResultFilesRoseDiagram     (void) {}
            virtual void   CloseResultFiles                (void);
  };
}
