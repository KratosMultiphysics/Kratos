//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#ifndef RVE_UTILITIES_H_INCLUDED
#define	RVE_UTILITIES_H_INCLUDED

// System includes

// External includes
#include "includes/model_part.h"

// Project includes
#include "custom_conditions/dem_wall.h"

namespace Kratos
{
  class KRATOS_API(DEM_APPLICATION) RVEUtilities
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(RVEUtilities);

      // Constructor / destructor methods
      RVEUtilities();
      ~RVEUtilities();

      // Public methods
      void ExecuteInitialize           (ModelPart& rDEMModelPart, ModelPart& rFEMModelPart);
      void ExecuteFinalizeSolutionStep (ModelPart& rDEMModelPart, ModelPart& rFEMModelPart);
      void AssembleWallVectors         (ModelPart& rFEMModelPart);
      void ComputeTotalVolume          (void);

      // Public properties
      bool   mCompressionStage;  // Flag for compression state of RVE generation
      int    mFrequency;         // Frequency (in steps) of RVE homogenization
      double mVolume;            // Total volume
      double mPorosity;          // Porosity (discounting overlaps)
      double mVoidRatio;         // Void ratio (discounting overlaps)

      std::vector<DEMWall*> mWallXMin; // Vector of wall elements in negative X direction
      std::vector<DEMWall*> mWallXMax; // Vector of wall elements in positive X direction
      std::vector<DEMWall*> mWallYMin; // Vector of wall elements in negative Y direction
      std::vector<DEMWall*> mWallYMax; // Vector of wall elements in positive Y direction
      std::vector<DEMWall*> mWallZMin; // Vector of wall elements in negative Z direction
      std::vector<DEMWall*> mWallZMax; // Vector of wall elements in positive Z direction

    private:

      // Assignment operator
      RVEUtilities& operator=(RVEUtilities const& rOther);

  }; // Class RVEUtilities
} // namespace Kratos

#endif // RVE_UTILITIES_H_INCLUDED
