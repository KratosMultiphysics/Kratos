//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes

// External includes
#include "includes/model_part.h"

#ifndef TETLIBRARY
#define TETLIBRARY
#endif

#if !defined(KRATOS_TETGEN_EXTERNAL_H_INCLUDED)
#define KRATOS_TETGEN_EXTERNAL_H_INCLUDED
#include "tetgen.h"
#endif

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) TesselationUtilities3D
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(TesselationUtilities3D);

      // Constructor / Destructor
      TesselationUtilities3D();
      ~TesselationUtilities3D();

      // Public methods
      void ExecuteInitialize             (ModelPart& rModelPart, bool update_voronoi, bool update_porosity);
      void ExecuteInitializeSolutionStep (ModelPart& rModelPart);

    protected:

      // Protected methods
      void Tetrahedralization (ModelPart& rModelPart);
      void UpdateVoronoi      (ModelPart& rModelPart, struct tetgenio& rOut);
      void UpdatePorosity     (ModelPart& rModelPart, struct tetgenio& rOut);
      void ComputeAlphaRadius (ModelPart& rModelPart, struct tetgenio& rOut);
      bool AlphaShape         (std::vector<double>& coords);
      void AddParticleVolume  (ModelPart& rModelPart, std::vector<int>& addedParticle, double& particle_volume, const int id);

      // Protected attributes
      bool        mUpdateVoronoi;
      bool        mUpdatePorosiy;
      double      mAlphaRadius;
      std::string mSwitches;

    private:

      // Assignment operator
      TesselationUtilities3D& operator=(TesselationUtilities3D const& rOther);

  }; // Class TesselationUtilities3D
} // namespace Kratos
