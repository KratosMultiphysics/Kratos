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
#ifndef TETLIBRARY
#define TETLIBRARY
#endif

#if !defined(KRATOS_TETGEN_EXTERNAL_H_INCLUDED)
#define KRATOS_TETGEN_EXTERNAL_H_INCLUDED
#include "tetgen.h"
#endif

#ifndef TRILIBRARY
#define TRILIBRARY
#endif

#include "triangle.h"
#include "includes/model_part.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  extern "C" {
    void triangulate(char*, struct triangulateio*, struct triangulateio*, struct triangulateio*);
    void trifree(void*);
  }

  class KRATOS_API(THERMAL_DEM_APPLICATION) TesselationUtilities2D
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(TesselationUtilities2D);

      // Constructor / Destructor
      TesselationUtilities2D();
      ~TesselationUtilities2D();

      // Public methods
      void ExecuteInitialize             (ModelPart& rModelPart, bool update_voronoi, bool update_porosity);
      void ExecuteInitializeSolutionStep (ModelPart& rModelPart);

    protected:

      // Protected methods
      void Triangulation      (ModelPart& rModelPart);
      void UpdateVoronoi      (ModelPart& rModelPart, struct triangulateio& rOut, struct triangulateio& rVorOut);
      void UpdatePorosity     (ModelPart& rModelPart, struct triangulateio& rOut, struct triangulateio& rVorOut);
      void ComputeAlphaRadius (ModelPart& rModelPart, struct triangulateio& rOut);
      bool AlphaShape         (std::vector<double>& coords);
      void AddParticleArea    (ModelPart& rModelPart, std::vector<int>& addedParticle, double& particle_area,   const int id);
      void ClearTriangle      (struct triangulateio& rTr);
      void FreeTriangle       (struct triangulateio& rTr);

      // Protected attributes
      bool        mUpdateVoronoi;
      bool        mUpdatePorosiy;
      double      mAlphaRadius;
      std::string mSwitches;

    private:

      // Assignment operator
      TesselationUtilities2D& operator=(TesselationUtilities2D const& rOther);

  }; // Class TesselationUtilities2D
} // namespace Kratos
