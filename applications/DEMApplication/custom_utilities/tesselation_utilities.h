//
// Author:  Rafael Rangel, rrangel@cimne.upc.edu
// Date:    November 2021
//

#ifndef TESSELATION_UTILITIES_H
#define	TESSELATION_UTILITIES_H

// System includes

// Project includes
#include "GeometryFunctions.h"
#include "includes/model_part.h"
#include "custom_elements/thermal_spheric_particle.h"

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

namespace Kratos {

  extern "C" {
    void triangulate(char*, struct triangulateio*, struct triangulateio*, struct triangulateio*);
    void trifree(void*);
  }

  class KRATOS_API(DEM_APPLICATION) TesselationUtilities {

  public:

    KRATOS_CLASS_POINTER_DEFINITION(TesselationUtilities);

    // Constructor / destructor methods
    TesselationUtilities();
    ~TesselationUtilities();

    // Public methods
    void ExecuteInitialize(ModelPart & r_modelpart, bool update_voronoi, bool update_porosity);
    void ExecuteInitializeSolutionStep(ModelPart& r_modelpart);

  protected:
    // Protected methods
    void Triangulation         (ModelPart& r_modelpart);
    void Tetrahedralization    (ModelPart& r_modelpart);
    void UpdateVoronoi2D       (ModelPart& r_modelpart, struct triangulateio& out, struct triangulateio& vorout);
    void UpdateVoronoi3D       (ModelPart& r_modelpart, struct tetgenio& out);
    void UpdatePorosity2D      (ModelPart& r_modelpart, struct triangulateio& out, struct triangulateio& vorout);
    void UpdatePorosity3D      (ModelPart& r_modelpart, struct tetgenio& out);
    void ComputeAlphaRadius2D  (ModelPart& r_modelpart, struct triangulateio& out);
    void ComputeAlphaRadius3D  (ModelPart& r_modelpart, struct tetgenio& out);
    bool AlphaShape2D          (std::vector<double>& coords);
    bool AlphaShape3D          (std::vector<double>& coords);
    void AddParticleArea       (ModelPart& r_modelpart, std::vector<int>& addedParticle, double& particle_area,   int id);
    void AddParticleVolume     (ModelPart& r_modelpart, std::vector<int>& addedParticle, double& particle_volume, int id);
    void ClearTriangle         (struct triangulateio& tr);
    void FreeTriangle          (struct triangulateio& tr);

    // Protected attributes
    bool   mUpdateVoronoi;
    bool   mUpdatePorosiy;
    double mAlphaRadius;
    char*  mSwitches;

  private:
    // Assignment operator
    TesselationUtilities& operator=(TesselationUtilities const& rOther);
  };

} // namespace Kratos

#endif  // TESSELATION_UTILITIES_H