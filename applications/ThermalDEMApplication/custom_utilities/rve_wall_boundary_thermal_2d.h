//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
// 
#pragma once
#include "custom_utilities/rve_wall_boundary_2d.h"
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
    class KRATOS_API(THERMAL_DEM_APPLICATION) RVEWallBoundaryThermal2D : public RVEWallBoundary2D
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(RVEWallBoundaryThermal2D);

            // Public attributes
            Matrix mConductivityTensor;      // Thermal conductivity tensor (all volume)
            Matrix mConductivityTensorInner; // Thermal conductivity tensor (inner volume)

            // Public methods
            RVEWallBoundaryThermal2D() {}
            RVEWallBoundaryThermal2D(int eval_freq, int write_freq, const std::string& consolidation_criterion, double consolidation_limit, double inner_vol_offset):
            RVEWallBoundary2D(eval_freq, write_freq, consolidation_criterion, consolidation_limit, inner_vol_offset) {}
            ~RVEWallBoundaryThermal2D() {}

        protected:
            // Protected methods
            void PreProcessGlobalResults       (void) override;
            void ProcessGlobalResults          (void) override;
            bool Homogenize                    (void) override;
            void WriteFileHeadersTensorResults (void) override;
            void WriteResultFilesTensorResults (void) override;

        private:
            // Private methods
            double ComputePipeAreaInner   (void);
            void   HomogenizeConductivity (void);
    };
}
