//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
//  RVE boundary nomenclature:
//
//               Y+ (Ymax)
//   Vertex 4 --------------- Vertex 3
//           |               |
//           |               |
// X- (Xmin) |      RVE      | X+ (Xmax)
//           |               |
//           |               | 
//   Vertex 1 --------------- Vertex 2
//               Y- (Ymin)
// 
#pragma once
#include "rve_utilities.h"

namespace Kratos
{
    class KRATOS_API(DEM_APPLICATION) RVEWallBoundary2D : public RVEUtilities
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(RVEWallBoundary2D);

            // Public methods
            RVEWallBoundary2D() {}
            RVEWallBoundary2D(int eval_freq, int write_freq, const std::string& consolidation_criterion, double consolidation_limit, double inner_vol_offset):
            RVEUtilities(eval_freq, write_freq, consolidation_criterion, consolidation_limit, inner_vol_offset) {}
            ~RVEWallBoundary2D() {}

        protected:
            // Protected methods
            void   AssembleWallElementVectors (void) override;
            void   SetVertexCoordinates       (void) override;
            void   SetVertexCoordinatesInner  (void) override;
            void   ProcessGlobalResults       (void) override;
            double ComputeBranchLengthInner   (std::vector<double>& coords1, std::vector<double>& coords2) override;
            double ComputeVolumeParticle      (SphericParticle& particle) override;
            double ComputeVolumeParticleInner (SphericParticle& particle) override;
            double ComputeVolumeRVE           (void) override;
            double ComputeVolumeRVEInner      (void) override;
            double ComputeSurfaceArea         (void) override;
            double ComputeFabricIndex         (Matrix fabric_tensor);
            void   AddContactToRoseDiagram    (std::vector<int>& rose_diagram, std::vector<double>& normal);
            void   EvaluateRoseUniformity     (void) override;
            void   WriteFileHeaders           (void) override;
            void   WriteResultFiles           (void) override;

        private:
            // Private methods
            void   SortVerticesCounterClockwise     (double& x1, double& x2, double& x3, double& x4, double& y1, double& y2, double& y3, double& y4);
            double ComputeAreaQuadrilateral         (double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
            double ComputeAreaTriangle              (double x1, double x2, double x3, double y1, double y2, double y3);
            double ComputeAreaCircleSector          (double r, double d);
            double ComputeDistancePointToSegment    (double x, double y, double xa, double ya, double xb, double yb);
            void   ComputeIntersectionCircleSegment (double r, double x, double y, double xa, double ya, double xb, double yb, double &xc, double& yc);
            bool   PointInsideQuadrilateral         (double x, double y, double vx1, double vy1, double vx2, double vy2, double vx3, double vy3, double vx4, double vy4);
            double ComputeRoseDiagramStdDev         (std::vector<int> rose_diagram);
    };
}
