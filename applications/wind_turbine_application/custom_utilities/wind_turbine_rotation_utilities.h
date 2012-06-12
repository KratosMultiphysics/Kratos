//
//   Project Name:        Kratos Wind Turbine Application Utility
//   Last Modified by:    $Author: efusto $
//   Date:                $Date: 2012-05-16 11:45:23 $
//   Revision:            $Revision: 0.1 $
//
//


#if !defined(KRATOS_WIND_TURBINE_ROTATION_UTILITIES_INCLUDED)
#define  KRATOS_WIND_TURBINE_ROTATION_UTILITIES_INCLUDED

extern "C" {
	#include "custom_external_libraries/triangle/triangle.h"
}

#include "custom_external_libraries/tetgen1.4.3/tetgen.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "wind_turbine_application.h"

namespace Kratos
{



class WindTurbineRegionMultiplicity
{
public:
	WindTurbineRegionMultiplicity();
	WindTurbineRegionMultiplicity(int);

	bool operator<(WindTurbineRegionMultiplicity) const;
	WindTurbineRegionMultiplicity operator++(int);

	int Region() { return mRegion; }
	int Multiplicity() { return mMultiplicity; }

private:
	int mRegion;
	int mMultiplicity;
};


class WindTurbineRotationUtilities
{
	enum {
            WIND_TURBINE_UNKNOWN_REGION = 0,
            WIND_TURBINE_INNER_INTERF_REGION = 1,
            WIND_TURBINE_OUTER_INTERF_REGION = 2,
            WIND_TURBINE_INNER_REGION = 3,
            WIND_TURBINE_OUTER_REGION = 4,
            WIND_TURBINE_CROWN_REGION = 5,
            WIND_TURBINE_INTERNAL_CYL_BASE = 6,   // cylindrical base region, used for nodes and tetrahedral bases (3D case only)
            WIND_TURBINE_EXTERNAL_CYL_BASE = 7,   //    " ...
            WIND_TURBINE_REGION_NUMBER = 8
	};

        enum {
            WIND_TURBINE_ECHOLEVEL_NONE = 0,
            WIND_TURBINE_ECHOLEVEL_INFO = 1,
            WIND_TURBINE_ECHOLEVEL_DEEPINFO = 2,
            WIND_TURBINE_ECHOLEVEL_WARNING = 3,
            WIND_TURBINE_ECHOLEVEL_SEVERITIES = 4,
            WIND_TURBINE_ECHOLEVEL_ALL = 5
        };

public:
	WindTurbineRotationUtilities(ModelPart&);

        void DoRotationAndRemesh(const int&, const double&, const double&);

        // debugging function
        void DoExtractFaceNodes(ModelPart&, const int&);

private:
	ModelPart& mrGlobalModelPart;
	ModelPart::ElementsContainerType mInnerElems;        // container of the elements inside the inner circumference
	ModelPart::ElementsContainerType mOuterElems;         // container of the elements out of the outer circumference
	ModelPart::ElementsContainerType mInnerInterfElems;  // container of the elements bordering the inner circumference
	ModelPart::ElementsContainerType mOuterInterfElems;  // container of the elements bordering the outer circumference
	ModelPart::ElementsContainerType mCrownElems;        // container of the elements inside the crown

	ModelPart::NodesContainerType mInterfaceNodes;
	ModelPart::NodesContainerType mInnerNodes;
        std::vector<int> mConstrainedBoundaryNodeAuxIndices;
	int mFirstOuterInterfaceNodeOffset;
        int mNumberOfBoundaryFaces;
        std::vector<double> mAngleHistory;
        unsigned int mEchoLevel;

	/// Percolate whole the ModelPart and put the elements and nodes in the proper containers
        void FillElementRegions();
	void FillNodeRegions();
        unsigned int DecideElementRegion(const unsigned int&, std::vector<WindTurbineRegionMultiplicity>&, unsigned int& edgeOppositeVertex, bool& warning) const;

	/// Remeshing primitives for the crown region
	void DestroyCrownElements();
        void DestroyCrownNodes();
	void RegenerateCrownElements2D();
	void RegenerateCrownElements3D();
        void RotateEntities(const double&, const double&, const double&);

        void InitTriangulationDataStructure( triangulateio& );
        void CleanTriangulationDataStructure( triangulateio& );

        double CalculateRotationVelocity(double NewRotAngle);
};

}

#endif
