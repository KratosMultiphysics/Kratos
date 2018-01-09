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
    #ifdef SINGLE
        #define REAL float
    #else /* not SINGLE */
        #define REAL double
    #endif /* not SINGLE */

    #include "triangle.h"
}

#include "tetgen.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
//#include "geometries/triangle_2d_3.h"
#include "wind_turbine_application.h"

// Extra includes for debug mesh printouts
#include "includes/gid_io.h"

namespace Kratos
{


class WindTurbineRegionMultiplicity
{
public:
	WindTurbineRegionMultiplicity();
	WindTurbineRegionMultiplicity(int);
        WindTurbineRegionMultiplicity(int, int);

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

        WindTurbineRotationUtilities(ModelPart&,double Zmin,double Zmax);

        void DoRotationAndRemesh(const int&, const double&, const double&);

        void GetModelPartForOutput(ModelPart&);

        // debugging functions
        void DoExtractFaceNodes(ModelPart&, const int&);
        int GetRemeshingRank() { return mRemeshingRank; }
        void PrintBaseBoundaries(GidIO<>& rIO);

private:
	ModelPart& mrGlobalModelPart;
	ModelPart::ElementsContainerType mInnerElems;        // container of the elements inside the inner circumference
	ModelPart::ElementsContainerType mOuterElems;        // container of the elements out of the outer circumference
	ModelPart::ElementsContainerType mInnerInterfElems;  // container of the elements bordering the inner circumference
	ModelPart::ElementsContainerType mOuterInterfElems;  // container of the elements bordering the outer circumference
	ModelPart::ElementsContainerType mCrownElems;        // container of the elements inside the crown

	ModelPart::NodesContainerType mInterfaceNodes;
	ModelPart::NodesContainerType mInnerNodes;
    
        std::vector<int> mConstrainedBoundaryNodeAuxIndices;
	int mFirstOuterInterfaceNodeOffset;
        int mNumberOfBoundaryFaces;
        std::vector<double> mAngleHistory;
        double mOmega;
        unsigned int mNumberOfNodesPerElement;
        unsigned int mEchoLevel;

        // Percolate the whole ModelPart and put the elements and nodes in the proper containers
        void InitializeUtility();
        void FillElementRegions();
	void FillNodeRegions();
        unsigned int DecideElementRegion(const unsigned int&, std::vector<WindTurbineRegionMultiplicity>&, unsigned int& edgeOppositeVertex, bool& warning) const;

        // Remeshing primitives for the crown region
	void DestroyCrownElements();
        void DestroyCrownNodes();
	void RegenerateCrownElements2D();
	void RegenerateCrownElements3D();
        void RotateEntities(const double&, const double&, const double&);

        void InitTriangulationDataStructure( triangulateio& );
        void CleanTriangulationDataStructure( triangulateio& );

        double CalculateRotationVelocity(double NewRotAngle);
        void SetMeshVelocity(double Rx, double Ry, array_1d<double,3>& rMeshVelocity);

        void RemoveLocalNodesWithNoElements();
        void InspectNodeContainerAndLogToFile(ModelPart::NodesContainerType&, const std::string&);

        // Tools to be used in case Tetgen decides to add nodes
        std::vector< Node<3>::Pointer > mCrownNodes;
        unsigned int mLastGlobalNodeId;

        void AddCrownNodesToModelPart( int FirstNewNodeOffset, tetgenio& TetgenOutput);
        void InterpolateNodalDataForNewNodes();
        void Parallel_FindLastNodeId();

        // Additional data and functions for the periodic case
        double mZmin;
        double mZmax;

        std::vector< Node<3>::Pointer > mInterfaceNodesAuxContainer;
        std::vector< Node<3>::Pointer > mMinSideInterfaceNodes;
        std::vector< Node<3>::Pointer > mMaxSideInterfaceNodes;

        std::vector<int> mMinSideBoundaryEdges;
        std::vector<int> mMaxSideBoundaryEdges;

        std::vector<int> mMinSideFacets;
        std::vector<int> mMaxSideFacets;

        static bool AuxIdComp(Node<3>::Pointer pThis, Node<3>::Pointer pThat);
        void FillBaseNodeRegions(std::vector< Node<3>::Pointer >& rAllNodes, std::vector< Node<3>::Pointer >& rNodeContainer, double Zpos);
        void FillBaseRegionEdges(std::vector< Node<3>::Pointer >& rInterfaceNodeList, std::vector<int>& rBoundaryElementList, std::vector<int>& rEdgeList, double Zpos);

        void CreateNewBaseFacets(std::vector< Node<3>::Pointer >& rBaseNodes, std::vector<int>& rEdgeList, std::vector<int>& rFacetList);

        // Nodes on interface containers are sorted, first nodes on the inner side, then these on the outer
        int mFirstOuterNodeOnMinSideOffset;
        int mFirstOuterNodeOnMaxSideOffset;

        // beginning of parallel stuff
        int mThisRank;
        int mRemeshingRank;
        int mNumberOfRanks;

        int mLastKratosGlobalElementId;

        std::vector<int> mRankInnerInterfaceOffsets;
        std::vector<int> mRankInitialInnerInterfNodes;
        std::vector<int> mRankInitialOuterInterfNodes;

        void Parallel_DecideRemeshingProcessor();
        void Parallel_MigrateQuantities();
        void Parallel_FindLastKratosGlobalElementId();
        template <class EntitiesContainer> void Parallel_MigrateEntities(const EntitiesContainer&);
        template <class EntitiesContainer> void Parallel_SerializerSave(EntitiesContainer&, std::string&);
        template <class EntitiesContainer> void Parallel_SerializerLoad(EntitiesContainer&, std::string&);
        void Parallel_SerializerLoadNodes(ModelPart::NodesContainerType&, const char*, const int&);
        template <class EntitiesContainer> void Parallel_SerializerLoad(EntitiesContainer&, const char*, const int&);

        void Parallel_TransferDataToRemeshingProcessor();
        void Parallel_GatherNodesAndFacets();
        static bool PointerIdComp(Node<3>::Pointer const* pThis, Node<3>::Pointer const* pThat);
        // end of parallel stuff
};

}

#endif
