//
// Authors: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

#if !defined(DEM_INLET_H)
#define DEM_INLET_H

// System includes
#include <string>
#include <iostream>

// External includes
// Project includes
#include "boost/smart_ptr.hpp"

// Project includes
//#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "containers/weak_pointer_vector.h"
#include "includes/constitutive_law.h"
//#include "custom_utilities/create_and_destroy.h"

#include "includes/condition.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_constitutive/DEM_continuum_constitutive_law.h"

namespace Kratos {

    class ParticleCreatorDestructor;
    
    class DEM_Inlet {
        
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;
        
    public:              
        
        /// Constructor:               
        DEM_Inlet(ModelPart& inlet_modelpart);

        /// Destructor.
        virtual ~DEM_Inlet();
        
        void InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const bool using_strategy_for_continuum);
        void DettachElements(ModelPart& r_modelpart, unsigned int& max_Id); 
        void DettachClusters(ModelPart& r_clusters_modelpart, unsigned int& max_Id);
        void CreateElementsFromInletMesh(ModelPart& r_modelpart, ModelPart& r_clusters_modelpart, ParticleCreatorDestructor& creator);

    private:
        Vector mPartialParticleToInsert; //array of doubles, must be resized in the constructor to the number of meshes
        Vector mLastInjectionTimes; //array of doubles, must be resized in the constructor to the number of meshes
        std::vector<int> mTotalNumberOfDetachedParticles;
        ModelPart& mInletModelPart; //The model part used to insert elements
        bool mFirstTime;
        boost::numeric::ublas::vector<bool> mLayerRemoved;
        bool mBallsModelPartHasSphericity;
        bool mBallsModelPartHasRotation;
        std::vector<PropertiesProxy> mFastProperties;
    };
}// namespace Kratos.

#endif // DEM_INLET_H defined
