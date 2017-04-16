// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#if !defined(DEM_FORCE_BASED_INLET_H)
#define DEM_FORCE_BASED_INLET_H

// Project includes
#include "inlet.h"

namespace Kratos {

    class DEM_Force_Based_Inlet: public DEM_Inlet {
        
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;
        
    public:              
        
        /// Constructor:               
        DEM_Force_Based_Inlet(ModelPart& inlet_modelpart, array_1d<double, 3> injection_force);

        /// Destructor.
        virtual ~DEM_Force_Based_Inlet(){}

    private:
        array_1d<double, 3> mInjectionForce;
        void FixInjectionConditions(Element* p_element) override;
        void RemoveInjectionConditions(Element &element) override;
        array_1d<double, 3> GetInjectionForce();
    };
}// namespace Kratos.

#endif // DEM_FORCE_BASED_INLET_H defined
