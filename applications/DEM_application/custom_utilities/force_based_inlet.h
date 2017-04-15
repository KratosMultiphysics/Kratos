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
        DEM_Force_Based_Inlet(ModelPart& inlet_modelpart);

        /// Destructor.
        virtual ~DEM_Force_Based_Inlet(){}

    private:
        void FixInjectionConditions(Element* p_element) override;
    };
}// namespace Kratos.

#endif // DEM_FORCE_BASED_INLET_H defined
