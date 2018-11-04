// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#if !defined(SWIMMING_DEM_BENTONITE_FORCE_BASED_INLET_H)
#define SWIMMING_DEM_BENTONITE_FORCE_BASED_INLET_H

// Project includes
#include "../../../DEM_application/custom_utilities/force_based_inlet.h"
#include "../../../DEM_application/custom_elements/nanoparticle.h"

namespace Kratos {
    class KRATOS_API(SWIMMING_DEM_APPLICATION) Bentonite_Force_Based_Inlet: public DEM_Force_Based_Inlet
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(Bentonite_Force_Based_Inlet);

        typedef NanoParticle* NanoParticlePointerType;
        typedef DEM_Force_Based_Inlet BaseClass;
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;                       
        typedef ModelPart::ElementsContainerType::iterator ElementIteratorType;

        /// Constructor:               
        Bentonite_Force_Based_Inlet(ModelPart& inlet_modelpart, array_1d<double, 3> injection_force);

        /// Destructor.
        virtual ~Bentonite_Force_Based_Inlet(){}

        void InitializeStep(ModelPart& r_receiver_model_part) override;

    private:
        double mCationConcentration;
        array_1d<double, 3> mInjectionForce;
        void FixInjectionConditions(Element* p_element, Element* p_injector_element) override;
        void FixInjectorConditions(Element* p_element) override;
        array_1d<double, 3> GetInjectionForce(Element* p_element) override;
        void UpdateInjectionForce(Element *p_element);
    };
}// namespace Kratos.

#endif // SWIMMING_DEM_BENTONITE_FORCE_BASED_INLET_H defined
