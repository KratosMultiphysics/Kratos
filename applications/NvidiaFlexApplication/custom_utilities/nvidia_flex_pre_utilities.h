#ifndef NVIDIA_FLEX_PRE_UTILITES_H
#define NVIDIA_FLEX_PRE_UTILITES_H

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/timer.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
//#include "cluster_information.h"
//#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos
{

class NvidiaFlexPreUtilities
    {
    public:

    typedef ModelPart::ElementsContainerType                         ElementsArrayType;
    typedef ModelPart::NodesContainerType::ContainerType             NodesContainerType;
    typedef WeakPointerVector<Element>                               ParticleWeakVectorType;
    typedef WeakPointerVector<Element>::iterator                     ParticleWeakIteratorType;

    KRATOS_CLASS_POINTER_DEFINITION(NvidiaFlexPreUtilities);

    /// Default constructor
    NvidiaFlexPreUtilities() {}

    NvidiaFlexPreUtilities(ModelPart& rModelPart)
    {
        //mInitialCenterOfMassAndMass = CalculateCenterOfMass(rModelPart);
        //mInitialMass                = CalculateTotalMass(rModelPart);
    }

    /// Destructor
    virtual ~NvidiaFlexPreUtilities() {}

    void RemoveSpheresInitiallyIndentedWithFEM(ModelPart& rSpheresModelPart) {

        ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel for
        for (int k = 0; k < (int)pElements.size(); k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Element* p_element = &(*it);
            SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(p_element);

            if (p_sphere->mNeighbourRigidFaces.size()) {
                p_sphere->Set(TO_ERASE);
                p_sphere->GetGeometry()[0].Set(TO_ERASE);
            }
        }
    }

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.

    virtual std::string Info() const
    {
        return "";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    vector<unsigned int>& GetElementPartition() {return (mElementPartition);};

    protected:

        vector<unsigned int> mElementPartition;

    private:

        array_1d<double, 3> mInitialCenterOfMassAndMass;
        double mInitialMass;

        /// Assignment operator
        NvidiaFlexPreUtilities & operator=(NvidiaFlexPreUtilities const& rOther);

    }; // Class NvidiaFlexPreUtilities

} // namespace Kratos

#endif // NVIDIA_FLEX_PRE_UTILITES_H
