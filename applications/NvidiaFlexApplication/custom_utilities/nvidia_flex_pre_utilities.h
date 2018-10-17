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

    bool CheckIfItsTimeToChangeGravity(ModelPart& rSpheresModelPart,
                                       double& time_of_last_gravity_shift,
                                       const double velocity_threshold_for_gravity_change,
                                       const double min_time_between_changes,
                                       const double max_time_between_changes) {

        const double& current_time = rSpheresModelPart.GetProcessInfo()[TIME];
        
        static double last_time_gravity_changed = 0.0;
        
        if (current_time < last_time_gravity_changed + min_time_between_changes) return false;
        if (current_time > last_time_gravity_changed + max_time_between_changes) {
            last_time_gravity_changed  = current_time;
            return true;
        }
        
        /*if (current_time < time_of_last_gravity_shift + min_time_between_changes) return false;
        if (current_time > time_of_last_gravity_shift + max_time_between_changes) {
            time_of_last_gravity_shift  = current_time;
            return true;
        }*/

        const size_t number_of_nodes = rSpheresModelPart.Nodes().size();
        double max_squared_velocity = 0.0;
        for (size_t i = 0; i < number_of_nodes; i++) {
            
            const auto node_it = rSpheresModelPart.Nodes().begin() + i;
            auto& vel = node_it->FastGetSolutionStepValue(VELOCITY);
            const double node_i_squared_velocity_module = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
            if (node_i_squared_velocity_module > max_squared_velocity) max_squared_velocity = node_i_squared_velocity_module;
        }

        if (max_squared_velocity < velocity_threshold_for_gravity_change * velocity_threshold_for_gravity_change) {
            //time_of_last_gravity_shift  = current_time;
            last_time_gravity_changed  = current_time;
            return true;
        } else return false;
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
