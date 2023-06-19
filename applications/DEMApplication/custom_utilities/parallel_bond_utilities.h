#ifndef PARALLEL_BOND_UTILITIES_H
#define PARALLEL_BOND_UTILITIES_H

#include "includes/define.h"
#include "custom_constitutive/DEM_KDEM_with_damage_parallel_bond_CL.h"
#include "custom_constitutive/DEM_parallel_bond_CL.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos {

    class ParallelBondUtilities {

    public:

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::NodesContainerType NodesContainerType;

        KRATOS_CLASS_POINTER_DEFINITION(ParallelBondUtilities);

        /// Default constructor.

        ParallelBondUtilities() {};

        /// Destructor.

        virtual ~ParallelBondUtilities() {};

        void SetCurrentIndentationAsAReferenceInParallelBonds(ModelPart& rModelPart)
        {
            ////WATCH OUT! This function respects the existing Id's!
            KRATOS_TRY;

            ElementsArrayType& rElements = rModelPart.GetCommunicator().LocalMesh().Elements();

            block_for_each(rElements, [&](ModelPart::ElementType& rElement) {

                SphericContinuumParticle& r_sphere = dynamic_cast<SphericContinuumParticle&>(rElement);

                for (int i = 0; i < (int) r_sphere.mContinuumInitialNeighborsSize; i++) {

                    if (!r_sphere.mNeighbourElements[i]) {
                        continue;
                    }

                    const auto& other_particle = r_sphere.mNeighbourElements[i];
                    array_1d<double, 3> other_to_me_vector;
                    noalias(other_to_me_vector) = r_sphere.GetGeometry()[0].Coordinates() - other_particle->GetGeometry()[0].Coordinates();
                    const double& this_radius = r_sphere.GetRadius();
                    const double& other_radius = other_particle->GetRadius();
                    const double distance = DEM_MODULUS_3(other_to_me_vector);
                    const double radius_sum = this_radius + other_radius;
                    const double& initial_delta = r_sphere.GetInitialDelta(i);
                    const double initial_dist = radius_sum - initial_delta;
                    const double indentation = initial_dist - distance;
                    Kratos::DEMContinuumConstitutiveLaw::Pointer ccl = r_sphere.mContinuumConstitutiveLawArray[i];
                    DEM_KDEM_with_damage_parallel_bond* p_parallel_bond_ccl = dynamic_cast<DEM_KDEM_with_damage_parallel_bond*>(&*ccl);
                    p_parallel_bond_ccl->mInitialIndentationForBondedPart = indentation;
                }
            });

            KRATOS_CATCH("");
        }

        void SetCurrentIndentationAsAReferenceInParallelBondsForPBM(ModelPart& rModelPart)
        {
            ////WATCH OUT! This function respects the existing Id's!
            KRATOS_TRY;

            ElementsArrayType& rElements = rModelPart.GetCommunicator().LocalMesh().Elements();

            block_for_each(rElements, [&](ModelPart::ElementType& rElement) {

                SphericContinuumParticle& r_sphere = dynamic_cast<SphericContinuumParticle&>(rElement);

                for (int i = 0; i < (int) r_sphere.mContinuumInitialNeighborsSize; i++) {

                    if (!r_sphere.mNeighbourElements[i]) {
                        continue;
                    }

                    const auto& other_particle = r_sphere.mNeighbourElements[i];
                    array_1d<double, 3> other_to_me_vector;
                    noalias(other_to_me_vector) = r_sphere.GetGeometry()[0].Coordinates() - other_particle->GetGeometry()[0].Coordinates();
                    const double& this_radius = r_sphere.GetRadius();
                    const double& other_radius = other_particle->GetRadius();
                    const double distance = DEM_MODULUS_3(other_to_me_vector);
                    const double radius_sum = this_radius + other_radius;
                    const double& initial_delta = r_sphere.GetInitialDelta(i);
                    const double initial_dist = radius_sum - initial_delta;
                    const double indentation = initial_dist - distance;
                    Kratos::DEMContinuumConstitutiveLaw::Pointer ccl = r_sphere.mContinuumConstitutiveLawArray[i];
                    DEM_parallel_bond* p_parallel_bond_ccl = dynamic_cast<DEM_parallel_bond*>(&*ccl);
                    p_parallel_bond_ccl->mInitialIndentationForBondedPart = indentation;
                }
            });

            KRATOS_CATCH("");
        }

    }; // Class ParallelBondUtilities

} // namespace Kratos.

#endif // PARALLEL_BOND_UTILITIES_H
