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

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "containers/global_pointers_vector.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../applications/DEMApplication/custom_utilities/properties_proxies.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    class ParticleCreatorDestructor;

    class  KRATOS_API(DEM_APPLICATION) DEM_Inlet
    {
    public:

        typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
        typedef GlobalPointersVector<Element> ParticleWeakVectorType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_Inlet);

        /// Constructor:
        DEM_Inlet(ModelPart& inlet_modelpart);

        /// Destructor.
        virtual ~DEM_Inlet(){}

        template<class TDataType> void CheckIfSubModelPartHasVariable(ModelPart& smp, const Variable<TDataType>& rThisVariable) {
            if(!smp.Has(rThisVariable)) KRATOS_ERROR<<"The SubModelPart '"<<smp.Name()<<"' does not have the variable '"<<rThisVariable.Name()<<"'";
        }
        virtual void CheckSubModelPart(ModelPart& smp);
        virtual void InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const bool using_strategy_for_continuum = false);
        virtual void InitializeStep(ModelPart& r_modelpart);
        void DettachElements(ModelPart& r_modelpart, unsigned int& max_Id);
        void DettachClusters(ModelPart& r_clusters_modelpart, unsigned int& max_Id);
        bool OneNeighbourInjectorIsInjecting(const Element::Pointer& element);
        virtual void CreateElementsFromInletMesh(ModelPart& r_modelpart, ModelPart& r_clusters_modelpart, ParticleCreatorDestructor& creator);
        ModelPart& GetInletModelPart();
        void SetNormalizedMaxIndentationForRelease(const double value);
        void SetNormalizedMaxIndentationForNewParticleCreation(const double value);
        int GetPartialNumberOfParticlesInjectedSoFar(const int i);
        int GetTotalNumberOfParticlesInjectedSoFar();
        double GetPartialMassInjectedSoFar(const int i);
        double GetTotalMassInjectedSoFar();
    protected:
        virtual void AddRandomPerpendicularComponentToGivenVector(array_1d<double, 3 >& vector, const double angle_in_degrees);
        virtual void AddRandomPerpendicularComponentToGivenVector2D(array_1d<double, 3 >& vector, const double angle_in_degrees);

    private:
        void UpdateInjectedParticleVelocity(Element &particle, Element &injector_element);
        virtual void FixInjectorConditions(Element* p_element);
        virtual void FixInjectionConditions(Element* p_element, Element* p_injector_element);
        virtual void RemoveInjectionConditions(Element &element, const int dimension);
        virtual void UpdateTotalThroughput(SphericParticle& r_spheric_particle);
        virtual void UpdateTotalThroughput(Cluster3D& r_cluster);
        virtual void UpdatePartialThroughput(SphericParticle& r_spheric_particle, const int i);
        virtual void UpdatePartialThroughput(Cluster3D& r_cluster, const int i);
        double GetInputNumberOfParticles(const ModelPart& mp);
        virtual void CheckDistanceAndSetFlag(ModelPart& r_modelpart);


        std::vector<double> mPartialParticleToInsert; //array of doubles, must be resized in the constructor to the number of meshes
        std::vector<double> mLastInjectionTimes; //array of doubles, must be resized in the constructor to the number of meshes

        bool mFirstInjectionIsDone;
        bool mBallsModelPartHasSphericity;
        bool mBallsModelPartHasRotation;
        bool mStrategyForContinuum;
        int  mTotalNumberOfParticlesInjected;
        std::vector<int> mNumberOfParticlesInjected;
        std::map<int, std::string> mOriginInletSubmodelPartIndexes;
        double mTotalMassInjected;
        std::vector<double> mMassInjected;
        // The following two ratios mark the limit indentation (normalized by the radius) for releasing a particle
        // and for allowing a new one to be injected. admissible_indentation_ratio_for_release should be smaller
        // (more strict), since we want to make sure that the particle is taken far enough to avoid interferences
        // with the newly created ones to come. Otherwise, an initial huge indentation could easily happen.
        double mNormalizedMaxIndentationForRelease;
        double mNormalizedMaxIndentationForNewParticleCreation;
        std::vector<PropertiesProxy> mFastProperties;
        std::vector<bool> mLayerRemoved;
        //std::vector<int> mTotalNumberOfDetachedParticles;
        ModelPart& mInletModelPart; //The model part used to insert elements

        bool mWarningTooSmallInlet;
        bool mWarningTooSmallInletForMassFlow;
        void ThrowWarningTooSmallInlet(const ModelPart& mp);
        void ThrowWarningTooSmallInletForMassFlow(const ModelPart& mp);
    };
}// namespace Kratos.

#endif // DEM_INLET_H defined

