//
// Authors:
// Miguel Angel Celigueta maceli@cimne.upc.edu
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#ifndef CREATE_AND_DESTROY_H
#define CREATE_AND_DESTROY_H


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "../DEM_application_variables.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "utilities/quaternion.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/spheric_particle.h"
#include "../custom_utilities/discrete_particle_configure.h"
#include "analytic_tools/analytic_watcher.h"


namespace Kratos {

class KRATOS_API(DEM_APPLICATION) ParticleCreatorDestructor {
friend class ExplicitSolverStrategy;

public:

        static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
        typedef DiscreteParticleConfigure<space_dim>        Configure;
        typedef Configure::ContainerType                    ParticlePointerVector;
        typedef ParticlePointerVector::iterator             ParticlePointerIterator;
        typedef Configure::IteratorType                     ParticleIterator;
        typedef PointerVectorSet<Element, IndexedObject>    ElementsContainerType;
        typedef ModelPart::ElementsContainerType            ElementsArrayType;
        unsigned int mMaxNodeId;

    KRATOS_CLASS_POINTER_DEFINITION(ParticleCreatorDestructor);

    /// Default constructor
    ParticleCreatorDestructor();

    ParticleCreatorDestructor(AnalyticWatcher::Pointer p_watcher);

    /// Destructor
    virtual ~ParticleCreatorDestructor();

    int FindMaxNodeIdInModelPart(ModelPart& r_modelpart);
    void FindAndSaveMaxNodeIdInModelPart(ModelPart& r_modelpart);
    int FindMaxElementIdInModelPart(ModelPart& r_modelpart);
    int FindMaxConditionIdInModelPart(ModelPart& r_modelpart);
    void RenumberElementIdsFromGivenValue(ModelPart& r_modelpart, const int initial_id);

    void NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                           Node < 3 > ::Pointer& pnew_node,
                                           int aId,
                                           Node < 3 > ::Pointer & reference_node,
                                           double radius,
                                           Properties& params,
                                           ModelPart& r_sub_model_part_with_parameters,
                                           bool has_sphericity,
                                           bool has_rotation,
                                           bool initial);

    void NodeForClustersCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                      Node < 3 > ::Pointer& pnew_node,
                                                      int aId,
                                                      Node < 3 > ::Pointer& reference_node,
                                                      Properties& params,
                                                      ModelPart& r_sub_model_part_with_parameters,
                                                      bool has_sphericity,
                                                      bool has_rotation,
                                                      bool initial);

    SphericParticle* ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                          int r_Elem_Id,
                                                          Node < 3 > ::Pointer reference_node,
                                                          Element::Pointer injector_element,
                                                          Properties::Pointer r_params,
                                                          ModelPart& r_sub_model_part_with_parameters,
                                                          const Element& r_reference_element,
                                                          PropertiesProxy* p_fast_properties,
                                                          bool has_sphericity,
                                                          bool has_rotation,
                                                          bool initial,
                                                          ElementsContainerType& array_of_injector_elements);

    SphericParticle* AddInitialDataToNewlyCreatedElementAndNode(ModelPart& r_modelpart,
                                                                Properties::Pointer r_params,
                                                                const double radius, Node<3>::Pointer& pnew_node,
                                                                Element::Pointer& p_particle);


    SphericParticle* CreateSphericParticleRaw(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const Element& r_reference_element);

    SphericParticle* CreateSphericParticleRaw(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const Element& r_reference_element);

    SphericParticle* CreateSphericParticleRaw(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    SphericParticle* CreateSphericParticleRaw(ModelPart& r_modelpart,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    SphericParticle* CreateSphericParticleRaw(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    SphericParticle* CreateSphericParticleRaw(ModelPart& r_modelpart,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    Element::Pointer CreateSphericParticle(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const Element& r_reference_element);

    Element::Pointer CreateSphericParticle(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const Element& r_reference_element);

    Element::Pointer CreateSphericParticle(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    Element::Pointer CreateSphericParticle(ModelPart& r_modelpart,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    Element::Pointer CreateSphericParticle(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);

    Element::Pointer CreateSphericParticle(ModelPart& r_modelpart,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_type);


    Cluster3D* ClusterCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                    ModelPart& r_clusters_modelpart,
                                                    int r_Elem_Id,
                                                    Node < 3 > ::Pointer reference_node,
                                                    Element::Pointer injector_element,
                                                    Properties::Pointer r_params,
                                                    ModelPart& r_sub_model_part_with_parameters,
                                                    const Element& r_reference_element,
                                                    PropertiesProxy* p_fast_properties,
                                                    bool has_sphericity,
                                                    bool has_rotation,
                                                    ElementsContainerType& array_of_injector_elements,
                                                    int& number_of_added_spheres,
                                                    const bool mStrategyForContinuum,
                                                    std::vector<SphericParticle*>& new_component_spheres);


    void NodeCreatorForClusters(ModelPart& r_modelpart,
                                Node < 3 > ::Pointer& pnew_node,
                                int aId,
                                array_1d<double, 3 >& reference_coordinates,
                                double radius,
                                Properties& params);

    void CentroidCreatorForRigidBodyElements(ModelPart& r_modelpart,
                                            Node<3>::Pointer& pnew_node,
                                            int aId,
                                            array_1d<double, 3>& reference_coordinates);

    SphericParticle* SphereCreatorForClusters(ModelPart& r_modelpart,
                                              Node < 3 > ::Pointer& pnew_node,
                                              int r_Elem_Id,
                                              double radius,
                                              array_1d<double, 3 >& reference_coordinates,
                                              double cluster_mass,
                                              Properties::Pointer r_params,
                                              const Element& r_reference_element,
                                              const int cluster_id,
                                              PropertiesProxy* p_fast_properties);

    SphericParticle* SphereCreatorForBreakableClusters(ModelPart& r_modelpart,
                                                       Node < 3 > ::Pointer& pnew_node,
                                                       int r_Elem_Id,
                                                       double radius,
                                                       array_1d<double, 3>& reference_coordinates,
                                                       Properties::Pointer r_params,
                                                       const Element& r_reference_element,
                                                       const int cluster_id,
                                                       PropertiesProxy* p_fast_properties);

    void CalculateSurroundingBoundingBox(ModelPart& r_balls_model_part,
                                         ModelPart& r_clusters_model_part,
                                         ModelPart& r_rigid_faces_model_part,
                                         ModelPart& r_dem_inlet_model_part,
                                         double scale_factor,
                                         bool automatic);

    void DestroyParticles(ModelPart& r_model_part);
    void DestroyParticleElements(ModelPart& r_model_part, Flags flag_for_destruction);
    void DestroyParticles(ModelPart::MeshType& rMesh);
    void DestroyContactElements(ModelPart& r_model_part);
    void MarkInitialNeighboursThatAreBeingRemoved(ModelPart& r_model_part);
    void RemoveUnusedNodesOfTheClustersModelPart(ModelPart& r_clusters_modelpart);
    void MarkDistantParticlesForErasing(ModelPart& r_model_part);
    void MarkParticlesForErasingGivenScalarVariableValue(ModelPart& r_model_part, const Variable<double>& rVariable, double value, double tol);
    void MarkParticlesForErasingGivenVectorVariableModulus(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& rVariable, double value, double tol);
    void MarkParticlesForErasingGivenBoundingBox(ModelPart& r_model_part, array_1d<double, 3> low_point, array_1d<double, 3> high_point);
    void MarkContactElementsForErasing(ModelPart& r_model_part, ModelPart& mcontacts_model_part);
    void DestroyParticlesOutsideBoundingBox(ModelPart& r_model_part);
    void MoveParticlesOutsideBoundingBoxBackInside(ModelPart& r_model_part);
    void DestroyContactElementsOutsideBoundingBox(ModelPart& r_model_part, ModelPart& mcontacts_model_part);
    Element::Pointer GetAnalyticReplacement(const Element& sample_element, Geometry<Node<3> >::PointsArrayType nodelist, Element::Pointer p_elem_to_be_replaced, ModelPart& spheres_model_part);
    static double rand_normal(const double mean, const double stddev, const double max_radius, const double min_radius);
    static double rand_lognormal(const double mean, const double stddev, const double max_radius, const double min_radius);

    array_1d<double, 3> GetHighNode();
    array_1d<double, 3> GetLowNode();
    array_1d<double, 3> GetStrictHighNode();
    array_1d<double, 3> GetStrictLowNode();

    double GetDiameter();
    double GetStrictDiameter();
    void SetHighNode(array_1d<double, 3> node);
    void SetLowNode(array_1d<double, 3> node);
    unsigned int GetCurrentMaxNodeId();
    unsigned int* pGetCurrentMaxNodeId();
    void SetMaxNodeId(unsigned int id);

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

protected:


private:

    void SetDoSearchNeighbourElements(bool true_or_false){mDoSearchNeighbourElements = true_or_false;}

    array_1d<double, 3 > mHighPoint;
    array_1d<double, 3 > mLowPoint;
    array_1d<double, 3 > mStrictHighPoint;
    array_1d<double, 3 > mStrictLowPoint;
    double mDiameter;
    double mStrictDiameter;
    double mScaleFactor;
    int mGreatestParticleId;
    bool mDoSearchNeighbourElements;
    AnalyticWatcher::Pointer mpAnalyticWatcher;
    void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size);
    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable);
    inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable);

    /// Assignment operator.
    ParticleCreatorDestructor & operator=(ParticleCreatorDestructor const& rOther);

}; // Class ParticleCreatorDestructor

} // namespace Kratos

#endif // CREATE_AND_DESTROY_H defined


