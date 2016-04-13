//
// Authors: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#ifndef CREATE_AND_DESTROY_H
#define CREATE_AND_DESTROY_H

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"

//Database includes
#include "../custom_utilities/discrete_particle_configure.h"

// Project includes
#include "includes/define.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/spheric_particle.h"
#include "../custom_elements/spheric_continuum_particle.h"
#include "../custom_utilities/GeometryFunctions.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../DEM_application_variables.h"


namespace Kratos {
        
class ParticleCreatorDestructor {
        
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

    /// Destructor
    virtual ~ParticleCreatorDestructor();

    int FindMaxNodeIdInModelPart(ModelPart& r_modelpart);
    void FindAndSaveMaxNodeIdInModelPart(ModelPart& r_modelpart);
    int FindMaxElementIdInModelPart(ModelPart& r_modelpart);
        
    void NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                           Node < 3 > ::Pointer& pnew_node,
                                           int aId,
                                           Node < 3 > ::Pointer & reference_node,
                                           double radius,
                                           Properties& params,
                                           bool has_sphericity,
                                           bool has_rotation,
                                           bool initial);
        
    Kratos::Element* ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node < 3 > ::Pointer reference_node, 
                                              Element::Pointer injector_element,
                                              Properties::Pointer r_params,
                                              const Element& r_reference_element,
                                              PropertiesProxy* p_fast_properties,
                                              bool has_sphericity,
                                              bool has_rotation,
                                              bool initial,
                                              ElementsContainerType& array_of_injector_elements);  
    
    void ClusterCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                            ModelPart& r_clusters_modelpart,
                                            int r_Elem_Id,
                                            Node < 3 > ::Pointer reference_node,
                                            Element::Pointer injector_element,
                                            Properties::Pointer r_params,
                                            const Element& r_reference_element,
                                            PropertiesProxy* p_fast_properties,
                                            bool has_sphericity,
                                            bool has_rotation,
                                            ElementsContainerType& array_of_injector_elements,
                                            int& number_of_added_spheres);
    
    
    void NodeCreatorForClusters(ModelPart& r_modelpart, 
                                Node < 3 > ::Pointer& pnew_node,
                                int aId,
                                array_1d<double, 3 >& reference_coordinates,                                  
                                double radius, 
                                Properties& params);
    
    Kratos::SphericParticle* SphereCreatorForClusters( ModelPart& r_modelpart, 
                                    int r_Elem_Id, 
                                    double radius,
                                    array_1d<double, 3 >& reference_coordinates, 
                                    double cluster_mass,
                                    Properties::Pointer r_params, 
                                    const Element& r_reference_element,
                                    const int cluster_id,
                                    PropertiesProxy* p_fast_properties);
    
    Kratos::SphericParticle* SphereCreatorForBreakableClusters(ModelPart& r_modelpart,
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
                                         double scale_factor,
                                         bool automatic);
    
    void DestroyParticles(ModelPart& r_model_part);
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
    
    static double rand_normal(const double mean, const double stddev, const double max_radius, const double min_radius) {

        if (!stddev) return mean;

        double return_value;

        do {
            double x, y, r;

            do {
                x = 2.0 * rand() / RAND_MAX - 1;
                y = 2.0 * rand() / RAND_MAX - 1;
                r = x * x + y*y;
            } while (r == 0.0 || r > 1.0);

            double d = sqrt(-2.0 * log(r) / r);
            return_value = x * d * stddev + mean;

        } while (return_value < min_radius || return_value > max_radius);

        return return_value;
    }

    static double rand_lognormal(const double mean, const double stddev, const double max_radius, const double min_radius) {
        
        const double normal_mean = log(mean * mean / sqrt(stddev * stddev + mean * mean));
        const double normal_stddev = sqrt(log(1 + stddev * stddev / (mean * mean)));
        double normally_distributed_value = rand_normal(normal_mean, normal_stddev, log(max_radius), log(min_radius));

        return exp(normally_distributed_value);
    }

    static void AddRandomPerpendicularVelocityToGivenVelocity(array_1d<double, 3 >& velocity, const double angle_in_degrees) {
        
        double velocity_modulus = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
        array_1d<double, 3 > unitary_velocity;
        noalias(unitary_velocity) = velocity / velocity_modulus;
        array_1d<double, 3 > normal_1;
        array_1d<double, 3 > normal_2;
        
        if (fabs(unitary_velocity[0])>=0.577) {
            normal_1[0]= - unitary_velocity[1];
            normal_1[1]= unitary_velocity[0];
            normal_1[2]= 0.0;
        }
        else if (fabs(unitary_velocity[1])>=0.577) {
            normal_1[0]= 0.0;
            normal_1[1]= - unitary_velocity[2];
            normal_1[2]= unitary_velocity[1];
        }        
        else {                   
            normal_1[0]= unitary_velocity[2];
            normal_1[1]= 0.0;
            normal_1[2]= - unitary_velocity[0];
        }        
        
        //normalize(normal_1);
        double distance0 = DEM_MODULUS_3(normal_1);
        double inv_distance0 = (distance0 != 0.0) ? 1.0 / distance0 : 0.00;
        normal_1[0] = normal_1[0] * inv_distance0;
        normal_1[1] = normal_1[1] * inv_distance0;
        normal_1[2] = normal_1[2] * inv_distance0;
        
        //CrossProduct(NormalDirection,Vector0,Vector1);
        normal_2[0] = unitary_velocity[1]*normal_1[2] - unitary_velocity[2]*normal_1[1];
        normal_2[1] = unitary_velocity[2]*normal_1[0] - unitary_velocity[0]*normal_1[2];
        normal_2[2] = unitary_velocity[0]*normal_1[1] - unitary_velocity[1]*normal_1[0];     
        
        double angle_in_radians = angle_in_degrees * KRATOS_M_PI / 180;
        double radius = tan(angle_in_radians) * velocity_modulus;
        double radius_square = radius * radius;
        double local_added_velocity_modulus_square = radius_square + 1.0; //just greater than the radius, to get at least one iteration of the while
        array_1d<double, 3> local_added_velocity; local_added_velocity[0] = local_added_velocity[1] = local_added_velocity[2] = 0.0;
        
        while (local_added_velocity_modulus_square > radius_square) {
            //Random in a range: (max - min) * ( (double)rand() / (double)RAND_MAX ) + min
            local_added_velocity[0] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_velocity[1] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_velocity_modulus_square = local_added_velocity[0]*local_added_velocity[0] + local_added_velocity[1]*local_added_velocity[1];
        }
        
        noalias(velocity) += local_added_velocity[0] * normal_1 + local_added_velocity[1] * normal_2;        
    }
    
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

    array_1d<double, 3 > mHighPoint;
    array_1d<double, 3 > mLowPoint;
    array_1d<double, 3 > mStrictHighPoint;
    array_1d<double, 3 > mStrictLowPoint;
    double mDiameter;
    double mStrictDiameter;
    double mScaleFactor;
    int mGreatestParticleId;

    void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size);
    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable);
    inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable);

    /// Assignment operator.
    ParticleCreatorDestructor & operator=(ParticleCreatorDestructor const& rOther);

}; // Class ParticleCreatorDestructor

} // namespace Kratos

#endif // CREATE_AND_DESTROY_H defined


