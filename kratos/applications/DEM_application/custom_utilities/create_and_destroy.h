//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_CREATE_AND_DESTROY )
#define  KRATOS_CREATE_AND_DESTROY

// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"

// Project includes
#include "includes/define.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/spheric_particle.h"
#include "../custom_elements/spheric_continuum_particle.h"
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
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
        unsigned int mMaxNodeId;
		
    KRATOS_CLASS_POINTER_DEFINITION(ParticleCreatorDestructor);

    /// Default constructor 
    ParticleCreatorDestructor();

    //Particle_Creator_Destructor() {};

    /// Destructor

    virtual ~ParticleCreatorDestructor();

    int FindMaxNodeIdInModelPart(ModelPart& r_modelpart);
    
    void NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                           Node < 3 > ::Pointer& pnew_node,
                                           int aId,
                                           Node < 3 > ::Pointer & reference_node,
                                           double radius,
                                           Properties& params,
                                           bool has_sphericity,
                                           bool has_rotation,
                                           bool initial);
        
    void ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart,
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
                                            bool initial,
                                            ElementsContainerType& array_of_injector_elements,
                                            int& number_of_added_spheres);
    
    
    void NodeCreatorForClusters(ModelPart& r_modelpart, 
                                Node < 3 > ::Pointer& pnew_node,
                                int aId,
                                array_1d<double, 3 >& reference_coordinates,                                  
                                double radius, 
                                Properties& params);
    
    Kratos::SphericParticle* ElementCreatorForClusters( ModelPart& r_modelpart, 
                                    int r_Elem_Id, 
                                    double radius,
                                    array_1d<double, 3 >& reference_coordinates, 
                                    double cluster_mass,
                                    Properties::Pointer r_params, 
                                    const Element& r_reference_element,
                                    const int cluster_id);

    void CalculateSurroundingBoundingBox(ModelPart& r_balls_model_part,
                                         ModelPart& r_clusters_model_part,
                                         ModelPart& r_rigid_faces_model_part,
                                         double scale_factor,
                                         bool automatic);

    void DestroyParticles(ModelPart& r_model_part);

    void DestroyContactElements(ModelPart& r_model_part);

    void MarkInitialNeighboursThatAreBeingRemoved(ModelPart& r_model_part);

    void MarkDistantParticlesForErasing(ModelPart& r_model_part);

    void MarkParticlesForErasingGivenScalarVariableValue(ModelPart& r_model_part, const Variable<double>& rVariable, double value, double tol);

    void MarkParticlesForErasingGivenVectorVariableModulus(ModelPart& r_model_part,
                                                           const Variable<array_1d<double, 3> >& rVariable,
                                                           double value,
                                                           double tol);
        
    void MarkParticlesForErasingGivenBoundingBox(ModelPart& r_model_part,
                                                 array_1d<double, 3> low_point,
                                                 array_1d<double, 3> high_point);

    void MarkContactElementsForErasing(ModelPart& r_model_part, ModelPart& mcontacts_model_part);

    void DestroyParticlesOutsideBoundingBox(ModelPart& r_model_part);

    void DestroyContinuumParticlesOutsideBoundingBox(ModelPart& r_model_part);

    void DestroyContactElementsOutsideBoundingBox(ModelPart& r_model_part, ModelPart& mcontacts_model_part);

    array_1d<double, 3> GetHighNode();

    array_1d<double, 3> GetLowNode();

    array_1d<double, 3> GetStrictHighNode();

    array_1d<double, 3> GetStrictLowNode();

    double GetDiameter();

    double GetStrictDiameter();

    void SetHighNode(array_1d<double, 3> node);

    void SetLowNode(array_1d<double, 3> node);

    unsigned int GetCurrentMaxNodeId();

    void SetMaxNodeId(unsigned int id);

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.

    virtual std::string Info() const;

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const;

protected:

    ///@}
    ///@name Protected member rVariables
    ///@{ template<class T, std::size_t dim>

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

/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}
///@}


} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined


