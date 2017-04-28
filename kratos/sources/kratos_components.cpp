//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//












// System includes

// External includes


// Project includes
#include "includes/kratos_components.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/flags.h"
#include "utilities/quaternion.h"
#include "geometries/point.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"
#include "containers/periodic_variables_container.h"


namespace Kratos
{
 
    
    typedef array_1d<double, 3> Array3;

    REGISTER_COMPONENT( Variable<bool> )
    REGISTER_COMPONENT( Variable<int> )
    REGISTER_COMPONENT( Variable<unsigned int> )
    REGISTER_COMPONENT( Variable<double> )
    REGISTER_COMPONENT( Variable<Array3 > )
    REGISTER_COMPONENT( Variable<Quaternion<double> > )
    REGISTER_COMPONENT( Variable<Vector> )
    REGISTER_COMPONENT( Variable<Matrix> )
    REGISTER_COMPONENT( Variable<std::string> )
    REGISTER_COMPONENT( VariableComponent<VectorComponentAdaptor<Array3>> )    
    REGISTER_COMPONENT( Variable<Flags> )
    REGISTER_COMPONENT( Variable<Element> )  
    REGISTER_COMPONENT( Variable<Element::Pointer> ) 
    REGISTER_COMPONENT( Variable<Condition> )
    REGISTER_COMPONENT( Variable<ConstitutiveLaw> )  
    REGISTER_COMPONENT( Variable<ConstitutiveLaw::Pointer> )  
    REGISTER_COMPONENT( Variable<vector<int> > )  
    REGISTER_COMPONENT( Variable<vector<Array3> > )  
    REGISTER_COMPONENT( Variable<WeakPointerVector<Node<3>>> )
    REGISTER_COMPONENT( Variable<WeakPointerVector<Element>> ) 
    REGISTER_COMPONENT( Variable<WeakPointerVector<Condition>> ) 
    REGISTER_COMPONENT( Variable<WeakPointerVector<GeometricalObject>>)
    REGISTER_COMPONENT( Variable<ConvectionDiffusionSettings::Pointer> ) 
    REGISTER_COMPONENT( Variable<RadiationSettings::Pointer> ) 
    REGISTER_COMPONENT( Variable<PeriodicVariablesContainer> )
    REGISTER_COMPONENT( Element )  
    
    REGISTER_COMPONENT( Condition )  
    REGISTER_COMPONENT( Flags )  
    
    REGISTER_COMPONENT( VariableData )  
    

    



}  // namespace Kratos.
