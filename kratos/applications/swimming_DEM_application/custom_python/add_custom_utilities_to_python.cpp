/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: G.Casas$
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.2 $
//
//

// System includes 

// External includes 
#include <boost/python.hpp>

// Project includes

#include "add_custom_utilities_to_python.h"
#include "custom_utilities/custom_functions.h"
#include "custom_utilities/binbased_DEM_fluid_coupled_mapping.h"
#include "custom_utilities/volume_averaging_tool.h"
#include "custom_utilities/embedded_volume_tool.h"
#include "custom_utilities/real_functions.h"
#include "custom_utilities/real_field.h"
#include "custom_utilities/real_field_linear_time_dependant_coeff.h"
#include "custom_utilities/time_dependant_porosity_field.h"
#include "custom_utilities/space_time_rule.h"
#include "custom_utilities/space_time_set.h"
#include "custom_utilities/field_utility.h"

namespace Kratos{

namespace Python{

typedef ModelPart::NodesContainerType::iterator PointIterator;
typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;

template<class TDataType>

  void AddNodalSolutionStepVariable(ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
  {
      rModelPart.AddNodalSolutionStepVariable(rThisVariable);
  }

template<int TDim>
void AddDEMCouplingVariable(BinBasedDEMFluidCoupledMapping<TDim,SphericParticle>& rProjectionModule, const VariableData& rThisVariable)
{
    rProjectionModule.AddDEMCouplingVariable(rThisVariable);
}

template<int TDim>
void AddFluidCouplingVariable(BinBasedDEMFluidCoupledMapping<TDim,SphericParticle>& rProjectionModule, const VariableData& rThisVariable)
{
    rProjectionModule.AddFluidCouplingVariable(rThisVariable);
}

void  AddCustomUtilitiesToPython(){
using namespace boost::python;

    class_<RealFunction> ("RealFunction", init<const double, const double>())
        .def("Evaluate", &RealFunction::Evaluate)
        .def("CalculateDerivative", &RealFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &RealFunction::CalculateSecondDerivative)
        ;

    class_<LinearFunction, bases<RealFunction> > ("LinearFunction", init<const double, const double>())
        .def("Evaluate", &LinearFunction::Evaluate)
        .def("CalculateDerivative", &LinearFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &LinearFunction::CalculateSecondDerivative)
        ;

    class_<PowerFunction, bases<RealFunction> > ("PowerFunction", init<const double, const double, const double>())
        .def("Evaluate", &PowerFunction::Evaluate)
        .def("CalculateDerivative", &PowerFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &PowerFunction::CalculateSecondDerivative)
        ;

    class_<AdditionFunction, bases<RealFunction> > ("AdditionFunction", init<const double, RealFunction&, RealFunction&>())
        .def("Evaluate", &AdditionFunction::Evaluate)
        .def("CalculateDerivative", &AdditionFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &AdditionFunction::CalculateSecondDerivative)
        ;

    class_<CompositionFunction, bases<RealFunction> > ("CompositionFunction", init<const double, RealFunction&, RealFunction&>())
        .def("Evaluate", &CompositionFunction::Evaluate)
        .def("CalculateDerivative", &CompositionFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &CompositionFunction::CalculateSecondDerivative)
        ;

    class_<RealField > ("RealField", boost::python::no_init)
        ;

    class_<VectorField<2> > ("VectorField2D", boost::python::no_init)
        ;

    class_<VectorField<3> > ("VectorField3D", boost::python::no_init)
        ;

    class_<LinearRealField, bases<RealField> > ("LinearRealField", init<const double&, const double&, const double&, RealFunction&, RealFunction&, RealFunction&>())
        .def("Evaluate", &LinearRealField::Evaluate)
        .def("CalculateTimeDerivative", &LinearRealField::CalculateTimeDerivative)
        ;

    class_<TimeDependantPorosityField, bases<RealField> > ("TimeDependantPorosityField", init<const double&>())
        .def("Evaluate", &TimeDependantPorosityField::Evaluate)
        .def("CalculateTimeDerivative", &TimeDependantPorosityField::CalculateTimeDerivative)
        .def("CalculateGradient", &TimeDependantPorosityField::CalculateGradient)
        .def("CalculateLaplacian", &TimeDependantPorosityField::CalculateLaplacian)
        ;

    class_<TimeDependantForceField, bases<VectorField<3> > > ("TimeDependantForceField", init<const double&>())
        .def("Evaluate", &TimeDependantForceField::Evaluate)
        .def("GetPorosityField", &TimeDependantForceField::GetPorosityField)
        ;

    class_<SpaceTimeRule> ("SpaceTimeRule", boost::python::no_init)
        ;

    class_<BoundingBoxRule, bases<SpaceTimeRule> > ("BoundingBoxRule", init<>())
        .def("SetTimeBoundingInterval", &BoundingBoxRule::SetTimeBoundingInterval)
        .def("SetXBoundingInterval", &BoundingBoxRule::SetXBoundingInterval)
        .def("SetYBoundingInterval", &BoundingBoxRule::SetYBoundingInterval)
        .def("SetZBoundingInterval", &BoundingBoxRule::SetZBoundingInterval)
        .def("SetSpaceTimeBoundingBox", &BoundingBoxRule::SetSpaceTimeBoundingBox)
        .def("CheckIfRuleIsMet", &BoundingBoxRule::CheckIfRuleIsMet)
        ;

    class_<MoreThanRule, bases<SpaceTimeRule> > ("MoreThanRule", init<RealField::Pointer, const double>())
        .def(init<const double, RealField::Pointer>())
        .def( init<RealField::Pointer, RealField::Pointer>())
        .def("CheckIfRuleIsMet", &MoreThanRule::CheckIfRuleIsMet)
        ;

    class_<SpaceTimeSet> ("SpaceTimeSet", init<>())
        .def("AddAndRule", &SpaceTimeSet::AddAndRule)
        .def("AddOrRule", &SpaceTimeSet::AddOrRule)
        .def("AddAndRules", &SpaceTimeSet::AddAndRules)
        .def("AddOrRules", &SpaceTimeSet::AddOrRules)
        .def("GetLowTime", &SpaceTimeSet::GetLowTime)
        .def("GetHighTime", &SpaceTimeSet::GetHighTime)
        .def("GetLowX", &SpaceTimeSet::GetLowX)
        .def("GetHighX", &SpaceTimeSet::GetHighX)
        .def("GetLowY", &SpaceTimeSet::GetLowY)
        .def("GetHighY", &SpaceTimeSet::GetHighY)
        .def("GetLowZ", &SpaceTimeSet::GetLowZ)
        .def("GetHighZ", &SpaceTimeSet::GetHighZ)
        .def("GetRules", &SpaceTimeSet::GetRules)
        ;

    // next we do what is needed to define the overloaded function 'FieldUtility::ImposeFieldOnNodes'

    typedef double (FieldUtility::*EvaluateDoubleFieldAtPoint)(const double&, const array_1d<double, 3>&, RealField::Pointer);
    typedef array_1d<double, 3> (FieldUtility::*EvaluateVectorFieldAtPoint)(const double&, const array_1d<double, 3>&, VectorField<3>::Pointer);

    EvaluateDoubleFieldAtPoint EvaluateDoubleField = &FieldUtility::EvaluateFieldAtPoint;
    EvaluateVectorFieldAtPoint EvaluateVectorField = &FieldUtility::EvaluateFieldAtPoint;

    // next we do what is needed to define the overloaded function 'FieldUtility::EvaluateFieldAtPoint'

    typedef void (FieldUtility::*ImposeDoubleFieldOnNodes)(Variable<double>&, const double, RealField::Pointer, ModelPart&, const ProcessInfo&, const bool);
    typedef void (FieldUtility::*ImposeVectorFieldOnNodes)(Variable<array_1d<double, 3> >&, const array_1d<double, 3>, VectorField<3>::Pointer, ModelPart&, const ProcessInfo&, const bool);

    ImposeDoubleFieldOnNodes ImposeDoubleField = &FieldUtility::ImposeFieldOnNodes;
    ImposeVectorFieldOnNodes ImposeVectorField = &FieldUtility::ImposeFieldOnNodes;

    class_<FieldUtility> ("FieldUtility", init<SpaceTimeSet::Pointer>())
        .def("MarkNodesInside", &FieldUtility::MarkNodesInside)
        .def("EvaluateFieldAtPoint", EvaluateDoubleField)
        .def("EvaluateFieldAtPoint", EvaluateVectorField)
        .def("ImposeFieldOnNodes", ImposeDoubleField)
        .def("ImposeFieldOnNodes", ImposeVectorField)
        ;

    class_<CustomFunctionsCalculator <2> > ("CustomFunctionsCalculator2D", init<>())
        .def("CalculatePressureGradient", &CustomFunctionsCalculator <2>::CalculatePressureGradient)
        .def("CalculateGradient", &CustomFunctionsCalculator <2>::CalculateGradient)
        .def("CalculateVelocityLaplacian", &CustomFunctionsCalculator <2>::CalculateVelocityLaplacian)
        .def("CalculateVelocityLaplacianRate", &CustomFunctionsCalculator <2>::CalculateVelocityLaplacianRate)
        .def("AssessStationarity", &CustomFunctionsCalculator <2>::AssessStationarity)
        .def("CalculateDomainVolume", &CustomFunctionsCalculator <2>::CalculateDomainVolume)
        .def("CalculateTotalHydrodynamicForceOnParticles", &CustomFunctionsCalculator <2>::CalculateTotalHydrodynamicForceOnParticles)
        .def("CalculateTotalHydrodynamicForceOnFluid", &CustomFunctionsCalculator <2>::CalculateTotalHydrodynamicForceOnFluid)
        .def("CalculateGlobalFluidVolume", &CustomFunctionsCalculator <2>::CalculateGlobalFluidVolume)
        ; 

    //**********************************************************************************************************************************************
    // WARNING!!: function RecoverSuperconvergentGradient uses an algorithm under a GPL 3.0 licence which CANNOT be included in comercial products.
    class_<CustomFunctionsCalculator <3> > ("CustomFunctionsCalculator3D", init<>())
        .def("CalculatePressureGradient", &CustomFunctionsCalculator <3>::CalculatePressureGradient)
        .def("CalculateGradient", &CustomFunctionsCalculator <3>::CalculateGradient)
        .def("CalculateVelocityLaplacian", &CustomFunctionsCalculator <3>::CalculateVelocityLaplacian)
        .def("RecoverSuperconvergentGradient", &CustomFunctionsCalculator <3>::RecoverSuperconvergentGradient)
        .def("CalculateVelocityLaplacianRate", &CustomFunctionsCalculator <2>::CalculateVelocityLaplacianRate)
        .def("AssessStationarity", &CustomFunctionsCalculator <3>::AssessStationarity)
        .def("CalculateDomainVolume", &CustomFunctionsCalculator <3>::CalculateDomainVolume)
        .def("CalculateTotalHydrodynamicForceOnParticles", &CustomFunctionsCalculator <3>::CalculateTotalHydrodynamicForceOnParticles)
        .def("CalculateTotalHydrodynamicForceOnFluid", &CustomFunctionsCalculator <3>::CalculateTotalHydrodynamicForceOnFluid)
        .def("CalculateGlobalFluidVolume", &CustomFunctionsCalculator <3>::CalculateGlobalFluidVolume)
        ;

    class_<BinBasedDEMFluidCoupledMapping <2, SphericParticle> >
            ("BinBasedDEMFluidCoupledMapping2D", init<double, int, int, int>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::InterpolateFromFluidMesh)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddDEMCouplingVariable)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddFluidCouplingVariable)
        ;
    //**********************************************************************************************************************************************

    class_<BinBasedDEMFluidCoupledMapping <2, NanoParticle> >
            ("BinBasedNanoDEMFluidCoupledMapping2D", init<double, int, int, int>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::InterpolateFromFluidMesh)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddDEMCouplingVariable)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddFluidCouplingVariable)
        ;

    class_<BinBasedDEMFluidCoupledMapping <3,SphericParticle> >
            ("BinBasedDEMFluidCoupledMapping3D", init<double, int, int>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::InterpolateFromFluidMesh)
        .def("InterpolateFromNewestFluidMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::InterpolateFromNewestFluidMesh)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddDEMCouplingVariable)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddFluidCouplingVariable)
        ;
    
    class_<BinBasedDEMFluidCoupledMapping <3,NanoParticle> >
            ("BinBasedNanoDEMFluidCoupledMapping3D", init<double, int, int>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::InterpolateFromFluidMesh)
        .def("InterpolateFromNewestFluidMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::InterpolateFromNewestFluidMesh)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddDEMCouplingVariable)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddFluidCouplingVariable)
        ;

    class_<EmbeddedVolumeTool <3> >("EmbeddedVolumeTool", init<>())
        .def("CalculateNegativeDistanceVolume", &EmbeddedVolumeTool <3> ::CalculateNegativeDistanceVolume)
        ;
    }
}  // namespace Python.

} // Namespace Kratos
