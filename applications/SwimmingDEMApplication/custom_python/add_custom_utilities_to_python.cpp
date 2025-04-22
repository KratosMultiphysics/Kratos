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

// Project includes

#include "add_custom_utilities_to_python.h"
#include "custom_utilities/custom_functions.h"
#include "custom_utilities/binbased_DEM_fluid_coupled_mapping.h"
#include "custom_utilities/volume_averaging_tool.h"
#include "custom_utilities/embedded_volume_tool.h"
#include "custom_utilities/fields/real_functions.h"
#include "custom_utilities/fields/real_field.h"
#include "custom_utilities/fields/real_field_linear_time_dependant_coeff.h"
#include "custom_utilities/fields/time_dependant_porosity_field.h"
#include "custom_utilities/fields/sets/space_time_rule.h"
#include "custom_utilities/fields/sets/space_time_set.h"
#include "custom_utilities/fields/field_utility.h"
#include "custom_utilities/fields/fluid_field_utility.h"
#include "custom_utilities/fields/vector_field.h"
#include "custom_utilities/fields/velocity_field.h"
#include "custom_utilities/fields/constant_velocity_field.h"
#include "custom_utilities/fields/cellular_flow_field.h"
#include "custom_utilities/fields/ethier_flow_field.h"
#include "custom_utilities/fields/product_of_sines_field.h"
#include "custom_utilities/fields/pouliot_flow_field.h"
#include "custom_utilities/fields/pouliot_flow_field_2D.h"
#include "custom_utilities/fields/shear_flow_1D_with_exponential_viscosity_field.h"
#include "custom_utilities/basset_force_tools.h"
#include "custom_utilities/statistics/sampling_tool.h"
#include "custom_utilities/derivative_recovery_meshing_tools.h"
#include "custom_utilities/inlets/bentonite_force_based_inlet.h"
#include "custom_utilities/swimming_dem_in_pfem_utils.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/mesh_rotation_utility.h"
#include "custom_utilities/renumbering_nodes_utility.h"
#include "custom_utilities/stationarity_check.h"
#include "custom_utilities/error_norm_calculator_utility.h"
#include "custom_utilities/error_norm_torus.h"
#include "custom_utilities/error_norm_ethier_field.h"
#include "custom_utilities/averaging_variables_utility.h"

namespace Kratos{

namespace Python{

typedef ModelPart::NodesContainerType::iterator PointIterator;
typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;

class VariableChecker{
public:
    VariableChecker(){}
    virtual ~VariableChecker(){}

    template<class TDataType>
    bool ModelPartHasNodalVariableOrNot(ModelPart& r_model_part, const Variable<TDataType>& rThisVariable)
    {
        return (r_model_part.GetNodalSolutionStepVariablesList()).Has(rThisVariable);
    }
};

template<class TDataType>
bool ModelPartHasNodalVariableOrNot(VariableChecker& rChecker, ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
{
    return rChecker.ModelPartHasNodalVariableOrNot(rModelPart, rThisVariable);
}

namespace py = pybind11;

// template<std::size_t TDim, class TParticle>
// void DeclareBinbasedCoupler(py::module &m, std::string sphere_type)
// {
// 	using Class = BinBasedDEMFluidCoupledMapping<TDim,TParticle>;
// 	std::string pyclass_name = std::string("BinBased")
// 							 + sphere_type
// 							 + std::string("FluidCoupledMapping")
// 							 + std::to_string(TDim)
// 							 + std::string("D");
// 	py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
// 	.def(py::init<Parameters&>())
// 	.def("InterpolateVelocityOnAuxVelocity", &Class::InterpolateVelocityOnAuxVelocity)
// 	.def("ImposeVelocityOnDEMFromFieldToAuxVelocity", &Class::ImposeVelocityOnDEMFromFieldToAuxVelocity)
// 	.def("InterpolateFromFluidMesh", &Class::InterpolateFromFluidMesh)
// 	.def("ImposeFlowOnDEMFromField", &Class::ImposeFlowOnDEMFromField)
// 	.def("InterpolateFromDEMMesh", &Class::InterpolateFromDEMMesh)
// 	.def("HomogenizeFromDEMMesh", &Class::HomogenizeFromDEMMesh)
// 	.def("ComputePostProcessResults", &Class::ComputePostProcessResults)
// 	.def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping<3,SphericParticle>::AddDEMCouplingVariable<double>)
// 	.def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddDEMCouplingVariable<array_1d<double, 3>>)
// 	.def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddFluidCouplingVariable<double>)
// 	.def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddFluidCouplingVariable<array_1d<double, 3>>)
// 	.def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddDEMVariablesToImpose<double>)
// 	.def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddDEMVariablesToImpose<array_1d<double, 3>>)
// 	.def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddFluidVariableToBeTimeFiltered<double>)
// 	.def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMapping <3,SphericParticle>::AddFluidVariableToBeTimeFiltered<array_1d<double, 3>>);
// }


void  AddCustomUtilitiesToPython(pybind11::module& m){

	// DeclareBinbasedCoupler<3, SphericParticle>(m, "DEM");

    py::class_<VariableChecker> (m, "VariableChecker").def(py::init<>())
        .def("ModelPartHasNodalVariableOrNot", ModelPartHasNodalVariableOrNot<double>)
        .def("ModelPartHasNodalVariableOrNot", ModelPartHasNodalVariableOrNot<array_1d<double, 3> >)
        ;

    py::class_<RealFunction> (m, "RealFunction").def(py::init<const double, const double>())
        .def("Evaluate", &RealFunction::Evaluate)
        .def("CalculateDerivative", &RealFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &RealFunction::CalculateSecondDerivative)
        ;

    py::class_<LinearFunction, RealFunction> (m, "LinearFunction")
        .def(py::init<const double, const double>())
        .def("Evaluate", &LinearFunction::Evaluate)
        .def("CalculateDerivative", &LinearFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &LinearFunction::CalculateSecondDerivative)
        ;

    py::class_<PowerFunction, RealFunction> (m, "PowerFunction")
        .def(py::init<const double, const double, const double>())
        .def("Evaluate", &PowerFunction::Evaluate)
        .def("CalculateDerivative", &PowerFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &PowerFunction::CalculateSecondDerivative)
        ;

    py::class_<AdditionFunction, RealFunction> (m, "AdditionFunction")
        .def(py::init<const double, RealFunction&, RealFunction&>())
        .def("Evaluate", &AdditionFunction::Evaluate)
        .def("CalculateDerivative", &AdditionFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &AdditionFunction::CalculateSecondDerivative)
        ;

    py::class_<CompositionFunction, RealFunction> (m, "CompositionFunction")
        .def(py::init<const double, RealFunction&, RealFunction&>())
        .def("Evaluate", &CompositionFunction::Evaluate)
        .def("CalculateDerivative", &CompositionFunction::CalculateDerivative)
        .def("CalculateSecondDerivative", &CompositionFunction::CalculateSecondDerivative)
        ;

    py::class_<RealField, RealField::Pointer> (m, "RealField").def(py::init<>())
        ;

    py::class_<VectorField<2>, VectorField<2>::Pointer> (m, "VectorField2D").def(py::init<>())
        ;

    py::class_<VectorField<3>, VectorField<3>::Pointer> (m, "VectorField3D").def(py::init<>())
        ;

    //typedef void (VelocityField::*Evaluate)(const double, const DenseVector<double>&, DenseVector<double>&, const int);
    //Evaluate EvaluateVector = &VelocityField::Evaluate;

    typedef void (VelocityField::*CalculateTimeDerivative)(const double, const DenseVector<double>&, DenseVector<double>&, const int);
    CalculateTimeDerivative CalculateTimeDerivativeVector = &VelocityField::CalculateTimeDerivative;

    typedef void (VelocityField::*CalculateGradient)(const double,
                                                     const array_1d<double, 3>&,
                                                     DenseVector< double>&,
                                                     DenseVector< double>&,
                                                     DenseVector< double>&,
                                                     const int);
    CalculateGradient CalculateGradientVector = &VelocityField::CalculateGradient;

    typedef double (VelocityField::*CalculateDivergence)(const double, const DenseVector<double>&, const int);
    CalculateDivergence CalculateDivergenceVector = &VelocityField::CalculateDivergence;

    typedef void (VelocityField::*CalculateRotational)(const double, const DenseVector<double>&, DenseVector<double>&, const int);
    CalculateRotational CalculateRotationalVector = &VelocityField::CalculateRotational;

    typedef void (VelocityField::*CalculateLaplacian)(const double, const DenseVector<double>&, DenseVector<double>&, const int);
    CalculateLaplacian CalculateLaplacianVector = &VelocityField::CalculateLaplacian;

    typedef void (VelocityField::*CalculateMaterialAcceleration)(const double, const DenseVector<double>&, DenseVector<double>&, const int);
    CalculateMaterialAcceleration CalculateMaterialAccelerationVector = &VelocityField::CalculateMaterialAcceleration;

    typedef void (VelocityField::*Evaluate)(const double, const DenseVector<double>&, DenseVector<double>&, const int);
    Evaluate EvaluateFieldVector = &VelocityField::Evaluate;

    py::class_<VelocityField, VelocityField::Pointer, VectorField<3>> (m, "VelocityField")
        .def(py::init<>())
        .def("CalculateTimeDerivative", CalculateTimeDerivativeVector)
        .def("CalculateGradient", CalculateGradientVector)
        .def("CalculateDivergence", CalculateDivergenceVector)
        .def("CalculateRotational", CalculateRotationalVector)
        .def("CalculateLaplacian", CalculateLaplacianVector)
        .def("CalculateMaterialAcceleration", CalculateMaterialAccelerationVector)
        .def("Evaluate", EvaluateFieldVector)
        ;

    py::class_<ConstantVelocityField, ConstantVelocityField::Pointer, VelocityField > (m, "ConstantVelocityField")
        .def(py::init<const double, const double, const double>())
        ;

    py::class_<ShearFlow1DWithExponentialViscosityField, ShearFlow1DWithExponentialViscosityField::Pointer, VelocityField > (m, "ShearFlow1DWithExponentialViscosityField")
        .def(py::init<const double, const double, const double>())
        .def("SetRimZoneThickness", &ShearFlow1DWithExponentialViscosityField::SetRimZoneThickness)
        ;

    py::class_<CellularFlowField, CellularFlowField::Pointer, VelocityField > (m, "CellularFlowField")
        .def(py::init<const double, const double, const double, const double>())
        ;

    py::class_<EthierFlowField, EthierFlowField::Pointer, VelocityField > (m, "EthierFlowField")
        .def(py::init<const double, const double>())
        ;

    py::class_<PouliotFlowField, PouliotFlowField::Pointer, VelocityField > (m, "PouliotFlowField")
        .def(py::init<>())
        ;

    py::class_<PouliotFlowField2D, PouliotFlowField2D::Pointer, VelocityField > (m, "PouliotFlowField2D")
        .def(py::init<>())
        ;

    py::class_<LinearRealField, LinearRealField::Pointer, RealField > (m, "LinearRealField")
        .def(py::init<const double&, const double&, const double&, RealFunction&, RealFunction&, RealFunction&>())
        .def("Evaluate", &LinearRealField::Evaluate)
        .def("CalculateTimeDerivative", &LinearRealField::CalculateTimeDerivative)
        ;

    py::class_<TimeDependantPorosityField, TimeDependantPorosityField::Pointer, RealField > (m, "TimeDependantPorosityField")
        .def(py::init<const double&>())
        .def("Evaluate", &TimeDependantPorosityField::Evaluate)
        .def("CalculateTimeDerivative", &TimeDependantPorosityField::CalculateTimeDerivative)
        .def("CalculateGradient", &TimeDependantPorosityField::CalculateGradient)
        .def("CalculateLaplacian", &TimeDependantPorosityField::CalculateLaplacian)
        ;

    py::class_<TimeDependantForceField, TimeDependantForceField::Pointer, VectorField<3>> (m, "TimeDependantForceField")
        .def(py::init<const double&>())
        .def("Evaluate", &TimeDependantForceField::Evaluate)
        .def("GetPorosityField", &TimeDependantForceField::GetPorosityField)
        ;

    py::class_<SpaceTimeRule> (m, "SpaceTimeRule")
        .def(py::init<>())
        ;

    py::class_<BoundingBoxRule, SpaceTimeRule > (m, "BoundingBoxRule")
        .def(py::init<>())
        .def(py::init<const double, const double, const double, const double, const double, const double, const double, const double>())
        .def("SetTimeBoundingInterval", &BoundingBoxRule::SetTimeBoundingInterval)
        .def("SetXBoundingInterval", &BoundingBoxRule::SetXBoundingInterval)
        .def("SetYBoundingInterval", &BoundingBoxRule::SetYBoundingInterval)
        .def("SetZBoundingInterval", &BoundingBoxRule::SetZBoundingInterval)
        .def("SetSpaceTimeBoundingBox", &BoundingBoxRule::SetSpaceTimeBoundingBox)
        .def("CheckIfRuleIsMet", &BoundingBoxRule::CheckIfRuleIsMet)
        .def("Info", &BoundingBoxRule::Info)
        ;

    py::class_<MoreThanRule, SpaceTimeRule > (m, "MoreThanRule")
        .def(py::init<const double, RealField::Pointer>())
        .def(py::init<RealField::Pointer, RealField::Pointer>())
        .def("CheckIfRuleIsMet", &MoreThanRule::CheckIfRuleIsMet)
        ;

    py::class_<SpaceTimeSet, SpaceTimeSet::Pointer> (m, "SpaceTimeSet")
        .def(py::init<>())
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


    py::class_<AveragingVariablesUtility, AveragingVariablesUtility::Pointer>(m, "AveragingVariablesUtility")
        .def(py::init<>())
        .def("AverageVariables", &AveragingVariablesUtility::AverageVariables)
        ;

    // next we do what is needed to define the overloaded function 'FieldUtility::ImposeFieldOnNodes'

    typedef double (FieldUtility::*EvaluateDoubleFieldAtPoint)(const double&, const array_1d<double, 3>&, RealField::Pointer);
    typedef array_1d<double, 3> (FieldUtility::*EvaluateVectorFieldAtPoint)(const double&, const array_1d<double, 3>&, VectorField<3>::Pointer);

    EvaluateDoubleFieldAtPoint EvaluateDoubleField = &FieldUtility::EvaluateFieldAtPoint;
    EvaluateVectorFieldAtPoint EvaluateVectorField = &FieldUtility::EvaluateFieldAtPoint;

    typedef void (FieldUtility::*ImposeDoubleFieldOnNodes)(Variable<double>&, const double, RealField::Pointer, ModelPart&, const ProcessInfo&, const bool);
    typedef void (FieldUtility::*ImposeVectorFieldOnNodes)(Variable<array_1d<double, 3> >&, const array_1d<double, 3>, VectorField<3>::Pointer, ModelPart&, const ProcessInfo&, const bool);
    typedef void (FieldUtility::*ImposeVelocityFieldOnNodes)(ModelPart&, const VariablesList&);
    typedef void (FieldUtility::*ImposeFieldOnNodes)(ModelPart&, const Variable<array_1d<double, 3> >&);

    ImposeDoubleFieldOnNodes ImposeDoubleField = &FieldUtility::ImposeFieldOnNodes;
    ImposeVectorFieldOnNodes ImposeVectorField = &FieldUtility::ImposeFieldOnNodes;
    ImposeVelocityFieldOnNodes ImposeVelocityField = &FieldUtility::ImposeFieldOnNodes;
    ImposeFieldOnNodes ImposeField = &FieldUtility::ImposeFieldOnNodes;


    py::class_<FieldUtility, FieldUtility::Pointer> (m, "FieldUtility")
        .def(py::init<SpaceTimeSet::Pointer, VectorField<3>::Pointer >())
        .def("EvaluateFieldAtPoint", EvaluateDoubleField)
        .def("EvaluateFieldAtPoint", EvaluateVectorField)
        .def("ImposeFieldOnNodes", ImposeDoubleField)
        .def("ImposeFieldOnNodes", ImposeVectorField)
        .def("ImposeFieldOnNodes", ImposeVelocityField)
        .def("ImposeFieldOnNodes", ImposeField)
        ;

    py::class_<FluidFieldUtility, FluidFieldUtility::Pointer, FieldUtility> (m, "FluidFieldUtility")
        .def(py::init<SpaceTimeSet::Pointer, VelocityField::Pointer, const double, const double >())
        ;


    typedef void (CustomFunctionsCalculator<3>::*CopyValuesScalar)(ModelPart&, const Variable<double>&, const Variable<double>&);
    typedef void (CustomFunctionsCalculator<3>::*CopyValuesVector)(ModelPart&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&);
    typedef void (CustomFunctionsCalculator<3>::*SetValuesScalar)(ModelPart&, const double&, const Variable<double>&);
    typedef void (CustomFunctionsCalculator<3>::*SetValuesVector)(ModelPart&, const array_1d<double, 3>&, const Variable<array_1d<double, 3>>&);

    CopyValuesScalar CopyValuesFromFirstToSecondScalar = &CustomFunctionsCalculator<3>::CopyValuesFromFirstToSecond;
    CopyValuesVector CopyValuesFromFirstToSecondVector = &CustomFunctionsCalculator<3>::CopyValuesFromFirstToSecond;
    SetValuesScalar SetValueOfAllNotesScalar = &CustomFunctionsCalculator<3>::SetValueOfAllNotes;
    SetValuesVector SetValueOfAllNotesVector = &CustomFunctionsCalculator<3>::SetValueOfAllNotes;


    py::class_<CustomFunctionsCalculator <2> > (m, "CustomFunctionsCalculator2D")
        .def(py::init<>())
        .def("CalculatePressureGradient", &CustomFunctionsCalculator <2>::CalculatePressureGradient)
        .def("AssessStationarity", &CustomFunctionsCalculator <2>::AssessStationarity)
        .def("CalculateDomainVolume", &CustomFunctionsCalculator <2>::CalculateDomainVolume)
        .def("CalculateTotalHydrodynamicForceOnParticles", &CustomFunctionsCalculator <2>::CalculateTotalHydrodynamicForceOnParticles)
        .def("CalculateTotalHydrodynamicForceOnFluid", &CustomFunctionsCalculator <2>::CalculateTotalHydrodynamicForceOnFluid)
        .def("CalculateGlobalFluidVolume", &CustomFunctionsCalculator <2>::CalculateGlobalFluidVolume)
        ;

    py::class_<CustomFunctionsCalculator <3> > (m, "CustomFunctionsCalculator3D")
        .def(py::init<>())
        .def("CalculatePressureGradient", &CustomFunctionsCalculator <3>::CalculatePressureGradient)
        .def("AssessStationarity", &CustomFunctionsCalculator <3>::AssessStationarity)
        .def("CalculateDomainVolume", &CustomFunctionsCalculator <3>::CalculateDomainVolume)
        .def("CalculateTotalHydrodynamicForceOnParticles", &CustomFunctionsCalculator <3>::CalculateTotalHydrodynamicForceOnParticles)
        .def("CalculateTotalHydrodynamicForceOnFluid", &CustomFunctionsCalculator <3>::CalculateTotalHydrodynamicForceOnFluid)
        .def("CalculateGlobalFluidVolume", &CustomFunctionsCalculator <3>::CalculateGlobalFluidVolume)
        .def("CopyValuesFromFirstToSecond", CopyValuesFromFirstToSecondScalar)
        .def("CopyValuesFromFirstToSecond", CopyValuesFromFirstToSecondVector)
        .def("SetValueOfAllNotes", SetValueOfAllNotesScalar)
        .def("SetValueOfAllNotes", SetValueOfAllNotesVector)
        ;

    //**********************************************************************************************************************************************
    // WARNING!!: function RecoverSuperconvergentGradient uses an algorithm under a GPL 3.0 licence which CANNOT be included in comercial products.

//    typedef void (DerivativeRecovery<3>::*RecoverGradientScalar)(ModelPart&, Variable<double>&, Variable<array_1d<double, 3> >&);
//    typedef void (DerivativeRecovery<3>::*RecoverGradientComponent)(ModelPart&, Variable<double>&, Variable<double>&);
//    RecoverGradientScalar RecoverSuperconvergentGradientScalar = &DerivativeRecovery<3>::RecoverSuperconvergentGradient<std::size_t TDim, class TScalarVariable>;
//    RecoverGradientComponent RecoverSuperconvergentGradientComponent = &DerivativeRecovery<3>::RecoverSuperconvergentGradient<std::size_t TDim, class TScalarVariable>;

    py::class_<FlowStationarityCheck> (m, "FlowStationarityCheck")
        .def(py::init<ModelPart& , const double>())
        .def("AssessStationarity", &FlowStationarityCheck::AssessStationarity)
        .def("GetCharacteristicPressureDerivative", &FlowStationarityCheck::GetCharacteristicPressureDerivative)
        .def("GetCurrentPressureDerivative", &FlowStationarityCheck::GetCurrentPressureDerivative)
        .def("GetTolerance", &FlowStationarityCheck::GetTolerance)
        .def("GetTransienceMeasure", &FlowStationarityCheck::GetTransienceMeasure)
        ;

    py::class_<DerivativeRecovery <3> > (m, "DerivativeRecoveryTool3D")
        .def(py::init<ModelPart&, Parameters&>())
        .def("AddTimeDerivativeComponent", &DerivativeRecovery <3>::AddTimeDerivativeComponent)
        .def("RecoverSuperconvergentGradient", &DerivativeRecovery <3>::RecoverSuperconvergentGradient< Variable<double> >)
        .def("RecoverSuperconvergentLaplacian", &DerivativeRecovery <3>::RecoverSuperconvergentLaplacian)
        .def("RecoverSuperconvergentVelocityLaplacianFromGradient", &DerivativeRecovery <3>::RecoverSuperconvergentVelocityLaplacianFromGradient)
        .def("RecoverSuperconvergentMatDeriv", &DerivativeRecovery <3>::RecoverSuperconvergentMatDeriv)
        .def("CalculateVectorMaterialDerivativeFromGradient", &DerivativeRecovery <3>::CalculateVectorMaterialDerivativeFromGradient)
        .def("CalculateVectorMaterialDerivativeComponent", &DerivativeRecovery <3>::CalculateVectorMaterialDerivativeComponent)
        .def("CalculateVorticityFromGradient", &DerivativeRecovery <3>::CalculateVorticityFromGradient)
        .def("CalculateVorticityContributionOfTheGradientOfAComponent", &DerivativeRecovery <3>::CalculateVorticityContributionOfTheGradientOfAComponent)
        .def("RecoverSuperconvergentMatDerivAndLaplacian", &DerivativeRecovery <3>::RecoverSuperconvergentMatDerivAndLaplacian)
        .def("CalculateGradient", &DerivativeRecovery <3>::CalculateGradient< Variable<double> >)
        .def("SmoothVectorField", &DerivativeRecovery <3>::SmoothVectorField)
        .def("CalculateVectorMaterialDerivative", &DerivativeRecovery <3>::CalculateVectorMaterialDerivative)
        .def("CalculateVectorMaterialDerivativeExactL2", &DerivativeRecovery <3>::CalculateVectorMaterialDerivativeExactL2)
        .def("CalculateVectorMaterialDerivativeExactL2Parallel", &DerivativeRecovery <3>::CalculateVectorMaterialDerivativeExactL2Parallel)
        .def("RecoverLagrangianAcceleration", &DerivativeRecovery <3>::RecoverLagrangianAcceleration)
        .def("CalculateVectorLaplacian", &DerivativeRecovery <3>::CalculateVectorLaplacian)
        .def("CalculateVelocityLaplacianRate", &DerivativeRecovery <3>::CalculateVelocityLaplacianRate)
        ;

    py::class_<DerivativeRecovery <2> > (m, "DerivativeRecoveryTool2D")
        .def(py::init<ModelPart&, Parameters&>())
        .def("AddTimeDerivativeComponent", &DerivativeRecovery <2>::AddTimeDerivativeComponent)
        .def("RecoverSuperconvergentGradient", &DerivativeRecovery <2>::RecoverSuperconvergentGradient< Variable<double> >)
        .def("RecoverSuperconvergentLaplacian", &DerivativeRecovery <2>::RecoverSuperconvergentLaplacian)
        .def("RecoverSuperconvergentVelocityLaplacianFromGradient", &DerivativeRecovery <2>::RecoverSuperconvergentVelocityLaplacianFromGradient)
        .def("RecoverSuperconvergentMatDeriv", &DerivativeRecovery <2>::RecoverSuperconvergentMatDeriv)
        .def("CalculateVectorMaterialDerivativeFromGradient", &DerivativeRecovery <2>::CalculateVectorMaterialDerivativeFromGradient)
        .def("CalculateVectorMaterialDerivativeComponent", &DerivativeRecovery <2>::CalculateVectorMaterialDerivativeComponent)
        .def("CalculateVorticityFromGradient", &DerivativeRecovery <2>::CalculateVorticityFromGradient)
        .def("CalculateVorticityContributionOfTheGradientOfAComponent", &DerivativeRecovery <2>::CalculateVorticityContributionOfTheGradientOfAComponent)
        .def("RecoverSuperconvergentMatDerivAndLaplacian", &DerivativeRecovery <2>::RecoverSuperconvergentMatDerivAndLaplacian)
        .def("CalculateGradient", &DerivativeRecovery <2>::CalculateGradient< Variable<double> >)
        .def("SmoothVectorField", &DerivativeRecovery <2>::SmoothVectorField)
        .def("CalculateVectorMaterialDerivative", &DerivativeRecovery <2>::CalculateVectorMaterialDerivative)
        .def("CalculateVectorMaterialDerivativeExactL2", &DerivativeRecovery <2>::CalculateVectorMaterialDerivativeExactL2)
        .def("CalculateVectorMaterialDerivativeExactL2Parallel", &DerivativeRecovery <2>::CalculateVectorMaterialDerivativeExactL2Parallel)
        .def("RecoverLagrangianAcceleration", &DerivativeRecovery <2>::RecoverLagrangianAcceleration)
        .def("CalculateVectorLaplacian", &DerivativeRecovery <2>::CalculateVectorLaplacian)
        .def("CalculateVelocityLaplacianRate", &DerivativeRecovery <2>::CalculateVelocityLaplacianRate)
        ;
    //**********************************************************************************************************************************************


    py::class_<BassetForceTools> (m, "BassetForceTools")
        .def(py::init<Parameters&>())
        .def("FillDaitcheVectors", &BassetForceTools::FillDaitcheVectors)
        .def("FillHinsbergVectors", &BassetForceTools::FillHinsbergVectors)
        .def("AppendIntegrands", &BassetForceTools::AppendIntegrands)
        .def("AppendIntegrandsWindow", &BassetForceTools::AppendIntegrandsWindow)
        .def("AppendIntegrandsImplicit", &BassetForceTools::AppendIntegrandsImplicit)
        ;

    py::class_<BinBasedDEMFluidCoupledMapping <2, SphericParticle> > (m, "BinBasedDEMFluidCoupledMapping2D")
        .def(py::init<Parameters&>())
        .def(py::init<Parameters&, SpatialSearch::Pointer>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::InterpolateFromFluidMesh)
        .def("ImposeFlowOnDEMFromField", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::ImposeFlowOnDEMFromField)
        .def("ImposeVelocityOnDEMFromFieldToAuxVelocity", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::ImposeVelocityOnDEMFromFieldToAuxVelocity)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddDEMCouplingVariable<double>)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddDEMCouplingVariable<array_1d<double, 3>>)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddFluidCouplingVariable<double>)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddFluidCouplingVariable<array_1d<double, 3>>)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddDEMVariablesToImpose<double>)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <2,SphericParticle> ::AddDEMVariablesToImpose<array_1d<double, 3>>)
        ;

    py::class_<BinBasedDEMFluidCoupledMapping <2, NanoParticle> > (m, "BinBasedNanoDEMFluidCoupledMapping2D")
        .def(py::init<Parameters&>())
        .def(py::init<Parameters&, SpatialSearch::Pointer>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::InterpolateFromFluidMesh)
        .def("ImposeFlowOnDEMFromField", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::ImposeFlowOnDEMFromField)
        .def("ImposeVelocityOnDEMFromFieldToAuxVelocity", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::ImposeVelocityOnDEMFromFieldToAuxVelocity)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddDEMCouplingVariable<double>)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddDEMCouplingVariable<array_1d<double, 3>>)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddFluidCouplingVariable<double>)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddFluidCouplingVariable<array_1d<double, 3>>)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddDEMVariablesToImpose<double>)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <2,NanoParticle> ::AddDEMVariablesToImpose<array_1d<double, 3>>)
        ;

	py::class_<BinBasedDEMFluidCoupledMapping <3, SphericParticle> > (m, "BinBasedDEMFluidCoupledMapping3D")
	    .def(py::init<Parameters&>())
        .def(py::init<Parameters&, SpatialSearch::Pointer>())
	    .def("InterpolateVelocityOnAuxVelocity", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::InterpolateVelocityOnAuxVelocity)
	    .def("ImposeVelocityOnDEMFromFieldToAuxVelocity", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::ImposeVelocityOnDEMFromFieldToAuxVelocity)
	    .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::InterpolateFromFluidMesh)
	    .def("ImposeFlowOnDEMFromField", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::ImposeFlowOnDEMFromField)
	    .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::InterpolateFromDEMMesh)
	    .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::HomogenizeFromDEMMesh)
	    .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::ComputePostProcessResults)
	    .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddDEMCouplingVariable<double>)
	    .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddDEMCouplingVariable<array_1d<double, 3>>)
	    .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddFluidCouplingVariable<double>)
	    .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddFluidCouplingVariable<array_1d<double, 3>>)
	    .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddDEMVariablesToImpose<double>)
	    .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddDEMVariablesToImpose<array_1d<double, 3>>)
	    .def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddFluidVariableToBeTimeFiltered<double>)
	    .def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMapping <3,SphericParticle> ::AddFluidVariableToBeTimeFiltered<array_1d<double, 3>>)
	    ;

    py::class_<BinBasedDEMFluidCoupledMapping <3, NanoParticle> > (m, "BinBasedNanoDEMFluidCoupledMapping3D")
        .def(py::init<Parameters&>())
        .def(py::init<Parameters&, SpatialSearch::Pointer>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::InterpolateFromFluidMesh)
        .def("ImposeFlowOnDEMFromField", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::ImposeFlowOnDEMFromField)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddDEMCouplingVariable<double>)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddDEMCouplingVariable<array_1d<double, 3>>)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddFluidCouplingVariable<double>)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddFluidCouplingVariable<array_1d<double, 3>>)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddDEMVariablesToImpose<double>)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddDEMVariablesToImpose<array_1d<double, 3>>)
        .def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddFluidVariableToBeTimeFiltered<double>)
        .def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMapping <3,NanoParticle> ::AddFluidVariableToBeTimeFiltered<array_1d<double, 3>>)
        ;

    py::class_<DerivativeRecoveryMeshingTools <2> > (m, "DerivativeRecoveryMeshingTools2D")
        .def(py::init<>())
        .def("FillUpEdgesModelPartFromSimplicesModelPart", &DerivativeRecoveryMeshingTools<2>::FillUpEdgesModelPartFromSimplicesModelPart)
        ;
    py::class_<DerivativeRecoveryMeshingTools <3> > (m, "DerivativeRecoveryMeshingTools3D")
        .def(py::init<>())
        .def("FillUpEdgesModelPartFromSimplicesModelPart", &DerivativeRecoveryMeshingTools<3>::FillUpEdgesModelPartFromSimplicesModelPart)
        ;

    py::class_<EmbeddedVolumeTool <3> > (m, "EmbeddedVolumeTool")
        .def(py::init<>())
        .def("CalculateNegativeDistanceVolume", &EmbeddedVolumeTool <3> ::CalculateNegativeDistanceVolume)
        ;

    py::class_<Bentonite_Force_Based_Inlet, Bentonite_Force_Based_Inlet::Pointer, DEM_Force_Based_Inlet > (m, "Bentonite_Force_Based_Inlet")
        .def(py::init<ModelPart&, array_1d<double, 3>, const int >())
        .def(py::init<ModelPart&, array_1d<double, 3>>())
        ;

    py::class_<SwimmingDemInPfemUtils> (m, "SwimmingDemInPfemUtils")
        .def(py::init<>())
        .def("TransferWalls", &SwimmingDemInPfemUtils::TransferWalls)
        ;

    py::class_<ErrorNormCalculator> (m, "ErrorNormCalculator")
        .def(py::init<>())
        .def("ComputeDofsErrors", &ErrorNormCalculator::ComputeDofsErrors)
        .def("GetL2VectorErrorNorm", &ErrorNormCalculator::GetL2VectorErrorNorm)
        .def("GetL2MaterialAccelNormUsingNodeDerivatives", &ErrorNormCalculator::GetL2MaterialAccelNormUsingNodeDerivatives)
        .def("GetL2ScalarErrorNorm", &ErrorNormCalculator::GetL2ScalarErrorNorm)
        .def("GetH1ScalarErrorSemiNorm", &ErrorNormCalculator::GetH1ScalarErrorSemiNorm)
        .def("GetH1VectorErrorSemiNorm", &ErrorNormCalculator::GetH1VectorErrorSemiNorm)
        ;

    py::class_<ErrorNormTorusCalculator> (m, "ErrorNormTorusCalculator")
        .def(py::init<const double, const double, const double>())
        .def("getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues", &ErrorNormTorusCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues)
        .def("getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues", &ErrorNormTorusCalculator::getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues)
        .def("getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues", &ErrorNormTorusCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues)
        .def("getL2NormFluidAccelWithRecoveryUsingGaussExactValues", &ErrorNormTorusCalculator::getL2NormFluidAccelWithRecoveryUsingGaussExactValues)
        ;

    py::class_<ErrorNormEthierFieldCalculator> (m, "ErrorNormEthierFieldCalculator")
        .def(py::init<const double, const double>())
        .def(py::init<const double, const double, const bool>())
        .def("getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues", &ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues)
        .def("getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues", &ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues)
        .def("getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues", &ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues)
        .def("getL2NormFluidAccelWithRecoveryUsingGaussExactValues", &ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithRecoveryUsingGaussExactValues)
        ;

    py::class_<MeshRotationUtility> (m, "MeshRotationUtility")
        .def(py::init<Parameters&>())
        .def("RotateMesh", &MeshRotationUtility::RotateMesh)
        .def("RotateDEMMesh", &MeshRotationUtility::RotateDEMMesh)
        .def("SetStationaryField", &MeshRotationUtility::SetStationaryField)
        .def("RotateFluidVelocities", &MeshRotationUtility::RotateFluidVelocities)
        ;
    py::class_<RenumberingNodesUtility> (m, "RenumberingNodesUtility")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&,ModelPart&,ModelPart&>())
        .def("Renumber", &RenumberingNodesUtility::Renumber)
        .def("UndoRenumber", &RenumberingNodesUtility::UndoRenumber)
        ;
    }

}  // namespace Python.

} // Namespace Kratos
