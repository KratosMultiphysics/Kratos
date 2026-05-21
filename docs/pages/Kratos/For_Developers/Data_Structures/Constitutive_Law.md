---
title: Constitutive Law
keywords: 
tags: [Constitutive Law Tutorial C++]
sidebar: kratos_for_developers
summary: 
---

## Introduction

The constitutive law behaviour is dealt with in kratos by the use of the class "ConstitutiveLaw", 
with a public interface defined in the file
   kratos/kratos/includes/constitutive_law.h
which also provides some rather extensive inline documentation (in the form of comments in the code).

By design such file aims to provide a very flexible interface to constitutive law modelling, with the specific goal of 
'''maximizing the flexibility in the implementation of complex constitutive behaviours'''. While such approach provide obvious advantages, it 
also implies that the API is more complex than what would be strictly needed for very simple constitutive laws.

The objective of current HowTo is to provide a brief introduction to the interface

###  Conventions
Through the whole section, the following convenctions will be employed:

''voigt notation:''

- 3D case:
  STRAIN Voigt Notation:  e00 e11 e22 2*e01 2*e12 2*e02
  STRESS Voigt Notation:  s00 s11 s22   s01   s12   s02

- 2D plane strain/axisymmetric case (4 stress components)
  STRAIN Voigt Notation:  e00 e11 e22 2*e01 
  STRESS Voigt Notation:  s00 s11 s22   s01   

- 2D plane stress/strain
  STRAIN Voigt Notation:  e00 e11 2*e01 
  STRESS Voigt Notation:  s00 s11   s01

The constitutive law works on the basis of the '''total deformation gradient F''', defined as

   F  := D(X) / D(X0) 

that is, as the deformation gradient connecting the original and deformed configuration

where the initial position X0 is the one obtained by 

```cpp 
const array_1d<double,3>& X0 = node->GetInitialPosition()
```
{: data-lang="C++"}

and the deformed one by

```cpp 
const array_1d<double,3>& X = node->Coordinates() 
//must coincide with      X = node->GetInitialPosition() + node.FastGetSolutionStepValue(DISPLACEMENT);
```
{: data-lang="C++"}

The ConstitutiveLaw '''always returns the total stress'''. Formulations expressed in terms of strain increments shall store internally the strain stresses from which the increment
shall be computed

## Usage API
The constitutive law API is based on the use of an auxiliary "Parameters" data structure, designed to encapsulate the data to be passed to the CL and received from it.
The parameters data structure should be initialized using the following constructor:

```cpp 
Parameters(
    const GeometryType& rElementGeometry,
    const Properties& rMaterialProperties,
    const ProcessInfo& rCurrentProcessInfo
)
```
{: data-lang="C++"}

Thus allowing to encapsulate the pointer to the elemental properties, to the element geometry and to the process info.

The data structure '''does not contain any internal storage''' and should be initialized with pointers to memory ''owned by the caller element''. 
Full documentation of the code can be found in the file constitutive_law.h (https://github.com/KratosMultiphysics/Kratos/blob/feature-clean-claw-flags/kratos/includes/constitutive_law.h).
For ease, the getter interface, ''returning a reference to the encapsulated data'', is reported here

```cpp 
GetOptions() // Returns a reference to a flag container, to be employed in passing options to the CL
GetDeterminantF()   
GetDeformationGradientF() 
GetShapeFunctionsValues() 
GetShapeFunctionsDerivatives() 
GetStrainVector() // INPUT/OUTPUT -- note that F will be used preferentially instead of the input strain
GetStressVector() 
GetConstitutiveMatrix() 
GetProcessInfo() 
GetMaterialProperties() 
GetElementGeometry()
```
{: data-lang="C++"}

The "Options" flag represents the fundamental tool in steering the control of the constitutive law behaviour. The interface provides a number of boolean flags that can be passed to the constitutive law:

'''fundamental flags''': these are flags that most users will need to employ in everydays usage

```cpp 
USE_ELEMENT_PROVIDED_STRAIN // (only valid in small strain or for rate form 
                            //  note that when employed large strain CLs can not be used)
COMPUTE_STRESS
COMPUTE_CONSTITUTIVE_TENSOR
```
{: data-lang="C++"}

The results of the computations done employing such flags is returned through the data, and can be calling the functions `GetStressVector()` and `GetConstitutiveMatrix()` respectively.

'''optional flags''':
these are flags that provide fine tuning in the expected behaviour of the CL. However they are '''strictly optional''', in the sense that '''many constitutive laws may not implement them.'''

```cpp 
ISOCHORIC_TENSOR_ONLY      
VOLUMETRIC_TENSOR_ONLY 
MECHANICAL_RESPONSE_ONLY
THERMAL_RESPONSE_ONLY   
```
{: data-lang="C++"}

The element should check if the feature is available (that is, if the CL provides the feature) in the Check function of the element

'''internal flags'''
These are flags that are useful for computation internally to the CL. '''they should NOT be used by the element'''
```cpp 
INITIALIZE_MATERIAL_RESPONSE
FINALIZE_MATERIAL_RESPONSE
```
{: data-lang="C++"}

## STRESS/STRAIN MEASURES

A fundamental feature of the constitutive law is to implement internally the transformations between different stress measures. 
When a user writes an element he should decide what
stress measure is desired. The list of available options is provided in the enum:

```cpp 
enum StressMeasure
{
    StressMeasure_PK1,            //stress related to reference configuration non-symmetric
    StressMeasure_PK2,            //stress related to reference configuration
    StressMeasure_Kirchhoff,      //stress related to current   configuration
    StressMeasure_Cauchy          //stress related to current   configuration
}
```
{: data-lang="C++"}

at the moment of writing the element the developer should hence query to the constitutive law for the desired stress measure. This is achieved by picking
one of the functions:

```cpp 
CalculateMaterialResponsePK1( parameters )
CalculateMaterialResponsePK2( parameters )
CalculateMaterialResponseKirchhoff( parameters )
CalculateMaterialResponseCauchy( parameters )
```
{: data-lang="C++"}

the elasticity tensor and the corresponding stresses will be stored in the "parameters" and shall be accessed by the Getter functions described above.

At the end of a given solution step (typically in the FinalizeSolutionStep function of the element), internal variables should be updated by calling the function

```cpp 
FinalizeMaterialResponsePK1( parameters )
FinalizeMaterialResponsePK2( parameters )
FinalizeMaterialResponseKirchhoff( parameters )
FinalizeMaterialResponseCauchy( parameters )
```
{: data-lang="C++"}

An example of usage of the ConstitutiveLaw from within a total lagrangian element could be as follows:

```cpp 
Parameters parameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

// Here we essentially set the input parameters
parameters.SetDeterminantF(detF) //assuming the determinant is computed somewhere else
parameters.SetDeformationGradientF(F) //F computed somewhere else

// Here we set the space on which the results shall be written
Matrix ConstitutiveTensor(strain_size,strain_size); //note that C is allocated in the element
Vector stress(strain_size);
parameters.SetConstitutiveMatrix(ConstitutiveTensor) //assuming the determinant is computed somewhere else
parameters.SetStressVector(stress) //F computed somewhere else

// Instruct the constitutive law that both stress and the constitutivetensor are needed
parameters.GetOptions().Set(COMPUTE_STRESS, True)
parameters.GetOptions().Set(COMPUTE_CONSTITUTIVE_TENSOR, True)

// Actually do the computations in the ConstitutiveLaw    
constitutivelaw->CalculateMaterialResponsePK2(parameters); //here the calculations are actually done 

// Here stress and C are already updated since we passed them to the CL
```
{: data-lang="C++"}

The constitutive law also provides helper functions to allow pushing/pulling stresses from one configuration to another.

## Providing the Deformation gradient F or the strain
Within Kratos the strain can be either provided to the constitutive law as an input parameter, OR obtained internally within the constitutive law
depending on the kinematics chosen. 

### Option 1 -- Provide strain as input parameter
This option is chosen by selecting

```cpp 
flags.Set(USE_ELEMENT_PROVIDED_STRAIN, true)
```
{: data-lang="C++"}

in this case the constitutive law reads the value of "StrainVector" and uses it as a basis for the computations of the stresses. Note that this will be implemented only by some of the constitutive laws. The user should hence check if this feature is implemented

### Option 2 -- Provide F as input parameter
The idea in this case is to provide a DeformationGradient '''F''' as an input parameter.

In the case of small deformation kinematics, is is not customary to define the deformation gradient. For this case (small deformations) one one can construct the deformation gradient given the strain vector as

```cpp 
F(0,0) = 1+strain(0);     F(0,1) = 0.5*strain(3); F(0,2) = 0.5*strain(5)
F(1,0) = 0.5*strain(3);   F(1,1) = 1+strain(1);   F(1,2) = 0.5*strain(4)
F(2,0) = 0.5*strain(5);   F(2,1) = 0.5*strain(4); F(2,2) = 1+strain(2)
```
{: data-lang="C++"}

for the 3D case.
The interest of this, is to allow using a large deformation law within a small deformation element.

## Compatibility Check Features
The Constitutive Law provides extensive facilities to check for compatibility between the element and the CL.

The element can ask to the CL for available features (one of the local flags of the CL). for example

```cpp 
//getting the struct containing the features
ConstitutiveLaw::Features LawFeatures = cl->GetLawFeatures(LawFeatures);

//check if the law is plane strain (for example)
LawFeatures.mOptions.Is(ConstitutiveLaw::PLANE_STRAIN_LAW)

//here we get the strain size and the space in which the CL lives
unsigned int expected_strain_size = LawFeatures.GetStrainSize();
unsigned int working_space_dimension = LawFeatures.GetSpaceDimension();
```
{: data-lang="C++"}

this relies on the CL being called implementing the function "GetLawFeatures", to make an example

```cpp   
void MyLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    // Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    // Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}
```
{: data-lang="C++"}

additionally each constitutive law provides a '''"Check"''' method, which verifies if the parameters needed are correctly set and if they are in the correct range (for example it may check that the poisson ratio does not exceed 0.5). An example of this could be
     
```cpp  
int HyperElasticPlastic3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;

    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;
}
```
{: data-lang="C++"}

## Conversion between different stress measures
The CL implements pull_back and push_forward mechanism for both the output stress and for the constitutive tensor.
An exampe of use is the following: let's consider one had a PK2 stress stored in the value PK2_stress_value, and wants to transform it to Cauchy

```cpp   
Matrix stress_value = //... here let's suppose that as input we have PK2 stress
stress_value = cl->TransformStresses(
    stress_value, 
    F, 
    detF, 
    ConstitutiveLawType::StressMeasure::StressMeasure_PK2,
    ConstitutiveLawType::StressMeasure::StressMeasure_Cauchy
);
```
{: data-lang="C++"}

the function implemented in the base class will do internally the necessary transformations so that on output "stress_value" we will get the corresponding Cauchy stress

The base class CL also implements the operations needed to do a pull back and push forward of the constitutive tensor, through the functions 

* PullBackConstitutiveMatrix
* PushForwardConstitutiveMatrix 

     ... here it would be good an example of usage ...


## Input Output features
The CL provides a GetValue/SetValue interface which allows setting providing or getting output from the CL.
It should be noted however that in order to keep the CL object as lightweight an automatic databas was not provided, meaning that the user needs to provide its own storage. To make an example let's consider a damage law which stores the variable "mDamage" as a class variable. It is possible to give access both in reading and in writing to this variable by implementing

```cpp   
double& GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    if(rThisVariable == DAMAGE)
    {
        rValue = mDamage;
        return rValue;
    }
    else
        KRATOS_ERROR << "sorry attempting to do a GetValue for the variable " << rThisVariable << " which is not available in the CL"; 
}
```
{: data-lang="C++"}

or

```cpp   
void SetValue(const Variable<double>& rVariable,
            const double& rValue,
            const ProcessInfo& rCurrentProcessInfo);
{
    if(rThisVariable == DAMAGE)
    {
        mDamage = rValue;
    }
    else
        KRATOS_ERROR << "sorry attempting to do a SetValue for the variable " << rThisVariable << " which is not available in the CL"; 
}
```
{: data-lang="C++"}

the CL also provides a function "CalculateValue" to allow accessing to variables that need to be calculated on the flight. For example to calculate and get the  STRAIN_ENERGY one may do

```cpp      
// ... here set the options ... 
double energy = 0.0
energy = CalculateValue(rParameterValues, STRAIN_ENERGY, energy);
```
{: data-lang="C++"}

## Storing State Variables at convgerence
By design, and with the goal of not storing unnecessary variables, the CLs **should not** store intermediate values of the state variables.
The moment at which InternalVariables should be update is the call to "FinalizeMaterialResponseXXX" which needs to be guaranteed to happen **stricly once** per time step, after convergence is reached. The function FinalizeMaterialResponse can use internally the private flag FINALIZE_MATERIAL_RESPONSE

## Using CLs in an explicit context
The options passed to the "CalculateMaterialResponse" allow the user to avoid the calculation of the ElasticityTensor. It is crucial for performance that in the implementation such flag is respected and that the calculation of the stresses does not imply computing the elasticity tensor

## Treatment of Initial Strains
TBD

## Interpolate State variables
**WARNING!: this is only a proposal and it is NOT currently respected**
The Kratos variable INTERNAL_VARIABLES (of type Vector<Variable>) **has a special meaning**. 
It is designed to contain consecutively all of the converged internal variables.
The function "Calculate" should compute this variable as needed. the function SetValue should allow setting such variables. 
Note that the load and serialize function would benefit making use of this function.
