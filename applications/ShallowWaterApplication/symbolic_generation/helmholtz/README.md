# Helmholtz Wave Equation Element

## ELEMENT DESCRIPTION:
Current directory contains the documentation for the symbolic derivation of the _"helmholtz"_ element. This element includes a formulation of _Helmoltz_ element for **2D** case.

## SYMBOLIC GENERATOR SETTINGS:
* This symbolic generator is under construction.

* Dimension to compute: This symbolic generator is valid for both **2D** and **3D**  cases. Since the element scope is to solve the SWE, only **2D** has a physical sense.
* Linearisation settings: _FullNR_ considers the convective velocity as _"v-vmesh"_, hence _v_ is taken into account in the derivation of the **LHS** and **RHS**. Picard (_a.k.a. QuasiNR_) considers the convective velocity as "a", thus it is considered as a constant in the derivation of the **LHS** and **RHS**.
* Artificial compressiblity: If set to true, the time derivative of the density is introduced in the mass conservation equation together with the state equation $dp/d\rho=c^2$ (being c the sound velocity). These assumptions add some extra terms to the usual **Navier-Stokes** equations that are intended to act  as a soft artificial compressibility, which is controlled by the value of "_c_", meaning that if this value is large enough the artificial  compressibility terms vanish.

## INSTRUCTIONS
Run:
~~~py
python generate_helmholtz_element.py
~~~
Then  file "_helmholtz.cpp_" is generated automatically. Such file should be copied within the "_custom_elements_" folder of the
**ShallowWaterApplication**. The corresponding header file "helmholtz.h", which implements the element is already stored in the custom_elements folder.
