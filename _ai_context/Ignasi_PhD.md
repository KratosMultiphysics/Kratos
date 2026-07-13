**Development of new computational methods for fluid-structure interaction analysis of multi-fractured media** 

## **Doctoral Thesis** 

## **Ignasi de Pouplana Sardà** 

**Supervisor: Eugenio Oñate Ibáñez de Navarra** 

**PhD in Structural Analysis** 

**Departament d’Enginyeria Civil i Ambiental (DECA) Universitat Politècnica de Catalunya** 

**January 2018** 

## **Acta de calificación de tesis doctoral** 

**Curso académico:** 

Nombre y apellidos Programa de doctorado Unidad estructural responsable del programa 

## **Resolución del Tribunal** 

Reunido el Tribunal designado a tal efecto, el doctorando / la doctoranda expone el tema de la su tesis doctoral titulada ____________________________________________________________________________________ 

__________________________________________________________________________________________. 

Acabada la lectura y después de dar respuesta a las cuestiones formuladas por los miembros titulares del tribunal, éste otorga la calificación: 

**==> picture [433 x 169] intentionally omitted <==**

**----- Start of picture text -----**<br>
NO APTO  APROBADO  NOTABLE  SOBRESALIENTE<br>(Nombre, apellidos y firma)  (Nombre, apellidos y firma)<br>Presidente/a  Secretario/a<br>(Nombre, apellidos y firma)  (Nombre, apellidos y firma)  (Nombre, apellidos y firma)<br>Vocal  Vocal  Vocal<br>______________________, _______ de __________________ de _______________<br>**----- End of picture text -----**<br>


El resultado del escrutinio de los votos emitidos por los miembros titulares del tribunal, efectuado por la Escuela de Doctorado, a instancia de la Comisión de Doctorado de la UPC, otorga la MENCIÓN CUM LAUDE: 

**==> picture [456 x 99] intentionally omitted <==**

**----- Start of picture text -----**<br>
SÍ  NO<br>(Nombre, apellidos y firma)  (Nombre, apellidos y firma)<br>Presidente de la Comisión Permanente de la Escuela de  Secretario de la Comisión Permanente de la Escuela de<br>Doctorado  Doctorado<br>**----- End of picture text -----**<br>


Barcelona a _______ de ____________________ de __________ 

## **Abstract** 

The objective of this thesis is the derivation and implementation of a robust Finite Element formulation for the solution of solid-pore fluid coupled problems in multi-fractured porous media. 

A coupled displacement-pore pressure FEM formulation for solving solid-pore fluid interaction problems is first introduced. The interaction between both components is governed by two equations: the balance of momentum for the mixture solid-fluid and the mass balance for the pore fluid. 

Under nearly undrained-incompressible conditions, such formulation suffers from instability problems because of the violation of Babuska-Brezzi conditions. In order to work with elements of equal order interpolation for the displacement and pore pressure, the formulation is stabilized by means of the Finite Increment Calculus method (FIC). The FIC-stabilized formulation is tested against stable elements with a higher order interpolation for the displacement field in 2D and 3D examples. 

Continuum damage mechanics is the basis of the crack growth strategy for the proposed fracture propagation technique. The strain softening models used for quasi-brittle materials favour spurious strain localization and ill-posedness of the boundary value problem if the damage variable only depends on the strain state at the point under consideration. 

An integral-type non-local damage model associated to a characteristic length parameter is presented as a method to control the size of the fracture process zone and fully regularize the problem. Two examples are solved assessing the robustness of the model in front of changes in the spatial discretization. 

Quasi-zero-thickness interface elements are formulated to represent discontinuities in the porous domain. A bilinear cohesive fracture model is used to describe its mechanical 

Abstract 

vi 

behaviour, and a formulation derived from the cubic law models the fluid flow through the crack. 

Finally, a new methodology for the simulation of fracture propagation processes in saturated porous media is presented. The non-local damage model is used in conjunction with the interface elements to predict the degradation pattern of the domain and insert new fractures followed by remeshing. Fluid-driven fracture propagation examples in 2D and 3D are presented to illustrate the accuracy of the proposed technique. 

## **Resumen** 

El objetivo de esta tesis es la derivación e implementación de una formulación robusta de Elementos Finitos para la solución de problemas acoplados de sólido-fluido de poro en medios porosos multi-fracturados. 

Una formulación del MEF acoplada desplazamiento-presión de poro para resolver problemas de interacción solido-fluido de poro es primeramente introducida. La interacción entre ambos componentes es gobernada por dos ecuaciones: el balance de momento para la mezcla sólido-fluido y el balance de masa para el fluido de poro. 

Bajo condiciones de impermeabilidad e incompresibilidad, esta formulación sufre problemas de inestabilidad debido a la violación de las condiciones Babuska-Brezzi. Para poder trabajar con elementos de igual orden de interpolación para los desplazamientos y la presión de poro, la formulación es estabilizada mediante el método de Finite Increment Calculus (FIC). La formulación estabilizada con FIC es testeada contra elementos estables de mayor orden de interpolación para el campo de desplazamientos en ejemplos 2D y 3D. 

La mecánica del daño continua es la base de la estrategia de crecimiento de fisura para la técnica de propagación de fracturas propuesta. Los modelos de deformación con reblandecimiento utilizados para materiales cuasi-frágiles favorecen la localización espuria de las deformaciones y el mal condicionamiento del problema de valores en el contorno si la variable de daño depende únicamente del estado de deformación en el punto considerado. 

Un modelo de daño no-local de tipo integral asociado a un parámetro de longitud carac- 

Resumen 

viii 

terística es presentado como un método para controlar el tamaño de la zona de fractura y regularizar totalmente el problema. Dos ejemplos son resueltos para evaluar la robustez del modelo frente a cambios en la discretización espacial. 

Elementos de interface de espesor cuasi-cero son formulados para representar discontinuidades en el dominio poroso. Un modelo de fractura cohesiva bilineal es utilizado para describir su comportamiento mecánico, y una formulación derivada de la ley cúbica modela el flujo de fluido a través de la fisura. 

Finalmente, una nueva metodología para la simulación de procesos de propagación de fractura en medios porosos saturados es presentada. El modelo de daño no-local es empleado juntamente con los elementos de interface para predecir el mapa de degradación del dominio e insertar nuevas fracturas seguido de un remallado. Ejemplos de fractura por fluido en 2D y 3D son presentados para ilustrar la precisión de la técnica propuesta. 

_A la meva família._ 

## **Acknowledgements** 

After three years of hard work I can say that most of the objectives set at the beginning of the PhD have been achieved. The time spent in the research and elaboration of this thesis has been truly rewarding for the acquired knowledge and for all the gained friends. Indeed, this project would have been impossible without all the people around me that have helped me in some way. 

First, I would like to express my sincere gratitude to Prof. Eugenio Oñate for offering me the opportunity of working at CIMNE in this doctoral thesis. As my supervisor he always gave me wise advice and helped me work better, but he also gave me freedom to think for myself and let me grow as a researcher, which I really appreciate. I also want to thank him along with Sonia Sagristà for encouraging me to go on a three-month secondment at Tsinghua University (Beijing) in the framework of the TCAINMAND project. I can say that it has been one of the best experiences in my life. 

I specially want to thank all the officemates at CIMNE that have been close to me kindly clarifying all my doubts, and with whom I enjoyed several meals during these years: Salva Latorre, Pablo Agustín Becker, Ferran Arrufat, Mercè López, Guillermo Casas, Lorenzo Gracia, Miguel Ángel Celigueta, Miquel Santasusana and Roberto Flores. 

I am deeply grateful to Prof. Zhuo Zhuang as the supervisor of my secondment in Tsinghua University. Although I was very far from home, he always made me feel comfortable. I really appreciate all the interesting discussions we had every Thursday during the group meetings. In the same way, I would like to thank Prof. Zhanli Liu and all the officemates, specially: Tao Wang, Qinglei Zeng, Yue Gao, Liu Fengxian, Liyuan Wang, Ramiro Mena and Jan Rojek. I am also very thankful to my room mate in Tsinghua Jaime Santirso, and to Duan Wenjie. 

Acknowledgements 

xii 

I would also like to thank other colleagues at CIMNE that have helped me in various matters: Fernando Salazar, Enrique Escolano, Miguel Pasenau, Javi Gárate, Josep Maria Carbonell, Pooyan Dadvand, Riccardo Rossi, Jordi Cotela, Carlos Roig, Rubén Zorrilla and Vicente Mataix. 

I also acknowledge the Agencia de Gestió d’Ajuts Universitaris i de Recerca (AGAUR) for the financial support. 

Finally, this thesis is dedicated to my family: Lluís, Marta, Oriol, M[a] Àngels and Gabriel. They are the most important people in my life and represent a fundamental and unique support. 

## **Contents** 

|**Abstract**|**Abstract**|**Abstract**|**v**|
|---|---|---|---|
|**Acknowledgements**|||**xi**|
|**1**|**Introduction**||**1**|
||1.1|Motivation and presentation . . . . . . . . . . . . . . . . . . . . . . . . .|1|
||1.2|Objectives . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|5|
||1.3|State of the art . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|5|
|||1.3.1<br>Coupled solid-pore fuid formulation . . . . . . . . . . . . . . . . .|5|
|||1.3.2<br>Fracture modelling . . . . . . . . . . . . . . . . . . . . . . . . . .|7|
||1.4|Organization<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|10|
||1.5|Related publications<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . .|11|
|**2**|**Coupled Solid-Pore Fluid Interaction Problem**||**13**|
||2.1|Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|13|
||2.2|Finite element formulation . . . . . . . . . . . . . . . . . . . . . . . . . .|16|
||2.3|Stabilization using Finite Increment Calculus (FIC) . . . . . . . . . . . .|22|



Contents 

xiv 

|||2.3.1<br>Introduction to the FIC-stabilization procedure<br>. . . . . . . . . .|22|
|---|---|---|---|
|||2.3.2<br>Modifed FIC-stabilized form of the mass balance equation . . . .|24|
|||2.3.3<br>Variational form of the FIC-stabilized mass balance equation . . .|27|
|||2.3.4<br>Discretized form of the momentum and FIC-stabilized mass bal-||
|||ance equations<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . .|30|
||2.4|Examples<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|32|
|||2.4.1<br>Elastic soil column subjected to surface loading<br>. . . . . . . . . .|32|
|||2.4.2<br>Elastic soil foundation subjected to surface loading<br>. . . . . . . .|39|
|**3**|**Continuum Damage Mechanics**||**49**|
||3.1|Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|49|
||3.2|Isotropic damage theory<br>. . . . . . . . . . . . . . . . . . . . . . . . . . .|51|
|||3.2.1<br>Elastic-damage tangent constitutive tensor . . . . . . . . . . . . .|58|
|||3.2.2<br>Equivalent strain . . . . . . . . . . . . . . . . . . . . . . . . . . .|61|
|||3.2.3<br>Damage evolution law<br>. . . . . . . . . . . . . . . . . . . . . . . .|64|
||3.3|Local and Non-local Damage Models . . . . . . . . . . . . . . . . . . . .|67|
|||3.3.1<br>Strain localization phenomenon . . . . . . . . . . . . . . . . . . .|68|
|||3.3.2<br>Regularization of the problem . . . . . . . . . . . . . . . . . . . .|73|
||3.4|Examples<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|84|
|||3.4.1<br>Three-Point Bending Test . . . . . . . . . . . . . . . . . . . . . .|85|
|||3.4.2<br>Four-Point Shear Test<br>. . . . . . . . . . . . . . . . . . . . . . . .|92|
|**4**|**Modelling Discontinuities in Porous Media**||**99**|
||4.1|Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|99|
||4.2|Quasi-zero-thickness interface elements . . . . . . . . . . . . . . . . . . .|102|



Contents 

xv 

|||4.2.1<br>Mechanical behaviour of the fracture . . . . . . . . . . . . . . . . 107|
|---|---|---|
|||4.2.2<br>Fluid fow in the fracture . . . . . . . . . . . . . . . . . . . . . . . 112|
||4.3|Fracture propagation approach . . . . . . . . . . . . . . . . . . . . . . . . 118|
|||4.3.1<br>Crack path estimation<br>. . . . . . . . . . . . . . . . . . . . . . . . 119|
|||4.3.2<br>Interface elements insertion<br>. . . . . . . . . . . . . . . . . . . . . 122|
||4.4|Examples<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 126|
|||4.4.1<br>Fluid-driven fracture propagation test . . . . . . . . . . . . . . . . 126|
|||4.4.2<br>Crack tracking test . . . . . . . . . . . . . . . . . . . . . . . . . . 133|
|||4.4.3<br>Parallel fracture propagation test . . . . . . . . . . . . . . . . . . 139|
|**5**|**Conclusions and Future Lines of Research**<br>**145**||
||5.1|Conclusions and contributions . . . . . . . . . . . . . . . . . . . . . . . . 145|
||5.2|Lines of future work<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . 149|
|**Bibliography**<br>**157**|||



## **List of Figures** 

|1.1|Examples of problems with important interaction between a solid skeleton||
|---|---|---|
||and a fuid phase. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|2|
|1.2|Kratos’ “core and applications approach”. . . . . . . . . . . . . . . . . . .|4|
|1.3|Two-dimensional zero-thickness interface elements. . . . . . . . . . . . . .|8|
|2.1|Infnitesimal volume of porous medium as a combination of two phases:||
||a solid skeleton and a fuid medium. . . . . . . . . . . . . . . . . . . . . .|14|
|2.2|Some 2D elements that fulfl the Babuska-Brezzi conditions.<br>. . . . . . .|23|
|2.3|Computation of the term _∂_˙_p_<br>_∂n_ at the side _ij_ of a triangle adjacent to the||
||external boundary Γ. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|29|
|2.4|Elastic soil column subjected to surface loading. Geometry and boundary||
||conditions. Dimensions in _m_.<br>. . . . . . . . . . . . . . . . . . . . . . . .|33|
|2.5|Surface step loading applied in the elastic soil column problem. . . . . . .|35|
|2.6|Normalized pore pressure along the soil column (_k_ = 1_·_10_−_14 _m_2, ∆_t_ =||
||0_._02 _s_, _t_= 2 _s_). . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|36|
|2.7|Surface cyclic loading applied in the elastic soil column problem. . . . . .|37|
|2.8|Evolution of the pore pressure with time at a node near the surface (_Q_=||
||1_·_1015_N/m_2, ∆_t_= 0_._05 _s_).<br>. . . . . . . . . . . . . . . . . . . . . . . . .|38|



List of Figures 

xviii 

|2.9|Time evolution of the pore pressure. Comparison between 20 Q4P4-FIC|||
|---|---|---|---|
||elements and 60 stable Q9P4 elements (_Q_= 1_·_1015_N/m_2, _k_ = 10_−_14_m_2,|||
||∆_t_= 0_._05 _s_).<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|.|39|
|2.10|Elastic soil foundation subjected to surface loading. Geometry and bound-|||
||ary conditions. Dimensions in _m_. . . . . . . . . . . . . . . . . . . . . .|.|40|
|2.11|Surface step loading applied in the elastic soil foundation problem. . . .|.|41|
|2.12|Spatial discretizations used for the elastic soil foundation problem. . . .|.|42|
|2.13|Time evolution of the maximum pore pressure under undrained-incompressible|||
||conditions (_Q_= 1_·_1015 _N/m_2, _k_ = 1_·_10_−_14 _m_2, ∆_t_= 0_._02 _s_). . . . . .|.|43|
|2.14|Pore pressure distribution under undrained-incompressible conditions (_Q_=|||
||1_·_1015 _N/m_2, _k_ = 1_·_10_−_14 _m_2, ∆_t_= 0_._02 _s_, _t_= 1 _s_). . . . . . . . . . .|.|45|
|2.15|Time evolution of the maximum pore pressure for partially compressible|||
||and drained conditions (_Q_= 1_·_1010 _N/m_2, _k_ = 1_·_10_−_11 _m_2, ∆_t_= 0_._02 _s_).||46|
|2.16|Pore pressure distribution for partially compressible and drained condi-|||
||tions (_Q_= 1_·_1010 _N/m_2, _k_ = 1_·_10_−_11 _m_2, ∆_t_= 0_._02 _s_, _t_= 1 _s_).<br>. . .|.|47|
|3.1|Generic damaged section. . . . . . . . . . . . . . . . . . . . . . . . . . .|.|50|
|3.2|Idealized material for the description of the uniaxial damage theory. . .|.|52|
|3.3|Scheme of a uniaxial damage model through a monotonic loading process.||53|
|3.4|Scheme of a uniaxial damage model through a non-monotonic loading|||
||process.<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|.|54|
|3.5|Unidimensional stress-strain diagram throughout a damage process. . .|.|59|
|3.6|3D elastic domain for a generic equivalent strain.<br>. . . . . . . . . . . .|.|61|
|3.7|Damage surfaces in the 2D principal stress space.<br>. . . . . . . . . . . .|.|62|
|3.8|Generic unidimensional stress-strain curves for diferent softening laws.|.|65|
|3.9|Bar under uniaxial tension. . . . . . . . . . . . . . . . . . . . . . . . . .|.|68|
|3.10|Stress-strain diagram of the linear softening law. . . . . . . . . . . . . .|.|69|



List of Figures 

xix 

|3.11|Axial force acting along the bar. . . . . . . . . . . . . . . . . . . . . . . .|69|
|---|---|---|
|3.12|Possible strain values corresponding to the same stress level. . . . . . . .|70|
|3.13|Piecewise constant strain profle along the bar. . . . . . . . . . . . . . . .|70|
|3.14|Fan of possible post-peak branches of the load-displacement diagram. . .|71|
|3.15|Mesh dependence of the numerical results. . . . . . . . . . . . . . . . . .|72|
|3.16|Characteristic length of a 2D fnite element.<br>. . . . . . . . . . . . . . . .|74|
|3.17|Domain of infuence protruding through the boundary of a body. . . . . .|77|
|3.18|Grid-based non-local search. . . . . . . . . . . . . . . . . . . . . . . . . .|80|
|3.19|Non-local assembly process.<br>. . . . . . . . . . . . . . . . . . . . . . . . .|83|
|3.20|General scheme for the stress evaluation and stifness matrix assembly.<br>.|84|
|3.21|TPB test. Geometry and boundary conditions. Dimensions in _mm_. . . .|85|
|3.22|TPB test. Spatial discretizations. . . . . . . . . . . . . . . . . . . . . . .|87|
|3.23|TPB test. Force-vertical defection diagrams. . . . . . . . . . . . . . . . .|89|
|3.24|TPB test. Damage initiation in the local model. . . . . . . . . . . . . . .|90|
|3.25|TPB test. Damage initiation in the non-local model.<br>. . . . . . . . . . .|91|
|3.26|TPB test. Zoom of damage growing in the local model. . . . . . . . . . .|91|
|3.27|TPB test. Evolution of damage propagation for _le_ = 3 _mm_. . . . . . . . .|92|
|3.28|FPS test. Geometry and boundary conditions. Dimensions in mm . . . .|92|
|3.29|FPS test. Spatial discretizations.<br>. . . . . . . . . . . . . . . . . . . . . .|94|
|3.30|FPS test. Force-Crack Mouth Sliding Displacement curves. . . . . . . . .|95|
|3.31|FPS test. Damage progression in the local model. . . . . . . . . . . . . .|96|
|3.32|FPS test. Damage progression in the non-local model. . . . . . . . . . . .|97|
|3.33|FPS test. Evolution of damage propagation for _le_ = 3 _mm_. . . . . . . . .|98|



List of Figures 

xx 

|4.1|Diferent methods to represent discontinuities. . . . . . . . . . . . . . .|. 100|
|---|---|---|
|4.2|Fluid fow inside a crack with bifurcation.<br>. . . . . . . . . . . . . . . .|. 103|
|4.3|Consolidation problem with a vertical fault.<br>. . . . . . . . . . . . . . .|. 103|
|4.4|Quasi-zero-thickness interface elements. . . . . . . . . . . . . . . . . . .|. 104|
|4.5|Scheme of a generic hexahedral interface element. Global and local coor-||
||dinate system. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|. 105|
|4.6|Bilinear cohesive fracture model. . . . . . . . . . . . . . . . . . . . . . .|. 108|
|4.7|Contact between two beams connected by interface elements. . . . . . .|. 110|
|4.8|Scheme of a contact with interface elements. . . . . . . . . . . . . . . .|. 111|
|4.9|Velocity profles between two smooth parallel plates. . . . . . . . . . . .|. 113|
|4.10|Natural coordinates at the mid plane of an hexahedral interface element.|114|
|4.11|Water fow in a pre-existing fractures network. . . . . . . . . . . . . . .|. 117|
|4.12|Finite element mesh with interface elements at all possible edges.<br>. . .|. 118|
|4.13|Crack tip neighbours search. . . . . . . . . . . . . . . . . . . . . . . . .|. 120|
|4.14|Scheme of the implemented crack tips.<br>. . . . . . . . . . . . . . . . . .|. 123|
|4.15|Flow chart of the fracture propagation technique.<br>. . . . . . . . . . . .|. 124|
|4.16|Mapping of variables between meshes. . . . . . . . . . . . . . . . . . . .|. 125|
|4.17|Fluid-driven fracture propagation test. Geometry and boundary condi-||
||tions. Dimensions in _m_.<br>. . . . . . . . . . . . . . . . . . . . . . . . . .|. 127|
|4.18|Fluid-driven fracture propagation test. Initial meshes. . . . . . . . . . .|. 129|
|4.19|Fluid-driven fracture propagation test. Time evolution of the pressure.|. 130|
|4.20|Fluid-driven fracture propagation test. Time evolution of the crack length.131||
|4.21|Fluid-driven fracture propagation test. Time evolution of the crack width. 132||



List of Figures 

xxi 

|4.22|Crack tracking test. Geometry and boundary conditions. Dimensions in|
|---|---|
||_m_.<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 134|
|4.23|Crack tracking test. Material incursions in front of the crack tip. . . . . . 136|
|4.24|Crack tracking test. Pore pressure and damage at _t_= 2 _s_.<br>. . . . . . . . 137|
|4.25|Crack tracking test. Fracture passing between moon shape incursions. . . 138|
|4.26|Parallel fracture propagation test. Geometry and boundary conditions.|
||Dimensions in _m_. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 139|
|4.27|Parallel fracture propagation test. Pore pressure and damage at _t_= 1_._5 _s_. 143|
|5.1|Scheme of the geothermal energy production. . . . . . . . . . . . . . . . . 150|
|5.2|MPI domain partitioning in 4 sections. . . . . . . . . . . . . . . . . . . . 152|
|5.3|Opening of the joints of Baserca’s arch dam during winter. Deformation|
||is exaggerated 200 times. . . . . . . . . . . . . . . . . . . . . . . . . . . . 153|
|5.4|Multi-delamination of a 3-layered beam.<br>. . . . . . . . . . . . . . . . . . 154|



## **List of Tables** 

|2.1|Elastic soil column problem. Material properties.<br>. . . . . . . . . . . . .|34|
|---|---|---|
|2.2|Elastic soil foundation problem. Material properties.<br>. . . . . . . . . . .|41|
|3.1|Computation of the stresses in the isotropic damage model. . . . . . . . .|57|
|3.2|Computation of tangent constitutive tensor in the isotropic damage model.|60|
|3.3|Computation of the stresses in the integral-type non-local damage model.|78|
|3.4|TPB test. Material properties for the Simo-Ju local model. . . . . . . . .|86|
|3.5|TPB test. Material properties for the Mazars non-local model. . . . . . .|86|
|3.6|FPS test. Material properties for the Simo-Ju local model. . . . . . . . .|93|
|3.7|FPS test. Material properties for the modifed von Mises non-local model.|93|
|4.1|Fluid-driven fracture propagation test. Material properties of the porous||
||domain.<br>. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|127|
|4.2|Fluid-driven fracture propagation test. Material properties of the crack. .|128|
|4.3|Crack tracking test. Material properties of the incursion. . . . . . . . . .|135|
|4.4|Parallel fracture propagation test. Material properties of the porous do-||
||main. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .|140|



List of Tables 

xxiv 

- 4.5 Parallel fracture propagation test. Material properties of the cracks. . . . 141 

- 4.6 Parallel fracture propagation test. Crack length and crack growth velocity.142 

## **List of Symbols** 

- _α_ Biot’s 

- _**δ**_ relative displacements of the interface 

- _**λ**_ rotation matrix 

- _**σ**_ total stress vector (poromechanics)/nominal stress vector (damage mechanics) 

- _**σ**[′]_ effective stress vector (poromechanics) 

- _**σ** e_ effective/undamaged stress vector (damage mechanics) 

- _**ε**_ strain vector 

- _**b**_ body force per unit mass 

- _**C**_ compressibility matrix 

- _**D**_ tangent constitutive tensor 

- _**H**_ permeability matrix 

- _**k**_ intrinsic permeability matrix 

- _**M**_ mass matrix 

- _**N** p_ matrix of shape functions for the pressure field 

- _**N** u_ matrix of shape functions for the displacement field 

List of Symbols 

xxvi 

- _**N** u,I_ matrix of shape functions for the displacements at the interface 

- _**Q**_ coupling matrix 

- _**u**_ displacement vector 

- _ε_ ˘ _eq_ non-local equivalent strain 

- _**δ**_ ˇ relative displacement of the interface in the local coordinate system 

- _δc_ critical relative displacement of the interface 

- _ϵ_ volumetric strain of the solid skeleton 

- _κ_ compressive to tensile strength ratio 

- d damage 

- d _p_ fracture propagation damage 

- _D_ distance 

- _µ_ dynamic viscosity 

- _µF_ friction coefficient of the interface 

- _ν_ Poisson’s ratio 

- _φ_ porosity 

- _ρ_ density of the solid-fluid mixture 

- _ρf_ fluid density 

- _ρs_ solid density 

- _σy_ yield stress of the interface 

- _σyc_ compressive strength 

- _σyt_ tensile strength 

- _τ_ stabilization parameter 

- ˜ _**t**_ prescribed surface tractions 

List of Symbols 

xxvii 

- _q_ ˜ _n_ prescribed normal flow rate 

- _εeq,I_ equivalent strain of the interface 

- _εeq_ local equivalent strain 

- _ϱ_ internal historical variable of the interface 

- _ϱy_ damage threshold of the interface 

- _ζ_ fluid mass content per unit volume 

- _E_ Young’s modulus 

- _G_ shear modulus 

- _Gf_ fracture energy per unit area 

- _K_ bulk modulus of the porous medium 

- _Kf_ fluid bulk modulus 

- _kn_ transversal permeability of the interface 

- _Ks_ solid bulk modulus 

- _lc_ non-local characteristic length 

- _lp_ fracture propagation length 

- _p_ pore pressure 

- _Q_ Biot’s modulus 

- _R_ residual strength 

- _r_ internal historical variable 

- _ry_ damage threshold 

- _S_ softening slope 

- _tI_ joint width 

- _tmin_ minimum joint width 

xxviii 

## List of Symbols 

_w_ integration coefficient 

_Z_ non-local weighting function 

## 1 Chapter 

## **Introduction** 

## **1.1 Motivation and presentation** 

Problems involving the flow of a fluid through a porous medium have always been the object of special interest in geomechanics. Many geomaterials contain pores and other cavities that are filled with fluid in saturated or unsaturated conditions. The pressure of such fluid and the deformation of the solid medium are closely related in a challenging problem that must be analysed in a coupled manner. 

Further complexity is added when one takes into consideration fluid-driven fracture propagation in the solid domain. Indeed, depending on the fluid pressure distribution and the external efforts acting over the domain, not only pre-existing discontinuities can open and close, but they can also propagate and conform new fractures, introducing an additional degree of coupling. 

From the classical contributions of Terzaghi [139] and Biot [18] important effort has been made in the last decades to develop competent numerical methods in poromechancis that allow analysing and understanding the complexity of the mechanisms involved in problems with a coupling between a solid skeleton and a fluid phase. Such techniques are becoming essential tools of for different fields and applications. 

In hydrogeology, for instance, these numerical methods are important to represent the water flow through porous soil when characterizing aquifers performance (storage ca- 

Introduction 

2 

pacity, water transmissivity, water content, etc.), and also to study the dispersion and transport of dissolved contaminants [15]. 

For soil mechanics, these tools are specially interesting for solving problems in which the seepage flow and pressures interact with the dynamic behaviour of the solid skeleton (Figure 1.1a). 

**==> picture [392 x 183] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Seepage through earth dam.<br>(b) Simplified scheme of an hydraulic frac-<br>turing process.<br>**----- End of picture text -----**<br>


Figure 1.1: Examples of problems with important interaction between a solid skeleton and a fluid phase. 

An important field of application is in petroleum engineering. Here, the oil-gas-soil interaction takes the leading role along with the hydraulic fracturing (Figure 1.1b) as a common technique to enhance reservoir permeability and well efficiency [3, 31]. 

Other applications for numerical methods in poromechanics have recently emerged for a variety of fields, including the underground storage of carbon dioxide [73], the geothermal energy production [159], the analysis of interstitial flow in bone tissue [44], and even the study of fractures in epithelial cell sheets [32, 86]. 

There are clearly different applications for this kind of problem, but in all of them one can always distinguish, at least, two distinct components: a solid skeleton and a fluid phase. The solid medium is usually composed by particles of different shapes and sizes enclosing little cavities that can be totally or partially filled by the fluid. 

Motivation and presentation 

3 

It seems natural that a proper analysis of these kind of problems should require taking into account both parts, i.e. solid and fluid, for the derivation of the governing equations and the definition of the boundary conditions. 

Thereby, when solving this coupled problem, it will be necessary to capture the seepage of the fluid through the porous medium, and the deformations of the solid skeleton at any point of the domain. From the numerical perspective, one of the most widespread approaches is the _displacement-pressure_ formulation, in which the main unknowns are the displacements of the solid skeleton _**u**_ and the pore pressure _p_ . 

However, in the limit of nearly incompressible pore fluid and small permeability, coupled poromechanics formulations suffer from instability problems. Finite elements exhibit locking in the pressure field and spurious oscillations appear due to the violation of the so-called Babuska-Brezzi conditions [6, 23]. The oscillations can be overcome by locally refining the mesh or by using shape functions of the displacement field one order higher than those of the pressure field. In practical applications this is, however, not feasible because of the increment in the computational cost. In this sense, stabilization methods have been found to be powerful tools to circumvent the Babuska-Brezzi conditions violation without compromising the efficiency of the computation. 

Discontinuities in porous materials such as concrete, soil and rock have a noticeable influence in the mechanical and hydraulic behaviour of a multi-phase system. They represent a jump in the displacement field as well as in the stress and strain fields. The presence of these fracture planes modifies the global material response in terms of strength and stiffness, and introduces directional preferences in the deformation of the solid and in the fluid flow. It is thereby essential to properly describe the initiation and propagation of these discontinuities with a suitable fracture mechanics approach. 

In the context of the Finite Element Method (FEM), modelling of discontinuities can be done by means of interface elements. These special elements, which can be used to define pre-existing or propagating cracks, act as joints between the discontinuities in the displacement field and can be used to represent the enhanced conductivity of the medium. 

Thereby, in the modelling of the fractured domain, it is necessary to implement a constitutive model that allows determining the initiation and propagation of new failure zones, and a model characterizing the mechanical behaviour of the fractures. In this 

Introduction 

4 

regard, the combination of a “smeared crack” model with a “discrete crack” approach permits working with pre-defined discontinuities and, at the same time, insert new joint elements when necessary. 

Concerning the propagation of fractures, suitable regularization methods and remeshing techniques become essential to minimize the mesh dependency when determining the direction of the crack growth. 

## **Kratos framework** 

This thesis has been implemented in Kratos Multiphysics [2]. Kratos is an Open-Source framework for building parallel multi-disciplinary simulation software. Modularity, extensibility and High Performance Computing (HPC) are some of its main objectives. Kratos has BSD licence and is written in C++ with extensive Python interface. 

It features a “core and applications approach”, where standard tools (databases, linear algebra, solvers, etc.) come as part of the “core” and are available as building blocks in the development of “applications”, which derive from the core and focus on the solution of particular problems of interest (see Figure 1.2). We note here that all the implementations of this thesis gave as a result the “Poromechanics Application” in Kratos. 

Figure 1.2: Kratos’ “core and applications approach”. 

Objectives 

5 

## **1.2 Objectives** 

The main objective of this thesis is to advance in the development of finite element methods for the analysis of coupled solid-pore fluid problems in which the initiation and propagation of cracks represents a key issue and is strongly influenced by the interaction between the different media. To this end, the following goals were set for the three-years PhD: 

- Derivation and testing of a stable FEM formulation for coupled solid-pore fluid problems under fully saturated conditions for 2D and 3D analysis. 

- Implementation and validation of a robust continuum damage model for the analysis of degradation of quasi-brittle materials. 

- Development of special interface elements compatible with the implemented formulation, and valid for the definition of pre-existing fractures in the porous media as well as propagating cracks. 

- Development of a crack tracking technique based on the analysis of the damage pattern around the crack tip. 

- Implementation of an efficient technique to insert new interface elements in the domain and remesh the model every time the crack propagates. 

- Application of the previous implementations for solving a problem of interest in engineering. An example of a problem where all the components are required is the analysis of hydraulic fracture processes. 

## **1.3 State of the art** 

## **1.3.1 Coupled solid-pore fluid formulation** 

The first work regarding the coupling between a fluid and a porous medium can be attributed to Terzaghi, who proposed a one-dimensional theory of consolidation in 1943 

Introduction 

6 

[139]. The generalization to three-dimensional cases was presented by Biot who also described methods for determining the elastic coefficients of the theory [18–20]. Although these first works were restricted to linear elastic materials, they have been the basis for a lot of subsequent works in geophysics, soil and rock mechanics. 

The transition between Biot’s theoretical formulation and the first successful numerical results is found in Ghaboussi and Wilson’s work, in which the assumption of incompressibility of fluid and solid particles was relaxed [61]. Zienkiewicz and Bettes [149] proposed a formulation based only on the displacement of solid and the pore pressure, which is suitable for the low-frequency situations such as earthquake problems. 

De Boer and Kowalski [46] developed a plasticity theory for saturated porous media and, due to the increasing interest in non-linear applications, Zienkiewicz and co-workers expanded the theory to a generalized incremental form for non-linear material behaviour and large strain effects for liquefaction analysis of soil structures [148, 156]. 

In order to overcome the violation of Babuska-Brezzi conditions, several stabilization techniques can be found in the literature in the context of fluid mechanics and solid mechanics, and have been also extended to poromechanics problems. 

Some of these techniques are based on the pioneer work of Chorin [42] and use special time stepping schemes to generate stabilization terms [69, 111, 152]. 

Other types of approaches have also achieved good stabilized results by modifying the fluid mass balance equation with the addition of a characteristic term that needs to be calibrated [24, 40, 70, 153]. Here one can also include the methods based on the concept of Polynomial-Pressure-Projections, in which the stabilizing terms use element-local projections of the pressure field to counteract the inherent instabilities [41, 135, 136, 142]. 

There is an additional group of stabilization techniques that make use of the Finite Increment Calculus method (FIC) developed by Oñate and co-workers. In this approach, stabilization terms are derived by expanding the residual of the governing equations in a Taylor series. Previously stabilized FEM formulations for quasi and fully incompressible fluids and solids were based on the first order form of the FIC balance equation in space [102, 103, 106, 107], while recent works are based on the second-order FIC form of the mass balance equation [104, 105]. 

State of the art 

7 

## **1.3.2 Fracture modelling** 

In the last decades various authors have developed numerical methods for the initiation and propagation of fractures in the context of porous media. 

One approach that has obtained notable attention in the last years is the extended finite element method (XFEM) developed for crack analysis [14, 81, 83, 93, 120, 146]. The XFEM handles discontinuities in an approximation space with jump functions, i.e. local enrichment functions, and so there is no need to represent cracks in a mesh. This considerably reduces the cost of remeshing during crack growth, but in return, the method demands a higher computational cost in terms of number of degrees of freedom and numerical integration. 

The other obvious approach that has shown significant impact includes the group of numerical techniques purely based on the finite element method (FEM). In this group, numerous methods can be found in the literature, but probably two main subgroups can be distinguished: the “smeared crack” approaches, continuum based methods in which the influence of developing fractures is incorporated into the constitutive stress-strain law [28, 66, 85, 119], and the “discrete crack” models, in which each single discontinuity is represented explicitly [34, 56, 60, 67]. 

For the analysis of material degradation and crack initiation in porous media, most of the techniques developed for continuum damage mechanics are applicable. 

The concept of damage was introduced by Kachanov in 1958 in the context of creep rupture in metals [75]. Rabotnov [118] gave the problem physical meaning by suggesting that the reduction of the sectional area was measured by means of the damage parameter. 

The thermodynamic formalism involved in the irreversible process of damage was developed by Lemaitre and Chaboche [82], and other important contributions to our knowledge were given by: Mazars and Pijaudier-Cabot [91], Simo and Ju [130], Oller et al. [101], Oliver et al. [99, 100], and Cervera et al. [35, 36], to name but a few. 

If the damage parameter depends only on the strain state at the point under consideration, numerical simulations exhibit a pathological mesh dependence and the energy consumed by the system tends to zero as the mesh is refined. Classical constitutive models require an extension in the form of a characteristic length to properly model the 

Introduction 

8 

thickness of localized zones. 

Such extension can be introduced by means of second gradient models [38], micro-polar [95], strain gradient [145], viscous [84] and non-local terms [74]. 

Non-local elasticity was developed by Eringen [54] who later extended it to non-local elasto-plasticity [53]. Subsequently, it was found that certain non-local formulations could act as efficient localization limiters with a regularizing effect on problems with strain localization [114]. 

Modelling of discontinuities between two bodies in contact has been performed via contact elements in solid mechanics [16, 62, 125], and the idea has also been used in fracture mechanics to model cohesive forces by means of cohesive elements [109, 144]. 

Since Goodman _et al._ proposed the “zero-thickness” interface element to describe the mechanical behaviour of pre-existing joints in rock masses [63], many authors have developed strategies to adapt that element for the solution of fracture processes in solidpore fluid coupled problems. 

Similarly as in the formulation of the continuum solid-pore fluid mixture, two equilibrium equations usually govern the coupled behaviour of the fractures. One equation deals with the mechanical behaviour of the crack, whereas the other equation describes the balance of mass within the fracture. 

**==> picture [434 x 83] intentionally omitted <==**

**----- Start of picture text -----**<br>
4 3<br>4 3<br>5 6<br>1 2<br>1 2<br>(a) Single-noded zero- 1 2<br>thickness interface element. (b) Double-noded zero- (c) Triple-noded zero-<br>thickness interface element.<br>thickness interface element.<br>**----- End of picture text -----**<br>


Figure 1.3: Two-dimensional zero-thickness interface elements. 

Concerning the modelling of the fluid flow through the discontinuity, three different types of zero-thickness interface elements can be found in the literature: single, double and triple noded elements (Figure 1.3). The single-noded elements are the simplest ones and represent pipe elements that are superimposed onto the standard continuum element edges [5]. With this elements the hydraulic head is uniform across the discontinuity and thus only the longitudinal conductivity is taken into account. The triple-noded elements 

State of the art 

9 

allow including the influence of a transversal conductivity through the discontinuity [64]. The two nodes at each side of the interface represent the potentials in the pore pressure between the two sides of the discontinuity. The third node, placed at the middle of the interface, stores the average potential of the longitudinal fluid through the fracture. Finally, the double-noded elements take into account both types of conductivity but the middle node variable is eliminated in terms of the variables of the external nodes. Ng and Small [98] used this double-noded zero-thickness interface element to model flow problems with pre-existing discontinuities but with no hydraulic potential drop between the two interface walls. Segura and Carol [128] introduced the transversal conductivity for zero-thickness elements to account for the exchange of fluid between the discontinuity and the porous material. 

For analysing the mechanical behaviour of fractures, we find essentially two different approaches: those based on linear elastic fracture mechanics (LEFM), and those based on non-linear fracture mechanics (NLFM). LEFM techniques were first proposed to solve fracture propagation problems by means of remeshing without considering a fracture process zone (FPZ) before the crack tip. This approach is applicable in large structures where the size of the FPZ is negligible. However, for quasi-brittle analyses, the consideration of a non-linear fracture process zone where the energy is dissipated before it completely fails was found to be essential. In those cases the NLFM technique is usually applied and a softening law relates the cohesive stress to the crack opening in the FPZ. The first procedure based on the cohesive fracture model was originally introduced by Barenblatt [8, 9] for brittle materials and by Dugdale [50] for plastic materials. Hillerborg _et al._ [68] developed the first fictitious crack model for Mode I fracture. It was extended later for the mixed mode fracture, from which Camacho and Ortiz [27] proposed a suitable fracture criterion that is widely used in the literature. 

A crucial part in the modelling of fracture propagation is the criterion for determining the direction of the crack growth. Some criteria are based on the local evaluation of the stress field at the crack tip, e.g. the maximum circumferential stress criterion [52]. Others are based on the energy distribution at the fractured zone, such as the minimum strain energy density criterion [129] or the maximal strain energy release rate criterion [72]. There is also another group of crack growth criteria based on continuum damage mechanics [140]. 

Schrefler _et al._ [126] and Secchi _et al._ [127] implemented the double-noded interface ele- 

Introduction 

10 

ment with cohesive fracture model, working together with an adaptive mesh refinement technique in porous materials. Also, Khoei _et al._ implemented a model of cohesive crack growth using an adaptive mesh refinement procedure for two-dimensional problems [79]. In that work, a bilinear cohesive law originally proposed by Espinosa and Zavattieri [55] is adopted, and a modified superconvergent patch recovery (SPR) technique is used for the error estimation. 

More recently, Moës _et al._ presented an alternative technique combining fracture mechanics with level set procedures [33, 94]. 

## **1.4 Organization** 

The thesis is organized as follows. 

In Chapter 2 the displacement-pore pressure formulation is described, starting from the governing equations of the problem and the fully coupled system of equations to be solved. Following that, the instability produced under undrained-incompressible conditions is presented and the stable formulation based on the second-order FIC form of the mass balance equation in space is thoroughly detailed. Two academic examples are solved at the end of the chapter comparing the FIC-FEM elements of equal order interpolation to stable elements with a higher order interpolation for the displacement 

Chapter 3 starts introducing the essential concepts on continuum damage mechanics. Details are given on the fundamental concepts of the isotropic damage theory, including the equivalent strain forms and damage evolution laws that have been implemented in this work. Next, the regularization techniques developed to overcome the problems associated to strain localization are illustrated. Special attention is given to the integraltype non-local damage model, highlighting the most relevant aspects of its numerical implementation. Two examples are shown at the end, where we assess the robustness of the implemented damage models against changes in the spatial discretization. 

In Chapter 4 the formulation of the developed quasi-zero-thickness interface elements is detailed, explaining the mechanical behaviour of the fracture and the model governing the fluid flowing in it. After that, we present the proposed fracture propagation technique combining the non-local damage model with the interface elements. Finally, two 

Related publications 

11 

plane-strain examples are solved to test the accuracy of the crack propagation methodology, and one additional case is included to show the performance of the generalized 3D formulation. 

Chapter 5 summarizes the most relevant conclusions of this work, pointing out its main contributions. In the end, the new perspectives to future works are briefly commented. 

## **1.5 Related publications** 

The following publications have emanated from the work carried out in this doctoral thesis. The list is ordered chronologically. 

- I. de Pouplana and E. Oñate. Combination of a non-local damage model for quasibrittle materials with a mesh-adaptive finite element technique. _Finite Element in Analysis and Design_ , vol. 112, pp. 26-39, 2016. DOI: 10.1016/j.finel.2015.12.011 

- I. de Pouplana and E. Oñate. A FIC-based stabilized mixed finite element method with equal order interpolation for solid-pore fluid interaction problems. _International Journal for Numerical and Analytical Methods in Geomechanics_ , vol. 41, pp. 110-134, 2016. DOI: 10.1002/nag.2550 

- I. de Pouplana, L. Gracia, F. Salazar and E. Oñate. Cracking of a concrete arch dam due to seasonal temperature variations. _Theme A -_ 14 _[th] International Benchmark Workshop on the Numerical Analysis of Dams_ , 2017. 

- L. Gracia, I. de Pouplana, F. Salazar and E. Oñate. Static and seismic analysis of an arch-gravity dam. _Theme B -_ 14 _[th] International Benchmark Workshop on the Numerical Analysis of Dams_ , 2017. 

- I. de Pouplana and E. Oñate. Finite element modelling of fracture propagation in saturated media using quasi-zero-thickness interface elements. Article sent to _Computers and Geotechnics_ , 2017. 

- E. Oñate and I. de Pouplana. A methodology for analysis of delamination in composites using joint elements. _In preparation_ . 

## 2 Chapter 

## **Coupled Solid-Pore Fluid Interaction Problem** 

## **2.1 Introduction** 

Sometimes two or more physical systems interact with each other in such a way that it is not possible to obtain a solution independent of any of the systems without considering all the systems at the same time. It is said that those systems belong to a coupled problem, with the coupling being stronger or weaker depending on the level of interaction between each system. 

Our problem of interest belongs to the type of coupled problems in which the different domains of the physical systems overlap. In such cases, the coupling occurs through the governing differential equations describing the different physical phenomena. In particular, we are going to study solid-pore fluid interaction problems in which a unique mesh represents a porous medium. 

The problem falls in the field of poromechanics [43], a branch of physics that studies porous materials whose mechanical behaviour is significantly influenced by the pore 

A porous medium is composed of a matrix, containing solid components of different shape and size, and a porous space, that can be totally or partially filled by a fluid. 

Coupled Solid-Pore Fluid Interaction Problem 

14 

**==> picture [446 x 219] intentionally omitted <==**

**----- Start of picture text -----**<br>
Thereby, any infinitesimal volume of a porous medium can be treated as the superim-<br>position of two material particles: the skeleton particle and the fluid particle (Figure<br>2.1).<br>rK Xs rX ¥<br>=<br>Py 6.9! y 6,1 +<br>SESS aa<br>Porous Medium Solid Skeleton Fluid Medium<br>Figure 2.1: Infinitesimal volume of porous medium as a combination of two phases: a<br>solid skeleton and a medium.<br>**----- End of picture text -----**<br>


A continuous description of such an heterogeneous medium requires taking into account the mechanical coupling between the solid skeleton and the fluid phase when formulating the governing equations of the problem. 

In the last decades important effort has been made in the development of powerful numerical methods for poromechanics that allow analysing porous media considering the mechanical coupling between the solid and the fluid phase. 

A one-dimensional theory of consolidation was first proposed by Terzaghi in 1943 [139] and then Biot extended it to three-dimensional cases [18–20]. Although these first works were restricted to linear elastic materials, they have been the basis for a lot of subsequent works in geophysics, soil and rock mechanics. 

The first numerical solution of Biot’s formulation was obtained by Ghaboussi and Wilson [61] and the work was further developed by Zienkiewicz _et al._ [149, 151]. Later on, due to the increasing interest in non-linear applications, Zienkiewicz and co-workers expanded the theory to a generalized incremental form for non-linear materials and large deformation problems [148, 156]. 

The mathematical formulation of solid skeleton and pore fluid interaction presented here is based on the model proposed by Zienkiewicz _et al._ [151]. The problem was originally formulated for fully saturated conditions in terms of the solid matrix displacement _ui_ , the mean fluid velocity relative to the solid phase _wi_ and the fluid pore pressure _p_ , but 

Introduction 

15 

in many geo-mechanical problems with no high-frequency phenomena involved, the fluid relative velocity _wi_ can be neglected and so the equations can be simplified to work with only two main variables: _ui_ and _p_ [150]. 

In the limit of nearly incompressible pore fluid and small permeability, the coupled poromechanics formulations suffer from instability problems. Finite elements exhibit locking in the pressure field and spurious oscillations in the numerical solution for the pressure appear due to the violation of the so-called Babuska-Brezzi conditions [6, 23]. The oscillations can be overcome by locally refining the mesh and by using shape functions of the displacement field one order higher than those of the pressure field. In practical applications this is, however, not the best approach because of the increment in the computational cost. In this sense, stabilization methods have been found to be powerful tools to circumvent the Babuska-Brezzi conditions violation without compromising the efficiency of the computation. 

Several stabilization techniques have been investigated in the past years in the context of computational fluid dynamics [24, 42, 70, 102, 103, 152] and (incompressible) solid mechanics [40, 106, 107, 110, 155], and have been also extended to poromechanics problems [69, 88, 111, 142, 153]. 

Although each stabilization approach has its differential aspects, they can be classified in three main categories. 

The first category comprises those techniques in which special time stepping schemes are applied in order to generate stabilization terms. Probably the earliest work in this category is due to Chorin [42] who proposed a technique to deal with incompressible fluid problems which is now referred as the fractional step method or the operatorsplitting method. Such a staggered time stepping algorithm has been found to provide stabilization in the steady-state when used in a convenient form [42, 69, 111, 152]. 

The second type of techniques are more direct stabilization methods based on the perturbation of the fluid mass conservation equation. Instead of using special time-stepping algorithms that give additional terms in the steady-state approximation, the fluid balance equation is modified by adding a stabilization term multiplied by a parameter that needs to be carefully calibrated [24, 40, 70, 153]. This group also includes the methods based on the concept of Polynomial-Pressure-Projections, in which the additional stabilizing terms use element-local projections of the pressure field to counteract the inherent 

Coupled Solid-Pore Fluid Interaction Problem 

16 

instabilities [41, 135, 136, 142]. 

The third category considered here is the Finite Increment Calculus method (FIC) [102–107], which is the approach also adopted in the present work. The FIC technique is based on expressing the equations of balance of mass and momentum in a space/time domain of finite size, and retaining higher order terms in the Taylor series expansion typically used for expressing the change in the transported variables within the balance domain. Apart from the standard terms of infinitesimal theory, the FIC form of the balance equations contains derivatives of the classical differential equations multiplied by characteristic distances in space and/or time. 

Previously stabilized FIC-FEM formulations were based on the first-order form of the FIC balance equation in space and can be found in the literature for quasi and fully incompressible fluids and solids [102, 103, 106, 107], and even for 1D and 2D poromechancis problems [116]. In this work, the more recent second-order FIC form of the mass balance equation [104, 105] has been adopted as the basis for deriving the stabilized FIC-FEM formulation for 2D and 3D poromechanics problems. This new formulation takes advantage of the second order derivatives terms to provide a simpler procedure for obtaining a residual-based stable form of the mass balance equation suitable for finite element analysis. Proof of the good results of this formulation regarding convergence and mass conservation are given in [104]. 

This chapter is organized as follows. First, the finite element formulation is introduced, starting from the governing equations of the problem and the fully coupled system of equations to be solved. After that, the instability problem is stated and the stabilization by means of the second-order FIC form of the mass balance equation in space is thoroughly detailed. Finally, two examples are solved in order to test the implemented FIC-stabilized element and compare it to stable elements with a higher order interpolation for the displacement field. 

## **2.2 Finite element formulation** 

In the theory of porous media, the effective stress is an essential concept for defining the stress state because it is responsible for the major deformations and rupture of the solid skeleton. Let _σij_ be the total stress (positive in tension) acting on the total area of the 

Finite element formulation 

17 

soil and the pores, and _p_ the pressure of the fluid in the pores (positive in compression), the stress is as: 

**==> picture [265 x 15] intentionally omitted <==**

where _δij_ is the Kronecker delta and _α_ is Biot’s coefficient [21]: 

**==> picture [264 x 27] intentionally omitted <==**

with _Ks_ being the bulk modulus of the solid phase and _K_ the bulk modulus of the porous medium: 

**==> picture [263 x 29] intentionally omitted <==**

where _E_ is the Young’s modulus and _ν_ is the Poisson’s ratio. 

In order to account for the coupling between the solid and fluid phases, the behaviour of a saturated porous medium is governed by the combination of two equations: the balance of momentum for the mixture solid-fluid, and the mass balance for the pore 

_Balance of momentum for the mixture solid-fluid_ 

**==> picture [264 x 30] intentionally omitted <==**

where _bi_ is the body force per unit mass, _u_ ¨ _i_ is the acceleration of the solid skeleton and _ρ_ is the density of the solid-fluid mixture _ρ_ = _φρf_ + (1 _− φ_ ) _ρs_ , being _φ_ the porosity of the soil, _ρf_ the density of the fluid and _ρs_ the density of the solid. 

Let us introduce the general definition for the effective stress (3.1) into the balance of momentum equation (2.4). This gives 

**==> picture [307 x 29] intentionally omitted <==**

_Mass balance for the pore fluid_ 

**==> picture [258 x 28] intentionally omitted <==**

Coupled Solid-Pore Fluid Interaction Problem 

18 

where _ζ_ is the fluid mass content per unit volume and _qi_ represents the flow rate. It should be noted that Equation (2.6) is in a linearized form of the general mass balance expression as the fluid density variation effect has been ignored. 

The mass balance equation is modified as follows. Let us first consider the constitutive equation presented in [18] relating the fluid mass content per unit volume _ζ_ with the volumetric strain of the solid skeleton _ϵ_ and the pore pressure _p_ . 

**==> picture [254 x 25] intentionally omitted <==**

where _Q_ is a combined compressibility of the fluid-solid phases, also called Biot’s modulus [21] 

**==> picture [268 x 29] intentionally omitted <==**

being _Kf_ the bulk modulus of the fluid phase. 

Let us now make use of the generalized form of Darcy’s law 

**==> picture [290 x 30] intentionally omitted <==**

where _µ_ is the dynamic viscosity of the fluid and _kij_ is the intrinsic permeability matrix of the porous medium. 

Taking into account the relations from equations (2.7) and (2.9), the mass balance equation (2.6) can be rewritten as 

**==> picture [331 x 30] intentionally omitted <==**

Equations (2.5) and (2.10) have to be supplemented by a constitutive law for the solid skeleton. In general, for any non-linear material we may consider an incremental form with a tangent constitutive tensor _Dijkl_ as 

with 

**==> picture [279 x 73] intentionally omitted <==**

Finite element formulation 

19 

The boundary conditions for this problem are specified as: 

**==> picture [238 x 13] intentionally omitted <==**

**==> picture [304 x 16] intentionally omitted <==**

where _u_ ˜ _i_ and _t_[˜] _i_ are the prescribed displacements and surface tractions, respectively. 

**==> picture [201 x 14] intentionally omitted <==**

**==> picture [332 x 30] intentionally omitted <==**

where _p_ ˜ and _q_ ˜ _n_ are the prescribed pore pressure and the normal flow rate, respectively. 

In order to express the resultant system of equations in a more compact manner, we rewrite equations (2.5) and (2.10) in matrix form using Voigt notation as follows. 

- Balance of momentum 

**==> picture [287 x 14] intentionally omitted <==**

_•_ Mass balance 

**==> picture [321 x 30] intentionally omitted <==**

For a general 3D case, we have 

**==> picture [332 x 109] intentionally omitted <==**

Coupled Solid-Pore Fluid Interaction Problem 

20 

**==> picture [328 x 25] intentionally omitted <==**

**==> picture [297 x 25] intentionally omitted <==**

The above system of partial differential equations (2.17) and (2.18) can be discretized by interpolating the displacement and pressure fields as: _**u**_ = _**N** u_ _**u**_ ¯ and _p_ = _**N** p_ _**p**_ ¯ where ¯( _·_ ) denotes nodal values and 

**==> picture [363 x 51] intentionally omitted <==**

**==> picture [280 x 23] intentionally omitted <==**

with _N_[ˆ] _i_ and _Ni_ being, respectively, the shape functions of the displacement and pressure interpolations, which do not necessarily need to coincide. 

Following the standard Galerkin technique [17, 58], we left multiply Equation (2.17) by test functions _**N**[T] u_[and][Equation][(2.18)][by][test][functions] _**[N]**[ T] p_[and][integrate][them][over] the spatial domain Ω to obtain the following set of ordinary differential equations 

- Balance of momentum 

**==> picture [315 x 28] intentionally omitted <==**

- Mass balance 

**==> picture [295 x 16] intentionally omitted <==**

where _**r** u_ and _**r** p_ are the residual vectors for the momentum and the mass balance equation, and 

**==> picture [346 x 28] intentionally omitted <==**

Finite element formulation 

21 

**==> picture [366 x 112] intentionally omitted <==**

with _**B**_ = _**SN** u_ . 

The time derivatives of _**u**_ and _p_ are approximated using the Generalized Newmark scheme. Thus, for a new time step _n_ + 1, we use the GN22 scheme for displacements [97]: 

**==> picture [352 x 30] intentionally omitted <==**

**==> picture [319 x 13] intentionally omitted <==**

and the GN11 scheme for pressure: 

**==> picture [446 x 43] intentionally omitted <==**

The residual vectors for the momentum and mass balance equations can be rewritten as 

**==> picture [413 x 77] intentionally omitted <==**

where _**u**_ ¯ _n_ +1 and _**p**_ ¯ _n_ +1 are the nodal unknowns at time _n_ + 1, and _**u**_ ¯[ˆ] _n_ , _**u**_ ¯[ˆ˙] _n_ and _**p**_ ¯[ˆ˙] _n_ stand for values that are computed from known parameters at time _n_ as: 

**==> picture [314 x 30] intentionally omitted <==**

**==> picture [311 x 25] intentionally omitted <==**

Coupled Solid-Pore Fluid Interaction Problem 

22 

**==> picture [283 x 26] intentionally omitted <==**

The Newton-Raphson method is used to solve the problem iteratively. Thus, at a time step _n_ + 1 and the iteration _i_ + 1 we have: 

**==> picture [368 x 127] intentionally omitted <==**

where _**K**_ is the tangent stiffness matrix: 

**==> picture [322 x 28] intentionally omitted <==**

Note that the system of equations (2.38) can be made symmetric by multiplying Equation (2.33) by _−β_ ∆ _t/γ_ . 

## **2.3 Stabilization using Finite Increment Calculus (FIC)** 

## **2.3.1 Introduction to the FIC-stabilization procedure** 

For a quasi-static problem the term involving the mass matrix is omitted. Also, in the undrained-incompressible limit, i.e. when _**k** →_ **0** and _Q →∞_ , the matrices _**C**_ and _**H**_ vanish and the system of equations (2.38) becomes 

**==> picture [325 x 40] intentionally omitted <==**

The resultant system matrix is almost identical to that frequently encountered in the solution of incompressible elasticity or incompressible fluid mechanics problems [70, 106]. In such cases, the spaces used to approximate the displacement (or the velocity) and 

Stabilization using Finite Increment Calculus (FIC) 

23 

pressure fields must fulfill the Babuska-Brezzi conditions [6, 23] or the Zienkiewicz-Taylor patch test [154, 157] in order to avoid spurious oscillations and locking in the pressure 

With the equations presented so far, these conditions can be accomplished using shape functions for the displacement field one order higher than those for the pressure field. Some examples of “stable” 2D elements are depicted in Figure 2.2. 

**==> picture [429 x 115] intentionally omitted <==**

**----- Start of picture text -----**<br>
u p<br>u p<br>(a) Quadratic displacement and linear (b) Biquadratic displacement and bilinear<br>pressure.<br>pressure.<br>**----- End of picture text -----**<br>


Figure 2.2: Some 2D elements that fulfil the Babuska-Brezzi conditions. 

However, low order finite elements with equal order interpolation for the displacements and the pressure are very attractive due to their simplicity and efficiency, and so stabilization techniques must be applied if one aims at solving large scale 3D computations. 

The stabilization approach implemented in this work relies on the Finite Increment Calculus (FIC) method [102–107] and affects only the continuity equation with the balance of momentum remaining unchanged. 

For the sake of clarity, let us redefine the mass balance equation (2.10) as: 

**==> picture [348 x 30] intentionally omitted <==**

As stated in the introduction, this technique is based on the second-order FIC form of the mass balance equation in space for a quasi-incompressible fluid. This can be written as [104]: 

**==> picture [321 x 30] intentionally omitted <==**

where _ns_ are the number of space dimensions (i.e. _ns_ = 3 for 3D problems). 

In a 2D case Equation (2.42) is obtained by expressing the balance of mass in a rectangular domain of finite size with dimensions _h_ 1 _× h_ 2, where _hi_ are arbitrary distances, 

Coupled Solid-Pore Fluid Interaction Problem 

24 

and retaining up to third-order terms in the Taylor series expansions used for expressing the change of mass within the balance domain [105]. The derivation of Equation (2.42) for a 1D problem is shown in [104]. 

The FIC term in Equation (2.42) plays the role of a space stabilization term, in which _hi_ are space dimensions related to the characteristic element dimensions. Note that for _hi →_ 0 the standard form of the mass balance equation (2.41), as given by the infinitesimal theory, is recovered. 

Equation (2.42) can be interpreted as a particular class of residual-based stabilization methods applied to the strong form of the mass conservation equation. This ensures the consistency of the stabilization method in the discretized problem. 

## **2.3.2 Modified FIC-stabilized form of the mass balance equation** 

Let us first expand Equation (2.42) using Equation (2.41) as 

**==> picture [385 x 63] intentionally omitted <==**

Next we split the stress tensor into its deviatoric and volumetric components as 

**==> picture [262 x 14] intentionally omitted <==**

where _sij_ is the deviatoric stress tensor and _σ_ is the hydrostatic pressure defined as _σ_ = _σii/_ 3. 

In the same way, we split the strain tensor into its deviatoric and volumetric parts as 

**==> picture [265 x 26] intentionally omitted <==**

with _eij_ being the deviatoric strain tensor and _ϵ_ = _εii_ the volumetric dilation. 

Now let us write the isotropic linear elastic constitutive equations: 

**==> picture [252 x 13] intentionally omitted <==**

Stabilization using Finite Increment Calculus (FIC) 

25 

where _G_ is the shear modulus, 

**==> picture [254 x 13] intentionally omitted <==**

where _σ[′]_ = _σii[′][/]_[3][is][the][mean][effective][stress][and] _[α]_[is][the][Biot’s][coefficient][defined][in] Equation (2.2). 

Combining Equations (2.44), (2.45), (2.46) and (2.47), we express the stress tensor as 

**==> picture [315 x 26] intentionally omitted <==**

At this point, we substitute the above expression in the balance of momentum (2.4). This gives 

**==> picture [376 x 30] intentionally omitted <==**

From this point forward, in the derivation of equations (2.43), (2.49) and in the following, we neglect the space and time changes of _α_ , _Q_ , _G_ and _ρ_ in the derivatives. 

If we now apply the divergence operator to both sides of Equation (2.49), differentiate it with respect to time, and rearrange the terms we can obtain 

**==> picture [382 x 37] intentionally omitted <==**

In the previous equation the term _ρ ∂x[∂][b]_[˙] _[i] i_[can][be][neglected][and][the][term][involving][the] third derivative of the displacements with respect to time can be obtained from the mass balance equation. Thus if we differentiate Equation (2.10) twice with respect to time assuming that the time derivative of the body force per unit mass is negligible, and we take into account the definition of the volumetric strain _ϵ_ = _[∂u] ∂xi[i]_[,][we][obtain] 

**==> picture [335 x 30] intentionally omitted <==**

Introducing (2.51) into (2.50) gives 

**==> picture [327 x 63] intentionally omitted <==**

Coupled Solid-Pore Fluid Interaction Problem 

26 

At this point, we can already substitute the above relation in Equation (2.43) giving 

**==> picture [395 x 96] intentionally omitted <==**

The above equation can be simplified if we take into account that the last term involves the fourth order spatial derivative of the pore pressure and the third order spatial derivative of the body force per unit mass. In practical problems these terms are either zero (for constant body forces) or much smaller than the others. Hence, these terms will be omitted for the rest of this work. Furthermore, numerical results have shown that the terms _∂x∂ j_ � _µ_ 1 _[k][jk] ∂x∂p_ ¨ _k_ � and _Q_ 1 _DDtp_ ¨[can][be][neglected][without][loss][of][accuracy.][Thereby,] Equation (2.53) is written in the simpler form as 

**==> picture [372 x 63] intentionally omitted <==**

In the following we will assume _hi_ = _h_ where _h_ is a characteristic length related to a typical average dimension of each element in the mesh. After rearranging terms, Equation (2.54) can be rewritten as 

**==> picture [348 x 63] intentionally omitted <==**

where _τ_ is a stabilization parameter given by 

**==> picture [245 x 27] intentionally omitted <==**

The form of the stabilization parameter in Equation (2.56) resembles that typically used in other stabilized procedures. We note that this term has naturally formed from the FIC formulation. 

Stabilization using Finite Increment Calculus (FIC) 

27 

In the examples solved in this work we have chosen _h_ = _l[e]_ , where _l[e]_ is a characteristic element length that, for 2D problems, is taken as the diameter of a circle with the area of the element, while for 3D problems, it is the diameter of a sphere with the volume of the element. In essence 

**==> picture [272 x 30] intentionally omitted <==**

**==> picture [272 x 29] intentionally omitted <==**

where _A[e]_ and _V[e]_ represent the area and the volume of the element, respectively. 

The presence of the characteristic element length _l[e]_ in Equation (2.56) helps us control the diffusion of the stabilized solution, avoiding over-diffusive numerical results provided that fine enough discretizations are used. This fact is clearly shown in the example of Section 2.4.2. 

## **2.3.3 Variational form of the FIC-stabilized mass balance equation** 

The variational expression of the FIC-stabilized mass balance equation is obtained by multiplying Equation (2.55) by arbitrary test functions _v_ and integrating over the domain Ω to give 

**==> picture [383 x 63] intentionally omitted <==**

Integrating by parts the last two integrals of Equation (2.59) and applying the divergence 

Coupled Solid-Pore Fluid Interaction Problem 

28 

theorem yields 

**==> picture [379 x 128] intentionally omitted <==**

where _ni_ are the components of the unit normal vector to the external boundary Γ of Ω. 

Introducing the boundary condition (2.16) and rearranging terms we have 

**==> picture [413 x 125] intentionally omitted <==**

Using the stress decomposition in (2.48), we can rewrite the boundary condition (2.14) as 

**==> picture [315 x 26] intentionally omitted <==**

Next if we apply the divergence operator to both sides of Equation (2.62), differentiate it with respect to time and rearrange terms we have 

**==> picture [344 x 34] intentionally omitted <==**

In all typical problems the divergence of the traction vector is zero and so it is the term _∂t_[˙˜] _i_[Moreover,][numerical][tests][have][shown][that][the][results][are][not][affected][by][the][term] _∂xi_[.] 

Stabilization using Finite Increment Calculus (FIC) 

29 

2 _G ∂ϵ_ ˙[Consequently,][the][stabilized][mass][balance][equation][is][rewritten][as] 3 _∂xi[n][i]_[.] 

**==> picture [362 x 125] intentionally omitted <==**

where _∂n[∂][p]_[˙][is][the][derivative][of] _[p]_[˙][in][the][direction][of][the][normal][to][the][external][boundary.] This derivative can be approximated as shown in Figure 2.3. 

**==> picture [293 x 133] intentionally omitted <==**

Figure 2.3: Computation of the term _∂n∂p_ ˙[at][the][side] _[ij]_[of][a][triangle][adjacent][to][the] external boundary Γ. 

Using the previous argument, the stabilized FIC form for the mass balance equation 

Coupled Solid-Pore Fluid Interaction Problem 

30 

- (2.64) finally becomes 

**==> picture [368 x 131] intentionally omitted <==**

where _hn_ is an arbitrary distance in the normal direction to the boundary. In our work _hn_ has been taken as the characteristic length _l[e]_ of the element adjacent to the external boundary. 

## **2.3.4 Discretized form of the momentum and FIC-stabilized mass balance equations** 

Combining the discretized form of Equation (2.65) along with the residual obtained from the original balance of momentum equation (2.24), the system of equations to be solved reads 

_•_ Balance of momentum 

**==> picture [314 x 28] intentionally omitted <==**

_•_ Mass balance 

**==> picture [354 x 15] intentionally omitted <==**

with the stabilizing terms 

**==> picture [386 x 77] intentionally omitted <==**

Stabilization using Finite Increment Calculus (FIC) 

31 

**==> picture [291 x 28] intentionally omitted <==**

where _**S**_[ˆ] = _**m∇**[T]_ and the matrix _**I** v_ appears due to the difference between the strain tensor and the strain vector expressed in Voigt notation. For a general 3D case 

**==> picture [293 x 101] intentionally omitted <==**

Note that since the test functions _v_ are zero on the Dirichlet boundaries Γ _p_ , matrix _**T** b_ must be only computed on those external boundaries where no pressure Dirichlet conditions are imposed. 

Using the stabilized residual vector (2.67), the system of equations for a linear elastic material, quasi-static and undrained-incompressible limit case becomes 

with 

**==> picture [376 x 96] intentionally omitted <==**

**==> picture [305 x 28] intentionally omitted <==**

Note that the diagonal of the iteration matrix is now different from zero. Matrices _**R**_ and _**L**_ emerging from the stabilized FIC formulation are essential for ensuring the consistency of the residual-based stabilization method used. 

The consistency of residual does not directly affect the diffusion of the numerical solution, but it helps improving the condition number of the system matrix and facilitates the convergence of the Newton-Raphson scheme. 

We highlight that the FIC-FEM formulation presented is applicable to equal order interpolation approximations of the displacements and the pressure of any degree. In this 

Coupled Solid-Pore Fluid Interaction Problem 

32 

work we have implemented the FIC-FEM formulation using 3-noded triangles and 4- noded quadrilaterals for 2D problems, and 4-noded tetrahedra and 8-noded hexahedra for 3D problems, with equal order interpolation for the displacements and pressure. 

## **2.4 Examples** 

This section includes two academic examples designed to assess the performance of the stabilized FIC-FEM formulation in situations near the undrained-incompressible limit. The first example is an elastic soil column subjected to a surface load. The second problem is an elastic soil foundation also subjected to a surface load. 

The soil column problem is analysed in a 2D framework under plane strain conditions, whereas the soil foundation case is tackled as a 3D problem. The problems have been solved with simple elements with equal order interpolation for the displacements and the pore pressure, using the direct (non-stabilized) formulation and the FIC-stabilized formulation presented in the previous sections. Stable mixed elements with higher order interpolation for the displacements have been used to validate the FIC-stabilized element. 

In both problems, the porous medium is under totally saturated conditions with isotropic permeability. The effect of gravity is not considered. 

## **2.4.1 Elastic soil column subjected to surface loading** 

This example consists on a 1 _×_ 30 _m_ column of saturated soil subjected to a surface loading that lies on a rigid rock bed. The geometry and boundary conditions of the problem are shown in Figure 2.4. 

The material properties of the soil column are summarized in Table 2.1. It must be noticed that the indices 1 and 2 in the solid bulk modulus _Ks_ , the fluid bulk modulus _Kf_ and the intrinsic permeability _k_ denote the different cases that have been considered here. Index 1 corresponds to the nearly undrained-incompressible limit, and index 2 corresponds to a more relaxed condition. 

This problem is solved in a 2D configuration under plane strain conditions. The geom- 

Examples 

33 

**==> picture [151 x 299] intentionally omitted <==**

**----- Start of picture text -----**<br>
1<br>q<br>H=30<br>y<br>x<br>**----- End of picture text -----**<br>


Figure 2.4: Elastic soil column subjected to surface loading. Geometry and boundary conditions. Dimensions in _m_ . 

etry is discretized with a structured mesh of 20 quadrilateral elements. We have solved the problem using 4-noded quadrilateral elements with bilinear shape functions for both the pressure and the displacements (Q4P4), and 9-noded quadrilateral elements with biquadratic shape functions for the displacements and bilinear shape functions for the pressure (Q9P4). 

Two different load cases have been considered: a surface step loading (Figure 2.5) and a surface cyclic loading (Figure 2.7). 

Coupled Solid-Pore Fluid Interaction Problem 

34 

|**Property**|**Value**|
|---|---|
|Young modulus (_E_)|2_._5_·_107 _N/m_2|
|Poisson’s ratio (_ν_)|0.2|
|Solid density (_ρs_)|2000_Kg/m_3|
|Fluid density (_ρf_)|1000_Kg/m_3|
|Porosity (_φ_)|0.3|
|Dynamic viscosity (_µ_)|0_._001_N/m_2 _· s_|
|Solid bulk modulus 1 (_Ks,_1)|1_._5_·_1017 _N/m_2|
|Fluid bulk modulus 1 (_Kf,_1)|3_·_1014 _N/m_2|
|Intrinsic Permeability 1 (_k_1)|1_·_10_−_14 _m_2|
|Solid bulk modulus 2 (_Ks,_2)|1_._5_·_1012 _N/m_2|
|Fluid bulk modulus 2 (_Kf,_2)|3_·_109 _N/m_2|
|Intrinsic Permeability 2 (_k_2)|1_·_10_−_11 _m_2|



Table 2.1: Elastic soil column problem. Material properties. 

## **Step loading case** 

The main purpose of this first case is to verify that the FIC-stabilized formulation captures properly the pressure distribution along the soil column under undrainedincompressible conditions. 

Figure 2.6 shows the normalized pore pressure along the normalized height of the soil column for nearly undrained conditions at a time _t_ = 2 _s_ . 

Understanding the results of Figure 2.6 requires taking into account the fact that gravity is not considered. Indeed, since the only force acting over the soil column is the surface load, the pressure applied by the surface load is completely transferred to the pore pressure throughout the column height and thus the ratio _p/q_ should equal 1. 

The results obtained with the non-stabilized formulation and the 4-noded quadrilaterals (hereafter called Q4P4-Direct) are not able to properly capture the pressure distribution along the column. The element locks and the pressure oscillates with arbitrary values. In contrast, by using the Q4P4 element stabilized with the FIC approach presented in this work (here called Q4P4-FIC), the correct pressure distribution is obtained. Similar 

Examples 

35 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
1200<br>1000<br>800<br>600<br>400<br>200<br>0<br>0 0.1 0.5 1 1.5 2<br>Time (s)<br>Surface Load (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.5: Surface step loading applied in the elastic soil column problem. 

good results are obtained with the Q9P4 element, as expected. 

Comparing the graphs in Figures 2.6a and 2.6b one can see the effect of the compressibility of the materials on the sought solution. Certainly, as we decrease the value of Biot’s modulus _Q_ , or equivalently, the bulk modulus of the solid and fluid phases _Ks_ and _Kf_ , the amplitude of the pressure oscillations obtained with the Q4P4-Direct element decreases near the rock bed. On the other hand, the results obtained with the Q4P4-FIC and Q9P4 elements remain unaltered. 

Coupled Solid-Pore Fluid Interaction Problem 

36 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
1<br>0.9<br>0.8<br>0.7<br>0.6<br>0.5<br>0.4<br>0.3<br>0.2<br>Q4P4-Direct<br>0.1 Q9P4<br>Q4P4-FIC<br>0<br>0 0.5 1 1.5 2<br>p/q<br>y/H<br>**----- End of picture text -----**<br>


(a) Quasi-incompressible limit ( _Q_ = 10[15] _N/m_[2] ). 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
1<br>0.9<br>0.8<br>0.7<br>0.6<br>0.5<br>0.4<br>0.3<br>0.2 Q4P4-Direct<br>Q9P4<br>0.1<br>Q4P4-FIC<br>0<br>0 0.5 1 1.5 2<br>p/q<br>y/H<br>**----- End of picture text -----**<br>


(b) Relaxed incompressibility ( _Q_ = 10[10] _N/m_[2] ). 

Figure 2.6: Normalized pore pressure along the soil column ( _k_ = 1 _·_ 10 _[−]_[14] _m_[2] , ∆ _t_ = 0 _._ 02 _s_ , _t_ = 2 _s_ ). 

Examples 

37 

## **Cyclic loading case** 

In this second case the objective is to analyse the evolution of the pore pressure in time, and assess the effect of the soil permeability on the dissipation of the pore pressure. 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
2000<br>q = 1000 + 500 sin(2t)<br>1500<br>1000<br>500<br>0<br>0 5 10 15 20 25 30<br>Time (s)<br>Surface Load (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.7: Surface cyclic loading applied in the elastic soil column problem. 

To that purpose, Figure 2.8 shows the temporal evolution of the pore pressure at a node located 1.5 m below the surface. 

Like in the previous case, the pressure obtained here should reflect the load applied on the surface. From the results in Figure 2.8, it is clear that the three tested elements show sinusoidal pressure evolutions in time. However, from the nearly undrained case of Figure 2.8a, it appears that while the Q4P4-FIC element shows a pressure that correctly oscillates from 500 to 1500 _N/m_[2] , the Q4P4-Direct element presents pressure values that vary from 1000 to 3000 _N/m_[2] , making no sense. One must notice that although the Q9P4 element shows higher pressures than the stabilized element, after refining the mesh up to 60 elements, the Q9P4 results are in good agreement with the solution obtained with 20 Q4P4-FIC elements (Figure 2.9). This evidences the good accuracy of the stabilized FIC-FEM Q4P4 element. 

Looking at the graph depicted in Figure 2.8b one can see that, by increasing the intrinsic 

Coupled Solid-Pore Fluid Interaction Problem 

38 

**==> picture [268 x 532] intentionally omitted <==**

**----- Start of picture text -----**<br>
3500<br>Q4P4-Direct<br>Q9P4<br>3000 Q4P4-FIC<br>2500<br>2000<br>1500<br>1000<br>500<br>0<br>0 5 10 15 20 25 30<br>Time (s)<br>(a) Quasi-undrained limit ( k = 10 [−] [14] m [2] ).<br>2000<br>Q4P4-Direct<br>Q9P4<br>1500<br>Q4P4-FIC<br>1000<br>500<br>0<br>-500<br>0 5 10 15 20 25 30<br>Time (s)<br>(b) Partially drained ( k = 10 [−] [11] m [2] ).<br>Pore Pressure (N/m2)<br>Pore Pressure (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.8: Evolution of the pore pressure with time at a node near the surface ( _Q_ = 1 _·_ 10[15] _N/m_[2] , ∆ _t_ = 0 _._ 05 _s_ ). 

Examples 

39 

permeability _k_ , the problem is no longer as bad conditioned as before, and so the gap in the amplitude of the three elements diminishes. Moreover, a greater permeability implies that the water can flow towards the surface due to the deformation of the column. Thereby, a dissipation in the pore pressure occurs as the computation advances in time. 

Finally, it is also interesting to notice that the pressure becomes negative after the soil has drained a certain amount of water. This is easily understood if one thinks about the elastic behaviour considered for the porous medium. Since the applied load is oscillating, the deformation of the elastic soil experiments a similar behaviour. In consequence, when the soil recovers its initial position as the load reduces, a suction appears that makes the pore pressure negative at some points near the surface. 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
1800<br>1600<br>1400<br>1200<br>1000<br>800<br>600<br>400<br>60 Q9P4<br>200 20 Q4P4-FIC<br>0<br>0 5 10 15 20 25 30<br>Time (s)<br>Pore Pressure (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.9: Time evolution of the pore pressure. Comparison between 20 Q4P4-FIC elements and 60 stable Q9P4 elements ( _Q_ = 1 _·_ 10[15] _N/m_[2] , _k_ = 10 _[−]_[14] _m_[2] , ∆ _t_ = 0 _._ 05 _s_ ). 

## **2.4.2 Elastic soil foundation subjected to surface loading** 

This second example consists on a pillar foundation on a saturated soil stratum lying over a rigid rock bed. The geometry and boundary conditions are sketched in Figure 

Coupled Solid-Pore Fluid Interaction Problem 

40 

2.10. 

**==> picture [284 x 241] intentionally omitted <==**

**----- Start of picture text -----**<br>
q<br>3 3<br>30<br>z<br>y x<br>30 30<br>**----- End of picture text -----**<br>


Figure 2.10: Elastic soil foundation subjected to surface loading. Geometry and boundary conditions. Dimensions in _m_ . 

The material properties for the soil are the same as for the previous example (Table 2.2). 

In this case the soil is subjected to a surface step loading of 1 _·_ 10[4] _N/m_[2] applied in a linear increment during 0 _._ 1 _s_ (Figure 2.11). 

As stated before, this problem is solved by means of a 3D analysis. We use 4-noded tetrahedral elements with linear interpolation for the displacements and the pressure in the non-stabilized form (T4P4-Direct) and the FIC-stabilized formulation (T4P4FIC). Ten-noded stable tetrahedral elements with quadratic shape functions for the displacement field and linear shape functions for the pressure field (T10P4) have also been used in the analysis for comparison purposes. 

This example has two objectives. First, it should help analysing the effect of the spatial discretization on the solution obtained with the FIC-stabilized formulation. Second, it may allow us to assess whether the three element types considered here converge to the expected solution when the mesh is refined. 

Examples 

41 

|**Property**|**Value**|
|---|---|
|Young modulus (_E_)|2_._5_·_107 _N/m_2|
|Poisson’s ratio (_ν_)|0.2|
|Solid density (_ρs_)|2000_Kg/m_3|
|Fluid density (_ρf_)|1000_Kg/m_3|
|Porosity (_φ_)|0.3|
|Dynamic viscosity (_µ_)|0_._001_N/m_2 _· s_|
|Solid bulk modulus 1 (_Ks,_1)|1_._5_·_1017 _N/m_2|
|Fluid bulk modulus 1 (_Kf,_1)|3_·_1014 _N/m_2|
|Intrinsic Permeability 1 (_k_1)|1_·_10_−_14 _m_2|
|Solid bulk modulus 2 (_Ks,_2)|1_._5_·_1012 _N/m_2|
|Fluid bulk modulus 2 (_Kf,_2)|3_·_109 _N/m_2|
|Intrinsic Permeability 2 (_k_2)|1_·_10_−_11 _m_2|



Table 2.2: Elastic soil foundation problem. Material properties. 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
12000<br>10000<br>8000<br>6000<br>4000<br>2000<br>0<br>0 0.1 1 2 3 4 5<br>Time (s)<br>Surface Load (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.11: Surface step loading applied in the elastic soil foundation problem. 

In order to do so, two different unstructured spatial discretizations have been used: a 

Coupled Solid-Pore Fluid Interaction Problem 

42 

coarse uniform mesh (Figure 2.12a) and a refined non-uniform mesh (Figure 2.12b). 

**==> picture [421 x 67] intentionally omitted <==**

**----- Start of picture text -----**<br>
z<br>z y x<br>y x<br>1 URS RSSANA LU Rare<br>(b) Refined non-uniform mesh: 13532 ele-<br>Coarse uniform mesh: 9985 elements.<br>ments.<br>**----- End of picture text -----**<br>


- (a) Coarse uniform mesh: 9985 elements. 

Figure 2.12: Spatial discretizations used for the elastic soil foundation problem. 

In Figure 2.13 we represent the evolution of the maximum pore pressure with time under nearly undrained-incompressible conditions for both the coarse and refined meshes using ∆ _t_ = 0 _._ 02 _s_ . The impermeability of the case makes that the pore pressure remains constant after the first 0 _._ 1 _s_ of application of the load for the two stable elements, the T10P4 and the T4P4-FIC. We see, though, a different behaviour in the unstable T4P4-Direct element, showing a certain dissipation for the finer mesh. 

Comparing the T10P4 and the T4P4-FIC elements we notice that, while the first one surpasses the threshold of 10000 _N/m_[2] for both meshes, the latter is always below the value of the applied load. In any case, both elements approach the expected solution as the mesh is refined, which indicates a consistent response. 

If we look at the contour lines of the pore pressure distribution obtained at time _t_ = 1 _s_ in Figure 2.14, we can understand the abnormal behaviour of the T4P4-Direct element. Indeed, the direct mixed formulation with equal order interpolation elements locks under undrained-incompressible conditions which leads to the incoherent pressure distribution obtained. 

Next we evaluate the capabilities of the FIC stabilized formulation to reproduce a solu- 

Examples 

43 

**==> picture [268 x 532] intentionally omitted <==**

**----- Start of picture text -----**<br>
18000<br>16000<br>14000<br>12000<br>10000<br>T4P4-Direct<br>8000 T10P4<br>T4P4-FIC<br>6000<br>4000<br>2000<br>0<br>0 0.2 0.4 0.6 0.8 1<br>Time (s)<br>(a) Coarse mesh.<br>20000<br>18000<br>16000<br>14000<br>12000<br>10000<br>8000<br>T4P4-Direct<br>6000<br>T10P4<br>4000 T4P4-FIC<br>2000<br>0<br>0 0.2 0.4 0.6 0.8 1<br>Time (s)<br>(b) Refined mesh.<br>Maximum Pore Pressure (N/m2)<br>Maximum Pore Pressure (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.13: Time evolution of the maximum pore pressure under undrainedincompressible conditions ( _Q_ = 1 _·_ 10[15] _N/m_[2] , _k_ = 1 _·_ 10 _[−]_[14] _m_[2] , ∆ _t_ = 0 _._ 02 _s_ ). 

Coupled Solid-Pore Fluid Interaction Problem 

44 

tion where stabilization is not needed. For this purpose the previous problem has been solved for partially drained and compressible conditions up to _t_ = 5 s (Figure 2.15). Results are again shown at _t_ = 1 s for ∆ _t_ = 0 _._ 02 s (Figure 2.16). Under these relaxed conditions, the equal order interpolation mixed element with the direct formulation is able to capture acceptable pressure distributions as expected (Figures 2.16a and 2.16b). 

The higher permeability of this second case also permits capturing the so-called MandelCryer effect, i.e. the increase of pore pressure due to the application of the load is immediate, but the dissipation due to the outflow of the pore fluid is delayed by the permeability of the porous medium [141]. Figures 2.15a and 2.15b show that, after the first 0.1 s of loading, the fluid starts draining due to the consolidation of the soil. 

In order to properly understand the low peak pore pressure obtained with the T4P4-FIC element in Figure 2.15a, it is useful to recall the definition of the stabilization parameter in Equation (2.56) and observe carefully the mesh in Figure 2.12a. It can be noticed that the area of application of the load is covered by only two triangular faces. This implies that the characteristic element length around that particular zone is relatively large, and so it is the stabilization parameter. As a result, one obtains an over-diffusive solution that leads to an underestimation of the maximum pore pressure. 

However, it is interesting to see that the larger peak pore pressure obtained with the T4P4-Direct and the T10P4 elements, caused by an initial locking of the pressure field, converges to the solution obtained with the T4P4-FIC element after the soil consolidates. 

On the other hand, the results for the maximum pressure obtained with the refined nonuniform mesh (Figures 2.15b, 2.16b, 2.16d and 2.16f) converge all to an almost identical value. This evidences that the FIC stabilization does not alter negatively the solution obtained with the original non-stabilized mixed formulation. 

Examples 

45 

(a) T4P4-Direct elem. in coarse mesh ( _pmax_ = 16282 _N/m_[2] ). 

(b) T4P4-Direct elem. in refined mesh ( _pmax_ = 15851 _N/m_[2] ). 

(c) T10P4 elem. in coarse mesh ( _pmax_ = 15197 _N/m_[2] ). 

(d) T10P4 elem. in refined mesh ( _pmax_ = 12283 _N/m_[2] ). 

(e) T4P4-FIC elem. in coarse mesh ( _pmax_ = 4390 _._ 8 _N/m_[2] ). 

(f) T4P4-FIC elem. in refined mesh ( _pmax_ = 9609 _._ 7 _N/m_[2] ). 

Figure 2.14: Pore pressure distribution under undrained-incompressible conditions ( _Q_ = 1 _·_ 10[15] _N/m_[2] , _k_ = 1 _·_ 10 _[−]_[14] _m_[2] , ∆ _t_ = 0 _._ 02 _s_ , _t_ = 1 _s_ ). 

Coupled Solid-Pore Fluid Interaction Problem 

46 

**==> picture [268 x 532] intentionally omitted <==**

**----- Start of picture text -----**<br>
16000<br>T4P4-Direct<br>14000<br>T10P4<br>12000 T4P4-FIC<br>10000<br>8000<br>6000<br>4000<br>2000<br>0<br>0 1 2 3 4 5<br>Time (s)<br>(a) Coarse mesh.<br>12000<br>T4P4-Direct<br>10000 T10P4<br>T4P4-FIC<br>8000<br>6000<br>4000<br>2000<br>0<br>0 1 2 3 4 5<br>Time (s)<br>(b) Refined mesh.<br>Maximum Pore Pressure (N/m2)<br>Maximum Pore Pressure (N/m2)<br>**----- End of picture text -----**<br>


Figure 2.15: Time evolution of the maximum pore pressure for partially compressible and drained conditions ( _Q_ = 1 _·_ 10[10] _N/m_[2] , _k_ = 1 _·_ 10 _[−]_[11] _m_[2] , ∆ _t_ = 0 _._ 02 _s_ ). 

Examples 

47 

(a) T4P4-Direct elem. in coarse mesh ( _pmax_ = 7268 _._ 1 _N/m_[2] ). 

(b) T4P4-Direct elem. in refined mesh ( _pmax_ = 8104 _._ 1 _N/m_[2] ). 

(c) T10P4 elem. in coarse mesh ( _pmax_ = 10427 _N/m_[2] ). 

(d) T10P4 elem. in refined mesh ( _pmax_ = 8173 _._ 8 _N/m_[2] ). 

(e) T4P4-FIC elem. in coarse mesh ( _pmax_ = 4089 _N/m_[2] ). 

**==> picture [153 x 28] intentionally omitted <==**

**----- Start of picture text -----**<br>
(f) T4P4-FIC elem. in refined<br>mesh ( pmax = 8086 . 8 N/m [2] ).<br>**----- End of picture text -----**<br>


Figure 2.16: Pore pressure distribution for partially compressible and drained conditions ( _Q_ = 1 _·_ 10[10] _N/m_[2] , _k_ = 1 _·_ 10 _[−]_[11] _m_[2] , ∆ _t_ = 0 _._ 02 _s_ , _t_ = 1 _s_ ). 

## 3 Chapter 

## **Continuum Damage Mechanics** 

## **3.1 Introduction** 

Many engineering materials subjected to unfavourable mechanical and environmental conditions undergo micro-structural changes that decrease their strength. Since these changes impair the mechanical properties of these materials, the term damage is generally used. 

This effect is particularly relevant and difficult to predict in brittle or quasi-brittle materials such as concrete, rocks, mortar or other geo-materials. For instance, concrete is a highly heterogeneous, anisotropic, brittle material with a very complex non-linear mechanical behaviour due to the occurrence of the localization of deformation. This localization of deformation can appear as cracks, if cohesive properties are dominant, or as shear zones, if frictional properties prevail. As a result of strain localization, material softening takes place and a significant reduction of the material stiffness occurs during cyclic loading. A good understanding of the mechanism of the formation of such localization is of crucial importance because it acts as a precursor of fracture and failure. 

Continuum damage mechanics is a branch of continuum mechanics that describes the progressive loss of material integrity due to the propagation and coalescence of microcracks, micro-voids, and similar defects (Figure 3.1). These changes in the microstructure lead to an irreversible material degradation, characterized by a loss of stiffness 

Continuum Damage Mechanics 

50 

that can be observed on the macro-scale. 

**==> picture [103 x 108] intentionally omitted <==**

Figure 3.1: Generic damaged section. 

The main objective of standard continuum damage mechanics is to propose a continuummechanics based framework that allows to characterize, represent and model at the macroscopic scale the effects of distributed defects and their growth on the material behaviour. 

The term “continuum damage mechanics” was first used by Hult [71], but the concept of damage was introduced by Kachanov in 1958 in the context of creep rupture [75]. In that work Kachanov introduced the concept of effective stress, and by using continuum damage he solved problems related to creep in metals. Rabotnov [118] gave the problem physical meaning by suggesting that the reduction of the sectional area was measured by means of the damage parameter. 

The thermodynamic formalism involved in the irreversible process of damage was developed by Lemaitre and Chaboche [82]. Other important contributions on damage mechanics include: Mazars and Pijaudier-Cabot [90, 91], Simo and Ju [130], Oller _et al._ [101], Oliver _et al._ [99, 100], and Cervera _et al._ [35–37]. 

If the damage parameter depends only on the strain state at the point under consideration, and non enriched kinematics are adopted to regularize the problem, numerical simulations exhibit a pathological mesh dependence and the energy consumed by the fracture process tends to zero as the mesh is refined. This is the typical behaviour of the so-called local damage models, which are not able to properly describe both the thickness of localization and the distance between damaged zones. They suffer from mesh sensitivity (for size and alignment) and produce unreliable results. Strains concentrate in one element wide zones and the computed force-displacement curves are 

Isotropic damage theory 

51 

mesh-dependent. The reason behind these misbehaviours is that the differential equations of motion change their type, from elliptic to hyperbolic in static problems, and the boundary value problem becomes ill-posed [11]. 

Classical constitutive models require an extension in the form of a characteristic length to properly model the thickness of localized zones. Such extension can be done by means of second gradient models [38], micro-polar [95], strain gradient [145], viscous [84] and non-local terms [74]. In this thesis we have worked with the latter approach using a weighted spatial averaging of the internal variables. In this manner the kinematic and equilibrium equations remain standard, and the notions of stress and strain keep their usual meaning. 

The first non-local models of this type, proposed in the 1960s, aimed at improving the description of elastic wave dispersions in crystals. Non-local elasticity was further developed by Eringen [54] who later extended it to non-local elasto-plasticity [53]. Subsequently, it was found that certain non-local formulations could act as efficient localization limiters with a regularizing effect on problems with strain localization [114]. 

The chapter is organized as follows. First, the basic concepts on continuum damage mechanics are introduced. Details are given on the essential components of the isotropic damage theory, including the equivalent strain forms and damage evolution laws implemented for this work. Following that, the regularization techniques developed to overcome the difficulties associated to strain localization are illustrated. Special attention is given to the integral-type non-local damage model, pointing out the most relevant aspects of its numerical implementation. Two examples are shown at the end, in order to assess the robustness of the implemented damage models against changes in the spatial discretization. 

## **3.2 Isotropic damage theory** 

Damage models work with certain internal variables that characterize the density and orientation of micro-defects. To introduce its concepts and understand the fundamental theory of damage mechanics, let us explain the stress evolution in a simple uniaxial tensile case (Figure 3.2). 

Continuum Damage Mechanics 

52 

Figure 3.2: Idealized material for the description of the uniaxial damage theory. 

For a better understanding, the material is idealized as a bundle of fibres parallel to the loading direction (Figure 3.2). Initially, all the fibres respond elastically, and the stress is carried by the total cross section of all fibres _A_ (Figure 3.3.1). As the applied strain is increased some fibres start breaking (Figure 3.3.2). Each fibre is assumed to be perfectly brittle, meaning that the stress in the fibre drops down to zero immediately after a critical strain level is reached. However, since the critical strain can differ from one fibre to another, the effective area _Ae_ (the area of unbroken fibres that can still carry stress) decreases gradually from _Ae_ = _A_ to _Ae_ = 0. Of course, when the applied force diminishes (Figure 3.4.2), the affected fibres remain broken and so the damaged area of the specimen is irrecoverable. 

It is important to make a distinction between the nominal stress _σ_ , defined as the force per unit initial area of the cross section, and the effective stress _σe_ , defined as the force per unit effective area[1] . The nominal stress enters the Cauchy equations of equilibrium on the macroscopic level, while the effective stress is the “actual” stress acting in the material micro-structure. Since the external applied force is the same regardless of the definition of the stress, i.e. _F_ = _σA_ = _σeAe_ , we can write: 

**==> picture [248 x 26] intentionally omitted <==**

The ratio of the effective area to the total area, _Ae/A_ , is a scalar characterizing the integrity of the material. Within the classical approach, a very simple measure of the damage amplitude in a given plane is obtained by measuring the area of the intersection of all defects with that plane. Thereby, we can define the damage variable as: 

**==> picture [296 x 26] intentionally omitted <==**

**==> picture [390 x 11] intentionally omitted <==**

Isotropic damage theory 

53 

**==> picture [68 x 35] intentionally omitted <==**

**==> picture [68 x 43] intentionally omitted <==**

Figure 3.3: Scheme of a uniaxial damage model through a monotonic loading process. 

where _Ad_ is the damaged part of the area. An undamaged material is characterized then by _Ae_ = _A_ , i.e. by d = 0. Due to propagation of micro-defects and their coalescence, the damage variable grows and at the late stages of degradation process it approaches asymptotically the limit value d = 1, corresponding to a complete damaged material with area reduced to zero. 

In the simplest version of the presented scheme, each fibre is supposed to remain linear elastic up to the strain level at which it breaks. Consequently, the 1D effective stress _σe_ is related to the elastic strain of the material _ε_ by the uniaxial Hooke’s law: 

**==> picture [244 x 13] intentionally omitted <==**

where E is the elastic modulus of the undamaged material. Combining (3.1), (3.2) and (3.3), it is straightforward to see that the constitutive law for the nominal stress _σ_ takes 

Continuum Damage Mechanics 

54 

**==> picture [65 x 49] intentionally omitted <==**

**==> picture [86 x 35] intentionally omitted <==**

Figure 3.4: Scheme of a uniaxial damage model through a non-monotonic loading process. 

the form: 

**==> picture [260 x 13] intentionally omitted <==**

For the uniaxial model formulation, Equation (3.4) must be complemented with the damage evolution law, which can be characterized by the dependence between the damage variable and the applied strain: 

**==> picture [280 x 13] intentionally omitted <==**

Function _g_ affects the shape of the stress-strain diagram and can be directly identified from a uniaxial test. 

Isotropic damage theory 

55 

The evolution of the effective stress, damage variable, and nominal stress in a material that remains elastic up to the peak stress is shown in Figure 3.3. This description is valid only for monotonic loading by an increasing applied strain _ε_ . 

When the material is first stretched up to a certain strain level _ε_ 1 that induces damage d1 = _g_ ( _ε_ 1) and then the strain decreases (Figure 3.4), the damaged area remains constant and the material responds as an elastic material with a reduced Young’s modulus _E_ 1 = (1 _−_ d1) _E_ . This means that, during unloading and partial reloading, the damage variable in (3.4) must be evaluated from the largest previously reached strain and not from the current strain _ε_ . This means that it is crucial to introduce an internal state variable _r_ characterizing the maximum strain level reached in the previous history of the material up to a given time _t_ : 

**==> picture [292 x 30] intentionally omitted <==**

The above expression implies that _r_ ( _t_ ) _≥ ry_ , where _ry_ is the so-called damage threshold, a material parameter that represents the value of strain at which damage starts. In this formula, _t_ is not necessarily the physical time (it can be any monotonically increasing parameter controlling the loading process). The damage evolution law (3.5) is then rewritten in the form: 

**==> picture [335 x 37] intentionally omitted <==**

which remains valid not only during monotonic loading but also during unloading and reloading. The evolution of the effective stress, damage variable, and nominal stress in a non-monotonic test is shown in Figure 3.4. Note that, upon a complete removal of the applied force, the strain returns to zero (due to elasticity of the unbroken fibres), i.e. the pure damage model does not take into account any permanent strains. Nevertheless, the material state is different from the initial virgin state because, according to (3.6) and (3.7), when the state variable _r_ becomes greater than _ry_ , the damage variable will not be zero again and thus the stiffness and strength mobilized in a new tensile loading process will be smaller than their initial values. Therefore, we can say that the loading history is always reflected by the value of the internal state variable _r_ . 

Instead of defining the variable _r_ like in (3.6), we can introduce a loading function _f_ ( _ε, r_ ) = _ε − r_ and postulate the loading-unloading conditions in the Kuhn-Tucker form: 

**==> picture [282 x 13] intentionally omitted <==**

Continuum Damage Mechanics 

56 

The first condition indicates that _r_ can never be smaller than _ε_ , while the second condition means that _r_ cannot decrease. Finally, according to the third condition, _r_ can grow only if the current values of _ε_ and _r_ are equal. 

At this point, we can already formulate the extension of the uniaxial damage theory to general multiaxial stress states. 

In this work we have chosen the simple isotropic damage model with a unique scalar variable. Isotropic damage models are based on the simplifying assumption that the stiffness degradation is isotropic, i.e. the stiffness moduli corresponding to different directions decrease proportionally, independently of the direction of loading. Since an isotropic elastic material is characterized by two independent elastic constrains, a general isotropic damage model should deal with two damage variables. The model with a single variable makes use of the additional assumption that the relative reduction of all the stiffness coefficients is the same, in other words, that the Poisson’s ratio is not affected by damage. In this regard, the stress-strain law in Voigt notation is postulated as: 

**==> picture [265 x 13] intentionally omitted <==**

where _**D** e_ is the elastic constitutive matrix of the intact material. One can clearly notice that (3.9) is a generalization of (3.4). Also, Equation (3.9) can alternatively be written as: 

**==> picture [260 x 13] intentionally omitted <==**

which is the multidimensional generalization of (3.1), and where _**σ** e_ is the effective stress vector as: 

**==> picture [249 x 13] intentionally omitted <==**

Similarly as in the uniaxial case, we introduce a loading function _f_ specifying the elastic domain and the states at which damage grows. The loading function now depends on the strain vector _**ε**_ , and on a variable _r_ that controls the evolution of the elastic domain. Physically. _r_ represents a scalar measure of the largest strain level ever reached in the history of the material. States for which _f_ ( _**ε** , r_ ) _<_ 0 are supposed to be below the current damage threshold. Damage can grow only if the current state reaches the boundary of the elastic domain, i.e. when _f_ ( _**ε** , r_ ) = 0. Essentially, we can postulate the damage criterion for a multiaxial isotropic damage model with: 

Isotropic damage theory 

57 

- The loading function: 

**==> picture [258 x 13] intentionally omitted <==**

- The loading-unloading conditions (3.8) 

_εeq_ in Eq (3.12) is the equivalent strain, a scalar magnitude of the strain level, and _r_ is the largest value of equivalent strain calculated in the previous deformation history of the material up to its current state. In this regard, (3.6) can now be generalized to: 

**==> picture [296 x 30] intentionally omitted <==**

An important advantage of an isotropic damage model is that the stress evaluation algorithm is usually explicit and there is no need to use an iterative solution for nonlinear equations. 

Thereby, for a particular strain increment, the corresponding stress is obtained by computing the current value of equivalent strain, updating the maximum previously reached equivalent strain and the damage variable, and reducing the effective stress according to (3.10). In essence, one must follow the scheme of Table 3.1. 

Step _n_ + 1 

1. Evaluate effective stress: _**σ** e,n_ +1 = _**D** e_ _**ε** n_ +1 2. Compute the new equivalent strain: _εeq,n_ +1 3. Update _r_ with (3.13): If _εeq,n_ +1 _> rn ⇒ rn_ +1 = _εeq,n_ +1 4. Update damage variable: d _n_ +1 = _g_ ( _rn_ +1) 5. Compute nominal stress: _**σ** n_ +1 = (1 _−_ d _n_ +1) _**σ** e,n_ +1 

Table 3.1: Computation of the stresses in the isotropic damage model. 

Continuum Damage Mechanics 

58 

## **3.2.1 Elastic-damage tangent constitutive tensor** 

The material non-linearity introduced by the damage model requires solving the system of equations by means of incremental-iterative techniques such as the Newton-Raphson’s method. If one aims at achieving reasonable convergence rates, it is crucial to properly compute the tangent constitutive tensor. 

To that purpose, it is essential to understand that the global stiffness of the structure changes through any damage process. In 1D problems one can easily distinguish the secant Young’s modulus _Esec_ from the tangent one _Etan_ and compute them as: 

**==> picture [288 x 26] intentionally omitted <==**

Figure 3.5 shows the stress-strain relation on a unidimensional damage process. The branch on the left side of the curve corresponds to the elastic loading case. The stress increases proportionally to the strain growth, and goes back following the same path when the strain diminishes. The secant and tangent elastic modulus coincide through all this stage. After a certain strain value is reached, i.e. _ε > ry_ , the relation between the stress and the strain becomes non-linear as the material starts degrading. In this situation, the secant and tangent Young’s moduli diverge and the internal state variable starts increasing ( _dr >_ 0). Finally, since no plastic deformation is considered, the strain is fully recovered upon unloading following a branch with a lower stiffness than the original. 

Understanding the previous concepts is not so straightforward in a general 3D case and so it is convenient to derive a general expression for the elastic-damage tangent constitutive tensor _**D** tan ≡_ _**D**_ . 

Let us remind first that the tangent constitutive tensor defines the tangent stiffness matrix _**K**_ (2.39) and must satisfy the following relation: 

**==> picture [251 x 13] intentionally omitted <==**

Considering the stress vector defined in (3.9) as a function of the strain vector and the damage variable, i.e. _**σ**_ = _**σ**_ ( _**ε** ,_ d), we compute its total derivative as: 

**==> picture [417 x 26] intentionally omitted <==**

Isotropic damage theory 

59 

**==> picture [85 x 31] intentionally omitted <==**

**==> picture [99 x 81] intentionally omitted <==**

Figure 3.5: Unidimensional stress-strain diagram throughout a damage process. 

Like in the 1D case (Figure 3.5), we distinguish two possible situations, bearing in mind that the Kuhn-Tucker conditions hold (Eq (3.8)): 

- Elastic loading or unloading process ( _f_ ( _**ε** , r_ ) _<_ 0 & _dr_ = 0) 

Whenever the structure is under an elastic loading or unloading process, there is no change in the internal historical variable _r_ (3.13). As a result, the damage variable d = _g_ ( _r_ ) also remains constant ( _d_ d = 0) and so the expression in (3.16) results: 

**==> picture [257 x 13] intentionally omitted <==**

Consequently, in this case the tangent constitutive tensor coincides with the secant constitutive tensor: 

**==> picture [269 x 13] intentionally omitted <==**

- Loading process with growing damage ( _f_ ( _**ε** , r_ ) = 0 & _dr >_ 0) 

If the structure is under a loading process with the equivalent strain exceeding the damage threshold _ry_ , the internal state variable increases with time, and so the damage variable ( _d_ d _>_ 0). 

From the definitions in (3.7), (3.12) and (3.13), we compute the total derivative of the damage variable as: 

**==> picture [298 x 26] intentionally omitted <==**

Continuum Damage Mechanics 

60 

Where the derivative of the equivalent strain with respect to the strain vector _dεeq/d_ _**ε**_ is estimated by means of the perturbation method [87]. This method is actually convenient as it allows us to compute the previous derivative regardless of the form of the equivalent strain: 

**==> picture [303 x 28] intentionally omitted <==**

with **∆** _**ε** i_ being the perturbation vector of the _i[th]_ strain component. 

**==> picture [285 x 25] intentionally omitted <==**

Now, introducing Equation (3.19) into (3.16) and rearranging terms gives 

**==> picture [348 x 96] intentionally omitted <==**

And so the expression defining the tangent constitutive tensor reads 

The typical scheme for the numerical implementation of the tangent constitutive tensor is summarized in Table (3.2). 

|||
|---|---|
|Step _n_+ 1||
|1. Compute secant constitutive tensor: **_D_**_sec,n_+1 =(1_−_d_n_+1)**_D_**_e_||
|If _εeq,n_+1 _≥rn_|2. Compute damage function derivative: (_dg/dr_)_n_+1<br>3. Calculate equivalent strain derivative: (_dεeq/d_**_ε_**)_n_+1<br>4. Compute tangent constitutive tensor:<br>**_D_**_n_+1 =**_D_**_sec,n_+1_−_(_dg_<br>_dr_)_n_+1**_σ_**_e,n_+1_×_ (_dεeq_<br>_d_**_ε_** )_n_+1|



Table 3.2: Computation of tangent constitutive tensor in the isotropic damage model. 

Isotropic damage theory 

61 

## **3.2.2 Equivalent strain** 

To some extent, the equivalent strain presented in (3.12) plays a role similar to the yield function in plasticity, because it directly affects the shape of the elastic domain (Figure 3.6). 

Figure 3.6: 3D elastic domain for a generic equivalent strain. 

There are numerous forms for the equivalent strain in the literature, but the simplest choice is to it as the Euclidean norm of the strain vector: 

**==> picture [293 x 15] intentionally omitted <==**

or as the energy norm: 

**==> picture [263 x 26] intentionally omitted <==**

The two norms of _**ε**_ introduced above, lead to symmetric elastic domain in tension and compression (Figure 3.7a). Nevertheless, several materials (rocks, concrete, ceramics, etc.) often show a non symmetric damage surface, i.e., the yield value in compression can be several times the value in tension. In order to overcome this limitation, it is necessary to adjust the definition of the equivalent strain. 

In this regard, one may work wit a definition of the equivalent strain based on the model proposed by Simo and Ju [130], using the energy norm of the strain and modifying it to 

Continuum Damage Mechanics 

62 

**==> picture [373 x 406] intentionally omitted <==**

**----- Start of picture text -----**<br>
Elastic domain<br>Elastic domain<br>(a) Symmetric energy norm. (b) Adapted Simo & Ju norm.<br>Elastic domain Elastic domain<br>(c) Mazars norm. (d) Modified von Mises norm.<br>**----- End of picture text -----**<br>


Figure 3.7: Damage surfaces in the 2D principal stress space. 

take into account the different degradation rate in tension and compression of concretelike materials (Figure 3.7b). Such definition takes the form: 

**==> picture [296 x 30] intentionally omitted <==**

where the variable _θ_ is a weighting factor computed from the eigenvalues of the effective stress tensor: 

**==> picture [259 x 33] intentionally omitted <==**

Isotropic damage theory 

63 

with _⟨.⟩_ denoting the “positive part” operator or McAuley brackets. For scalars _⟨x⟩_ = max(0 _, x_ ), i.e., 

**==> picture [278 x 37] intentionally omitted <==**

The parameter _κ_ in (3.26) is defined as the ratio between the compression elastic limit and the tension elastic limit, i.e. 

**==> picture [243 x 26] intentionally omitted <==**

Eq (3.26) is not the only possible form of equivalent strain valid for quasi-brittle materials. Micro-cracks in concrete grow mainly if the material is stretched, and so it is natural to take into account only strains that are positive (tensile strains) and neglect those that are negative (compressive strains). This leads to the so-called Mazars definition of equivalent strains [89]: 

**==> picture [264 x 44] intentionally omitted <==**

where _εi_ are the principal values of the strain tensor. 

Alternatively, the modified von Mises definition [49], reads: 

**==> picture [360 x 37] intentionally omitted <==**

where _I_ 1 is the first invariant of the strain tensor and _J_ 2 is the second invariant of the deviatoric strain tensor. Given a generic symmetric strain tensor: 

**==> picture [278 x 51] intentionally omitted <==**

The first invariant _I_ 1 is the trace of the strain tensor: 

**==> picture [296 x 13] intentionally omitted <==**

Also, one can always decompose the strain tensor into its volumetric and deviatoric parts [ _ε_ ] = [ _εv_ ] + [ _εd_ ]: 

**==> picture [273 x 51] intentionally omitted <==**

Continuum Damage Mechanics 

64 

**==> picture [318 x 51] intentionally omitted <==**

From the deviatoric strain tensor [ _εd_ ] we can calculate _J_ 1 and _J_ 2: 

**==> picture [380 x 26] intentionally omitted <==**

**==> picture [408 x 54] intentionally omitted <==**

The Mazars form of equivalent strain (3.30) and the modified von Mises equation (3.31) lead to different damage surfaces as shown in Figures 3.7c and 3.7d. However, in both cases the elastic limit in tension _σyt_ is considerably lower than the elastic limit in compression _σyc_ , and thus they are convenient for modelling damage progression of quasibrittle materials such as rocks or concrete. 

## **3.2.3 Damage evolution law** 

Once we have defined the energy norm in the strain space or equivalent strain _εeq_ , let us introduce the law _g_ ( _r_ ) governing the evolution of the damage variable (3.7). 

There are various damage governing laws that can be effectively used to model damage growth in geo-materials. One option, especially suited for simplified analyses, is the linear softening law: 

**==> picture [292 x 27] intentionally omitted <==**

Such a relation limits the range of the internal state variable _r_ between the damage threshold and a maximum admissible value _rmax_ , and leads to a linear softening branch in the stress-strain curve (Figure 3.8a). 

Two typical choices to describe the evolution of degradation in quasi-brittle materials are the polynomial law [115]: 

**==> picture [321 x 29] intentionally omitted <==**

Isotropic damage theory 

65 

**==> picture [300 x 213] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Linear softening law. (b) Polynomial softening law.<br>(c) Exponential softening law.<br>**----- End of picture text -----**<br>


Figure 3.8: Generic unidimensional stress-strain curves for different softening laws. 

and the exponential softening law [90]: 

**==> picture [333 x 26] intentionally omitted <==**

In the above expressions, the parameter _R_ is associated to the residual strength of the material, whereas the parameter _S_ controls the slope of the softening branch after the peak of the stress-strain curve. 

The polynomial equation (3.39) and the exponential softening law (3.40) mainly differ from the linear law (3.38) in that the material preserves some residual strength in the post-peak branch, as shown in the 1D stress-strain curves of Figures 3.8b and 3.8c. 

Continuum Damage Mechanics 

66 

An alternative exponential softening model was proposed in [99]: 

**==> picture [312 x 30] intentionally omitted <==**

The parameter _Af_ is obtained from the following expression: 

**==> picture [293 x 32] intentionally omitted <==**

where _Gf_ is the specific fracture energy per unit area, _lf_ is the characteristic length for the fractured domain, usually taken as the characteristic length of the finite elements, and _ry_ is the aforementioned damage threshold which, for quasi-brittle materials, can be estimated from the tensile strength _σyt_ : 

**==> picture [248 x 25] intentionally omitted <==**

Another popular damage evolution law specifically designed for concrete was proposed by Mazars [89, 90]. He divided the damage variable in two components: d _c_ and d _t_ , which are computed from the equivalent strain of Eq (3.30), but using two different damage laws: _gc_ and _gt_ . Function _gc_ is characterized by the uniaxial compressive test, while _gt_ corresponds to the uniaxial tensile test. Both functions were proposed in the same exponential form of (3.40): 

**==> picture [342 x 26] intentionally omitted <==**

**==> picture [340 x 26] intentionally omitted <==**

The damage evolution law _g_ ( _r_ ) results from the combination of the two parts: 

**==> picture [303 x 16] intentionally omitted <==**

The coefficient _αt_ ranges from 0 to 1 and takes into account the character of the stress state. It is evaluated from: 

**==> picture [263 x 36] intentionally omitted <==**

where _εti_ are the principal strains due to positive effective stresses (or tensile stresses). 

Local and Non-local Damage Models 

67 

The parameter _β_ in (3.46) was equal to 1 in the original version of the model [89]. For values higher than 1, _β_ allows to slow down the evolution of damage under shear loading (when principal stresses have different sign). 

The definition of _εti_ tells us that if all principal stresses are negative then _αt_ = 0 and d = d _c_ = _gc_ ( _r_ ), corresponding to a “purely compressive” stress state. On the other hand if all principal stresses are positive, i.e. in a “purely tensile” stress states, we have _αt_ = 1 and d = d _t_ = _gt_ ( _r_ ). 

## **3.3 Local and Non-local Damage Models** 

In the last section we have presented the fundamental theory of the isotropic damage models with a unique scalar variable d. Although these models are relatively simple, they are often used to model the failure of concrete, rocks and other quasi-brittle materials that show important strain localization. If the damage parameter depends only on the strain state at the point under consideration, the numerical simulations exhibit a pathological mesh dependence, and physically unrealistic results are obtained. 

This is the typical behaviour of the so-called local damage models, which are not able to properly describe both the thickness of localization and distance between them. They suffer from mesh sensitivity (for size and alignment) and produce unreliable results. The strains concentrate in one element wide zones and the computed force-displacement curves are mesh-dependent. The reason is that differential equations of motion change their type (from elliptic to hyperbolic in static problems) and the rate boundary value problem becomes ill-posed [11]. 

Classical constitutive models require an extension in the form of a characteristic length to properly model the thickness of localized zones. Such extension can be done with micro-polar [95, 138], strain gradient [92, 96, 145], viscous [51, 84, 117] and non-local terms [10, 74, 114, 123]. In this thesis we have developed the latter approach. 

First of all, we present the problems derived from the strain localization phenomenon, using a simple uniaxial case as example. Following that, the basic concepts of the implemented non-local damage model are introduced. 

Continuum Damage Mechanics 

68 

## **3.3.1 Strain localization phenomenon** 

The idea of modelling damaged concrete and other quasi-brittle materials as strainsoftening continua, was not immediately accepted by all the scientific community. Actually, most of the early analyses were not truly objective and, upon mesh refinement, their results would not converge to a robust solution. Let us explain the nature of the problem by means of a unidimensional example. 

Consider a straight bar with a constant cross section _A_ and a total length _L_ 0 under uniaxial tension (Figure 3.9). The material is assumed to obey a simple stress-strain law with linear elasticity up to the peak stress _σy_ , followed by linear softening (Figure 3.10). If the bar is loaded in tension by an applied displacement _u_ at its right extreme, the response remains linear elastic up to _uy_ = _L_ 0 _εy_ , instant at which the force transmitted by the bar (reaction at the support) attains its maximum value _Fy_ = _σyS_ . After that, the resistance of the bar starts decreasing until the strain reaches _εf_ and the transmitted stress completely disappears. 

Figure 3.9: Bar under uniaxial tension. 

The equilibrium equations at a section in the middle imply that the axial forces are constant along the bar and so the stress profile must remain also uniform (Figure 3.11). 

Local and Non-local Damage Models 

69 

Figure 3.10: Stress-strain diagram of the linear softening law. 

**==> picture [173 x 118] intentionally omitted <==**

Figure 3.11: Axial force acting along the bar. 

However, at a given stress level _σ_ 1 between 0 and _σy_ , there are actually two values of strain, _ε_ 1 _,e_ and _ε_ 1 _,s_ , that satisfy the constitutive equation (Figure 3.12). This is quite straightforward if one considers that any cross section can either be damaged, or undamaged. Indeed, an undamaged section is on the elastic branch with _σ_ 1 = _Eε_ 1 _,e_ , whereas a damaged one falls in the softening branch with _σ_ 1 = (1 _−_ d1) _Eε_ 1 _,s_ . 

Thereby, the strain profile along the bar does not have to be necessarily uniform. In fact, any piecewise constant strain distribution that jumps between the two possible strain values represents a valid solution. 

Continuum Damage Mechanics 

70 

Figure 3.12: Possible strain values corresponding to the same stress level. 

Figure 3.13: Piecewise constant strain profile along the bar. 

In Figure 3.13 we have denoted by _Ls_ the cumulative length of the softening regions and by _Le_ = _L_ 0 _− Ls_ the cumulative length of the elastic regions. When the stress is completely relaxed, the strain in the elastic region is _εe_ = 0 and the strain in the softening branch is _εs_ = _εf_ . Thus the total elongation of the bar in this case is _uf_ = _Leεe_ + _Lsεs_ = _Lsεf_ . At this point, although _εf_ is perfectly known from the linear softening law (Figure 3.10), _Ls_ is totally undetermined. This means that the problem has infinite possible solutions with its corresponding post-peak branches in the loaddisplacement diagram (Figure 3.14). 

This fan of post-peak branches is bounded on one side by the solution with uniform softening ( _uf_ = _L_ 0 _εf_ ) and on the other side by the solution with uniform elasticity 

Local and Non-local Damage Models 

71 

Figure 3.14: Fan of possible post-peak branches of the load-displacement diagram. 

( _uf_ = 0). The first limit corresponds to a totally damaged bar while the latter represents the case in which the bar is unloaded before any damage takes place. All the other solutions describe processes in which a part of the bar is damaged. However, it is not trivial to determine which of all these solutions is the one that reflects better the actual failure process. 

The ambiguity is removed if imperfections are taken into account. The material properties and sectional dimensions of a real bar are not perfectly uniform. Thereby, supposing that there is a small region where the strength is lower than in the remaining portion of the bar, when the applied stress reaches the peak of that weaker region, softening starts and the stress decreases. Consequently, the material outside the damaged region must unload elastically because its strength has not been exhausted. This leads to the conclusion that the size of the softening region is related to the size of the region with minimum strength. The problem is that such a region can be arbitrarily small and so the corresponding softening branch is arbitrarily close to the elastic unload branch. Therefore, the standard strain-softening continuum formulation leads to a solution with several pathological features: 

- The softening region is infinitely small. 

- The load-displacement diagram always shows snap-back, regardless of the size of the structure and the ductility of the material. 

- The total amount of energy dissipated during the failure process tends to zero. 

Continuum Damage Mechanics 

72 

From the mathematical point of view, these problematic features are related to the loss of ellipticity of the governing differential equation. As a result, the boundary value problem no longer has a unique solution with continuous dependence on the given data. 

From the numerical point of view, these inconveniences are manifested by a pathological sensitivity of the results to the size of the finite elements. 

For instance, let us assume that the bar is uniformly discretized by _n_ two-node elements with linear displacement interpolation and that the weakest region is located at the middle of the bar. The numerical algorithm will capture a very localized solution with a softening region extending over one only element. In consequence, the softening length will decrease as the number of elements increases ( _Ls_ = _L_ 0 _/n_ ) and thus the softening post-peak branch will completely depend on the number of elements of the mesh. 

Indeed, as it is shown in Figure 3.15a, for _n_ = 1 all the bar is damaged and the softening length is the total length of the bar _Ls_ = _L_ 0, whereas for _n >_ 1 the softening region is more localized with strains becoming especially important at the damaged element. In the limit case of _n →∞_ the softening branch would coincide with the elastic branch (see Figure 3.15b). 

**==> picture [50 x 38] intentionally omitted <==**

**==> picture [373 x 32] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) Load-displacement diagrams for<br>(a) Strain profiles for a prescribed im-<br>different number of elements.<br>posed displacement.<br>**----- End of picture text -----**<br>


Figure 3.15: Mesh dependence of the numerical results. 

In this section we presented a problem that commonly arises in the simulation of damage processes involving quasi-brittle materials: the strain localization phenomenon. Al- 

Local and Non-local Damage Models 

73 

though only a uniaxial case has been discussed, this problem is also present in multidimensional problems with similar effects on the numerical results. 

The simple one-dimensional example illustrates the essence of the problem with localization of inelastic strain into a zone of an arbitrarily small width. In uniaxial cases, localization occurs when the peak of the stress-strain diagram is reached, independently of the specific constitutive model used. In multiple dimensions, the analysis of the localization process is more complicated and the derivation of a criteria for the start of localization represents a challenging problem. 

Mesh refinement in multiple dimensions leads to a reduction of the total dissipated energy (area under the load-displacement curve) with a lower peak load and a more brittle response. Apart from this, upon further refinement, one can also face convergence difficulties due to the abrupt change of strain distribution, from a smoothly distributed to a highly localized one. These effects will be shown more clearly with the simulation of a bi-dimensional case later in this chapter (see Section 3.4). 

## **3.3.2 Regularization of the problem** 

When simulating fracture propagation processes, it is essential to ensure that the direction of crack growth is not affected by numerical conditioners such as the mesh size or the type of element. For this reason, in the present work two different continuum damage models for the porous medium were implemented in the research of a robust method that avoids the pathological sensitivity of the finite element results to the mesh size. 

## **Partially regularized local damage model** 

As a first attempt, we used a simple partially regularized local damage model, based on the crack band models [13]. We combined the equivalent strain form of Eq (3.26) and the damage evolution law in Eq (3.41). 

Essentially, the model is an isotropic damage model following the classical local damage theory, in which an energy based adjustment of the stress-strain diagram, depending on the size of the element, is introduced in the definition of the damage parameter. 

Continuum Damage Mechanics 

74 

Indeed, the damage evolution law in (3.41) depends on the characteristic length of the fractured domain _lf_ included in the definition of the variable _Af_ (3.42). This characteristic length of the fractured domain is what, in the end, partially regularizes the damage model. If the mesh size is modified, the energy dissipated by the structure is scaled according to the element characteristic length and, ideally, the energy consumption should remain unaltered by the changes on the spatial discretization. 

Thereby, before updating the damage variable, one just needs to compute the characteristic length of the element _l[e]_ , which can be estimated from the dimensions of the element like it was done in Equations (2.57) and (2.58). In a two-dimensional analysis, for instance, the characteristic length of the element can be defined as the diameter of the circle that contains the area of the element (see Figure 3.16). 

**==> picture [182 x 77] intentionally omitted <==**

Figure 3.16: Characteristic length of a 2D finite element. 

This approach is endowed with some, but not all of the properties of fully regularized damage models. It can ensure a correct energy dissipation in a localized damage band, but the width of the numerically resolved fracture process zone depends on the element size and tends to zero as the mesh is refined. This is why this approach cannot be considered as a true localization limiter. It provides only a partial regularization of the problem in the sense that global response characteristics do not exhibit spurious mesh sensitivity, but the mesh-induced directional bias is still present. 

Moreover, the scaling of the stress-strain diagram is straightforward only for models that explicitly control the evolution of inelastic strain, e.g. for softening plasticity [113] or smeared crack models [99]. In those cases, the desired scaling effect is achieved just by modifying the hardening modulus (derivative of stress with respect to inelastic strain). Nevertheless, in continuum damage mechanics, non-linearity and softening are controlled by the damage evolution law, and the reduction factor (1 _−_ d) multiplies the total strain. For this reason, it is not easy to scale only the post-peak part of strain localization while keeping the unloading part unaffected. 

Local and Non-local Damage Models 

75 

In some cases, diffuse softening damage patterns in certain parts of the structure can coexist with localized cracks in other parts, and they may even change during the loading process. In such cases it is virtually impossible to define a coherent rule for the adjustment of the stress-strain diagram according to the element size. 

## **Fully regularized non-local damage model** 

The introduction of a characteristic length into the constitutive model, and the formulation of a non-local strain-softening model, have been shown to prevent the spurious localization of strain-softening damage, to regularize the boundary value problem, and to ensure numerical convergence to physically meaningful solutions [53, 114, 123]. 

Fully regularized description of localized inelastic deformation is achieved by a proper generalization of the underlying continuum theory. Two different approaches are normally used: generalization of the kinematic relations, i.e. continua with micro-structure (Cosserat-type continua or strain gradient theories), and continua with non-local strain (non-local elasticity); and generalization of constitutive equations, i.e., material models with gradients of internal variables, and materials models with weighted spatial averages of internal variables. 

In this thesis we have worked with the second kind of generalization, because kinematic and equilibrium equations remain standard, and the notions of stress and strain keep their usual meaning. Also, in the generalization of constitutive equations through nonlocal models, we have focused on the integral-type methods. 

Integral-type non-local models abandon the classical assumption of locality and admit that stress a certain point depends, not only on the state variables at that point, but also on the distribution of the state variables over the whole body, or over a finite neighbourhood of the point under consideration. 

In a general manner, the non-local integral approach consists in replacing a certain local variable by its non-local counterpart, which is usually obtained by a weighted averaging over a spatial neighbourhood of each point under consideration. 

Let _f_ ( _**x**_ ) be some local field in a domain Ω. The corresponding non-local field is defined 

Continuum Damage Mechanics 

76 

as: 

**==> picture [314 x 28] intentionally omitted <==**

where _Z_ ( _**x** ,_ _**χ**_ ) is the non-local weighting function. 

In an infinite, isotropic and homogeneous medium, the weighting function depends only on the distance _Dxχ_ between the source point _**χ**_ , and the receiver point _**x**_ . Thereby, we usually write _Z_ ( _**x** ,_ _**χ**_ ) = _Z_ 0( _∥_ _**x** −_ _**χ** ∥_ ) = _Z_ 0( _Dxχ_ ), where _Z_ 0 is usually chosen as a non-negative function monotonically decreasing for _D ≥_ 0. 

One possible choice for _Z_ 0 is the Gauss distribution function: 

**==> picture [300 x 37] intentionally omitted <==**

where _lc_ is the characteristic length, a material parameter reflecting the internal length of the non-local continuum. 

If a bounded support is preferred, one can also truncate the previous function as follows: 

**==> picture [343 x 46] intentionally omitted <==**

where the interaction radius _Rin_ is a parameter related to the characteristic length _lc_ . In the present work, we have considered _Rin_ = _lc_ . 

In this thesis we worked with the previous Gauss distribution (3.50), but there are other alternatives that can be effectively used, e.g. the following truncated quartic polynomial function: 

**==> picture [338 x 51] intentionally omitted <==**

In essence, the interaction radius _Rin_ represents the smallest distance between points _**x**_ and _**χ**_ at which the interaction weight _Z_ 0( _Dxχ_ ) vanishes (for weighting functions with a bounded support) or becomes negligible (for weighting functions with an unbounded support). It must be carefully chosen as it controls the size of the softening region. 

The interval, circle, or sphere of radius _Rin_ , centered at _**x**_ , is called the domain of influence of point _**x**_ . In the vicinity of the boundary of a finite body, it is simply 

Local and Non-local Damage Models 

77 

assumed that the averaging is performed only on the part of the domain of influence that lies within the body (Figure 3.17). 

**==> picture [198 x 186] intentionally omitted <==**

**----- Start of picture text -----**<br>
Domain of influence<br>Receiver point<br>Source points<br>Solid domain<br>**----- End of picture text -----**<br>


Figure 3.17: Domain of influence protruding through the boundary of a body. 

Thereby, if a weighting function with bounded support is chosen, the non-local average of (3.48) is calculated as a weighted sum over the values at all the finite element integration points _**χ**_ lying within the non-local interaction radius _Rin_ . 

In the application to softening materials, it is often required that the non-local operator do not alter a uniform field (consistency of order 0), which means that the weighting function must satisfy the normalizing condition: 

**==> picture [293 x 28] intentionally omitted <==**

In order to satisfy (3.52), the weighting function is usually rescaled as: 

**==> picture [300 x 44] intentionally omitted <==**

In the present work, the damage variable is computed from the non-local equivalent strain, which corresponds to a suitable non-local damage formulation that restores wellposedness of the boundary value problem [12]. 

Continuum Damage Mechanics 

78 

Therefore, in the loading function (3.12), the local equivalent strain _εeq_ is replaced by its weighted spatial average: 

**==> picture [298 x 28] intentionally omitted <==**

And the internal state variable _r_ is then the largest previously reached value of the non-local equivalent strain: 

**==> picture [296 x 29] intentionally omitted <==**

It is important to note that while the damage variable is evaluated from the non-local equivalent strain _ε_ ˘ _eq_ , the strains _**ε**_ used in (3.9) to compute the stresses are considered as local. This way, during the elastic range, when the damage variable remains equal to zero, the stress-strain relation is completely local. The process for the evaluation of stresses with the non-local damage model is schematically represented in Table 3.3. 

## Step _n_ + 1 

1. Evaluate effective stress: _**σ** e,n_ +1 = _**D** e_ _**ε** n_ +1 2. Compute the new local equivalent strain: _εeq,n_ +1 3. Compute the new non-local equivalent strain (3.54): _ε_ ˘ _eq,n_ +1 ˘ ˘ 4. Update _r_ with (3.55): If _εeq,n_ +1 _> rn ⇒ rn_ +1 = _εeq,n_ +1 5. Update damage variable: d _n_ +1 = _g_ ( _rn_ +1) 6. Compute nominal stress: _**σ** n_ +1 = (1 _−_ d _n_ +1) _**σ** e,n_ +1 

Table 3.3: Computation of the stresses in the integral-type non-local damage model. 

As it can be deduced from Table 3.3, the numerical implementation of the non-local damage model based on averaging of the equivalent strain is relatively simple. The 

Local and Non-local Damage Models 

79 

evaluation of the stresses remains explicit, and no internal iteration is required. One just needs to implement the algorithm of weighted spatial averaging and, before damage is evaluated, replace the local equivalent strain by its non-local counterpart. 

The values of the non-local equivalent strain must be traced at the individual Gauss integration points of the finite element mesh, because these are the points at which stresses are computed. 

Thereby, the averaging integral in (3.54) is evaluated numerically as follows: 

**==> picture [278 x 28] intentionally omitted <==**

where _wχ_ is a coefficient containing the product of the determinant of the Jacobian and the integration weight of Gauss point _χ_ , and _Zx,χ_ is the aforementioned weight of non-local interaction between points _x_ and _χ_ , computed as: 

**==> picture [280 x 31] intentionally omitted <==**

In the previous two equations, subscript _x_ represents the receiver point under consideration, whereas indexes _χ_ and _ϕ_ correspond to source points. Since the chosen weighting function _Z_ 0 (3.50) has bounded support, _Zxχ_ vanishes if the distance between points _x_ and _χ_ is larger than the interaction radius _Rin_ . Therefore, the sums in (3.56) and (3.57) do not need to be taken over all Gauss points, but only over those that are located inside the domain of influence of point _x_ . 

Note that, at least around the damage process zone, one must always use an element size smaller than the interaction radius in order to account for the non-local interaction. Otherwise the damage model would become local. In this regard, an interaction radius relatively small can restrict the applicability of non-local models to small or medium size domains because of very fine meshes being necessary. 

Each Gauss point must have a non-local interaction table that gives access to its neighbours. This table must be constructed at the beginning of the problem and every time the model is remeshed. In this work, the search of neighbours is performed by means of a grid-based algorithm. A general rectangular grid is defined in the entire domain and all the integration points are positioned in the cells. This way, the neighbour search that must be performed for each Gauss point is restricted to a limited number of cells, 

Continuum Damage Mechanics 

80 

i.e. the ones that fall inside the domain of influence of the considered point (see Figure 3.18). 

**==> picture [278 x 266] intentionally omitted <==**

**----- Start of picture text -----**<br>
Solid domain<br>Search grid<br>Receiver Point Source Point<br>**----- End of picture text -----**<br>


Figure 3.18: Grid-based non-local search. 

The stress evaluation procedure (Table 3.3), repeatedly called during the incrementaliterative strategy, makes use of the non-local interaction tables when the non-local equivalent strain is computed. To obtain _ε_ ˘ _eq_ we first compute the local equivalent strains at all Gauss points, and then we calculate the non-local counterpart using (3.56). 

Finally, if one aims to obtain an iterative solver with quadratic convergence when working with a non-local damage model, it is necessary to construct the tangent stiffness matrix in a consistent manner, which slightly differs from the standard procedure already explained in Section 3.2.1. Here we derive the analytical expression of the tangent stiffness matrix, but we also note that some authors report complete consistency when obtaining the global tangent operators by finite-differencing the residuals [112]. 

Let us start expressing the tangent stiffness matrix (2.39) as a numerical integration 

Local and Non-local Damage Models 

81 

over the Gauss points of the model: 

**==> picture [279 x 32] intentionally omitted <==**

Using the stress-strain law in (3.9) and the classical strain-displacement relation _**ε**_ = _**Bu**_ ¯, we expand Eq (3.58) as follows: 

**==> picture [313 x 31] intentionally omitted <==**

Since the damage variable depends on the nodal displacements through the equivalent strain, we compute first the derivative of the damage variable with respect to the dis¯ ˘ placement vector _**u**_ . Taking into account that d = _g_ ( _r_ ) (3.7), _r_ depends on _εeq_ (3.55), and _ε_ ˘ _eq_ depends on _**u**_ ¯ through the interpolated strains, we use the chain rule to write the derivative of the damage variable for an integration point _x_ as follows: 

**==> picture [313 x 29] intentionally omitted <==**

where _gx[′]_[is the derivative of the damage evolution law with respect to the internal state] variable _r_ , and _lx_ is the loading-unloading factor that is 0 in an elastic loading or in an unloading regime, and 1 in a loading regime with growing damage: 

**==> picture [318 x 37] intentionally omitted <==**

Using expression (3.56), we can differentiate the non-local equivalent strain of a Gauss point _x_ with respect to the nodal displacements as 

**==> picture [419 x 33] intentionally omitted <==**

The derivative of the equivalent strain with respect to the strain vector _∂εeq,χ/∂_ _**ε** χ_ is a row matrix that depends on the chosen form of equivalent strain. 

At this point, we can already differentiate Eq (3.59), substitute (3.60) and rearrange terms to obtain the following expression: 

**==> picture [376 x 32] intentionally omitted <==**

Continuum Damage Mechanics 

82 

Note that the first term in (3.63) is the secant stiffness matrix, which coincides with the tangent matrix in an elastic loading or in an unloading regime ( _lx_ = 0). The second term in Eq (3.63) is the non-local part of the tangent stiffness matrix. Substituting (3.62) into (3.63) yields: 

**==> picture [357 x 70] intentionally omitted <==**

Defining for convenience the column matrix _**b**[c] x_[=] _**[ B]**[T] x_ _**[σ]**[e,x]_[, the row matrix] _**[ b]**[r] χ_[=] _[∂ε][e][q,χ]_ _**B** χ_ , _∂_ _**ε** χ_ and the coefficient _wxχ_ = _wxwχZxχ_ , Equation (3.64) can be rewritten as: 

**==> picture [298 x 28] intentionally omitted <==**

The double index of the sum, caused by the non-local interaction, implies that the term on the right part of Equation (3.65) cannot be assembled only from the standard element loop. Essentially, each pair of Gauss points _x_ and _χ_ contributes to the global stiffness matrix with a block of the same size as that of the classical element stiffness matrix. The difference is that the assembling routine differs from the usual one because in this case one needs to take into account the elements of both points _x_ and _χ_ (Figure 3.19). As a consequence, the global stiffness matrix is always non-symmetric and its bandwidth increases due to the non-local interaction. 

In order to avoid the additional non-zero entries that the non-local interaction introduces into the global stiffness matrix, one could neglect the non-local terms by using _wxχ_ = _δxχwχ_ , where _δxχ_ is the Kronecker delta. This way, Equation (3.65) reduces to 

**==> picture [304 x 27] intentionally omitted <==**

where the sum is performed over one index only. Note that the resulting local tangent matrix _**K** local_ is no longer consistent, and quadratic convergence is lost. 

Probably the most important issue caused by non-locality is the evolutionary character of the profile of the stiffness matrix. For the simulation of quasi-brittle materials like concrete or rock, the consistent stiffness matrix remains local through the elastic branch, 

Local and Non-local Damage Models 

83 

**==> picture [94 x 70] intentionally omitted <==**

**==> picture [144 x 71] intentionally omitted <==**

Figure 3.19: Non-local assembly process. 

and so the initial distribution of non-zero entries is the same as in the local case. However, when the damage threshold is exceeded and the damage zone starts propagating, new non-zero entries appear due to the non-local interaction between Gauss points belonging to different elements, and the profile of the stiffness matrix must be dynamically adapted. The number of additional non-zero entries depends on each particular case, but if a finer mesh is used in the expected softening zones, this number can be relatively high. 

Thereby, although the non-local damage model completely removes the pathological sensitivity to the mesh size and substantially alleviates the mesh-induced directional bias, it also increases the computational cost with respect to the local model, specially in the memory usage. 

Figure 3.20 compares the general scheme for the stress evaluation and stiffness matrix assembly between the local and non-local damage models. Apart from a greater memory usage, the non-local damage model requires performing more elements loops than the classical local model. 

Continuum Damage Mechanics 

84 

**==> picture [376 x 416] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Local damage model.<br>(b) Integral-type non-local damage model.<br>**----- End of picture text -----**<br>


**==> picture [42 x 131] intentionally omitted <==**

Figure 3.20: General scheme for the stress evaluation and stiffness matrix assembly. 

## **3.4 Examples** 

The fundamental concepts on damage mechanics theory and the most relevant aspects concerning the implementation of a non-local damage model within the Finite Element Method have already been introduced. 

This section is devoted to test and validate the presented damage models by solving two 

Examples 

85 

classical examples in the damage mechanics field: the three-point bending test, and the four-point shear test. 

For each one of the tests, we solve the problem with the two implemented damage models: the partially regularized local damage model, and the integral-type non-local damage model. The objective of this experiment is to analyse and compare both models, pointing out the strengths and limitations of the approaches, and assessing whether the non-local procedure is a valid and robust technique. 

Both examples are carried out under a 2D plane stress framework. The problems are solved as quasi-static step by step loading cases by means of the so-called modified Riks-Wempner arc-length strategy. Self weight is not taken into account. 

## **3.4.1 Three-Point Bending Test** 

This test is performed with a notched beam subjected to three-point bending (TPB). The beam has a square cross section of 40 _×_ 320 _mm_ , a span of 1280 _mm_ , and the notch is 3 _mm_ thick and extends over one tenth of the beam depth (see Figure 3.21). 

**==> picture [43 x 37] intentionally omitted <==**

Figure 3.21: TPB test. Geometry and boundary conditions. Dimensions in _mm_ . 

Plane stress conditions are assumed, and the geometry is meshed by means of bilinear 4-node quadrilaterals with 2 _×_ 2 integration points. 

We solve the problem using two different damage models: a partially regularized local damage model using the form of the equivalent strain proposed by Simo and Ju (3.26) 

Continuum Damage Mechanics 

86 

and the damage evolution law presented in (3.41), and a non-local damage approach, which is defined with the Mazars model regarding the equivalent strain (3.30) and the damage evolution law (3.46). The weighting function for the non-local interaction relies on a Gauss distribution function of bounded support (3.50). 

The material properties for the local damage model are summarized in Table 3.4 and have been obtained after some calibration, trying to fit the experimental results in [80]. For the non-local approach we use the same material parameters of [80] (see Table 3.5). 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|3_._85_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._24|
|Compressive strength (_σyc_)|4_._5_·_107 _N/m_2|
|Tensile strength (_σyt_)|3_._8_·_106 _N/m_2|
|Fracture energy (_Gf_)|100 _Nm/m_2|



Table 3.4: TPB test. Material properties for the Simo-Ju local model. 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|3_._85_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._24|
|Damage threshold (_ry_)|3_·_10_−_5|
|Residual strength in compression (_Rc_)|1_._25|
|Softening slope in compression (_Sc_)|1000|
|Residual strength in tension (_Rt_)|0_._95|
|Softening slope in tension (_St_)|9000|
|Characteristic length (_lc_)|0_._04 _m_|



Table 3.5: TPB test. Material properties for the Mazars non-local model. 

In order to assess the robustness of the models in terms of mesh sensitivity, we solve the problem for different spatial discretizations. Here we use three unstructured meshes of quadrilaterals with a refined area in the center (see Figure 3.22). The first model was obtained from a characteristic element size of _l[e]_ = 15 _mm_ , with a resultant mesh of 2024 elements and 2191 nodes. The second mesh, with 2679 elements and 2859 nodes, 

Examples 

87 

resulted from an element size of _l[e]_ = 7 _mm_ . The third model was obtained after defining a characteristic element size of _l[e]_ = 3 _mm_ and lead to a mesh of 6543 elements and 6772 nodes. 

(a) _l[e]_ = 15 _mm_ : 2024 elements. ~ eeeattt aiseiiaeEESaa Has eae BK O ee ee oO FF] | [>reas **Hee** SererECCTECR |[|eil **e** eediLy cannag **e** atitd eeee_a| fae = Ree (b) _l[e]_ = 7 _mm_ : 2679 elements. gsoeann:sate Saeee acooe posse ere care (c) _l[e]_ = 3 _mm_ : 6543 elements. Figure 3.22: TPB test. Spatial discretizations. Figure 3.23 shows the relation between the applied load and the vertical deflection of the beam. As one can see from the discontinued equilibrium curves, we had serious difficulties in tracing the response of a full test. The reason behind the convergence problems could be the use of a too global arc-length method. Indeed, to account for the localized nature of quasi-brittle failure, a more specific control parameter, like the Crack Mouth Opening Displacement (CMOD), could help improving the convergence near snap-back zones. 

That aside, if we look at the curves in Figure 3.23a we can see that, although there is 

Continuum Damage Mechanics 

88 

no relevant difference between the response obtained with the two coarser meshes, the peak force actually decreases with the finest mesh, and so the total dissipated energy. On the other hand, the response diagram in Figure 3.23b shows virtually the same peak load for the three meshes. 

Focusing on Figure 3.23b and comparing the computed results with the experimental solution [80], one sees that the non-local damage model properly captures the peak load of the expected solution for all meshes. 

Nonetheless, the post-peak branch of the numerical solution falls faster than in the reference solution, and even shows a certain amount of snap-back behaviour. Moreover, looking at the elastic branch of the responses, it seems that the stiffness degradation starts before in the numerical curves and, in consequence, the peak is slightly displaced to the right. 

Using an advancing technique purely based on the control of displacements instead of the mentioned arc-length method would reduce the differences in the post-peak response. 

From the depicted curves in Figure 3.23, it seems that the partially regularized local damage model is more sensitive to changes in the spatial discretization than the nonlocal model. In any case, though, the behaviour of the numerical model seems more brittle than the observed in experiments. 

In order to properly understand the response of both damage models, let us present some snapshots with the damage distributions. 

Figures 3.24 and 3.25 show an initial stage of damage progression in the local and nonlocal model, respectively. Just by looking at the size and shape of the damaged zone, we can clearly state the most differential trait of each model: localization on the one hand, and diffusion on the other. This concept is crucial to understand the observed behaviour in the load-deflection diagrams of Figure 3.23. 

Essentially, in the local approach, damage starts at the most stressed element and then “jumps” to the next one when the first is totally damaged. As a result, shape and direction of progression of damage strongly depends on the size and distribution of the elements of the mesh (see Figure 3.26). 

On the other side, in Figures 3.25a, 3.25b and 3.25c we see a very similar diffusive damaged area. Indeed, in the non-local approach, the damage size is controlled by 

Examples 

89 

**==> picture [268 x 242] intentionally omitted <==**

**----- Start of picture text -----**<br>
12<br>10<br>8<br>le=15mm<br>6<br>le=7mm<br>le=3mm<br>4 Experimental<br>2<br>0<br>0 0.05 0.1 0.15 0.2<br>Deflection (mm)<br>Force (kN)<br>**----- End of picture text -----**<br>


- (a) Partially regularized local damage model. 

**==> picture [268 x 262] intentionally omitted <==**

**----- Start of picture text -----**<br>
10<br>9<br>8<br>7<br>le=15mm<br>6 le=7mm<br>le=3mm<br>5<br>Experimental<br>4<br>3<br>2<br>1<br>0<br>0 0.05 0.1 0.15 0.2<br>Deflection (mm)<br>(b) Non-local damage model.<br>Force (kN)<br>**----- End of picture text -----**<br>


Figure 3.23: TPB test. Force-vertical deflection diagrams. 

Continuum Damage Mechanics 

90 

- (a) Model with _l[e]_ = 15 _mm_ . 

**==> picture [130 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) Model with l [e] = 7 mm .<br>**----- End of picture text -----**<br>


**==> picture [130 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(c) Model with l [e] = 3 mm .<br>**----- End of picture text -----**<br>


Figure 3.24: TPB test. Damage initiation in the local model. 

the interaction radius _Rin_ = _lc_ = 0 _._ 04 _m_ and thus even when reducing the size of the elements the damaged area is practically unaltered. 

The distinction between localization and diffusion can also be seen in the evolution of the damage variable depicted in Figure 3.27. 

While in the local damage model (Figures 3.27a, 3.27c and 3.27e) only a thin line of damage is propagated, i.e. a line of the size of the elements, in the non-local case (Figures 3.27b, 3.27d and 3.27f) damage occupies a wider area corresponding to the defined crack length _lc_ , and even some zones at each side of the notch are affected by the damage progression at the center, stressing the idea that the damage depends, not only on the state of the point under consideration, but also on the state of the neighbouring points. 

Examples 

91 

**==> picture [135 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Model with l [e] = 15 mm .<br>**----- End of picture text -----**<br>


**==> picture [131 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) Model with l [e] = 7 mm .<br>**----- End of picture text -----**<br>


**==> picture [129 x 11] intentionally omitted <==**

**----- Start of picture text -----**<br>
(c) Model with l [e] = 3 mm .<br>**----- End of picture text -----**<br>


Figure 3.25: TPB test. Damage initiation in the non-local model. 

**==> picture [282 x 15] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Model with l [e] = 15 mm . (b) Model with l [e] = 3 mm .<br>**----- End of picture text -----**<br>


Figure 3.26: TPB test. Zoom of damage growing in the local model. 

Continuum Damage Mechanics 

92 

**==> picture [401 x 178] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Local model. Initial stage. (b) Non-local model. Initial stage.<br>(c) Local model. Intermediate stage. (d) Non-local model. Intermediate stage.<br>(e) Local model. Advanced stage. (f) Non-local model. Advanced stage.<br>**----- End of picture text -----**<br>


Figure 3.27: TPB test. Evolution of damage propagation for _l[e]_ = 3 _mm_ . 

## **3.4.2 Four-Point Shear Test** 

In this second example a single-edge notched beam is subjected to four-point shear (FPS). The analysed beam has a square cross section of 100 _×_ 200 _mm_ , a span of 840 _mm_ , and the notch is 10 _mm_ thick and 40 _mm_ depth (see Figure 3.28). 

Figure 3.28: FPS test. Geometry and boundary conditions. Dimensions in mm 

Examples 

93 

Plane stress conditions have been again assumed, but in this case the geometry is meshed by means of linear 3-node triangular elements with one integration point. 

Again, we solve the problem using the two implemented damage models. The local approach is modelled with the same Simo and Ju model of the previous example, using the material parameters shown in Table 3.6. On the other hand, the non-local damage model is defined in this case using the equivalent strain form of the modified von Mises model (3.31), and the exponential damage evolution law that we presented in (3.40). The material parameters for the non-local approach have been obtained from [122] and are summarized in Table 3.7. 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|2_._8_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._1|
|Compressive strength (_σyc_)|3_._5_·_107 _N/m_2|
|Tensile strength (_σyt_)|3_._2_·_106 _N/m_2|
|Fracture energy (_Gf_)|140 _Nm/m_2|



Table 3.6: FPS test. Material properties for the Simo-Ju local model. 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|2_._8_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._1|
|Damage threshold (_ry_)|1_._5_·_10_−_4|
|Compressive to tensile strength ratio (_κ_)|10|
|Residual strength (_R_)|0_._8|
|Softening slope (_S_)|9000|
|Characteristic length (_lc_)|0_._01 _m_|



Table 3.7: FPS test. Material properties for the modified von Mises non-local model. 

As we did for the previous example, we solve the problem for three different unstructured meshes of triangles (Figure 3.29). The first model was obtained by refining the central area with a characteristic element size of _l[e]_ = 8 _mm_ , which resulted in a mesh of 2216 elements and 1264 nodes. The second model is represented by a mesh of 3502 elements 

Continuum Damage Mechanics 

94 

and 1923 nodes and was obtained from an element size of _l[e]_ = 5 _mm_ . The last spatial discretization is the result of defining a characteristic element size of _l[e]_ = 3 _mm_ and lead to a total of 7183 elements and 3796 nodes. 

**==> picture [146 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) l [e] = 8 mm : 2216 elements.<br>**----- End of picture text -----**<br>


**==> picture [147 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) l [e] = 5 mm : 3502 elements.<br>**----- End of picture text -----**<br>


**==> picture [145 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(c) l [e] = 3 mm : 7183 elements.<br>**----- End of picture text -----**<br>


Figure 3.29: FPS test. Spatial discretizations. 

In Figure 3.30 we represent the relation between the applied load and the Crack Mouth Sliding Displacement (CMSD) for each mesh. The curves in Figure 3.30a show that in this case the peak and the dissipated energy also decrease as the mesh is refinement. However, this reduction is not so clear as in the previous example. On the other hand, Figure 3.30b shows virtually the same peak load for the three meshes. Also, although one can notice some oscillations in the post-peak region for the mesh with _l[e]_ = 8 _mm_ and the mesh with _l[e]_ = 5 _mm_ , we can say that the residual force at the right part of the graph is actually similar for all cases. 

Examples 

95 

**==> picture [268 x 532] intentionally omitted <==**

**----- Start of picture text -----**<br>
90<br>80<br>70<br>60<br>50<br>40<br>30<br>le=8mm<br>le=5mm<br>20<br>le=3mm<br>10 Experimental<br>0<br>0 0.02 0.04 0.06 0.08<br>CMSD (mm)<br>(a) Partially regularized local damage model.<br>70<br>60<br>50<br>40<br>le=8mm<br>30<br>le=5mm<br>le=3mm<br>20<br>Experimental<br>10<br>0<br>0 0.02 0.04 0.06 0.08<br>CMSD (mm)<br>(b) Non-local damage model.<br>Force (kN)<br>Force (kN)<br>**----- End of picture text -----**<br>


Figure 3.30: FPS test. Force-Crack Mouth Sliding Displacement curves. 

Continuum Damage Mechanics 

96 

We note that the “bilinear-like response” obtained with the finest mesh of _l[e]_ = 3 _mm_ in Figure 3.30b is due to the low precision of the arch-length strategy used in this problems. Again, an advancing technique based on the control of displacements, like the CMSD, could be better suited in this case. 

Like in the previous case, we can compare the curves obtained from the non-local damage model with an experimental solution [30]. 

Thereby, focusing on Figure 3.30b we can see that the computed responses are very similar to the experimental solution, specially for the mesh with _l[e]_ = 3 _mm_ . Also, if one carefully looks at the elastic loading branches, it may be noticed that the numerical solutions exhibit a slightly stiffer behaviour than the experimental curve in the beginning of the degradation process. However, this difference is subtle and the post-peak branches converge to the results of the experiment. 

(a) Model with _l[e]_ = 8 _mm_ . 

**==> picture [196 x 164] intentionally omitted <==**

**----- Start of picture text -----**<br>
_<br>(b) Model with l [e] = 5 mm .<br>RREEEYa AEE<br>iL<br>aoeMae:<br>ndadad<br>(c) Model with l [e] = 3 mm .<br>**----- End of picture text -----**<br>


Figure 3.31: FPS test. Damage progression in the local model. 

Figures 3.31 and 3.32 show the damage progression in the peak-load region of the re- 

Examples 

97 

- (a) Model with _l[e]_ = 8 _mm_ . 

- (b) Model with _l[e]_ = 5 _mm_ . 

- (c) Model with _l[e]_ = 3 _mm_ . 

Figure 3.32: FPS test. Damage progression in the non-local model. 

sponse, for the local and non-local model, respectively. Like before, from the size and shape of the damaged zone, we can clearly see a more localized response in the local damage model, and a more diffusive response in the non-local case. 

Comparing the different damage patterns of the local model (Figures 3.31a, 3.31b and 3.31c), we notice an additional vertical damage line that appears only in the coarser meshes, stressing the idea that the local model suffers from mesh sensitivity. 

The damage patterns in the non-local damage model (Figures 3.32a, 3.32b and 3.32c) are virtually identical for all meshes. 

Finally, a comparative evolution of the damage variable with the mesh corresponding to _l[e]_ = 3 _mm_ is represented in Figure 3.33. We can see that, at the beginning of damage propagation (Figures 3.33a and 3.33b) both approaches show a very similar damage zone. The initially localized damage zone in the non-local model is mainly due 

Continuum Damage Mechanics 

98 

to the small interaction radius of this example _Rin_ = _lc_ = 0 _._ 01 _m_ . After that, the difference between both models is clear, with a wider damage mark of the order of the characteristic length in the non-local case. 

(a) Local model. Initial stage. (b) Non-local model. Initial stage. 

(c) Local model. Intermediate stage. (d) Non-local model. Intermediate stage. ae (e) Local model. Advanced stage. (f) Non-local ee model. Advanced stage. 

Figure 3.33: FPS test. Evolution of damage propagation for _l[e]_ = 3 _mm_ . 

## 4 Chapter 

## **Modelling Discontinuities in Porous Media** 

## **4.1 Introduction** 

Modelling the fluid flow in a multi-fractured porous domain implies taking into account that the cracks in the solid skeleton introduce preferential flow paths, apart from jumps in the displacement field. A proper understanding of discontinuities is crucial not only because they influence the behaviour of the local surroundings of the cracks, but also because they modify the global permeability and the mechanical response of the medium, specially whenever it undergoes a crack growth process. 

In the last decades, many efforts have been made to develop numerical models for the accurate analysis of discontinuities in solids and porous media. 

The extended finite element method (XFEM) has obtained notable attention in the past years [14, 83, 120, 146]. The essential idea of the XFEM is to enhance the solution by decomposing the displacement field into a continuous and a discontinuous part. The discontinuity is captured by means of enrichment functions that introduce the jumps in the displacement field. The most remarkable advantage of the method is that there is no need to explicitly represent cracks in a mesh provided that enriched nodes are considered (Figure 4.1a). This avoids the necessity of remeshing during crack growth, 

Modelling Discontinuities in Porous Media 

100 

but in return it demands a higher computational cost in terms of number of degrees of freedom and numerical integration. 

The present work focuses on the methods purely based on the finite element method (FEM). In this category, numerous methods can be found in the literature, but two main subgroups can be distinguished: the "smeared crack" approaches (Figure 4.1b), continuum based methods in which the influence of developing fractures is incorporated into the constitutive stress-strain law [28, 66, 85, 119], and the "discrete crack" models (Figure 4.1c), in which each single discontinuity is represented explicitly [34, 56, 60, 67]. 

**==> picture [67 x 77] intentionally omitted <==**

(a) XFEM. Enriched nodes have been coloured red. 

**==> picture [60 x 78] intentionally omitted <==**

**==> picture [115 x 112] intentionally omitted <==**

**==> picture [284 x 28] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) FEM: smeared crack ap- (c) FEM: discrete crack ap-<br>proach. proach.<br>**----- End of picture text -----**<br>


Figure 4.1: Different methods to represent discontinuities. 

Since Goodman _et al._ proposed the "zero-thickness" interface element to describe the mechanical behaviour of pre-existing joints in rock masses [63] many authors have developed strategies to adapt this element for the solution of fracture processes in coupled 

Introduction 

101 

## solid-pore fluid problems. 

Three different types of zero-thickness interface elements can be found in the literature concerning the way the fluid is modelled: single, double and triple noded elements. The single-noded element is the simplest one and only considers longitudinal conductivity with no pressure drop across the interface [5]. The triple-noded element was meant to include the effect of the transversal conductivity through the discontinuity [64]. The two nodes at each side of the interface represent the potentials in the pore pressure, while the third node in the middle stores the average potential of longitudinal fluid through the fracture. Finally, the double-noded elements take into account both types of conductivity but the external nodes variables substitute the middle node. Ng and Small used this double-noded zero-thickness interface element to model flow problems with pre-existing discontinuities, but did not consider hydraulic potential drop between the two interface walls [98]. Segura and Carol introduced the transversal conductivity in double-noded zero-thickness elements to account for the exchange of fluid between the discontinuity and the porous media [128]. 

Regarding the mechanical behaviour of fractures, two different approaches are typically used: the linear elastic fracture mechanics (LEFM), and the non-linear fracture mechanics (NLFM). LEFM was first proposed to solve fracture propagation problems by means of remeshing without considering a fracture process zone (FPZ) before the crack tip. This approach is applicable in large structures where the size of the FPZ is negligible. However, for quasi-brittle analyses, the consideration of a non-linear fracture process zone where the energy is dissipated before it completely fails was found to be essential. In those cases NLFM is usually applied and a softening law relates the cohesive stress to the crack opening in the FPZ. The first technique based on the cohesive fracture model was originally introduced by Barenblatt for brittle materials [8, 9] and by Dugdale for plastic materials [50]. Hillerborg _et al._ developed the first fictitious crack model for Mode I fracture [68]. It was extended later for the mixed mode fracture, from which Camacho and Ortiz proposed a suitable fracture criterion that is widely used in the literature [27]. 

One of the most important parts in the modelling of fracture propagation is the criterion for determining the direction of the crack growth. Some methodologies are based on the local evaluation of the stress field at the crack tip, such as the maximum circumferential stress [52] and the maximum principal stress criteria [22, 79]. Others measure the energy 

Modelling Discontinuities in Porous Media 

102 

distribution at the fractured zone, e.g. the minimum strain energy density criterion [129] or the maximum strain energy release rate criterion [72]. Finally, some authors have developed crack growth criteria based on continuum damage mechanics [140] and, more recently, combined with level sets [33, 94]. 

In order to minimize the mesh-induced directional bias, in this thesis we use the non-local damage model of Chapter 3 in combination with a discrete crack approach, in which discontinuities are represented by quasi-zero-thickness interface elements. A special utility allows us inserting new interface elements according to the computed damage map, finishing with a remeshing of the model to ensure a conformal spatial discretization. The low permeability and high compressibility of fluid-driven fracture propagation problems makes the pressure field oscillate spuriously if equal order interpolation elements are used without stabilization. Here we solve the solid-pore fluid interaction problem with the FIC-FEM stabilized formulation presented in Chapter 2. 

The present chapter is organized as follows. First, the formulation of the developed quasi-zero-thickness interface elements is detailed, explaining the mechanical behaviour of the fracture and the model governing the fluid flowing in it. After that, we present the proposed fracture propagation technique combining the non-local damage model with the interface elements. Finally, two plane-strain examples are solved to test the accuracy of the crack propagation methodology, and one additional case is included to show the performance of the generalized 3D formulation. 

## **4.2 Quasi-zero-thickness interface elements** 

Before going through the governing equations of the interface elements, let us clarify the “quasi-zero-thicknes” adjective used in this work. 

Essentially, the quasi-zero-thickness interface elements follow the same idea as the double-noded zero-thickness interface elements of [98] or [128] but, in fact, they are surface elements in 2D and volume elements in 3D and so they can be defined using a width larger than zero if necessary. 

This fact is particularly useful in case we need to define pre-existing joints with a certain opening, or when we want to capture the distribution of the fluid flow in the intersection between fractures (Figure 4.2). 

Quasi-zero-thickness interface elements 

103 

Figure 4.2: Fluid flow inside a crack with bifurcation. 

However, one can also define them with a zero thickness. In such case, the interface elements work as double lines in 2D or as a double surfaces in 3D, which is useful to represent very thin layers or closed faults that could eventually open (Figure 4.3). 

Figure 4.3: Consolidation problem with a vertical fault. 

Modelling Discontinuities in Porous Media 

104 

Three different types of geometries for joint elements have been implemented: a 4-node quadrilateral interface element with two Lobatto integration points for bi-dimensional problems (Figure 4.4a), a 6-node wedge interface element with three Lobatto points and an 8-node hexahedral interface element with four Lobatto points for three-dimensional cases (Figures 4.4b and 4.4c). The 6-node prismatic interface element is used in problems where the porous domain is meshed with tetrahedra, whereas the hexahedral interface element is meant for problems meshed with hexahedra in the rest of the domain. 

**==> picture [44 x 49] intentionally omitted <==**

**==> picture [45 x 50] intentionally omitted <==**

**==> picture [364 x 132] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Four-node quadrilateral.<br>(b) Six-node triangular prism. (c) Eight-node hexahedron.<br>**----- End of picture text -----**<br>


Figure 4.4: Quasi-zero-thickness interface elements. 

It is important to note that the choice of Lobatto integration over Gauss integration responds to the spurious traction oscillations that may appear in joint elements with a stiff cohesive zone [76, 125, 131]. As it has been reported in the literature, a common strategy to overcome this issue is to make use of reduced Lobatto integration along the mid plane of the element [47, 65, 137]. 

Furthermore, it is essential to understand that the quantities of interest in the interface elements are in a different coordinate system from the one used to solve the standard elements in the porous domain. 

Indeed, while for the porous domain we work in a unique coordinate system to obtain the global deformations of the structure, for each interface element we aim at obtaining the normal and tangential relative displacements at any point along the crack, and so we must work in the local coordinate system of the fracture. Moreover, as it is usually done for beams (in 2D) or shells (in 3D), this local system is located at the mid plane 

Quasi-zero-thickness interface elements 

105 

of the interface element, and it is necessary to distinguish the upper and lower faces of the joint to compute its relative displacements (Figure 4.5). 

**==> picture [45 x 27] intentionally omitted <==**

**----- Start of picture text -----**<br>
,<br>,<br>,<br>**----- End of picture text -----**<br>


Figure 4.5: Scheme of a generic hexahedral interface element. Global and local coordinate system. 

Let us first define the vector of relative displacements in a joint _**δ**_ as the difference between the displacements at the upper and lower faces: 

**==> picture [262 x 13] intentionally omitted <==**

The vector of nodal displacements in global axis _**u**_ **¯** is written as in Chapter 2. For a general 3D case it reads: 

**==> picture [317 x 25] intentionally omitted <==**

where _n_ is the number of nodes of the element. 

We also define the matrix of shape functions for the displacements at the interface element _**N** u,I_ as: 

**==> picture [346 x 52] intentionally omitted <==**

in which the shape functions _Ni_ are preceded by a negative sign for nodes located at the lower face of the joint and by a positive sign for nodes at the upper face. 

Modelling Discontinuities in Porous Media 

106 

Thereby, taking into account the definition in Eq (4.1), we can compute the vector of relative displacements at any point of the interface as: 

**==> picture [412 x 123] intentionally omitted <==**

Also, in order to obtain the normal and tangential relative displacements, we need to transform the computed relative displacements from the global to the local coordinate system. To do so, we need to compute a rotation matrix _**λ**_ from the direction cosines between the global ( _x_ , _y_ , _z_ ) and local ( _x_ ˇ, _y_ ˇ, _z_ ˇ) coordinates (Figure 4.5). In essence: 

**==> picture [321 x 51] intentionally omitted <==**

The transformation is written then as: 

**==> picture [372 x 51] intentionally omitted <==**

In the following paragraphs, we present the governing equations of the implemented quasi-zero-thickness interface elements. 

Similarly to the formulation presented for the continuum solid-pore fluid mixture, two equilibrium equations govern the behaviour of the interface elements. One equation deals with the mechanical behaviour of the crack, whereas the other equation describes the balance of mass within the fracture. 

It is important to stress here that the formulation of the interface elements is completely compatible with the displacement-pore pressure formulation presented in Chapter 2. Thereby, one can combine both type of elements in the same problem and solve them together in a unique coupled system of equations. 

Quasi-zero-thickness interface elements 

107 

## **4.2.1 Mechanical behaviour of the fracture** 

The starting point is the strong form of the balance of momentum equation for the solidfluid mixture (2.5) defined at the fracture plane. We discretize the relative displacements at the joint and integrate the equation by parts to obtain the following weak form: 

**==> picture [345 x 30] intentionally omitted <==**

Equation (4.7) shows the same structure as Equation (2.24), but the resulting matrices and vectors are defined in the coordinate system in which the interface element is solved, i.e.: 

**==> picture [368 x 29] intentionally omitted <==**

**==> picture [321 x 30] intentionally omitted <==**

where _**B** I_ is the deformation matrix of the interface element, defined as: 

**==> picture [257 x 27] intentionally omitted <==**

with _tI_ being the distance between the two faces of the interface (Figure 4.5). If the interface is closed, then we use a small value _tI_ = _tmin_ to avoid dividing by zero. 

It is important to note that the terms in (4.7) are integrated over a quasi-zero-thickness interface element and thus: 

**==> picture [255 x 13] intentionally omitted <==**

The effective stress vector in the fracture _**σ**[′] I_[must][be][computed][first][in][the][local][coor-] dinate system and then transformed into the global one as 

**==> picture [365 x 54] intentionally omitted <==**

_T_ Also, the vector _**m** I_ = 0 _,_ 0 _,_ 1 comes from the definition of the local effective � � stresses in the joint: 

**==> picture [316 x 51] intentionally omitted <==**

Modelling Discontinuities in Porous Media 

108 

The constitutive law governing the mechanical behaviour of the interface element is a bilinear cohesive fracture model based on the fracture criteria of Camacho and Ortiz [27] and Song _et al._ [133] (see Figure 4.6). 

The fundamental concept behind the cohesive fracture model is that, although the interface element represent discontinuities, the stress field can still be transferred through the crack sides of the cohesive zone (Figure 4.6a). From the micro-mechanical point of view, the cohesive zone is the local region ahead of the crack tip where the micro-voids or micro-cracks initiate, grow and coalesce with the main crack. 

**==> picture [365 x 200] intentionally omitted <==**

**----- Start of picture text -----**<br>
Real crack tip Cohesive Zone<br>Fictitious crack tip<br>(a) Cohesive fracture zone.<br>(b) Bilinear cohesive law.<br>**----- End of picture text -----**<br>


Figure 4.6: Bilinear cohesive fracture model. 

Figure 4.6b relates the normalized equivalent stress _ςeq,I_ with the internal state variable _ϱ_ and shows that the evolution of the cohesive zone is, in essence, an irreversible damage process. 

The variable _ϱ_ plays the same role as the internal scalar variable _r_ ( _t_ ) (Eq (3.13)) and thus it characterizes the maximum strain level reached in the previous history of the joint: 

**==> picture [326 x 30] intentionally omitted <==**

Here _ϱy_ is the damage threshold of the joint, a material parameter ranging from 0 to 1 that represents the value of strain at which damage starts. 

The chosen model is suitable for mixed mode fracture and thus the evolution of the 

Quasi-zero-thickness interface elements 

109 

damage propagation in the cohesive zone must depend on the simultaneous activation of tangential and normal relative displacements. We define the equivalent strain of the interface _εeq,I_ as 

**==> picture [280 x 31] intentionally omitted <==**

where _δc_ is the critical relative displacement, i.e. it is the relative displacement at which the cohesive zone stops transmitting forces. 

Expression (4.14) implies that _ϱy ≤ ϱ_ ( _t_ ) _≤_ 1, which shows that, in fact, the internal state variable _ϱ_ can also be seen as the damage variable of the cohesive fracture zone. 

Similarly as in Eq (4.15) we can define the normalized equivalent stress as: 

**==> picture [302 x 31] intentionally omitted <==**

with _σy_ being the yield stress at the cohesive zone, i.e. the stress at which the cohesive zone starts damaging. Each local component of the effective stress is computed as follows: 

**==> picture [308 x 29] intentionally omitted <==**

The definition of the equivalent strain in (4.15) and the expressions of the stresses in (4.17) are valid only for cases in which the two faces of the joint are separated by a certain gap ( _tI >_ 0). 

However, in any interface it is possible that both faces push each other whenever they get in contact (Figure 4.7). In such cases, it is necessary to redefine some variables to properly model the mechanical behaviour of the joint. 

Thereby, when there is contact between the two opposite bodies of an interface ( _tI <_ 0) the tangential effective stresses in Eq (4.17) are defined as the shear strength of the Mohr-Coulomb criterion: 

**==> picture [270 x 13] intentionally omitted <==**

**==> picture [275 x 14] intentionally omitted <==**

where _µF_ is the friction coefficient, and _cl_ and _cm_ are the tangential cohesive stresses. They are computed as: 

**==> picture [270 x 30] intentionally omitted <==**

Modelling Discontinuities in Porous Media 

110 

**==> picture [409 x 163] intentionally omitted <==**

**----- Start of picture text -----**<br>
Y-DISPLACEMENT<br>0<br>: -0.0043191<br>-0.0086383<br>-0.012957<br>-0.017277<br>y y<br>z x z x<br>eh eh<br>(a) Mesh in the original position. (b) Vertical displacement after contact.<br>**----- End of picture text -----**<br>


Figure 4.7: Contact between two beams connected by interface elements. 

**==> picture [275 x 29] intentionally omitted <==**

We note that the absolute values in Equations (4.18) and (4.19) are used to represent that, regardless of the sign of the stresses, the cohesion and the friction terms have always the same sign. 

The definition of the normal component of the effective stress _σn[′]_[also][changes][during] contact. It imitates the expression for the normal force between two neighbouring rigid particles in the Discrete Element Method (DEM) [108, 124]. 

Let us assume that the interface is composed by linear springs with a stiffness per unit area of _Kn_ = _E/_ ( _ϱyδc_ ) (Figure 4.8). 

The normal effective stress in contact is computed as: 

**==> picture [275 x 29] intentionally omitted <==**

where _II_ is the indentation between the two faces of the joint. 

Also, supposing that during contact _δn_ = 0, the equivalent strain of the interface is computed as: 

**==> picture [268 x 30] intentionally omitted <==**

Finally, in order to compute the stiffness matrix of the interface element _**K** I_ , one just 

Quasi-zero-thickness interface elements 

111 

**==> picture [32 x 43] intentionally omitted <==**

Figure 4.8: Scheme of a contact with interface elements. 

needs to derive the term �Ω _I_ _**[B]** I[T]_ _**[σ]**[′] I[d]_[Ω] _[I]_[in][Equation][(4.7)][with][respect][to][the][displace-] ment vector _**u**_ ¯: 

**==> picture [338 x 29] intentionally omitted <==**

where _**D** I_ is the tangent constitutive tensor of the interface in the global reference system. It is obtained from: 

**==> picture [264 x 15] intentionally omitted <==**

The matrix _**D**_[ˇ] _I_ in Eq (4.25) is the local tangent constitutive tensor per unit length, which can be obtained by deriving the stresses in the interface with respect to their relative displacements: 

**==> picture [324 x 51] intentionally omitted <==**

Considering that the definition of the stresses in Eq (4.17) are used, the expressions defining the local tangent constitutive tensor of the interface (4.26) can be analytically derived. Taking into account that _ϱ_ = _ϱ_ ( _δl, δm, δn_ ), we can distinguish two different scenarios: 

Modelling Discontinuities in Porous Media 

112 

- Elastic loading or unloading process ( _∂ϱ_ = 0): 

**==> picture [386 x 187] intentionally omitted <==**

- Loading process with growing damage ( _∂ϱ >_ 0): 

## **4.2.2 Fluid in the fracture** 

The fluid flow in the fracture is modelled by the following mass balance equation: 

**==> picture [320 x 16] intentionally omitted <==**

where 

**==> picture [284 x 29] intentionally omitted <==**

**==> picture [346 x 76] intentionally omitted <==**

¯ The pressure at the joint is interpolated as we did in Chapter 2, i.e. _p_ = _**N** p_ _**p**_ . Thereby, the compressibility matrix _**C** I_ (4.30) is virtually the same as the one used in standard elements (2.27), with the only difference being the domain of integration _δ_ Ω _I_ (4.11). 

Matrix _**H** I_ (4.31) introduces the enhanced permeability of the fractures in the porous medium. It is assumed that the longitudinal flow is coupled to the interface opening according to the so-called cubic law [132, 143]. 

Quasi-zero-thickness interface elements 

113 

The name “cubic law” comes from the expression of the flow in the joint as a function of the cube of the local aperture _tI_ . Such relation arises from the Reynolds equation [121]: 

**==> picture [373 x 30] intentionally omitted <==**

Equation (4.33) originates from the lubrication theory, which was formulated for fluid flows in narrow and smooth parallel plates. The terms on the left hand side correspond to the flow induced by a pressure gradient or Poiseuille flow (Figure 4.9a), whereas the terms on the right hand side are given by the shear flow or Couette flow, generated by the movement of the upper face with a velocity _V_ (Figure 4.9b). 

**==> picture [279 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) Poiseuille flow. (b) Couette flow.<br>**----- End of picture text -----**<br>


Figure 4.9: Velocity profiles between two smooth parallel plates. 

The cubic law is naturally incorporated in the matrix _**H** I_ of Equation (4.31) through the gradient of the matrix of shape functions ( _**∇N** p_ ) _I_ and the intrinsic permeability matrix _**k** I_ . 

The matrix ( _**∇N** p_ ) _I_ is defined as [78]: 

**==> picture [381 x 51] intentionally omitted <==**

The two first rows in (4.34) allow computing the derivatives of the pressure with respect to the tangential directions of the joint _l_ and _m_ , whereas the third row approximates 

Modelling Discontinuities in Porous Media 

114 

the derivative of the pressure in the normal direction _n_ by means of the pressure drop between the upper and lower faces of the interface: 

**==> picture [265 x 28] intentionally omitted <==**

In essence, the computation of the derivatives of the shape functions with respect to the three coordinates of the space _l_ , _m_ and _n_ requires obtaining first the complete Jacobian matrix and then compute its inverse. The problem is that, when working with quasizero-thickness interface elements, the Jacobian matrix can be singular if the joint is nearly closed. Because of that, the tangential derivatives of the pressure are computed in the mid plane of the joint, i.e. a line in 2D and a triangle or a quadrilateral in 3D (Figure 4.10). The normal derivative is computed through finite differences as in (4.35). 

**==> picture [182 x 139] intentionally omitted <==**

**----- Start of picture text -----**<br>
, ,<br>, ,<br>**----- End of picture text -----**<br>


Figure 4.10: Natural coordinates at the mid plane of an hexahedral interface element. 

If we define _ξ_ and _η_ as the two natural coordinates at the mid plane of an interface element (Figure 4.10), we can obtain the expression for the derivative of the shape functions with respect to the local directions _l_ and _m_ by deriving the shape functions with respect to the natural coordinates and applying the chain rule. For every node _i_ we have: 

**==> picture [288 x 72] intentionally omitted <==**

Quasi-zero-thickness interface elements 

115 

or equivalently: 

**==> picture [410 x 36] intentionally omitted <==**

And so the derivatives of the shape functions with respect to the local tangential directions are computed as: 

**==> picture [309 x 30] intentionally omitted <==**

The terms of the Jacobian matrix _**J** I_ in (4.37) are obtained by rotating the derivatives of the global coordinates with respect to the natural coordinates: 

**==> picture [356 x 51] intentionally omitted <==**

Regarding the intrinsic permeability matrix of the fracture _**k** I_ (4.31), we make a distinction between the permeability in the tangential directions and the permeability in the normal one. 

**==> picture [273 x 51] intentionally omitted <==**

where _kl_ and _km_ are defined as [132]: 

**==> picture [258 x 27] intentionally omitted <==**

The parameter _kn_ represents the transversal permeability of the fracture, which is usually given a value similar to the intrinsic permeability of the porous domain. 

Let us note that, although we write a 2 for the exponential of the joint width _tI_ in Eq (4.41), we actually have the expected 3 of the cubic law if we take into account the additional _tI_ appearing when integrating over the interface domain (4.11). 

The anisotropic permeability introduced in the permeability matrix of the joint (4.40) can become a problem whenever we have intersection of different fractures. Indeed, an 

Modelling Discontinuities in Porous Media 

116 

intersection between two cracks represents a void in the porous domain through which fluid can flow towards any direction. In this case there are not clear longitudinal and transversal directions for the flow because all directions can be preferential, and thus a variant of the conventional interface element with a modified intrinsic permeability matrix _**k**_[ˆ] _I_ is introduced: 

**==> picture [269 x 51] intentionally omitted <==**

Figure 4.11 exemplifies the effect that a set pre-defined interface elements introduce in a porous medium. Water is injected at a constant rate of 1 _m/s_ through the left side of the horizontal crack, and after only 0 _._ 001 _s_ , the displacement and pore pressure fields show noticeable jumps between each side of the joints (see Figures 4.11a, 4.11b, 4.11c and 4.11d). 

Figure 4.11 also displays an application case of the modified intrinsic permeability matrix (4.42). The cyan square element in Figure 4.11e has the same structure as the rest of interface elements in pink, and just by using matrix _**k**_[ˆ] _I_ , a proper flow distribution is captured at the intersection of each fracture (Figure 4.11f). 

Quasi-zero-thickness interface elements 

117 

**==> picture [195 x 134] intentionally omitted <==**

**----- Start of picture text -----**<br>
0.000003<br>0.0000025<br>0.000002<br>0.0000015<br>0.000001<br>0.0000005<br>0<br>0 5 10 15<br>-5E-07<br>-0.000001<br>Distance (m)<br>Displacement-X (m)<br>**----- End of picture text -----**<br>


- (a) Contour lines of displacement-X. 

- (b) Displacement-X along selected line. 

**==> picture [195 x 134] intentionally omitted <==**

**----- Start of picture text -----**<br>
80<br>70<br>60<br>50<br>40<br>30<br>20<br>10<br>0<br>0 5 10 15<br>Distance (m)<br>Water Pressure (N/m2)<br>**----- End of picture text -----**<br>


**==> picture [169 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(c) Contour lines of water pressure.<br>**----- End of picture text -----**<br>


**==> picture [183 x 11] intentionally omitted <==**

**----- Start of picture text -----**<br>
(d) Water pressure along selected line.<br>**----- End of picture text -----**<br>


**==> picture [201 x 26] intentionally omitted <==**

**----- Start of picture text -----**<br>
(e) Detail of a intersection between two<br>cracks.<br>**----- End of picture text -----**<br>


**==> picture [177 x 11] intentionally omitted <==**

**----- Start of picture text -----**<br>
(f) Fluid flow at a crack intersection.<br>**----- End of picture text -----**<br>


Figure 4.11: Water flow in a pre-existing fractures network. 

Modelling Discontinuities in Porous Media 

118 

## **4.3 Fracture propagation approach** 

After introducing the quasi-zero-thickness interface element as a means to represent discontinuities in the porous domain, the current section is devoted to describe the proposed method for propagating cracks through the solid skeleton. 

There are essentially two main strategies for solving fracture propagation problems by means of interface elements. The first technique has been widely used in the literature [25, 26, 29, 133] and is based on the definition of an initial finite element mesh in which the joint elements are introduced at every possible location between standard finite elements (Figure 4.12). 

**==> picture [156 x 156] intentionally omitted <==**

Figure 4.12: Finite element mesh with interface elements at all possible edges. 

Such an approach requires defining a number of interface elements much larger than the strictly necessary, which may alter the mechanical response of the solid skeleton and condition the permeability of the whole porous domain. In addition, when predicting the propagation of fractures, it is crucial that the type of the mesh does not affect the direction of propagation. With the mentioned method the influence of the mesh on the obtained crack path is virtually unavoidable. 

The alternative method that can be encountered in the literature is based on the insertion of new interface elements where the strength of the porous domain reaches a threshold value [7, 78, 147]. This method provides a realistic approach for the simulation of dynamic crack propagation, but it requires remeshing the potential fracture zone after the insertion of every new interface element. 

Fracture propagation approach 

119 

In this work, the fracture propagation strategy is based on the latter approach. 

## **4.3.1 Crack path estimation** 

Essentially, we combine the integral-type non-local damage model introduced in Chapter 3 with the presented quasi-zero-thickness interface elements (Section 4.2). The fundamental idea is to use the information of the damaged points around the crack tip to determine, as accurately as possible, the direction of the fracture propagation. 

The implemented method relies on two main assumptions: 

- Any crack growth process must start from a predefined crack tip. 

- Such propagation can either follow one direction, or bifurcate in two. 

The first assumption implies that the implemented method can not model crack initiation. However, this is not major inconvenience because one can always use the damage model to determine the initial damage pattern, and then insert crack tips to model fracture propagation. The second assumption means that the presented approach allows capturing crack branching for any propagation process. 

Next, we present the main steps configuring the crack path estimation technique. 

## **Crack tip neighbours search** 

The information of the damaged elements around the crack tip is used to estimate the direction of the crack path. In order to do so, a neighbour search must be performed at the beginning of each time step. 

Let _lp_ be the propagation length, a material parameter determining the domain of influence of the fracture Ω _f_ . For each fracture tip, we search the neighbouring Gauss points falling inside a circle (or a sphere in 3D) of radius _lp_ (Figure 4.13a). 

At this stage we also divide the domain of influence Ω _f_ into quarters and distribute the neighbour points between the Top-Front-Quarter (TFQ) and the Bottom-Front-Quarter (BFQ) (Figure 4.13b). 

Modelling Discontinuities in Porous Media 

120 

**==> picture [392 x 212] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) Division of Ω f into quarters.<br>(a) Grid-based search around the tip.<br>**----- End of picture text -----**<br>


Figure 4.13: Crack tip neighbours search. 

## **Crack growth** 

In order to determine the start of the crack growth, we estimate the damage at the crack tip d _t_ and compare it with a prefixed threshold value, called propagation damage d _p_ . 

Damage at the crack tip is defined here as a non-local measure of the damage inside the region of influence of the fracture Ω _f_ . This measure follows the same idea as the non-local average performed in Chapter 3 for the equivalent strain (Eq (3.56)). The damage at the tip is evaluated as: 

**==> picture [294 x 28] intentionally omitted <==**

where _wχ_ is the coefficient with the product of the determinant of the Jacobian and the integration weight of Gauss point _χ_ , and _Ztχ_ is the weight of non-local interaction between the tip and any other point _χ_ , defined as in Eq (3.57), i.e. 

**==> picture [279 x 32] intentionally omitted <==**

When the damage at the tip exceeds the threshold value (d _t >_ d _p_ ) the crack is supposed to grow, but it is still necessary to discern the direction of propagation. First, we com- 

Fracture propagation approach 

121 

pute the potential location for the new crack tip by weighted averaging the coordinates of all the points in front of the crack tip. 

Unlike the average performed in (4.43), the coordinates are weighted in terms of the damage of each Gauss point under consideration. In essence, the coordinate of the new crack tip _**x** t_ is computed as: 

**==> picture [323 x 27] intentionally omitted <==**

where _Z_[ˆ] _tχ_ is a modified version of the weighting function in (4.44) as 

**==> picture [264 x 31] intentionally omitted <==**

Once the coordinates for the new crack tip have been found, it is necessary to check whether it is at a feasible location of the domain. Basically, the new crack tip must fulfil the following two conditions: 

- It must fall inside an existing element in the mesh. 

- The average damage at the element must be larger than a minimum value d _[e] ≥_ d _min_ 

The first condition is obviously necessary to avoid propagating the fracture outside the limits of the model, whereas the second condition prevents propagating the fracture towards undamaged regions. 

In case one of the two conditions is not fulfilled, we recalculate two possible coordinates for the crack tip, distinguishing the elements of the top-front-quarter (TFQ) and the bottom-front-quarter (BFQ) as 

**==> picture [307 x 52] intentionally omitted <==**

Of course, the two coordinates in (4.47) must also fulfil the aforementioned conditions to make sure they are feasible crack tips. In this regard, three different scenarios are possible. If the two new tips are valid, then the fracture bifurcates in two separate ways. If only one of the computed coordinates is feasible, then the fracture will propagate only towards that valid tip. Finally, if neither of the tips is valid, the fracture does not grow. 

Modelling Discontinuities in Porous Media 

122 

## **4.3.2 Interface elements insertion** 

After estimating direction of the propagating fracture, one can actually insert new interface elements into the model. 

To do so, it is essential to make sure that the mesh after the insertion is conformal, and that the primary unknowns and internal variables are properly mapped from the old mesh to the new one. 

## **Insertion and remeshing** 

For the insertion of new interface elements and regeneration of the model, we take advantage of the meshing capacities of the pre and post-processor GiD [1]. 

Thereby, not only GiD meshes the geometry at the beginning of the numerical simulation, but it is also GiD which generates the new spatial discretization every time we need to adapt the model with new interface elements. 

The flexible design of GiD makes it possible to expand its basic features by means of customizable modules called “Problem-Types”. Such modules are quite versatile and allow, not only to create a graphical user interface (GUI) for defining the conditions and properties of the problem, but also to perform automatic GiD procedures by means of Tcl scripts. 

Essentially, the computation of any problem starts with a Tcl procedure from the Problem-Type that writes the necessary input files with the conditions defined in the GUI: 

- _ProjectName_ .mdpa: lists the material properties, the nodes and the elements. 

- ProjectParameters.json: saves the parameters for the project such as the dimension, the time step or the linear solver, and writes all the boundary conditions of the problem. 

- FracturesData.json: registers all the pre-defined fracture tips in the model. 

Fracture propagation approach 

123 

The “FracturesData.json” file must contain all the necessary information of the preexisting crack tips in the model so as to perform the crack path estimation procedure already explained. This information includes the coordinates of all the points defining every crack tip. 

In this regard, it is essential to have a clear image of how the crack tip is designed. In 2D any crack tip can be defined from just three points (Figure 4.14a), but in 3D this is not so straightforward. In the present thesis, the three-dimensional fracture tip is assumed to be conical, but it is approximated as a pyramid of quadrilateral base (Figure 4.14b). 

**==> picture [271 x 61] intentionally omitted <==**

**----- Start of picture text -----**<br>
(a) 2D tip.<br>(b) 3D tip.<br>**----- End of picture text -----**<br>


Figure 4.14: Scheme of the implemented crack tips. 

The input files are read by the main Python script of the program, and then the first step of the problem is solved in a C++ strategy. Next, another method is called to check whether any of the pre-existing fractures propagates. If none of the crack tips progresses, the calculation simply jumps to the next time step, and in case one of the fractures propagates, the new crack tips are estimated following the procedure presented in Section 4.3.1, and the information of these tips is saved in another file called “PropagationData.tcl”. 

Python calls then another Tcl procedure from the Problem-Type, which reads the “PropagationData.tcl” file, draws the new fractures in the geometry, and regenerates the model making sure that the new mesh is conformal. 

In this work, a uniform mesh is generated for simplicity, but one could also apply the method proposed in [48] for adapting the mesh of the whole domain in a non-structured 

Modelling Discontinuities in Porous Media 

124 

manner according to the stress state. The spatial discratization could be then finer around the cracks and coarser in those areas where no damage appears, thus making the computation cost more efficient. 

After the mesh is generated, the input files are updated with the new configuration of fractures, and the problem continues in the next time step after a mapping of variables between the old and new meshes is performed. 

In the end, the implemented fracture propagation technique is based on the combination of the modelling capabilities of GiD with the solving power of Kratos. It is schematically represented in Figure 4.15. 

**==> picture [211 x 232] intentionally omitted <==**

Figure 4.15: Flow chart of the fracture propagation technique. 

## **Mapping of variables** 

After the new model is generated, a final step must be performed: the mapping of primary unknowns and internal variables between the old and new meshes. 

In essence, if one aims at continuing the analysis from the current state, instead of restarting it from scratch every time the model is regenerated, it is necessary to apply some transfer algorithms so as not to lose the nodal unknowns, and the history variables 

Fracture propagation approach 

125 

at the Gauss points. 

The mapping of the primary unknowns, i.e. nodal displacements and nodal pore pressure, is achieved by using the shape function projections. To do so, we first place each new node inside an old element, by means of another grid-based search, and then we interpolate the displacements of the old nodes at the position of the new one (Figure 4.16a). 

On the other hand, the mapping of the internal state variables, i.e. _r_ and _ϱ_ respectively governing the damage evolution of the standard elements and the interface elements, is done through a weighted spatial averaging similar to the one used in the non-local damage model (3.48). The difference is that, in this case, the source points are the integration points of the old mesh, and the receiver points are the integration points of the new mesh (Figure 4.16b). Like before, an efficient grid-based search is performed in order to determine the Gauss points of the old mesh that fall inside the interaction radius _Rin_ = _lp_ of each new Gauss point. 

**==> picture [64 x 57] intentionally omitted <==**

(a) Nodal variables mapping. 

**==> picture [112 x 112] intentionally omitted <==**

**==> picture [92 x 59] intentionally omitted <==**

**==> picture [176 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(b) Internal state variables mapping.<br>**----- End of picture text -----**<br>


Figure 4.16: Mapping of variables between meshes. 

Modelling Discontinuities in Porous Media 

126 

## **4.4 Examples** 

This section includes three test cases devoted to assess the performance of the presented fracture propagation strategy and the behaviour of the implemented interface elements. The first example is meant to validate the fracture propagation approach against an analytical solution, whereas the second and third cases are used to test the proposed technique under complex conditions and evaluate its strengths and limitations. 

The three examples model fluid-driven fracture propagation problems in nearly undrainedincompressible porous media. FIC-stabilized elements of equal order interpolation for the displacements and the pore pressure have been used in that regard. 

In all cases, the material of the porous domain is considered to undergo isotropic degradation by means of the same non-local damage model as the one used in Section 3.4.2. In essence, the equivalent strain form is defined as the von Mises model (3.31), and the damage follows the exponential law presented in (3.40). 

The two first problems are analysed in a 2D framework under plane strain conditions, and the last test is solved with a 3D model. The porous medium is considered to have isotropic permeability and the effect of gravity is neglected. 

## **4.4.1 Fluid-driven fracture propagation test** 

This example consists on a semi-cylindrical portion of soil of 80 _m_ radius, laying on a rectangle 10 _m_ high that has a notch of 5 _×_ 50 _cm_ at the center. Water is injected at a constant rate of _qin_ = 0 _._ 2 _m/s_ through an incipient crack tip of 0 _._ 5 _mm_ thick and 125 _mm_ long, which results in a volumetric flow of 1 _·_ 10 _[−]_[4] _m_[3] _/s_ . The geometry and boundary conditions of the problem are shown in Figure 4.17. 

The material of the porous domain is defined as a non-local damage model combining the von Mises equivalent strain form of Eq (3.31) with the exponential damage law in (3.40). The material parameters are summarized in Table 4.1. 

Examples 

127 

**==> picture [90 x 97] intentionally omitted <==**

**==> picture [32 x 34] intentionally omitted <==**

Figure 4.17: Fluid-driven fracture propagation test. Geometry and boundary conditions. Dimensions in _m_ . 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|1_._596_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._33|
|Solid density (_ρs_)|2000 _Kg/m_3|
|Fluid density (_ρf_)|1000 _Kg/m_3|
|Porosity (_φ_)|0_._19|
|Solid bulk modulus (_Ks_)|3_._6_·_1010 _N/m_2|
|Fluid bulk modulus (_Kf_)|3_·_109 _N/m_2|
|Intrinsic permeability (_k_)|6_·_10_−_15 _m_2|
|Dynamic viscosity (_µ_)|0_._001 _N/m_2 _· s_|
|Damage threshold (_ry_)|1_._5_·_10_−_4|
|Compressive to tensile strength ratio (_κ_)|25|
|Residual strength (_R_)|0_._95|
|Softening slope (_S_)|1_._05_·_104|
|Characteristic length (_lc_)|0_._1 _m_|



Table 4.1: Fluid-driven fracture propagation test. Material properties of the porous domain. 

On the other hand, the predefined crack tip is represented by quasi-zero-thickness in- 

Modelling Discontinuities in Porous Media 

128 

terface elements with the properties in Table 4.2. 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|1_._596_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._33|
|Solid density (_ρs_)|2000 _Kg/m_3|
|Fluid density (_ρf_)|1000 _Kg/m_3|
|Porosity (_φ_)|0_._19|
|Solid bulk modulus (_Ks_)|3_._6_·_1010 _N/m_2|
|Fluid bulk modulus (_Kf_)|3_·_109 _N/m_2|
|Transversal permeability (_kn_)|6_·_10_−_15 _m_2|
|Dynamic viscosity (_µ_)|0_._001 _N/m_2 _· s_|
|Damage threshold (_ϱy_)|1_·_10_−_4|
|Minimum joint width (_tmin_)|1_·_10_−_4 _m_|
|Critical displacement (_δc_)|0_._01 _m_|
|Yield stress (_σy_)|7_·_106 _N/m_2|
|Friction coefcient (_µF_)|0_._4|
|Propagation length (_lp_)|0_._06 _m_|
|Propagation damage (d_p_)|0_._6|



Table 4.2: Fluid-driven fracture propagation test. Material properties of the crack. 

As commented above, this problem is approached in a 2D configuration under planestrain assumption. The porous domain is solved with the FIC-stabilized formulation, using 3-noded triangular elements with linear interpolation for the displacements and the pore pressure. The interface elements of the crack tip are 4-noded quadrilateral elements. 

The main objective of this test is to validate the implemented propagation technique against an analytical solution obtained by Spence and Sharp [134] and replicated by Khoei _et al._ [78]. The fluid pressure along the crack mouth, the crack length and the crack width have been monitored through 10 seconds of simulation. In order to assess the robustness of the numerical solution, a convergence analysis was performed in terms of the mesh size and the time step. 

The influence of the mesh size was analysed by fixing the time step at ∆ _t_ = 0 _._ 02 _s_ and 

**==> picture [458 x 432] intentionally omitted <==**

**----- Start of picture text -----**<br>
Examples 129<br>7as < | |D : RY<br>: ee |<br><eek EES He a x<br>ee : aDe Rt _<br>se rane<br>. co<br>a Vv Se [_] SERB Wa<br>fee<br>:eal x“A 7 in\a [LU][_] ao8 a<br>Ooate a C7 oT<br>- . |.LO [a.][.] -ieee [_] KeKea [we] reBOavisBE,OO a [VW][a] =oS8 [oS][ae][nae][ksh] WX|e [he] —-a = [|][a][Ly] ean<br>»Yara [oS] Otes ae<br>Wi<br>* a -Z\OO: [ae][i?][S] [e] [i] sue vAVAAg [RR]<br>Tera ne yz oo a x ee [oe][ae][ae] LXaS ae [ORY] a [a] ae [Se] SER 7 [r][ee] ce84 KXVeVan_we o [oe][ae][p] OOe [a] | v)/\e/\ Or)L\VW [u] :meee:VAN;iy”\ [A]  ZA:Vs [I] [i] [7][A][A][J] —aN yz a [a] = x =. _ _7 [7] - - -:;:RaKAA)Ze- [e][2] [=][a][=] i = 2eau= [oe] _ _ [os] oe S _:-<br>(a) l [e] = 4 cm : 15330 elements. (b) l [e] = 2 cm : 44708 elements.<br>|ee;: BKbsRK [>] Z| [|]<br>RS<br>~_ [|]<br>- ag<br>. vas(Ue. |@<br>_ TeWe| Bs ue [4]<br>: [a]<br>ol<br>/ i. [ll]<br>_ ~ [ay)] [e] [t][ ra][ aK]<br>8<br>y<br>Se<br>z x<br>| = =—<br>(c) l [e] = 1 cm : 153582 elements.<br>**----- End of picture text -----**<br>


Figure 4.18: Fluid-driven fracture propagation test. Initial meshes. 

solving the problem for three different spatial discretizations (Figure 4.18). The first model was obtained from a characteristic element size of _l[e]_ = 4 _cm_ , with a resultant mesh of 15330 elements and 7846 nodes. The second mesh, with 44708 elements and 22580 nodes, resulted from an element size of _l[e]_ = 2 _cm_ . The third model was obtained after defining a characteristic element size of _l[e]_ = 1 _cm_ and lead to a mesh of 153582 elements and 77095 nodes. 

In order to study the effect of the time step on the solution, we fixed the mesh size at _l[e]_ = 2 _cm_ and run the case for three different time steps: ∆ _t_ = 0 _._ 04 _s_ , ∆ _t_ = 0 _._ 02 _s_ and ∆ _t_ = 0 _._ 01 _s_ . 

Modelling Discontinuities in Porous Media 

130 

**==> picture [268 x 533] intentionally omitted <==**

**----- Start of picture text -----**<br>
3.5<br>3<br>le=4cm<br>le=2cm<br>2.5 le=1cm<br>Analytical<br>2<br>1.5<br>1<br>0.5<br>0<br>0 2 4 6 8 10<br>Time (s)<br>(a) Pressure for different mesh sizes with ∆ t  = 0 . 02 s .<br>3.5<br>3<br>Δt=0.04s<br>2.5<br>Δt=0.02s<br>2 Δt=0.01s<br>Analytical<br>1.5<br>1<br>0.5<br>0<br>0 2 4 6 8 10<br>Time (s)<br>Pressure (MPa)<br>Pressure (MPa)<br>**----- End of picture text -----**<br>


(b) Pressure for different time steps with _l[e]_ = 2 _cm_ . 

Figure 4.19: Fluid-driven fracture propagation test. Time evolution of the pressure. 

Examples 

131 

**==> picture [268 x 252] intentionally omitted <==**

**----- Start of picture text -----**<br>
4.5<br>4<br>3.5<br>3<br>2.5<br>le=4cm<br>2 le=2cm<br>le=1cm<br>1.5<br>Analytical<br>1<br>0.5<br>0<br>0 2 4 6 8 10<br>Time (s)<br>Crack length (m)<br>**----- End of picture text -----**<br>


(a) Crack length for different mesh sizes with ∆ _t_ = 0 _._ 02 

_s_ . 

**==> picture [268 x 252] intentionally omitted <==**

**----- Start of picture text -----**<br>
4.5<br>4<br>3.5<br>3<br>2.5<br>Δt=0.04s<br>2<br>Δt=0.02s<br>1.5<br>Δt=0.01s<br>1 Analytical<br>0.5<br>0<br>0 2 4 6 8 10<br>Time (s)<br>Crack length (m)<br>**----- End of picture text -----**<br>


(b) Crack length for different time steps with _l[e]_ = 2 _cm_ . 

Figure 4.20: Fluid-driven fracture propagation test. Time evolution of the crack length. 

Modelling Discontinuities in Porous Media 

132 

**==> picture [268 x 252] intentionally omitted <==**

**----- Start of picture text -----**<br>
3.5<br>3<br>2.5<br>2 le=4cm<br>le=2cm<br>1.5 le=1cm<br>Analytical<br>1<br>0.5<br>0<br>0 2 4 6 8 10<br>Time (s)<br>Crack width (mm)<br>**----- End of picture text -----**<br>


(a) Crack width for different mesh sizes with ∆ _t_ = 0 _._ 02 

_s_ . 

**==> picture [268 x 252] intentionally omitted <==**

**----- Start of picture text -----**<br>
3.5<br>3<br>2.5<br>Δt=0.04s<br>2<br>Δt=0.02s<br>Δt=0.01s<br>1.5<br>Analytical<br>1<br>0.5<br>0<br>0 2 4 6 8 10<br>Time (s)<br>Crack width (mm)<br>**----- End of picture text -----**<br>


(b) Crack width for different time steps with _l[e]_ = 2 _cm_ . 

Figure 4.21: Fluid-driven fracture propagation test. Time evolution of the crack width. 

Examples 

133 

Figures 4.19a, 4.20a and 4.21a show the pressure, crack length and crack width for different mesh size and ∆ _t_ = 0 _._ 02 _s_ . We highlight that the combination of an integral-type non-local damage model with the insertion of interface elements has no mesh dependence. 

The results displayed in Figures 4.19b, 4.20b and 4.21b evidence that the time step is an important parameter in this kind of problems. Indeed, despite integrating the time variable with an unconditionally stable Newmark scheme, a large time step together with the material non-linearity of the problem lead to inaccurate results. Such effect is particularly remarkable in Figure 4.21b, mainly because of the greater sensitivity of the crack width variable (in the range of millimetres) with respect to the crack length and the pore pressure variables (in the range of meters and mega pascals, respectively). In any case, using a small enough time step ensures obtaining consistent results. In this case good results have been obtained for ∆ _t ≤_ 0 _._ 02 _s_ . 

Regarding the validation of the present test, a general good agreement is observed between the numerical and the analytical curves. Looking at Figures 4.19a and 4.19b, it is evident that the pressure at the mouth falls faster in the beginning of the numerical solutions. Such pressure drop is related to a sudden lost of integrity of the soil, and thus an accurate calibration of the parameters could smooth this effect. However, after a few seconds all curves converge to the same value of 0 _._ 5 _MPa_ , which remains virtually constant due to the near undrained condition of the soil. 

In Figures 4.20a and 4.20b, one observes that the computed crack grows at a slightly slower pace than the expected solution. This fact is clearly linked to the fast drop in the pressure already commented. 

Finally, it is noticeable the similarity between the theoretical crack width and the numerical one. Except for the solution with ∆ _t_ = 0 _._ 04 _s_ , the slope of the curve and the maximum value reached in the rest of the cases are virtually identical to the analytical result. 

## **4.4.2 Crack tracking test** 

This second test models a rectangular saturated porous domain of 6 _×_ 5 _m_ . The same notch and crack tip of the previous example have been horizontally placed at the left 

Modelling Discontinuities in Porous Media 

134 

side of the model, with equivalent conditions at the inlet. In this case, an incursion of stronger soil is defined over 60 _cm_ apart from the initial crack tip, as represented in Figure 4.22. 

**==> picture [164 x 121] intentionally omitted <==**

**==> picture [37 x 41] intentionally omitted <==**

Figure 4.22: Crack tracking test. Geometry and boundary conditions. Dimensions in _m_ . 

The porous media and the crack tip are characterized with the same properties of the previous example (see Tables 4.1 and 4.2). The only difference is in the material incursion in front of the crack tip, the properties of which are listed in Table 4.3. 

Examples 

135 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|4_._0_·_1010 _N/m_2|
|Poisson’s ratio (_ν_)|0_._33|
|Solid density (_ρs_)|2000 _Kg/m_3|
|Fluid density (_ρf_)|1000 _Kg/m_3|
|Porosity (_φ_)|0_._19|
|Solid bulk modulus (_Ks_)|3_._6_·_1010 _N/m_2|
|Fluid bulk modulus (_Kf_)|3_·_109 _N/m_2|
|Intrinsic permeability (_k_)|6_·_10_−_15 _m_2|
|Dynamic viscosity (_µ_)|0_._001 _N/m_2 _· s_|
|Damage threshold (_ry_)|5_._5_·_10_−_3|
|Compressive to tensile strength ratio (_κ_)|25|
|Residual strength (_R_)|0_._95|
|Softening slope (_S_)|1_._05_·_104|
|Characteristic length (_lc_)|0_._1 _m_|



Table 4.3: Crack tracking test. Material properties of the incursion. 

Again, this case is solved as a plane-strain 2D problem, and FIC-stabilized 3-noded triangles are used to circumvent the violation of Babuska-Brezzi conditions due to the impermeability of the porous medium. Also, the characteristic element size chosen is _l[e]_ = 2 _cm_ and ∆ _t_ = 0 _._ 02 _s_ . 

The main purpose of the present test is to verify the “crack-tracking” capabilities of the implemented approach against anisotropic porous domains. In order to do so, three different “obstacles” have been placed in front of the propagating crack (Figure 4.23). The first two incursions (Figures 4.23a and 4.23b) are meant to simulate branching situations for different bifurcation angles, whereas the third obstacle introduces an irregular preferential path for the advancing fracture (Figure 4.23c). 

**==> picture [458 x 412] intentionally omitted <==**

**----- Start of picture text -----**<br>
136 Modelling Discontinuities in Porous Media<br>y y<br>z x z x<br>Pp<br>(a) Obstacle 1: 45-90-45 isosceles triangle. (b) Obstacle 2: 75-30-75 isosceles triangle.<br>y<br>z x<br>**----- End of picture text -----**<br>


**==> picture [186 x 12] intentionally omitted <==**

**----- Start of picture text -----**<br>
(c) Obstacle 3: moon shape incursions.<br>**----- End of picture text -----**<br>


Figure 4.23: Crack tracking test. Material incursions in front of the crack tip. 

The three crack paths, after 2 seconds of simulation, are displayed in Figure 4.24. The fluid pressure and the damage of the porous domain are also illustrated to give a deeper insight into the current analysis. 

Examples 

137 

**==> picture [364 x 174] intentionally omitted <==**

**----- Start of picture text -----**<br>
WATER_PRESSURE DAMAGE_VARIABLE<br>6.6419e+05 1<br>4.9814e+05 0.75<br>3.3209e+05 0.5<br>1.6605e+05 0.25<br>0 0<br>y y<br>z x z x<br>(a) Pressure with obstacle 1. (b) Damage with obstacle 1.<br>**----- End of picture text -----**<br>


**==> picture [367 x 181] intentionally omitted <==**

**----- Start of picture text -----**<br>
e WATER_PRESSURE e DAMAGE_VARIABLE<br>6.7931e+05 1<br>5.0948e+05 0.75<br>3.3966e+05 0.5<br>1.6983e+05 0.25<br>0 0<br>y y<br>z x z x<br>~ ae<br>(c) Pressure with obstacle 2. (d) Damage with obstacle 2.<br>**----- End of picture text -----**<br>


**==> picture [368 x 200] intentionally omitted <==**

**----- Start of picture text -----**<br>
(c) Pressure with obstacle 2. (d) Damage with obstacle 2.<br>a WATER_PRESSURE a 2ee e a DAMAGE_VARIABLE<br>7.3184e+05 1<br>5.4888e+05 0.75<br>3.6592e+05 0.5<br>1.8296e+05 0.25<br>0 0<br>IN \ \ ( { [ = ————<br>y y<br>z x z x<br>(e) Pressure with obstacle 3. (f) Damage with obstacle 3.<br>**----- End of picture text -----**<br>


Figure 4.24: Crack tracking test. Pore pressure and damage at _t_ = 2 _s_ . 

Modelling Discontinuities in Porous Media 

138 

Focusing first on Figures 4.24a, 4.24c and 4.24e, it is evident that the proposed technique captures the anisotropy introduced in the model, and the inserted interface elements properly represent the enhanced permeability of the porous medium. Moreover, the smooth contour lines of the pressure field show that the FIC-stabilized 3-noded triangles have an excellent behaviour near the undrained-incompressible conditions. It is interesting to note that the maximum water pressure is different for each case, with the highest value in Figure 4.24e and the lowest in Figure 4.24a. Although the difference is subtle, one may infer that the pressure dissipates faster in a branching case with a wider angle of bifurcation. 

**==> picture [382 x 71] intentionally omitted <==**

**----- Start of picture text -----**<br>
AVAVAVAVAVAVAVIN AV W AY eg KAVA AVA NY AVA AYA S VAVVay VAVAVAVAVANINVAVA9S VAVAVAVAVAVAVAA<br>VAVAVAVAVAVAV y  GINNee A VA NY AVNVAV ACNYAYAYAVAVAVAV AVAV NV NV V A VAVN A VAVAVVA A N CASTIAVAVAVAV y AVA VAVAVAV AVA SOS 5 NOAA AYA<br>SVAVAVAVAVAV AY NU NVAVAVAVAVAN IVAN AVAVAVAVAVAVAVAVAV IVAWAVAVAVAVAVAVAVAVAVAVAVAABETS COO<br>TAEES,ATAVAV z x  AY AYA AV AN AVAVAVAVAVAN AVA AN YAY AV AV AVA AN TAYAYAVAVANAVAVAYAYAVAVAVAVAVAYAVAUAO AVIGERAIS z x Pa 4s7AVAVANAVAYAV4\AVACAYAVAOO<br>(a) Mesh at initial stage. (b) Mesh at intermediate stage.<br>**----- End of picture text -----**<br>


**==> picture [179 x 73] intentionally omitted <==**

**----- Start of picture text -----**<br>
REREKESREEPERSE SRE SOAK DERERE RIS<br>REPS y  SRERER ERR ASR ERK IOL<br>GEEK ERPS SEE RE EESK ISX<br>PPKKEKE z 5 x E e cheeREREKPKER eeekheeRkee BREEPEKEREEK REEERE E KERRee]RPEKKSA<br>SRERESEERERSEREREREERDACERERD<br>(c) Mesh at advanced stage.<br>**----- End of picture text -----**<br>


Figure 4.25: Crack tracking test. Fracture passing between moon shape incursions. 

Examples 

139 

In Figures 4.24b, 4.24d and 4.24f, the non-locality of the damage model is clearly visible in the width of the damage mark around the crack. Indeed, with a characteristic length of _lc_ = 0 _._ 1 _m_ , the crack tip is surrounded by a damaged area with a diameter of over 0 _._ 2 _m_ . Such diffusive damage pattern is what leads the crack through the right path, avoiding mesh-induced directional bias. Thereby, in the first two cases (Figures 4.24b and 4.24d) the fracture eventually runs into the stronger incursion and, given the low degradation of the material in front of the tip, the crack bifurcates with an angle that depends on the shape of the incursion. For the moon shape obstacles (Figure 4.24f), the undamaged areas on both sides of the tip let the crack no choice but to pass between the two obstacles. A detailed view of the mesh is shown in Figure 4.25. 

## **4.4.3 Parallel fracture propagation test** 

This last example consists on a 4 _×_ 4 _×_ 6 _m_ block of soil with two parallel crack tips separated by a distance _s_ . Water is injected at a constant rate of _qin_ = 6 _m/s_ through each one of the cracks of 0 _._ 5 _mm_ thick and 125 _mm_ long, giving a total volumetric flow of around 1 _._ 5 _·_ 10 _[−]_[6] _m_[3] _/s_ . A scheme of the geometry and the boundary conditions is represented in Figure 4.26. 

**==> picture [91 x 65] intentionally omitted <==**

**==> picture [94 x 120] intentionally omitted <==**

Figure 4.26: Parallel fracture propagation test. Geometry and boundary conditions. Dimensions in _m_ . 

Modelling Discontinuities in Porous Media 

140 

In this case a softer material has been defined in order to reduce the maximum pore pressure and minimize the computational cost. The material properties of the porous domain and the two incipient cracks are given in Tables 4.4 and 4.5, respectively. 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|2_·_109 _N/m_2|
|Poisson’s ratio (_ν_)|0_._1|
|Solid density (_ρs_)|2000 _Kg/m_3|
|Fluid density (_ρf_)|1000 _Kg/m_3|
|Porosity (_φ_)|0_._3|
|Solid bulk modulus (_Ks_)|1_._5_·_1017 _N/m_2|
|Fluid bulk modulus (_Kf_)|3_·_1014 _N/m_2|
|Intrinsic permeability (_k_)|1_·_10_−_14 _m_2|
|Dynamic viscosity (_µ_)|0_._001 _N/m_2 _· s_|
|Damage threshold (_ry_)|7_._5_·_10_−_6|
|Compressive to tensile strength ratio (_κ_)|10|
|Residual strength (_R_)|0_._4|
|Softening slope (_S_)|1_._5_·_104|
|Characteristic length (_lc_)|0_._08 _m_|



Table 4.4: Parallel fracture propagation test. Material properties of the porous domain. 

Since the problem is approached in a 3D framework, the porous domain is represented by stabilized 4-noded tetrahedra with linear interpolation for the displacements and the pore pressure. The two crack tips are meshed with six-node wedge interface elements. In this case, the characteristic element size and the time step are larger than in previous examples with _l[e]_ = 5 _cm_ and ∆ _t_ = 0 _._ 05 _s_ . 

Examples 

141 

|**Property**|**Value**|
|---|---|
|Young’s modulus (_E_)|2_·_109 _N/m_2|
|Poisson’s ratio (_ν_)|0_._1|
|Solid density (_ρs_)|2000 _Kg/m_3|
|Fluid density (_ρf_)|1000 _Kg/m_3|
|Porosity (_φ_)|0_._3|
|Solid bulk modulus (_Ks_)|1_._5_·_1017 _N/m_2|
|Fluid bulk modulus (_Kf_)|3_·_1014 _N/m_2|
|Transversal permeability (_kn_)|1_·_10_−_14 _m_2|
|Dynamic viscosity (_µ_)|0_._001 _N/m_2 _· s_|
|Damage threshold (_ϱy_)|1_·_10_−_4|
|Minimum joint width (_tmin_)|1_·_10_−_4 _m_|
|Critical displacement (_δc_)|0_._01 _m_|
|Yield stress (_σy_)|3_._5_·_106 _N/m_2|
|Friction coefcient (_µF_)|0_._4|
|Propagation length (_lp_)|0_._05 _m_|
|Propagation damage (d_p_)|0_._4|



Table 4.5: Parallel fracture propagation test. Material properties of the cracks. 

Basically, the purpose of this problem is to analyse the influence of the distance between neighbouring cracks in a parallel fracture propagation case. Thereby, in order to give a more realistic insight and account for the Poisson’s effect in the transversal deformation of the soil, a 3D configuration has been chosen for this last test. 

As reported in [39], to ensure the efficiency of a hydraulic fracturing process, the spacing between parallel cracks should be of the order of 0 _._ 1 _m_ . Here three different spacings have been studied: _s_ = 0 _._ 075 _m_ , _s_ = 0 _._ 15 _m_ and _s_ = 0 _._ 3 _m_ , and the crack length in each case has been measured after 1 _._ 5 seconds of simulation. Table 4.6 summarizes the results along with the obtained crack growth velocity. 

The values in Table 4.6 show that the closer the initial cracks are defined, the faster they propagate, as expected. However, the relation between the spacing and the crack growth velocity seems to be not linear, which also makes sense if one takes into account the material non-linearity of this problem. 

Modelling Discontinuities in Porous Media 

142 

|Spacing|Crack Length at _t_= 1_._5 _s_<br>Crack growth velocity|
|---|---|
|0_._075 _m_<br>0_._4 _m_<br>0_._18 _m/s_<br>0_._15 _m_<br>0_._37 _m_<br>0_._16 _m/s_<br>0_._3 _m_<br>0_._28 _m_<br>0_._1 _m/s_||



Table 4.6: Parallel fracture propagation test. Crack length and crack growth velocity. 

In order to properly appreciate the previous results, the contour lines of the pressure field and the damage variable are represented on two orthogonal planes cutting both cracks (Figure 4.27). 

Looking first at the pressure field in Figures 4.27a, 4.27c and 4.27e, it is quite clear that a closer spacing between cracks makes the pore pressure grow, and consequently the propagation velocity also increases. Also, to understand the interaction between two parallel cracks, it is useful to observe Figures 4.27b, 4.27d and 4.27f. Indeed, the nonlocal damage model generates a diffusive damage mark around the cracks that extends up to the radius of influence _lc_ . In Figures 4.27b and 4.27d, in which the spacing is smaller than twice the characteristic length of the material, the damage patterns of the cracks seem to overlap each other. In Figure 4.27f, with _s >_ 2 _lc_ , the two damage marks are almost independent and hence the velocity of propagation noticeably decreases. 

This example shows that the characteristic length of the non-local damage model _lc_ is not only a mathematical parameter used to regularize the strain localization problem, but it actually has a physical meaning. 

Examples 

143 

**==> picture [422 x 364] intentionally omitted <==**

**----- Start of picture text -----**<br>
WATER_PRESSURE DAMAGE_VARIABLE<br>4.0044e+05 1<br>3.0033e+05 0.75<br>z 2.0022e+05 z 0.5<br>y x 1.0011e+05 y x 0.25<br>0 0<br>Le lL [ t. | |<br>(a) Pressure with s  = 0 . 075 m . (b) Damage with s  = 0 . 075 m .<br>WATER_PRESSURE DAMAGE_VARIABLE<br>4.0209e+05 1<br>3.0157e+05 0.75<br>z 2.0104e+05 z 0.5<br>y x 1.0052e+05 y x 0.25<br>0 0<br>(c) Pressure with s  = 0 . 15 m . (d) Damage with s  = 0 . 15 m .<br>WATER_PRESSURE DAMAGE_VARIABLE<br>3.7788e+05 1<br>2.8341e+05 0.75<br>z 1.8894e+05 z 0.5<br>y x 94471 y x 0.25<br>0 0<br>(e) Pressure with s  = 0 . 3 m . (f) Damage with with s  = 0 . 3 m .<br>**----- End of picture text -----**<br>


Figure 4.27: Parallel fracture propagation test. Pore pressure and damage at _t_ = 1 _._ 5 _s_ . 

5 Chapter 

## **Conclusions and Future Lines of Research** 

## **5.1 Conclusions and contributions** 

The main objective of this thesis was the derivation of a robust Finite Element formulation for the solution of solid-pore fluid coupled problems with multi-fracture propagation. The following conclusions can be highlighted. 

A displacement-pore pressure FEM formulation for solving coupled solid-pore fluid interaction problems has been successfully implemented for 2D and 3D analysis, and stabilized by means of the Finite Increment Calculus method (FIC). In this formulation, a unique mesh represents a porous medium composed by two phases: a solid skeleton and a fluid phase. The interaction between both components is governed by means of the coupling of two equations: the balance of momentum for the mixture solid-fluid and the mass balance for the pore fluid. 

The derived FIC-FEM formulation represents the first application of the second order FIC mass balance equation in space to the stabilization of a poromechanics problem. In the presented approach, the balance of momentum equation is unaltered and only the continuity equation is modified with the stabilization terms. 

Comparing this work with other classical stabilization methods, we notice that the FIC 

Conclusions and Future Lines of Research 

146 

stabilization allows solving the problem in a fully coupled manner without relying on time stepping algorithms. It also avoids the calibration of the stabilization parameter, which simply depends on the characteristic element size. 

The presented FIC-based approach does not require additional basis functions or element level condensation, but relies on higher-order derivatives to obtain the terms that ensure the consistency of the residual-based procedure. 

As presented in the Examples’ section in Chapter 2, the FIC-FEM formulation has been used to solve an elastic saturated soil column subjected to surface loading in a 2D problem under plane strain conditions. It has proven to avoid arbitrary oscillations along the column, and it has shown consistent results for both transient and cyclic loading cases, despite modifying the material compressibility and permeability. 

Additionally, the FIC-FEM element has been tested in a 3D problem of an elastic saturated soil foundation. The effect of the spatial discretization on the solution has been addressed in this case. This problem has shown that the residual-based character of the FIC-stabilization favours that the numerical solution converges to the expected solution as the mesh is refined. Moreover, when relaxing the undrained-incompressible conditions, the FIC stabilization does not alter negatively the response obtained in the original non-stabilized formulation, and yields an accurate solution. 

Since the proposed fracture propagation technique relies on continuum damage mechanics for the estimation of the crack path, Chapter 3 focuses on the seek of a robust damage model for the analysis of the degradation and failure of quasi-brittle materials. 

We saw that a local damage model, although partially regularized with an energy based adjustment, shows mesh sensitiveness and presents severe mesh-induced directional bias, which makes it virtually unusable for tracking the crack path. 

However, we must also notice that for certain cases and if carefully calibrated, a partially regularized local damage model can become a useful tool to obtain a first estimation of the damage pattern with a very low computational cost. 

Of course, in solid-pore fluid coupled problems, where cracks represent preferential paths for the fluid flow and have a critical influence on the global permeability of the porous medium, it is essential to work with a fully regularized damage model. In this regard, the implemented integral-type non-local damage model has shown a robust behaviour 

Conclusions and contributions 

147 

for different spatial discretizations, avoiding pathological mesh sensitivity and meshinduced directional bias. 

An important aspect to note about the non-local approach is that the neighbour search and the averaging performed to account for the non-local interaction, considerably increase the computational cost of the solution, as compared to the classical local approach, and also modify the traditional assembly process of the global tangent stiffness matrix. 

Both damage models have been tested in a 2D three-point bending test under plane stress assumption. Excluding the mesh-dependent local model, the non-local damage model based on the Mazars definition of the equivalent strain has proved to properly capture the peak load of the Force-Deflection curve for all meshes. Nonetheless, the postpeak branch of the numerical solution is more brittle than in the reference solution, and even shows a certain snap-back behaviour. The use of an advancing technique based on the control of displacements should smooth the behaviour in the post-peak region and improve the convergence of the solution. 

Furthermore, a plane stress four-point shear test has been solved to check that the non-local damage model accurately reproduces the Force-CMSD curve and converges to the expected solution as the mesh is refined. This last test also showed that the non-local form of the modified von Mises definition of the equivalent strain, along with the exponential damage evolution law actually represent an effective combination given its remarkable robustness and convergence, in spite of its simplicity, specially when compared to the tested alternatives. 

Chapter 4 presents the quasi-zero-thickness interface elements which, combined with the coupled formulation of Chapter 2 and the non-local damage model of Chapter 3, are used to develop a methodology for the 2D and 3D analysis of fracture propagation problems in saturated porous media. 

The quasi-zero-thickness interface elements are governed by the same two equations presented in Chapter 2, i.e. the balance of momentum for the mixture solid-fluid and the mass balance for the pore fluid, but they are conveniently adapted to account for the special behaviour of the joints. This fact is particularly useful as it allows combining the standard elements with the interface elements in the same model and solve a unique coupled system of equations. 

Conclusions and Future Lines of Research 

148 

The formulation of the interface elements, distinguishing its upper and lower faces, has shown to adequately introduce jumps in the displacement field, and capture the enhanced permeability of the porous medium. 

Three examples involving crack propagation have been analysed with different objectives. In the 2D fluid-driven fracture propagation test, the combination of the non-local damage model with the automatic insertion of interface elements has shown to be robust and effective. It has been validated against analytical values of the pressure at the crack mouth, the crack length, and the crack width with promising results. 

Furthermore, as seen in Chapter 3, the integral-type non-local damage model fully regularized the boundary value problem and thus the curves of the monitored variables during crack propagation were virtually unaffected by the changes in the spatial discretization. 

On the other hand, the time step appeared as an important parameter in problems involving material non-linearities and fluid diffusion. A small enough value must be used to ensure accurate results. 

In the crack tracking test, we checked that the search of damaged points around the crack tip, combined with the non-local damage model, is a robust tool for predicting the crack path in anisotropic media, including the possibility of branching when the proper conditions are met. 

The FIC-stabilized elements of equal order interpolation for the displacements and the pore pressure have also been proved to prevented the blocking of the pressure field during crack propagation. 

Finally, the proposed technique has been positively tested in a 3D configuration for the analysis of the complex interaction between parallel cracks. As expected, the crack growth velocity is inversely proportional to the spacing between fractures. 

Moreover, this example provided useful information to understand that the characteristic length of the non-local damage model _lc_ is not only a numerical parameter, but it actually has a physical meaning, as it plays a crucial role in the interaction between parallel cracks. 

We note that the memory usage in 3D fracture problems is remarkable due to the storage of the neighbours for the non-local damage model and thus the computational cost can 

Lines of future work 

149 

be intense. The optimization of the code is therefore essential before the characteristic element size can be reduced. 

In conclusion, the main contributions and innovative points of the present thesis can be summarized in the following list: 

- Derivation of a stable FIC-FEM formulation for the solution of coupled solid-pore fluid interaction problems, which is based on the second order FIC mass balance equation in space. 

- Evolution and assessment of a modular 2D and 3D integral-type non-local damage model, relying on a grid-based search of Gauss points neighbours. 

- Generalization of the interface elements presented in [78, 128], in the form of 2D and 3D quasi-zero-thickness interface elements, which are valid for the solution of the fully coupled solid-pore fluid problem. 

- Development of a crack tracking technique based on the non-local evaluation of the region of influence around the crack tip. 

- Implementation of a innovative utility for the automatic insertion of interface elements for fracture propagation and bifurcation using GiD. 

## **5.2 Lines of future work** 

This thesis opens new perspectives to feasible and attractive future works. The most relevant are outlined below. 

The porous domain has been considered as a two-phase medium composed by a solid skeleton and a fluid phase. In reality, though, this is not usually the case and the porous space can be partially filled by a fluid and partially filled by a gas. In this regard, the generalization to three-phase problems, such as those encountered in unsaturated soils or in oil-gas-soil interaction, are possible extensions of the numerical approach presented here. 

Conclusions and Future Lines of Research 

150 

Fredlund and Morgenstern considered an unsaturated soil as a four-phase system, with the air-water interface being the fourth phase, and proposed suitable stress state variables based on the equilibrium of forces for each phase [57]. A simple extension of two phase formulation to semi-saturated problems was proposed by Zienkiewicz _et al._ [158] based on the assumption that the air or gas present in the pores remains at atmospheric pressure. Alonso _et al._ [4] formulated a constitutive model within the framework of hardening plasticity for describing the stress-strain behaviour of partially saturated soils which are slightly expansive. 

Further development has been done by Gawin and Schrefler [59], and Khoei _et al._ [77] for solving geotechnical problems of unsaturated soils. A more recent application of a solid-fluid-gas formulation is found in the underground storage of carbon dioxide [73]. 

Another possible expansion of the current formulation is the inclusion of the temperature field as a variable of the problem. A particularly interesting application is in the modelling of geothermal energy production where the solid-pore fluid-thermal formulation can be combined with interface elements in a displacement-pore pressure-temperature system of equations [159]. 

**==> picture [144 x 134] intentionally omitted <==**

**----- Start of picture text -----**<br>
Generating<br>station<br>a<br>Bs<br>[Leepf =] BSy—\SE<br>“ .' we ‘\ vou<br>\ \ ?<br>ve _ at ‘a1\ vos<br>Cold water Steam and<br>pumped down hot water<br>**----- End of picture text -----**<br>


Figure 5.1: Scheme of the geothermal energy production. 

Figure 5.1 schematically represents the process for the generation of geothermal energy: 

Lines of future work 

151 

water is injected into deep formations and, because of the higher temperatures of the underground rocks, steam and hot water is returned to the surface, driving turbines and electricity generators. 

Regarding future enhancements for the current code, we highlight two main topics: the generalization of the presented fracture propagation technique, and the optimization of the code for improving its speed and efficiency. 

In the fracture propagation utility, many advances can be introduced. The most important are summarized below: 

- Automatically generate interface elements from a cloud of damaged points, so as to avoid pre-defining initial crack tips. 

- Introduce the technique presented in [48] to adapt the mesh according to the stress state after the insertion of new interface elements. 

- Implement alternative crack designs for 3D fracture propagation. 

- Add the possible scenario of two different cracks crossing each other. 

In terms of code optimization, there is still room for cleaning the implemented procedures to make them more efficient. However, there is one aspect that should be addressed first: the parallel programming. 

As commented in Chapter 1, the code has been implemented inside the Kratos programming framework based on C++ language. Kratos is prepared for parallel computing using OpenMP (Open Multiprocessing) and MPI (Message Passing Interface) (Figure 5.2). However, certain features of the code, such as the non-local damage model or the crack tracking utility, are currently designed to work only in OpenMP configuration, but not in MPI. In order to develop a computational technology that can compete with other commercial software, specially in large 3D problems, the extension of all the code to MPI can not be a secondary task. 

Conclusions and Future Lines of Research 

152 

Figure 5.2: MPI domain partitioning in 4 sections. 

Finally, this work is also related to alternative lines of research that, although differ from the main topic of the thesis, are supplied with features that derive from it and will be further developed in the near future. 

One of these lines of work is the numerical analysis of concrete dams through the finite element method. 

This thesis was involved in a convenient synergy between the ACOMBO project for the development of numerical tools for the thermo-mechanical analysis of arch dams, and Lorenzo Gracia’s PhD thesis “Unified finite element modelling of the life cycle of concrete dams”. 

The main objective of both the ACOMBO project and Lorenzo Gracia’s thesis is the development of numerical tools that help analyse the behaviour of a concrete dam from the design period to the end of its service life. Such analysis is divided in three different stages: construction of the dam, in which the evolution of the heat of hydration is studied to prevent uncontrolled cracking, operation period, where the solution of a thermomechanical problem considering the temperature variations of the year is used to detect possible anomalies when compared to monitor devices, and extraordinary events, which include non-linear analysis for the study of damage initiation and crack propagation and dynamic analysis for the prevention of catastrophic failure under seismic loading. 

Lines of future work 

153 

The present thesis took also part in that synergy as the main contributor in the development of damage models and interface elements. Although the quasi-zero-thickness interface elements were initially implemented for modelling fractures in saturated porous media, the project required them for the definition of joints in the body of the dam or in the rock formation, and for the modelling of the contact between the dam and the rock foundation. 

The ACOMBO project, which congregates Jesús Granell Engineering, Endesa, Universidad Politécnica de Madrid and CIMNE, gave a practical insight into the derived formulations in the thesis, and provided useful data to understand the applicability of the implemented numerical tools. 

An example of this synergy is displayed in Figure 5.3, where the joints of Baserca’s arch dam open because of the thermal contraction produced by the cold temperatures of winter. 

Figure 5.3: Opening of the joints of Baserca’s arch dam during winter. Deformation is exaggerated 200 times. 

Related to the numerical analysis of dams, two publications have already been generated for the ICOLD Benchmark Workshop 2017. They are listed in Chapter 1. 

Conclusions and Future Lines of Research 

154 

Additionally, the implemented interface elements open another alternative line of work: the analysis of interlaminar cracks in laminated composite materials. 

Essentially, the idea is to take advantage of the natural definition of the interface elements to perform a direct numerical simulation of composite delamination. This means that every layer of the composite is modelled as a solid material, and quasi-zero-thickness interface elements are placed in between them. 

As it is widely known, 3D finite element analysis is the more appropriate tool to accurately model plates and shells of laminated composite material. Combined with cohesive interface elements it can effectively reproduce delamination in composite beams, as it is shown in Figure 5.4. However, for composites with hundreds of piles, such model becomes considerably expensive, specially if we take into account that interface elements introduce duplicated nodes at every layer. 

Figure 5.4: Multi-delamination of a 3-layered beam. 

In order to partially alleviate the computational cost, an innovative strategy for the deactivation and activation of interface elements has already been implemented. 

In essence, all the interface elements between the layers of the laminate would be initially defined as inactive, so that no extra degrees of freedom are solved. As the beam is loaded, the interface elements with a stress value higher than a pre-fixed threshold would be dynamically activated along with its corresponding nodes. We remark that no remeshing is required here to activate the interface elements, because they are all defined at the 

Lines of future work 

155 

beginning of the problem. 

If this technique is combined with efficient MPI parallel computing, solving problems that had been prohibitively expensive in the past could be now feasible. 

A publication regarding the analysis of composite delamination using the mentioned approach is currently in preparation. 

## **Bibliography** 

- [1] GiD the personal pre and post processor. `https://www.gidhome.com/` , July 2017. 

- [2] Kratos multiphysics. `http://www.cimne.com/kratos/` , July 2017. 

- [3] J. Adachi, E. Siebrits, A. Peirce, and J. Desroches. Computer simulation of hydraulic fractures. _International Journal of Rock Mechanics and Mining Sciences_ , 44(5):739–757, 2007. 

- [4] E. E. Alonso, A. Gens, and A. Josa. A constitutive model for partially saturated soils. _Géotechnique_ , 40(3):405–430, 1990. 

- [5] J. Andersson and B. Dverstorp. Conditional simulations of fluid flow in three-dimensional networks of discrete fractures. _Water Resources Research_ , 23(10):1876–1886, 1987. 

- [6] I. Babuška. The finite element method with lagrangian multipliers. _Numerische Mathematik_ , 20(3):179–192, 1973. 

- [7] O. R. Barani, A. R. Khoei, and M. Mofid. Modeling of cohesive crack growth in partially saturated porous media; a study on the permeability of cohesive fracture. _International Journal of Fracture_ , 167(1):15–31, 2011. 

- [8] G. I. Barenblatt. The formation of equilibrium cracks during brittle fracture. general ideas and hypotheses. axially-symmetric cracks. _Journal of Applied Mathematics and Mechanics_ , 23(3):622–636, 1959. 

Bibliography 

158 

- [9] G. I. Barenblatt. The mathematical theory of equilibrium cracks in brittle fracture. _Advances in Applied Mechanics_ , 7(1):55–129, 1962. 

- [10] Z. P. Bazant. Mechanics of distributed cracking. _Applied Mechanics Reviews_ , 39(5):675–705, 1986. 

- [11] Z. P. Bazant, T. B. Belytschko, and T.-P. Chang. Continuum theory for strainsoftening. _Journal of Engineering Mechanics_ , 110(12):1666–1692, 1984. 

- [12] Z. P. Bazant and M. Jirásek. Nonlocal integral formulations of plasticity and damage: survey of progress. _Journal of Engineering Mechanics_ , 128(11):1119– 1149, 2002. 

- [13] Z. P. Bazant and B. H. Oh. Crack band theory for fracture of concrete. _Matériaux et construction_ , 16(3):155–177, 1983. 

- [14] T. Belytschko and T. Black. Elastic crack growth in finite elements with minimal remeshing. _International Journal for Numerical Methods in Engineering_ , 45(5):601–620, 1999. 

- [15] B. Berkowitz, J. Bear, and C. Braester. Continuum models for contaminant transport in fractured porous formations. _Water Resources Research_ , 24(8):1225–1236, 1988. 

- [16] G. Bfer. An isoparametric joint/interface element for finite element analysis. _International journal for numerical methods in engineering_ , 21(4):585–600, 1985. 

- [17] C. B. Biezeno and R. Grammel. _Technische Dynamik_ . Springer-Verlag, 1939. 

- [18] M. A. Biot. General theory of three-dimensional consolidation. _Journal of applied physics_ , 12(2):155–164, 1941. 

- [19] M. A. Biot. Theory of elasticity and consolidation for a porous anisotropic solid. _Journal of applied physics_ , 26(2):182–185, 1955. 

- [20] M. A. Biot. Theory of propagation of elastic waves in a fluid-saturated porous solid. I: Low-frequency range. II: Higher frequency range. _the Journal of the Acoustical Society of America_ , 28(2):168–191, 1956. 

- [21] M. A. Biot and D. G. Willis. The elastic coefficients of the theory of consolidation. _Journal of applied mechanics_ , 24:594–601, 1957. 

Bibliography 

159 

- [22] P. O. Bouchard, F. Bay, and Y. Chastel. Numerical modelling of crack propagation: automatic remeshing and comparison of different criteria. _Computer methods in applied mechanics and engineering_ , 192(35):3887–3908, 2003. 

- [23] F. Brezzi. On the existence, uniqueness and approximation of saddle-point problems arising from lagrangian multipliers. _Revue française d’automatique, informatique, recherche opérationnelle. Analyse numérique_ , 8(2):129–151, 1974. 

- [24] F. Brezzi and J. Pitkäranta. On the stabilization of finite element approximations of the Stokes equations. In W. Hackbusch, editor, _Efficient Solutions of Elliptic Systems_ , volume 10 of _Notes on Numerical Fluid Mechanics_ , pages 11–19. Vieweg, Wiesbaden, 1984. 

- [25] A. Caballero, C. M. López, and I. Carol. 3D meso-structural analysis of concrete specimens under uniaxial tension. _Computer Methods in Applied Mechanics and Engineering_ , 195(52):7182–7195, 2006. 

- [26] A. Caballero, K. J. Willam, and I. Carol. Consistent tangent formulation for 3D interface modeling of cracking/fracture in quasi-brittle materials. _Computer Methods in Applied Mechanics and Engineering_ , 197(33):2804–2822, 2008. 

- [27] G. T. Camacho and M. Ortiz. Computational modelling of impact damage in brittle materials. _International Journal of Solids and Structures_ , 33(20):2899– 2938, 1996. 

- [28] P. P. Camanho, M. A. Bessa, G. Catalanotti, M. Vogler, and R. Rolfes. Modeling – 

- the inelastic deformation and fracture of polymer composites Part II: smeared crack model. _Mechanics of Materials_ , 59:36–49, 2013. 

- [29] I. Carol, C. M. López, and O. Roa. Micromechanical analysis of quasi-brittle materials using fracture-based interface elements. _International Journal for Numerical Methods in Engineering_ , 52(1-2):193–215, 2001. 

- [30] A. Carpinteri, S. Valente, G. Ferrara, and G. Melchiorrl. Is mode II fracture energy a real material property? _Computers & structures_ , 48(3):397–413, 1993. 

- [31] B. Carrier and S. Granet. Numerical modeling of hydraulic fracture problem in permeable medium using cohesive zone model. _Engineering fracture mechanics_ , 79:312–328, 2012. 

Bibliography 

160 

- [32] L. Casares, R. Vincent, D. Zalvidea, N. Campillo, D. Navajas, M. Arroyo, and X. Trepat. Hydraulic fracture during epithelial stretching. _Nature materials_ , 14(3):343–351, 2015. 

- [33] F. Cazes and N. Moës. Comparison of a phase-field model and of a thick level set model for brittle and quasi-brittle fracture. _International Journal for Numerical Methods in Engineering_ , 103(2):114–143, 2015. 

- [34] J. Cervenka. _Discrete crack modeling in concrete structures_ . University of Colorado at Boulder, 1994. 

- [35] M. Cervera and M. Chiumenti. Mesh objective tensile cracking via a local continuum damage model and a crack tracking technique. _Computer Methods in Applied Mechanics and Engineering_ , 196(1):304–320, 2006. 

- [36] M. Cervera, M. Chiumenti, and C. A. de Saracibar. Shear band localization via local j 2 continuum damage mechanics. _Computer Methods in Applied Mechanics and Engineering_ , 193(9):849–880, 2004. 

- [37] M. Cervera, L. Pelà, R. Clemente, and P. Roca. A crack-tracking technique for localized damage in quasi-brittle materials. _Engineering Fracture Mechanics_ , 77(13):2431–2450, 2010. 

- [38] R. Chambon, D. Caillerie, and T. Matsuchima. Plastic continuum with microstructure, local second gradient theories for geomaterials: localization studies. _International Journal of Solids and Structures_ , 38(46):8503–8527, 2001. 

- [39] V. T. Chau, Z. P. Bažant, and Y. Su. Growth model for large branched threedimensional hydraulic crack system in gas or oil shale. _Phil. Trans. R. Soc. A_ , 374(2078):20150418, 2016. 

- [40] M. Chiumenti, Q. Valverde, C. A. De Saracibar, and M. Cervera. A stabilized formulation for incompressible elasticity using linear displacement and pressure interpolations. _Computer Methods in Applied Mechanics and Engineering_ , 191(46):5253–5264, 2002. 

- [41] J. Choo and R. I. Borja. Stabilized mixed finite elements for deformable porous media with double porosity. _Computer Methods in Applied Mechanics and Engineering_ , 293:131–154, 2015. 

Bibliography 

161 

- [42] A. J. Chorin. A numerical method for solving incompressible viscous flow problems. _Journal of Computational Physics_ , 2(1):12–26, 1967. 

- [43] O. Coussy. _Poromechanics_ . John Wiley & Sons, 2004. 

- [44] S. C. Cowin and L. Cardoso. Interstitial flow in the hierarchical pore space architecture of bone tissue. In _Poromechanics V_ , pages 1203–1212. ASCE, 2013. 

- [45] W. J. T. Daniel. Modal methods in finite element fluid–structure eigenvalue problems. _International Journal for Numerical Methods in Engineering_ , 15(8):1161– 1175, 1980. 

- [46] R. De Boer and S. J. Kowalski. A plasticity theory for fluid saturated porous media. _International Journal of Engineering Science_ , 21:11–16, 1983. 

- [47] R. De Borst. Numerical aspects of cohesive-zone models. _Engineering fracture mechanics_ , 70(14):1743–1757, 2003. 

- [48] I. de Pouplana and E. Oñate. Combination of a non-local damage model for quasibrittle materials with a mesh-adaptive finite element technique. _Finite Elements in Analysis and Design_ , 112:26–39, 2016. 

- [49] J. H. P. De Vree, W. A. M. Brekelmans, and M. A. J. Van Gils. Comparison of nonlocal approaches in continuum damage mechanics. _Computers & Structures_ , 55(4):581–588, 1995. 

- [50] D. S. Dugdale. Yielding of steel sheets containing slits. _Journal of the Mechanics and Physics of Solids_ , 8(2):100–104, 1960. 

- [51] W. Ehlers and T. Graf. Adaptive computation of localization phenomena in geotechnical applications. _Bifurcations and Instabilities in Geomechanics, Swets, Zeitlinger_ , pages 247–262, 2003. 

- [52] F. Erdogan and G. C. Sih. On the crack extension in plates under plane loading and transverse shear. _Journal of Fluids Engineering_ , 85(4):519–525, 1963. 

- [53] A. C. Eringen. On nonlocal plasticity. _International Journal of Engineering Science_ , 19(12):1461–1474, 1981. 

- [54] A. C. Eringen and D. Edelen. On nonlocal elasticity. _International Journal of Engineering Science_ , 10(3):233–248, 1972. 

Bibliography 

162 

- [55] H. D. Espinosa and P. D. Zavattieri. A grain level model for the study of failure initiation and evolution in polycrystalline brittle materials. Part I: Theory and numerical implementation. _Mechanics of Materials_ , 35(3):333–364, 2003. 

- [56] G. Etse, A. Caggiano, and S. Vrech. Multiscale failure analysis of fiber reinforced concrete based on a discrete crack model. _International journal of fracture_ , 178(12):131–146, 2012. 

- [57] D. G. Fredlund and N. R. Morgenstern. Stress state variables for unsaturated soils. _Journal of Geotechnical and Geoenvironmental Engineering_ , 103(ASCE 12919), 1977. 

- [58] B. G. Galerkin. Series solution of some problems of elastic equilibrium of rods and plates. _Vestnik Inzhenerov i Tekhnikov_ , 19(7):897–908, 1915. 

- [59] D. Gawin and B. A. Schrefler. Thermo-hydro-mechanical analysis of partially saturated porous materials. _Engineering Computations_ , 13(7):113–143, 1996. 

- [60] G. Geißler, C. Netzker, and M. Kaliske. Discrete crack path prediction by an adaptive cohesive crack model. _Engineering Fracture Mechanics_ , 77(18):3541– 3557, 2010. 

- [61] J. Ghaboussi and E. L. Wilson. Flow of compressible fluid in porous elastic media. _International Journal for Numerical Methods in Engineering_ , 5(3):419–442, 1973. 

- [62] J. Ghaboussi, E. L. Wilson, and J. Isenberg. Finite element for rock joints and interfaces. _Journal of the Soil Mechanics and Foundations Division_ , 99(10):849– 862, 1973. 

- [63] R. E. Goodman, R. L. Taylor, and T. L. Brekke. A model for the mechanics of jointed rock. _Journal of the Soil Mechanics and Foundations Division_ , 1968. 

- [64] C. Guiducci, A. Pellegrino, J.-P. Radu, F. Collin, and R. Charlier. Numerical modeling of hydro-mechanical fracture behavior. In _NUMOG VIII_ , pages 293– 299. Balkema, 2002. 

- [65] F. Hashagen and R. De Borst. An interface element for modelling the onset and growth of mixed-mode cracking in aluminium and fibre metal laminates. _Structural Engineering and Mechanics_ , 5(6):817–837, 1997. 

Bibliography 

163 

- [66] C. Heinrich and A. M. Waas. Investigation of progressive damage and fracture in laminated composites using the smeared crack approach. _CMC-Computers Mater Continua_ , 35:155–181, 2013. 

- [67] A. Hillerborg. Numerical methods to simulate softening and fracture of concrete. In _Fracture mechanics of concrete: Structural application and numerical calculation_ , pages 141–170. Springer, 1985. 

- [68] A. Hillerborg, M. Modéer, and P.-E. Petersson. Analysis of crack formation and crack growth in concrete by means of fracture mechanics and finite elements. _Cement and concrete research_ , 6(6):773–781, 1976. 

- [69] M. Huang, S. Wu, and O. C. Zienkiewicz. Incompressible or nearly incompressible soil dynamic behaviour–a new staggered algorithm to circumvent restrictions of mixed formulation. _Soil Dynamics and Earthquake Engineering_ , 21(2):169–179, 2001. 

- [70] T. J. R. Hughes, L. P. Franca, and M. Balestra. A new finite element formulation for computational fluid dynamics: V. Circumventing the Babuška-Brezzi condition: A stable Petrov-Galerkin formulation of the Stokes problem accommodating equal-order interpolations. _Computer Methods in Applied Mechanics and Engineering_ , 59:85–99, 1986. 

- [71] J. Hult. Creep in continua and structures. In _Topics in applied continuum mechanics_ , pages 137–155. Springer, 1974. 

- [72] M. A. Hussain, S. L. Pu, and J. Underwood. Strain energy release rate for a crack under combined mode I and mode II. _Fracture Analysis_ , 560(1), 1974. 

- [73] B. Jha and R. Juanes. Coupled modeling of multiphase flow and fault poromechanics during geologic CO2 storage. _Energy Procedia_ , 63:3313–3329, 2014. 

- [74] M. Jirásek and S. Rolshoven. Comparison of integral-type nonlocal plasticity models for strain-softening materials. _International Journal of Engineering Science_ , 41(13):1553–1602, 2003. 

- [75] L. M. Kachanov. On the time of fracture under conditions of creep. _Izv. AN SSSR, Otd. Tekh. Nauk_ , (8):26–35, 1958. 

Bibliography 

164 

- [76] V. N. Kaliakin and J. Li. Insight into deficiencies associated with commonly used zero-thickness interface elements. _Computers and Geotechnics_ , 17(2):225– 252, 1995. 

- [77] A. R. Khoei, A. R. Azami, and S. M. Haeri. Implementation of plasticity based models in dynamic analysis of earth and rockfill dams: A comparison of Pastor– Zienkiewicz and cap models. _Computers and Geotechnics_ , 31:385–410, 2004. 

- [78] A. R. Khoei, O. R. Barani, and M. Mofid. Modeling of dynamic cohesive fracture propagation in porous saturated media. _International Journal for Numerical and Analytical Methods in Geomechanics_ , 35(10):1160–1184, 2011. 

- [79] A. R. Khoei, H. Moslemi, K. M. Ardakany, O. R. Barani, and H. Azadi. Modeling of cohesive crack growth using an adaptive mesh refinement via the modified-SPR technique. _International Journal of Fracture_ , 159(1):21–41, 2009. 

- [80] C. Le Bellégo. _Couplages chimie-mécanique dans les structures en béton attaquées par l’eau: Étude expérimentale et analyse numérique_ . LMT-ENS Cachan, 2001. 

- [81] B. Lecampion. An extended finite element method for hydraulic fracture problems. _Communications in Numerical Methods in Engineering_ , 25(2):121–133, 2009. 

- [82] J. Lemaitre and J. L. Chaboche. Mécanique des matériaux solides, 1985. _Dunod, Paris_ . 

- [83] Z. Lin, Z. Zhuang, X. You, H. Wang, and D. Xu. Enriched goal-oriented error estimation applied to fracture mechanics problems solved by XFEM. _Acta Mechanica Solida Sinica_ , 25(4):393–403, 2012. 

- [84] B. Loret and J. H. Prevost. Dynamic strain localization in elasto-(visco-) plastic solids, part 1. general formulation and one-dimensional examples. _Computer Methods in Applied Mechanics and Engineering_ , 83(3):247–273, 1990. 

- [85] J. Lubliner, J. Oliver, S. Oller, and E. Oñate. A plastic-damage model for concrete. _International Journal of Solids and Structures_ , 25(3):299–326, 1989. 

- [86] A. Lucantonio, G. Noselli, X. Trepat, A. DeSimone, and M. Arroyo. Hydraulic fracture and toughening of a brittle layer bonded to a hydrogel. _Physical review letters_ , 115(18):188105, 2015. 

Bibliography 

165 

- [87] X. Martínez. Micro-mechanical simulation of composite materials using the serial/parallel mixing theory. _Unpublished doctoral dissertation) Universidad Politécnica de Cataluña. Spain_ , 2008. 

- [88] A. Masud and T. J. R. Hughes. A stabilized mixed finite element method for Darcy flow. _Computer Methods in Applied Mechanics and Engineering_ , 191(39):4341– 4370, 2002. 

- [89] J. Mazars. _Application de la mécanique de l’endommagement au comportement non linéaire et à la rupture du béton de structure_ . PhD thesis, 1984. 

- [90] J. Mazars. A description of micro-and macroscale damage of concrete structures. _Engineering Fracture Mechanics_ , 25(5):729–737, 1986. 

- [91] J. Mazars and G. Pijaudier-Cabot. From damage to fracture mechanics and conversely: a combined approach. _International Journal of Solids and Structures_ , 33(20):3327–3342, 1996. 

- [92] F. Meftah and J. M. Reynouard. A multilayered mixed beam element in gradient plasticity for the analysis of localized failure modes. _Mechanics of Cohesivefrictional Materials_ , 3(4):305–322, 1998. 

- [93] N. Moës, J. Dolbow, and T. Belytschko. A finite element method for crack growth without remeshing. _International Journal for Numerical Methods in Engineering_ , 46(1):131–150, 1999. 

- [94] N. Moës, C. Stolz, P. E. Bernard, and N. Chevaugeon. A level set based model for damage growth: the thick level set approach. _International Journal for Numerical Methods in Engineering_ , 86(3):358–380, 2011. 

- [95] H. B. Mühlhaus. Shearband analysis in antigranulocytes material by cosserattheory. _Ingenieur Archiv_ , 56(5):389–399, 1986. 

- [96] H. B. Mühlhaus and E. C. Alfantis. A variational principle for gradient plasticity. _International Journal of Solids and Structures_ , 28(7):845–857, 1991. 

- [97] N. M. Newmark. A method of computation for structural dynamics. _Journal of the Engineering Mechanics Division_ , 85(3):67–94, 1959. 

Bibliography 

166 

- [98] K. L. A. Ng and J. C. Small. Behavior of joints and interfaces subjected to water pressure. _Computers and Geotechnics_ , 20(1):71–93, 1997. 

- [99] J. Oliver, M. Cervera, S. Oller, and J. Lubliner. Isotropic damage models and smeared crack analysis of concrete. In _Second international conference on computer aided analysis and design of concrete structures_ , volume 2, pages 945–958. Pineridge Press, 1990. 

- [100] J. Oliver, M. Cervera, S. Oller, and J. Lubliner. A simple damage model for concrete, including long term effects. In _Second International Conference on Computer Aided Analysis And Design of Concrete Structures_ , volume 2, pages 945–958, 1990. 

- [101] S. Oller, E. Oñate, J. Oliver, and J. Lubliner. Finite element nonlinear analysis of concrete structures using a “plastic-damage model”. _Engineering Fracture Mechanics_ , 35(1):219–231, 1990. 

- [102] E. Oñate. Derivation of stabilized equations for numerical solution of advectivediffusive transport and fluid flow problems. _Computer Methods in Applied Mechanics and Engineering_ , 151(1):233–265, 1998. 

- [103] E. Oñate. A stabilized finite element method for incompressible viscous flows using a finite increment calculus formulation. _Computer Methods in Applied Mechanics and Engineering_ , 182(3):355–370, 2000. 

- [104] E. Oñate, A. Franci, and J. M. Carbonell. Lagrangian formulation for finite element analysis of quasi-incompressible fluids with reduced mass losses. _International Journal for Numerical Methods in Fluids_ , 74(10):699–731, 2014. 

- [105] E. Oñate, S. R. Idelsohn, and C. Felippa. Consistent pressure laplacian stabilization for incompressible continua via higher-order finite calculus. _International Journal for Numerical Methods in Engineering_ , 87(1-5):171–195, 2011. 

- [106] E. Oñate, J. Rojek, R. L. Taylor, and O. C. Zienkiewicz. Finite calculus formulation for incompressible solids using linear triangles and tetrahedra. _International Journal for Numerical Methods in Engineering_ , 59(11):1473–1500, 2004. 

- [107] E. Oñate, R. L. Taylor, O. C. Zienkiewicz, and J. Rojek. A residual correction method based on finite calculus. _Engineering Computations_ , 20(5/6):629–658, 2003. 

Bibliography 

167 

- [108] E. Oñate, F. Zárate, J. Miquel, M. Santasusana, M. A. Celigueta, F. Arrufat, R. Gandikota, K. Valiullin, and L. Ring. A local constitutive model for the discrete element method. application to geomaterials and concrete. _Computational particle mechanics_ , 2(2):139–160, 2015. 

- [109] M. Ortiz and S. Suresh. Statistical properties of residual stresses and intergranular fracture in ceramic materials. _Journal of Applied Mechanics_ , 60(1):77–84, 1993. 

- [110] M. Pastor, M. Quecedo, and O. C. Zienkiewicz. A mixed displacement-pressure formulation for numerical analysis of plastic failure. _Computers & Structures_ , 62(1):13–23, 1997. 

- [111] M. Pastor, O. C. Zienkiewicz, T. L. L. Xiaoqing, and M. Huang. Stabilized finite elements with equal order of interpolation for soil dynamics problems. _Archives of Computational Methods in Engineering_ , 6(1):3–33, 1999. 

- [112] A. Pérez-Foguet, A. Rodríguez-Ferran, and A. Huerta. Numerical differentiation for local and global tangent operators in computational plasticity. _Computer Methods in Applied Mechanics and Engineering_ , 189(1):277–296, 2000. 

- [113] S. T. Pietruszczak and Z. Mroz. Finite element analysis of deformation of strainsoftening materials. _International Journal for Numerical Methods in Engineering_ , 17(3):327–334, 1981. 

- [114] G. Pijaudier-Cabot and Z. P. Bazant. Nonlocal damage theory. _Journal of Engineering Mechanics_ , 113(10):1512–1533, 1987. 

- [115] G. Pijaudier-Cabot and A. Huerta. Finite element analysis of bifurcation in nonlocal strain softening solids. _Computer methods in applied mechanics and engineering_ , 90(1):905–919, 1991. 

- [116] M. Preisig and J. H. Prévost. Stabilization procedures in coupled poromechanics problems: A critical assessment. _International Journal for Numerical and Analytical Methods in Geomechanics_ , 35(11):1207–1225, 2011. 

- [117] J. H. Prevost and B. Loret. Dynamic strain localization in elasto-(visco-) plastic solids, part 2. plane strain examples. _Computer Methods in Applied Mechanics and Engineering_ , 83(3):275–294, 1990. 

Bibliography 

168 

- [118] Y. N. Rabotnov. Damage from creep. _Zhurn. Prikl. Mekh. Tekhn. Phys_ , 2:113–123, 1963. 

- [119] Y. R. Rashid. Ultimate strength analysis of prestressed concrete pressure vessels. _Nuclear Engineering and Design_ , 7(4):334–344, 1968. 

- [120] J. Réthoré, R. De Borst, and M.-A. Abellan. A two-scale approach for fluid flow in fractured porous media. _International Journal for Numerical Methods in Engineering_ , 71(7):780–800, 2007. 

- [121] O. Reynolds. On the Theory of Lubrication and its application to Mr. Beauchamp Tower’s experiments, including an experimental determination of the viscosity of olive oil. _Proceedings of the Royal Society of London_ , 40(242-245):191–203, 1886. 

- [122] A. Rodríguez Ferran, I. Arbós, A. Huerta, et al. Adaptive analysis based on error estimation for nonlocal damage models. 2001. 

- [123] A. Rodríguez-Ferran, I. Morata, and A. Huerta. A new damage model based on non-local displacements. _International Journal for Numerical and Analytical Methods in Geomechanics_ , 29(5):473–493, 2005. 

- [124] M. Santasusana, J. Irazábal, E. Oñate, and J. M. Carbonell. The Double Hierarchy Method. a parallel 3D contact method for the interaction of spherical particles with rigid FE boundaries using the DEM. _Computational Particle Mechanics_ , 3(3):407–428, 2016. 

- [125] J. C. J. Schellekens and R. De Borst. On the numerical integration of interface elements. _International Journal for Numerical Methods in Engineering_ , 36(1):43– 66, 1993. 

- [126] B. A. Schrefler, S. Secchi, and L. Simoni. On adaptive refinement techniques in multi-field problems including cohesive fracture. _Computer Methods in Applied Mechanics and Engineering_ , 195(4):444–461, 2006. 

- [127] S. Secchi, L. Simoni, and B. A. Schrefler. Mesh adaptation and transfer schemes for discrete fracture propagation in porous materials. _International Journal for Numerical and Analytical Methods in Geomechanics_ , 31(2):331, 2007. 

Bibliography 

169 

- [128] J. M. Segura and I. Carol. On zero-thickness interface elements for diffusion problems. _International journal for numerical and analytical methods in geomechanics_ , 28(9):947–962, 2004. 

- [129] G. C. Sih. Strain-energy-density factor applied to mixed mode crack problems. _International Journal of Fracture_ , 10(3):305–321, 1974. 

- [130] J. C. Simo and J. W. Ju. Strain-and stress-based continuum damage models-I. Formulation. _International journal of solids and structures_ , 23(7):821–840, 1987. 

- [131] A. Simone. Partition of unity-based discontinuous elements for interface phenomena: computational issues. _Communications in Numerical Methods in Engineering_ , 20(6):465–478, 2004. 

- [132] D. T. Snow. _A parallel plate model of fractured permeable media_ . PhD thesis, University of California, Berkeley, 1965. 

- [133] S. H. Song, G. H. Paulino, and W. G. Buttlar. A bilinear cohesive zone model tailored for fracture of asphalt concrete considering viscoelastic bulk material. _Engineering Fracture Mechanics_ , 73(18):2829–2848, 2006. 

- [134] D. A. Spence and P. Sharp. Self-similar solutions for elastohydrodynamic cavity flow. In _Proceedings of the Royal Society of London A: Mathematical, Physical and Engineering Sciences_ , volume 400, pages 289–313. The Royal Society, 1985. 

- [135] W. Sun. A stabilized finite element formulation for monolithic thermo-hydromechanical simulations at finite strain. _International Journal for Numerical Methods in Engineering_ , 103(11):798–839, 2015. 

- [136] W. Sun, J. T. Ostien, and A. G. Salinger. A stabilized assumed deformation gradient finite element formulation for strongly coupled poromechanical simulations at finite strain. _International Journal for Numerical and Analytical Methods in Geomechanics_ , 37(16):2755–2788, 2013. 

- [137] E. Svenning. A weak penalty formulation remedying traction oscillations in interface elements. _Computer Methods in Applied Mechanics and Engineering_ , 310:460– 474, 2016. 

- [138] J. Tejchman and W. Wu. Numerical study on patterning of shear bands in a cosserat continuum. _Acta Mechanica_ , 99(1-4):61–74, 1993. 

Bibliography 

170 

- [139] K. Terzaghi. _Theoretical soil mechanics_ . Wiley, New York, 1943. 

- [140] J. C. W. Van Vroonhoven and R. De Borst. Combination of fracture and damage mechanics for numerical failure analysis. _International Journal of Solids and Structures_ , 36(8):1169–1191, 1999. 

- [141] A. Verruijt. _Theory and problems of poroelasticity_ . Delft University of Technology, 2013. 

- [142] J. A. White and R. I. Borja. Stabilized low-order finite elements for coupled soliddeformation/fluid-diffusion and their application to fault zone transients. _Computer Methods in Applied Mechanics and Engineering_ , 197(49):4353–4366, 2008. 

- [143] P. A. Witherspoon, J. S. Y. Wang, K. Iwai, and J. E. Gale. Validity of Cubic Law for fluid flow in a deformable rock fracture. _Water Resources Research_ , 16(6):1016– 1024, 1980. 

- [144] X. P. Xu and A. Needleman. Numerical simulations of fast crack growth in brittle solids. _Journal of the Mechanics and Physics of Solids_ , 42(9):1397–1434, 1994. 

- [145] H. M. Zbib and E. C. Aifantis. A gradient-dependent flow theory of plasticity: application to metal and soil instabilities. _Applied Mechanics Reviews_ , 42(11S):S295– S304, 1989. 

- [146] Q. Zeng, Z. Liu, D. Xu, H. Wang, and Z. Zhuang. Modeling arbitrary crack propagation in coupled shell/solid structures with X-FEM. _International Journal for Numerical Methods in Engineering_ , 2015. 

- [147] F. Zhou and J.-F. Molinari. Dynamic crack propagation with cohesive elements: a methodology to address mesh dependency. _International Journal for Numerical Methods in Engineering_ , 59(1):1–24, 2004. 

- [148] O. C. Zienkiewicz. Basic formulation of static and dynamic behaviour of soil and other porous media. In _Numerical methods in geomechanics_ , pages 39–55. Springer, 1982. 

- [149] O. C. Zienkiewicz and P. Bettess. Soils and other saturated media under transient, dynamic conditions: general formulation and the validity of various simplifying assumptions, 1982. 

Bibliography 

171 

- [150] O. C. Zienkiewicz, A. H. C. Chan, M. Pastor, D. K. Paul, and T. Shiomi. Static and dynamic behaviour of soils: a rational approach to quantitative solutions. I. Fully saturated problems. In _Proceedings of the Royal Society of London_ , volume 429, pages 285–309. The Royal Society, 1990. 

- [151] O. C. Zienkiewicz, C. T. Chang, and P. Bettess. Drained, undrained, consolidating and dynamic behaviour assumptions in soils. _Géotechnique_ , 30(4):385–395, 1980. 

- [152] O. C. Zienkiewicz and R. Codina. A general algorithm for compressible and incompressible flow. Part I. The split, characteristic-based scheme. _International Journal for Numerical Methods in Fluids_ , 20(8-9):869–885, 1995. 

- 

- [153] O. C. Zienkiewicz, M. Huang, and M. Pastor. Computational soil dynamics A new algorithm for drained and undrained conditions. In H. J. Siriwardane and M. M. Zaman, editors, _Computer Methods and Advances in Geomechanics_ , pages 47–59. Balkema, Rotterdam, 1994. 

- [154] O. C. Zienkiewicz, S. Qu, R. L. Taylor, and S. Nakazawa. The patch test for mixed formulations. _International Journal for Numerical Methods in Engineering_ , 23(10):1873–1883, 1986. 

- [155] O. C. Zienkiewicz, J. Rojek, R. L. Taylor, and M. Pastor. Triangles and tetrahedra in explicit dynamic codes for solids. _International Journal for Numerical Methods in Engineering_ , 43(3):565–583, 1998. 

- [156] O. C. Zienkiewicz and T. Shiomi. Dynamic behaviour of saturated porous media; The generalized Biot formulation and its numerical solution. _International journal for numerical and analytical methods in geomechanics_ , 8(1):71–96, 1984. 

- [157] O. C. Zienkiewicz, R. L. Taylor, and J. Z. Zhu. _The Finite Element Method_ , volume 1. Butterworth-Heinemann, 6th edition, 2005. 

- [158] O. C. Zienkiewicz, Y. M. Xie, B. A. Schrefler, A. Ledesma, and N. Bicanic. Static and dynamic behaviour of soils: a rational approach to quantitative solutions. II. Semi-saturated problems. In _Proceedings of the Royal Society of London_ , volume 429, pages 311–321. The Royal Society, 1990. 

- [159] G. Zimmermann and A. Reinicke. Hydraulic stimulation of a deep sandstone reservoir to develop an Enhanced Geothermal System: laboratory and field experiments. _Geothermics_ , 39(1):70–77, 2010. 

