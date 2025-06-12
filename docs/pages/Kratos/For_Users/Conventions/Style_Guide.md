---
title: Style Guide
keywords: 
tags: [Style Guide]
sidebar: kratos_for_users
summary: 
---

## Repository

### Workflow

The [Feature Branch](https://www.atlassian.com/git/tutorials/comparing-workflows#feature-branch-workflow) workflow has been adopted for this project. So the standard procedure would be:
* Create a branch for any feature we want to develop or bug we want to fix
* Work on the branch and commit/push changes to the branch
* Make a pull request when it is time to merge
* Interact with other developers and comments
* After approval, merge and **delete** the branch.

### Branch 

Clear definition and correct management of the branches are very important in order to successfully apply the feature branch workflow. Some hints for good practice are:
* Choose descriptive but not too long names like `feature-io-write-element` but **not** `feature-writing-the-element-block-in-model-part-io`
* Use dash "-" for separating words
* use standard prefixes:
    * `issue` for changes related to a given issue (#123) usually for small tasks. (like `issue-123`)
    * `test`  for experimental branch will won’t be merged. (like `test-my-stupid-idea`)
    * `feature` for new developments or enhancements. Or when deals with different issues. (like `feature-necessary-and-requested`)

* **Delete** the branch if there is no specific reason to keep it

### Community tasks

The community members (regardless of being developer or not) are highly encouraged to actively participate in the community tasks in order to improve the overall experience of the project. Following actions are welcomed:
* Frequently check for issues and pull requests
    * Comments are welcomed 
    * Avoid unnecessary noise. If everything is going well just leave it going well!
* Participate in code review of pull requests
    * It is very educative
    * Be constructive!

* Use the issue system for report a bug or propose a new feature/enhancement
    * Self assigning an issue is a way of announcing your ongoing work and get feedback.

## Project Structure
The root directory of Kratos has following main folders:
* **applications** all the applications are here
* **cmake_build** The default build directory
* **kratos** The core of the Kratos
* **libs** The output dlls of the kratos 

It is very important to avoid adding new folders to the root.

## C++ Coding Convention
Imposing a standard usually is restrictive and annoying. But having a coding convention can help in understanding better and faster the code written by a team or even by only one programmer. The coding convection used in Kratos are described briefly in following sections.

### Character encoding standard

In general it is a good practice to write codes that only use [`ASCII`](https://en.wikipedia.org/wiki/ASCII) (not extended) characters, since some compilers may be likely to experience problems trying to compile those files, specially if the encoding of the source code is not supported.

### Classes

Class names must be in lower case except the first letter of each word in upper case as separator. This format is called **UpperCamelCase**. For example:

  `class MyClassNameExample`

Try to describe the type of the class with its name. Some typical cases are:

* Elements name must be finished with <tt>Element</tt> word:

  `class MyShellElement`

* Conditions name must be finished with <tt>Condition</tt> word:

  `class My2DLoadCondition`

* Linear solvers name must be finished with <tt>Solver</tt> word:

  `class GmresSolver`

* Application classes must be started <tt>Kratos</tt> and finished with <tt>Application</tt>:

  `class KratosStructuralApplication`

It can be seen that in all cases only the first letter of each word is in capital. It's important to mention that there is NO initial letter (like "C") to distinguish a class.

### Class Members
Class members are in lower case with first letter of each word in capital starting with following prefixes:

* lower case <tt>m</tt> to indicate that they are member variables:

    `int mNumberOfElements;`

* lower case <tt>mp</tt> to indicate that they are member pointers:

    `Node::Point mpFirstNode;`

* lower case <tt>mr</tt> to indicate that they are member references:

    `Matrix& mrLocalMatrix;`

### Class Methods
Is UpperCamelCase **without** any underline "_". It has the same convention as the class name.

### Macros 
Macros in Kratos are all uppercase with underline "_" as seperator and always starting with <tt>KRATOS_</tt> as prefix:

  `#define KRATOS_AN_EXAMPLE_OF_MACRO`

### Arguments
Write arguments with a UpperCamelCase format. First input arguments, then output arguments. If the argument is passed as a reference, add 'r' at the beginning of the name. If the argument is a pointer, add 'p' at the beginning of the name.

`void CopyModelPart(ModelPart& rSourceModelPart, ModelPart& rDestinationModelPart);`

`void AddNodeToModelPart(Node::Pointer pNodeToInsert, ModelPart& rModelPart){}`

### Local Variables 
Local variables in Kratos must be written in lowercase, with underscores between words if necessary.

 `double a = 0.0;`

 `const unsigned int vector_of_neighbours_size = mVectorOfNeighbours.size();`

## Python style conventions

The following style convetion, which (mostly) follows PEP8, is recommended:

### Classes
* Class names          :`UpperCamelCase`
* Member Variable names: `lower_case_with_underscores`
* "Public" functions   : `UpperCamelCase`
* "Protected" Functions: `_UpperCamelCaseWithOneLeadingUnderscore`
* "Private" Functions  : `__UpperCamelCaseWithTwoLeadingUnderscores`

### Function arguments
* Generally: `lower_case_with_underscores`

### Local Variables
* Generally: `lower_case_with_underscores`

### Imports
Kratos follows the [PEP8 Style Guide](https://www.python.org/dev/peps/pep-0008/#imports) for importing python modules (see also [the discussion in Kratos](https://github.com/KratosMultiphysics/Kratos/issues/5402)):

**Asterisk imports \*** are absolutely prohibited unless there is no alternative. 

**Absolute imports** are the recommended way to import python modules in Kratos. This means that the imports should be done by using the full path, e.g.:

```python
import KratosMultiphysics.FluidDynamicsApplication.vms_solver as vms_solver
from KratosMultiphysics.FluidDynamicsApplication import vms_solver
```

The other - **not promoted** - way for imports are **explicit relative imports**, e.g.:

```python
# in vms_sovler.py:
from KratosMultiphysics.FluidDynamicsApplication import fluid_solver # absolute import
# OR:
from . import fluid_solver # explicit relative import
```

## JSON Configuration File
The [JSON configuration file](How-to-write-a-JSON-configuration-file) must follow certain style conventions. These are based on the [Google coding style](https://google.github.io/styleguide/jsoncstyleguide.xml), adapted for Kratos. The following notation is used:

```json
{
  "property_name": "property_value"
}
```

* A property consists of name and value
* NO comments, use descriptive names
* The name of the property should be written in lowercase with underscores "_" in between words if necessary
* Use only double quotes for property names. Same applies if the property value is a string. Other property values such as integers or booleans should not be surrounded by quotes
* Indent Sub-Properties with 4 Spaces
* Singular vs Plural Property Names: Properties that are not array types (numbers, strings, booleans) should have singular names. Array types should have singular names if they define a physical vector maginitude or plural property names if they define a generic list of items.
* Enum values should be represented as strings. 

## Good Hints in programming
Here are some basic hints to keep in mind while developing inside Kratos:
* Performance is really important but modularity is a must.
* Avoid using raw pointers and C-arrays 
* Do not use the structure new/delete or even worst malloc/free. Do memory management via `Kratos::shared_ptr` & `Kratos::unique_ptr` (declared [here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/shared_pointers.h))
* Do not re-implement linear algebra! use ublas instead
* Use the Kratos data structures as much as possible. They will simplify your life.
* Use the const attribute as much as you can ... it will avoid many errors

## Configuring your IDE for easier edition of Kratos

### Visual Studio Code
Codacy checks the style of the code every time you open a PR. It usually complains about 'trailing whitespaces'. You can configure VSCode to remove them automatically when you save the file by adding `"files.trimTrailingWhitespace":true` as an option in your settings.json (File->Preferences->Settings).
 
## Doxygen 

Commenting for classes consist of a brief description and a detailed one as follow:
```cpp
 /// Brief description.
 /** Detailed description. 
 */
```

which must be placed just befor the class definition. for example in Geometry.h:
```cpp
 /// Geometry base class.
 /** 
  As a base class Geometry has all the common
  interface of Kratos' geometries. Also it contains array of
  pointers to its points, reference to shape functions values in
  all integrations points and also local gradients of shape
  functions evaluated in all integrations points.
  
  Geometry is a template class with just one template parameter: 
  - TPointType which reperesent the type of the point this geometry
  type contain and build on.
  
  \see Point
  \see Node
  \see Formulation
  \see GeometryAndFormulationElement
*/
template<class TPointType> 
class Geometry : public PointerVector<TPointType>
{
  public:
      ....
}
```

It is very recommended to use the @see command which appears as see also with a link to indicate some related classes.

For methods of the class the same structure can be used but with information about parameters and return value if there is:

```cpp
 /** 
  Jacobian in specific integration point of given integration
  method. This method calculate jacobian matrix in given
  integration point of given integration method.
  
  \param IntegrationPointIndex index of integration point which jacobians has to
  be calculated in it.
  
  \param ThisMethod integration method which jacobians has to
  be calculated in its integration points.
  
  \return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
  i \f$ is the given integration point index of given
  integration method.
  
  \see DeterminantOfJacobian
  \see InverseOfJacobian
*/
virtual Matrix& Jacobian(Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const
{
  ....
}
```

The following commands can be used in documentation of methods and functions:

* `\param` to indicate a parameter of the method
* `\return` to indicate the return value of the method
* `\see` to include a see also 
* `\f$ \f$` Create a Latex formula in the document
* `\f[ \f]` for creating a centered latex formula

Here is an example:

```cpp
/** 
  Calculates center of this geometry by a simple averaging algorithm.
  Each center point component calculated using:
  \f[
  c_i = \sum_j^n(x_i^j) / n
  \f]

  where \f$ c_i \f$ is component i of center point and \f$
  X_i^j \f$ is component i of j'th point of geometry and n is
  number of the points in this geometry.

  @return PointType which is the calculated center of this geometry.
*/
virtual PointType Center() const
{
  ...
}
```

A quick reference to the doxygen commands can be found [here](http://www.digilife.be/quickreferences/QRC/Doxygen%20Quick%20Reference.pdf)

In `kratos/kratos/template` directory there is a `header_template` file prepared with grouping commands for class methods in Kratos.
