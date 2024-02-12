---
title: Creating a Base Application
keywords: 
tags: [Creating Base Application]
sidebar: kratos_for_users
summary: 
---

## Quick start

For more advanced use cases in Kratos you will probably want to create a custom application that better adapts to your necessities. In this section we will review the process of creating a simple empty application and how to integrate it with Kratos so you can start your development. In following tutorials, you will find how to add content to this application and a more detailed view of its parts.

The process of creating such a basic application is very simple and it is done through the create application script. Simply navigate to `kratos/python_scripts/application_generator`.

You will find a couple of python files there. If you are the **TL;DR** kind of person you can execute right away the file `createApplication.py` with your application name and it will generate your application or you can jump to "In detail" to learn how to customize your app. It is important that your application name is in camel case format (first letter of every word in caps). For example:

```console
python createApplication.py MyExample
```
{: data-lang="console"}

That's it! Your application has been generated, You can check that it has been created in `application/MyExampleApplication`

In order to compile it just add:

```console
-DMY_EXAMPLE_APPLICATION=ON \
```
{: data-lang="Console"}

In your configure file.

You can now jump to the [Basic application tutorial](Stationary_Heat_Transfer) or read the details as well.

## In detail

The `createApplication.py` file is only a set of instructions that the application generator will use to create your app. While it cannot code for you, it can provide some of the most elemental classes you will need in your application.

### Imports
The first block we will find is the import block. You should not touch this unless you really know what you do or have used python before:

### The generator class
The second block is the generator itself. It will read the name of your app from the terminal and initialize an `ApplicationGenerator` object with it

```python
# Read the application name and generate Camel, Caps and Low
appNameCamel = sys.argv[1]

# Fetch the applications directory
debugApp = ApplicationGenerator(appNameCamel)
```
{: data-lang="Python"}

This is the very minimum you need to create the app, but in the example we also make use of some additional features:

#### Variables
```python
# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='1D_VARIABLE', vtype='double'),
    VariableCreator(name='3D_VARIABLE', vtype='double', is3D=True),
])
```
{: data-lang="python"}

The `AddVariables` function let you add an initial set of kratos variables. You only need to specify the **name** of the variable, its **type** and set the `is3D` to true if you want to be a variable with components X, Y and Z. The variable will be created, defined and registered in kratos automatically, and also exported into the python interface. We recommend you to define at least a couple of them so you can check how the code is organized once the application is created.

#### Elements
```python
# Add test element
debugApp.AddElements([
    ElementCreator('CustomTestElement')
    .AddDofs(['DOF_1', 'DOF_2'])
    .AddFlags(['FLAG_1', 'FLAG_2'])
    .AddClassMemberVariables([
        ClassMemberCreator(name='Variable1', vtype='double *', default='nullptr'),
        ClassMemberCreator(name='Variable2', vtype='int', default='0'),
        ClassMemberCreator(name='Variable3', vtype='std::string &', default='0'),
        # ClassMemberCreator(name='Warnvar1', vtype='std::string &',)
    ])
])
```
{: data-lang="Python"}

The `AddElements` function will tell the generator to also add an element when generating your application code. This is useful if you know that you will need one. The element (or elements) can be created using the `ElementCreator` class. This class also provide methods to add **dofs**, **flags** and even class members. Is not the aim of the tutorial to go through all this options but feel free to play a bit with them.

### Conditions
```python
debugApp.AddConditions([
    ConditionCreator('CustomTestCondition')
])
```
{: data-lang="Python"}

The `AddCondition` function is the analogue of `AddElements` but for conditions. It allows you to do the same but will generate a condition instead.

Finally you can create processes. There is no custom options for processes yet but its name and you will probably need at least one.

### Generating the app
```python
debugApp.Generate()
```
{: data-lang="Python"}

All its left is call the `Generate` function and the code will be created. If there has been any error it will be reported during the creation with a suggested solution when possible. It will also warn you if any of the variables, elements, etc.. you try to create already exist on Kratos or have a very similar name. The generator will also tell you if you application has been already added to Kratos (for example if you are not satisfied with the result and you delete it).

Now you can continue to [Basic application tutorial](www.google.com).

## Removing an application

It is normal that at the beginning you make tests or by some reason you are not satisfied with the way an application is going and you want to remove form kratos. If you want an application to be removed from kratos, you will have to take these steps:

1 - Remove it from the configure file
2 - Remove the folder in applications. For example `application/MyExampleApplication`
3 - Delete it from `application/CMakeLists.txt`. you will find in a couple of places, typically at the bottom of every block:

```CMake
message("MY_EXAMPLE_APPLICATION..................... ${MY_EXAMPLE_APPLICATION}")

# ...

if(${MY_EXAMPLE_APPLICATION} MATCHES ON)
  add_subdirectory(MyExampleApplication)
endif(${MY_EXAMPLE_APPLICATION} MATCHES ON)
```
{: data-lang="CMake"}

4 - Delete if from `application/applications_interface.py`. Again it will appear in several blocks:

```python
Import_MyExampleApplication = False

# ...

print("Import_MyExampleApplication: False")

# ...

print("Import_MyExampleApplication: " + str(Import_MyExampleApplication))

# ...

if(Import_MyExampleApplication):
    print("importing KratosMyExampleApplication ...")
    sys.path.append(applications_path + '/MyExample/python_scripts')
    sys.path.append(applications_path + '/MyExample/Linux')
    from KratosMyExampleApplication import *
    my_example_application = KratosMyExampleApplication()
    kernel.AddApplication(my_example_application)
    print("KratosMyExampleApplication Succesfully imported")

# ...
# And finally:

if(Import_MyExampleApplication):
    kernel.InitializeApplication(my_example_application)

```
{: data-lang="Python"}

## Troubleshot

Here is a list of common problems that can appear while using the generator:

### File exists error:

- **Symptom**:

```console
FileExistsError: [Errno 17] File exists: '/path/Kratos/kratos/python_scripts/application_generator
/utils/../../../../applications/MyExampleApplication'
```
{: data-lang="Console"}

- **Cause**:
The application name you used is used by another application

- **Solution**:
Use another name.
