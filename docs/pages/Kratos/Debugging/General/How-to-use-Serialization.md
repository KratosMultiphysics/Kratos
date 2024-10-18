---
title: Serialization
keywords: 
tags: [How-to-use-Serialization.md]
sidebar: kratos_debugging
summary: 
---
## Introduction
The serialization consists of storing the state of an object into a storage format like data file or memory buffer and also retrieving the object from such a media. The idea of serialization is based on saving all object's data consecutively in the file or buffer and then load it in the same order. In Kratos a serialization mechanism is used for creating the restart file. So for storing an object into a restart file and retrieving it afterward one must add the necessary components used by serialization. In the following section we will specify the necessary steps for making an object serializable in Kratos.
## Preparing a Class for Serialization
The following steps are necessary in order to prepare a class for serialization:

Adding a include for the serializer class:
```cpp
// Project includes
#include "includes/serializer.h"
```
{: data-lang="C++"}

Adding a doxygen statement for documentation of the serialization methods by adding the following comment in the **private** section of the class: (A good option is after the Member Variables section) 

```cpp
private:
///@} 
///@name Member Variables 
///@{ 

///@}
///@name Serialization
///@{

///@} 
///@name Private Operators
///@{ 
```
{: data-lang="C++"}

Making the Serializer a friend of the class by adding a friend statement in the **private** part of the class

```cpp
private:
///@} 
///@name Member Variables 
///@{ 
     
///@}
///@name Serialization
///@{

friend class Serializer;

///@} 
///@name Private Operators
///@{  
```
{: data-lang="C++"}

If the class **does not** provide a default constructor (an empty constructor) you have to provide one in the serialization part as follows:

```cpp
private:
///@}
///@name Serialization
///@{

friend class Serializer;

// A private default constructor necessary for serialization 
ClassName() : BaseClassName() { }

///@} 
///@name Private Operators
///@{
```
{: data-lang="C++"}

Adding the corresponding declarations to the the `own_application.h`:
```cpp
///@}
///@name Member Variables
///@{

// ...

const OwnConstitutive mOwnConstitutive;

//...
```
{: data-lang="C++"}

Finally include the following declarations to the the `own_application.cpp`: 
```cpp
KratosOwnApplication::KratosOwnApplication(): 

// ...

mOwnConstitutive();

// ...
```
{: data-lang="C++"}

And:
```cpp
void KratosOwnApplication::Register(){

// ...

Serializer::Register( "OwnConstitutive", mOwnConstitutive );

//...
```
{: data-lang="C++"}

## Serialization of a Non-polymorphic Class
A non-polymorphic class is a class that does not belong to any inheritance. In other words it means a class that is not derived from other classes and is not the base for other classes. Making non-polymorphic classes serializable  consists in the following steps: 

The first step is adding the serialization save and load methods as follow:
```cpp
///@}
///@name Serialization
///@{

friend class Serializer;

// A private default constructor necessary for serialization 
ClassName() : BaseClassName() { }

virtual void save(Serializer& rSerializer) const { }

virtual void load(Serializer& rSerializer) { }

///@} 
///@name Private Operators
///@{
```
{: data-lang="C++"}

Now the class is ready for adding the serialization statements in order to serialize its data. There are two overloaded methods of Serializer which are in charge of saving and loading the internal data of the class:

```cpp
Serializer::save(std::string rTag, TDataType const& Data) const
 
Serializer::load(std::string rTag, TDataType& Data)
```
{: data-lang="C++"}

Where the Tag string is the tag given to the saved or loaded data. This tag is only for text format and debugging purpose, and actually won't be written in binary format buffers and files. Here is an example of using these methods in a class with two member variables: 

```cpp
///@} 
///@name Member Variables 
///@{ 

int mMyInteger;

Vector mMyVector;
     
///@}
///@name Serialization
///@{

friend class Serializer;

// A private default constructor necessary for serialization 
ClassName() {}

virtual void save(Serializer& rSerializer) const
{
     rSerializer.save("My Integer", mMyInteger);
     rSerializer.save("My Vector", mMyVector);
}

virtual void load(Serializer& rSerializer)
{
     rSerializer.load("My Integer", mMyInteger);
     rSerializer.load("My Vector", mMyVector);
}

///@} 
///@name Private Operators
///@{
```
{: data-lang="C++"}

{% include note.html content="It is **VERY IMPORTANT** to have the **SAME SEQUENCE** of variables in the save and load method. Altering this order results in corrupted loaded object or even crash of the program!" %}

## Serialization of a Polymorphic Class
Serializing of a polymorphic class is similar to the non-polymorphic one with the following two additional steps:

{% include note.html content="Saving the name of the class **IS NOT NECESSARY** anymore. Please remove it from the old codes because it will prevent them from being serialized." %}

{% include note.html content="The default constructor mentioned in the non-Polymorphic part has to be added in **protected** section to be able to be accessed from inheriting classes." %}

In general it is necessary to call the base class save and load in the derived class. This can be done using the following defined macro:

```cpp
KRATOS_SERIALIZE_SAVE_BASE_CLASS(Serializer, BaseType)
KRATOS_SERIALIZE_LOAD_BASE_CLASS(Serializer, BaseType)
```
{: data-lang="C++"}

Here is an example of use for the `TotalLagrangian` element derived from `Element`:
```cpp
virtual void save(Serializer& rSerializer) const
{
     KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
}

virtual void load(Serializer& rSerializer)
{
     KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
}
```
{: data-lang="C++"}

The class must be registered in serializer using the Register() method. For elements and conditions the registeration is done via KRATOS_REGISTER_ELEMENT and KRATOS_REGISTER_CONDITION with the name of the class (NOT the given name respect to its geometry). Here is an example: 
```cpp
KRATOS_REGISTER_ELEMENT( "TotalLagrangian", mTotalLagrangian )
```
{: data-lang="C++"}

This is equivalent to write: 
```cpp
Serializer::Register("TotalLagrangian", mTotalLagrangian);
```
{: data-lang="C++"}

The above form can be used to register other polymorphic classes in Kratos if necessary. 

## Debugging
There is a trace option which helps debugging the serialization. The tracing option can be passed as an additional parameter to the constructor. Here is an example: 

```cpp
serializer = Serializer("filename", SERIALIZER_TRACE_ERROR)
```
{: data-lang="C++"}

The three options are:
* **SERIALIZER_NO_TRACE** disables the tracing. This is the default option and won't write any additional information to the restart file in order to minimize the size of the file. 
* **SERIALIZER_TRACE_ERROR** enables the tracing but only reports the errors. It will write the tags given in time of saving to the restart file and then in time of loading checks if the written one coincides with expected one. In case of error it will report it as follow: 

```lua
 File "cantilever3dstatic.py", line 96, in <module>
   serializer.Load("StructureModelPart", loaded_model_part);
ValueError: in bool Kratos::Serializer::load_trace_point(const std::string&) [ /home/kratos/kratos/includes/serializer.h , Line 561 ]

with subject    :  In ine 386 the trace tag is not the expected one:
   Tag found : Variables
   Tag given : Variables List
```
{: data-lang="Error Message"}

* **SERIALIZER_TRACE_ALL** also enables the tracing and saves tags to check them in loading. The only difference with **SERIALIZER_TRACE_ERROR** is the additional messages which remark the tags correctly loaded. VERY VERBOSE!!

## Examples
As an example we serialize the Isotropic3D class which is derived from ConstitutiveLaw: 

```cpp
class Isotropic3D : public ConstitutiveLaw
{
public:
     /**
     * Default constructor.
     */
     Isotropic3D();

     /**
     * Destructor.
     */
     virtual ~Isotropic3D();

     /**
     * Operations
     */

//....
//....

private:

     /**
     * Member Variables
     */
     double mE, mNU, mDE;
     Vector mInSituStress;
     Matrix mCtangent;
     Vector mCurrentStress;
     Vector mMaterialParameters;

//....
//....

}; // Class Isotropic3D
```
{: data-lang="C++"}

Now we add the serialization comment block and make Serializer a friend of this class:

```cpp
  class Isotropic3D : public ConstitutiveLaw
   {
   public:
       /**
        * Default constructor.
        */
       Isotropic3D();
 
      /**
        * Destructor.
        */
       virtual ~Isotropic3D();
 
       /**
        * Operations
        */
 
    //....
    //....
 
   private:
 
        /**
        * Member Variables
        */
       double mE, mNU, mDE;
       Vector mInSituStress;
       Matrix mCtangent;
       Vector mCurrentStress;
       Vector mMaterialParameters;
 
     ///@} 
     ///@name Serialization
     ///@{    
       friend class Serializer;
  
    //....
    //....
 
    }; // Class Isotropic3D
```
{: data-lang="C++"}

Then add the save and load methods with serializing the base classes

```cpp
  class Isotropic3D : public ConstitutiveLaw
   {
   public:
       /**
        * Default constructor.
        */
       Isotropic3D();
 
      /**
        * Destructor.
        */
       virtual ~Isotropic3D();
 
       /**
        * Operations
        */
 
    //....
    //....
 
   private:
 
        /**
        * Member Variables
        */
       double mE, mNU, mDE;
       Vector mInSituStress;
       Matrix mCtangent;
       Vector mCurrentStress;
       Vector mMaterialParameters;
 
     ///@} 
     ///@name Serialization
     ///@{    
       friend class Serializer;
 
       virtual void save(Serializer& rSerializer) const
       {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
       }
                  
       virtual void load(Serializer& rSerializer)
       {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
       }
 
    //....
    //....
 
    }; // Class Isotropic3D
```
{: data-lang="C++"}

And finally we load and save all the member variables of this class as follows:

```cpp
  class Isotropic3D : public ConstitutiveLaw
   {
   public:
       /**
        * Default constructor.
        */
       Isotropic3D();
 
      /**
        * Destructor.
        */
       virtual ~Isotropic3D();
 
       /**
        * Operations
        */
 
    //....
    //....
 
   private:
 
        /**
        * Member Variables
        */
       double mE, mNU, mDE;
       Vector mInSituStress;
       Matrix mCtangent;
       Vector mCurrentStress;
       Vector mMaterialParameters;
 
     ///@} 
     ///@name Serialization
     ///@{    
       friend class Serializer;
 
       virtual void save(Serializer& rSerializer) const
       {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
            rSerializer.save("E",mE);
            rSerializer.save("NU",mNU);
            rSerializer.save("DE",mDE);
            rSerializer.save("InSituStress",mInSituStress);
            rSerializer.save("Ctangent",mCtangent);
            rSerializer.save("CurrentStress",mCurrentStress);
            rSerializer.save("MaterialParameters",mMaterialParameters);
        }
                  
       virtual void load(Serializer& rSerializer)
       {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
            rSerializer.load("E",mE);
            rSerializer.load("NU",mNU);
            rSerializer.load("DE",mDE);
            rSerializer.load("InSituStress",mInSituStress);
            rSerializer.load("Ctangent",mCtangent);
            rSerializer.load("CurrentStress",mCurrentStress);
            rSerializer.load("MaterialParameters",mMaterialParameters);
      }
 
    //....
    //....
 
    }; // Class Isotropic3D
```
{: data-lang="C++"}