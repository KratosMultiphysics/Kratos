---
title: How to apply a function to a list of variables
keywords: 
tags: [How-to-apply-a-function-to-a-list-of-variables.md]
sidebar: kratos_for_developers
summary: 
---

Variables in Kratos can have different types. This makes them difficult to be used in a vector and passed to a function. The main difficulty is that the standard vector cannot hold entities with different types. On the other side, the VariablesList class can store variables with different types but you may access them as a VariableData which is the base class and does not have the type information needed to manipulate data.

In fact, the [Visitor Pattern](https://en.wikipedia.org/wiki/Visitor_pattern) is the standard way to overcome such situation. The main idea behind this is to let each variable accept a function class (A.K.A functional) and then call the visit function of it using its known type. 

This pattern is implemented in Kratos and variables accepting a functional which should be derived from the `VariablesListVisitorBase` and implementing the `Visit` methods for each type of the variables which wants to work with them. For the list of the supported types in the base class, one could verify directly the definition of the `VariablesListVisitorBase` class but the list includes `bool`, `int`, `double`, `array_1d<double,3>`, `Vector`, and `Matrix`. It is important to emphasize that if the visitor is applied to a variable that is not supported by the derived class, the base class will send a runtime error.

In order to see better the details, let's create a class that sets the nodal values in a model part to one. First of all I should create a class derived from the `VariablesListVisitorBase`:

```cpp
class SetModelPartVariableToOne : public VariablesListVisitorBase
{
    ModelPart& mrModelPart;
  public:
    SetModelPartVariableToOne(ModelPart& TheModelPart) : mrModelPart(TheModelPart){}
};
```
As you can see, we take the modelpart we want to act over in the constructor. Now we should implemente the `Visit` methods:

 ```cpp
class SetModelPartVariableToOneVisitor : public VariablesListVisitorBase
{
    ModelPart& mrModelPart;
  public:
    SetModelPartVariableToOneVisitor(ModelPart& TheModelPart) : mrModelPart(TheModelPart){}

    virtual void Visit(Variable<bool> const& TheVariable) override {
        SetScalarToOne(TheVariable);
    }
    virtual void Visit(Variable<int> const& TheVariable) override {
        SetScalarToOne(TheVariable);
    }
    virtual void Visit(Variable<unsigned int> const& TheVariable) override {
        SetScalarToOne(TheVariable);
    }
    virtual void Visit(Variable<double> const& TheVariable) override {
        SetScalarToOne(TheVariable);
    }
    virtual void Visit(Variable<array_1d<double,3>> const& TheVariable) override {
        SetVectorToOne(TheVariable);
    }
    virtual void Visit(Variable<Vector> const& TheVariable) override {
        SetVectorToOne(TheVariable);
    }
    virtual void Visit(Variable<Matrix> const& TheVariable) override {
        // Do nothing.
    }
  private:
    template<typename TVariableType> 
    void SetScalarToOne()(TVariableType const& TheVariable) const {
        for (auto& r_node : mrModelPart.Nodes()) {
            auto& r_value = r_node.FastGetSolutionStepValue(TheVariable);
            r_value = 1.0;    
        }                     
    }

    template<typename TVectorType> 
    void SetVectorToOne()(Variable<TVectorType> const& TheVariable) const {
        for (auto& r_node : mrModelPart.Nodes()) {
            auto& r_value = r_node.FastGetSolutionStepValue(TheVariable);
            r_value = ScalarVector(3,-1.0);    
        }                     
    }   
};
```
You may observe that the real implementation is in the two template private methods which take a variable and set its nodal value in the model part to one. It is important to mention that one can skip any type of the variable just by creating an empty `Visit` method for that type. (like the one for matrix). This would avoid a runtime error from the base class if the variable list contains a `Variable<Matrix>` 

Now to use this one may simply call the VariablesList::ApplyVisitor method of an existing VariablesList (like the one of the modelpart):

```cpp
SetModelPartVariableToOneVisitor set_to_one_visitor(model_part);
my_variables_list.ApplyVisitor(set_to_one_visitor);
```

or we may make a container of pointers to the VariableData on the fly:


```cpp
VariablesList::VariablesContainerType my_variables_list{&TEMPERATURE, &VELOCITY, &DISPLACEMENT_X, &DISPLACEMENT_Z};
SetModelPartVariableToOneVisitor set_to_one_visitor(model_part);
set_to_one_visitor.VisitVariables(my_variables_list);
```
Now the visitor would be applied to the nodal temperature,velocity, and the displacement x and z components (but not the y component)

Please note the `&` before the name of the variables. 


 



