---
title: Parallel Utilities
keywords: 
tags: [Parallel-Utilities.md]
sidebar: kratos_for_developers
summary: 
---

# Parallel Utilities
The landscape of parallelization within C++ is currently undergoing a change, as since C++11 the language is acquiring a set of standard extensions which aim to rival with OpenMP.
We propose a number of utilities which should simplify the parallelization of loops by using C++11 lambda syntax.
A number of advantages of the proposed functions are the following:
* direct iteration over containers with random access.
* allow defining custom reductions, a feature that we were currently prevented from using due to limitations of OpenMP support in MSVC.
* will (in a future) allow switching to pure-c++ implementation (no openmp), which will **allow handling exceptions in parallel!!**

## Basic, statically scheduled, for loop 
the implementation is divded, under the hood, in a partitioner class and a "for" function, both contained in the header "parallel_utilities.h"

to make an example, making a parallel for loop to set the value of a variable to a prescribed value would look in the user code as:

```cpp
block_for_each(model_part.Nodes(), [&](Node& rNode){
    noalias(rNode.FastGetSolutionStepValue(rVariable)) = value;
});
```

It is worth focusing on the use of the C++ lambda function. In the current code everything is captured by reference (the [&]) thus emulating OpenMP's default shared.
The user should refer to the C++11 documentation or to some tutorial to understand the syntax of lambdas. For example (https://www.cprogramming.com/c++11/c++11-lambda-closures.html).

A limitation of using C++11 is that the type of iterator must be specified explicitly in the lambda function. Once we eventually switch to C++14, the former will simplify to:

```cpp
block_for_each(model_part.Nodes(), [&](auto  rNode){
    noalias(rNode.FastGetSolutionStepValue(rVariable)) = value;
});
```

This type of loops is good for implementing simple cases, in which no `firstprivate` or `private` functions are needed. The point here is that passing a private value to the lambda can be trivially achieved 
by simply capturing the relevant call by argument. For example we could have made a local copy of value by doing (which would imply capturing by value the variable with the same name)
```cpp
block_for_each(model_part.Nodes(), [value](auto it) {
    noalias(it->FastGetSolutionStepValue(rVariable)) = value;
}
);
```
Note that **this usage has performance implications, since lambda capturing happens `_once` per execution of the `lambda_` and not once per thread as was the case for OpenMP**.
We plan in the future to give a more general version of this function to allow "per chunk" thread local storage

## Simple Reductions
The proposed functions also provide the essential tools to allow reductions:

The idea is that:
* reduction type is specified by the template parameter
* reduction type is specified by the return paramter of the lambda function

let's imagine that we want to find the maximum value of a given input "data_vector". This is achieved by
```cpp
double max_value = block_for_each<MaxReduction<double>>(data_vector, [](double& item){
    return item; //NOTE THAT THE RETURN THE VALUE TO BE REDUCED!
});
``` 

the same interface also allows multiple reductions to be combined at the same time. For example we could compute the max_value and sum_value at once by 

```cpp
//here we define how to combine more than one reductions at once
typedef CombinedReduction< SumReduction<double>, MaxReduction<double>> MultipleReduction; 

double sum_value,max_value;
std::tie(sum_value,max_value) =  block_for_each<MultipleReduction>(data_vector,  [&](unsigned int i){
    double to_sum = data_vector[i];
    double to_max = data_vector[i];
    return std::make_tuple( to_sum, to_max ); //note that these may have different types
});
``` 

## defining "custom" reductions
The approach described also allows users to eventually define their own custom reductions. This is as simple as defining a "Reducer" class which essnetially provides 3 functions:

* GetValue - returning the reduced value
* LocalReduce - reduction local to a thread, is the "serial" implementation of the reduce, without any need of threadsafety
* ThreadSafeReduce - threadsafe version of the reduction (will be called as few times as possible) to sync the reduced value between the threads

An example of implementation for computing at once the reduction of the max_value and of the max_absolute_value is shown in the following code snippet

```cpp
class CustomReducer{
    public:
        typedef std::tuple<double,double> value_type;
        double max_value = -std::numeric_limits<double>::max();
        double max_abs = 0.0;

        value_type GetValue()
        {
            value_type values;
            std::get<0>(values) = max_value;
            std::get<1>(values) = max_abs;
            return values;
        }

        void LocalReduce(double function_return_value){
            this->max_value = std::max(this->max_value,function_return_value);
            this->max_abs   = std::max(this->max_abs,std::abs(function_return_value));
        }
        void ThreadSafeReduce(CustomReducer& rOther){
            #pragma omp critical
            {
                this->max_value = std::max(this->max_value,rOther.max_value);
                this->max_abs   = std::max(this->max_abs,std::abs(rOther.max_abs));
            }
        }
};
```

## Index Based loops
The last relevant feature is the possibility of providing parallelization "by integer index". 
A very simple example, in which we elevate to the power 0.01 every entry in a given vector is 

```cpp
std::vector<double> data_vector(100000)
IndexPartition<unsigned int>(data_vector.size()) //here we partition the index span of data vector, and then use [i] access
    .for_each([&](unsigned int i){
        output[i] = std::pow(data_vector[i],0.01);
    });
```

This feature is important since it also opens the way to the definition of "parallel" sections (one could take the number of partitions to be equal to the number of threads to be used in total)

## TODO: 
some features that will (well, let's say "may", also depending on requests) be added in a future are
* Add a proxy for "atomic" operations
* Add version of for loops with explicit TLS storage => https://github.com/KratosMultiphysics/Kratos/pull/7261
* add possibility to merge multiple loops (emulating OpenMP's nowait)
* possibility to switch between OpenMP-based and c++11-based implementation depending on compilation macro => https://github.com/KratosMultiphysics/Kratos/pull/7262