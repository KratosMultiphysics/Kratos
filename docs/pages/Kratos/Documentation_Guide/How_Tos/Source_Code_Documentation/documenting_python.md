---
title: Documenting Python Code
keywords: documentation, doxygen, python
tags: [documenting_python.md]
sidebar: documentation_guide
summary: 
---
## Intro

This page is meant to provide a quick primer on how <a>Kratos</a> code should be documented in-source, using doxygen.

Doxygen is the most popular tool for generating documentation from annotated source files of major languages (such as C++ or Python). It performs static code analysis on the repository to deduce a lot useful information from namespace contents to class composition and inheritance hierarchies. Although Doxygen can do this automatically without any comments or annotations, the proper description of classes and functions is essential for a useful documentation, which in turn is the cornerstone of attracting new users and developers to the project (and remembering what you were thinking a month ago when you wrote that class template).

## Documenting Python Code

For the most part, annotating python code is identical to documenting C++, save for the different comment character ```#``` and different languages features (such as package scopes instead of ```namespace```s). However, python is dynamically typed so type information must be manually provided if possible. Thankfully, python already provides a built-in remedy that is useful during development as well: [type hinting](https://peps.python.org/pep-0484/).

### Type Hints
Kratos developers seem to be oblivious to this feature of the language, so here's quick intro:
```py
# Type hinting a variable:
an_integer: int = 0 # Tell the interpreter that 'an_integer' is supposed to be an int

# Type hinting a function:
# The arguments' types can be hinted just like regular variables'
# The name followed by the arrow '->' after the argument list denotes the return type of the function.
def GetTime(model_part: KratosMultiphysics.ModelPart) -> float:
    return model_part.ProcessInfo(KratosMultiphysics.TIME)

# Type hinting optional arguments in a function:
def NewModelPart(model: KratosMultiphysics.Model, name: str = "Main") -> KratosMultiphysics.ModelPart:
    return model.CreateModelPart(name)

# Compound types (available since python3.9):
list_of_processes: list[KratosMultiphysics.Process] = [] # the interpreter will assume this list holds processes
mixed_tuple: tuple[str, int, float] = ("some_string", 1, 1.0)
material_parameters: dict[str, float] = {"density" : 2700.0, "youngs_modulus" : 6.9e10}

# Use the built-in 'typing' module for more complex hinting (available since python3.8)
import typing
cosine: typing.Callable[[float], float] = lambda x: math.cos(x)

# If the hinted type is not available in the script, or you are using an old version
# of python that doesn't support array hints, you can hint the types in strings instead.
# However, this won't be too useful for your IDE or the documentation, so this should be
# a last resort.
def GetOrigin() -> "tuple[float]":
    return (0.0, 0.0, 0.0)
```

Type hinted python code has the added benefit that your IDE will be able to provide you with more support (such as auto-complete suggestions).

### Annotations

As for annotations, the main difference is that functions' and classes' documentation must be written **inside** their **docstrings**. Variables' annotations are the same as in C++ and must immediately precede their definition. Python has no exact equivalent feature to C++'s namespace, so excluding parts of the code must be done manually. The best option to do this is to enclose the undocumented region between ```@cond impl``` and ```@endcond``` tags (the name of the conditional region in this case is ```impl``` but you are free to choose something else).

*The docstrings must not contain blank lines, otherwise doxygen commands won't be parsed in it. This is probably a doxygen bug or a configuration error on my part (let me know if you find out which).*

Here's an example and its [generated documentation](data/example_3/html/index.html) to give you an idea what python annotations should look like in practice:
```py
import typing

##@addtogroup utilities
##@{

class Memoizer:
    """@brief Utility class for caching and recalling results of expensive calculations.
    @details The memoized function must be an instance of a class implementing the @a __call__
    method. Unfortunately, normal functions are immutable in Python and so cannot be replaced.
    """

    def __init__(self, function: typing.Callable[[int], int]):
        """@brief Construct a memoizer from a callable object.
        @param function: Callable object to be memoized.
        @warning @p function must not keep internal state. Otherwise it won't
        always reproduce the same result and memoizing it won't make any sense.
        """
        self.__function: typing.Callable[[int], int] = function
        self.__cache: dict[int, int] = {}
        self.__backup: typing.Callable[[int], int] = function.__call__

    def __enter__(self) -> "Memoizer":
        """@brief Begin rerouting invocations."""
        type(self.__function).__call__ = lambda instance, argument: Memoizer.__call__(self, argument)
        return self

    def __call__(self, argument: int) -> int:
        """@brief Perform a lookup in the internal cache and invoke the wrapped function if that fails.
        @details Calls to the original function are redirected here first.
        @param argument: argument of the wrapped @a function.
        @return identical to the return value of the wrapped @a function
        """
        # Check whether this call is cached
        value: typing.Optional[int] = self.__cache.get(argument, None)
        if value == None: # It isn't => dispatch the original function and cache the result!
            value = self.__backup(argument)
            self.__cache[argument] = value
        return value

    def __exit__(self, *args) -> None:
        """@brief Stop rerouting invocations."""
        # For the sake of brevity, no error handling is done in this example.
        self.__function.__call__ = self.__backup


def Memoize(function: typing.Callable[[int], int]) -> Memoizer:
    """@brief Construct a @ref Memoizer for handling memoized contexts.
    @param function: callable to be memoized in the scope of the context.
    """
    return Memoizer(function)
##@}


##@cond impl
class FibonacciImpl:
    def __call__(self, n: int) -> int:
        return 0 if n < 1 else 1 if n==1 else self(n - 2) + self(n - 1)
##@endcond


##@defgroup maths

def Fibonacci(n: int, memoizer: Memoizer = Memoizer(FibonacciImpl())) -> int:
    """@brief Compute the @a n-th term of the fibonacci series.
    @details @f[
        F(n) = F(n-1) + F(n-2) \\
        F(0) = 0 \\
        F(1) = 1
    @f]
    @param n: 0-based index of the term to be computed.
    @param memoizer: instance of a @ref Memoizer for avoiding runaway recursion.
    @warning Term indices begin with 0!
    @ingroup maths
    """
    # A class instance is needed to overwrite __call__
    fibonacci = FibonacciImpl()
    with memoizer as MemoizedFibonacci:
        return MemoizedFibonacci(n)
```