---
title: Clang Tidy to automatically correct code
keywords: 
tags: [How-to-use-Clang-Tidy-to-automatically-correct-code.md]
sidebar: kratos_for_developers
summary: 
---

# **Clang-Tidy: Automatically correct your code**

## **Introduction**

*Clang-Tidy* is a tool developed and maintained by the *Clang/LLVM* community. The official documentation can be found at [http://clang.llvm.org/extra/clang-tidy/](http://clang.llvm.org/extra/clang-tidy/). 

*Clang-Tidy* is a powerful tool for **automated code improvements**. It helps enforce modern C++ best practices, improve performance, and maintain high-quality code. By integrating it into your workflow, you can **automate code reviews and refactor large codebases efficiently**. 🚀

## **Install** 

When running *GNU/Linux*, *clang-tidy* is usually easy to get via your distribution’s package manager. On **Ubuntu/Linux**, Clang-Tidy can be installed using the package manager:

```sh
sudo apt-get install clang-tidy clang-tools
```
Alternatively, you can download it directly from the [**LLVM project page**](http://releases.llvm.org/download.html).

## **Using Clang-Tidy**
 
 A typical invocation of the command-line tool looks like this:
 
 ```console
 clang-tidy test.cpp -- -Imy_project/include -DMY_DEFINES ...
 ```
 
Executing it like this, the tool will print a bunch of warnings and notes (if applicable), in exactly the same way **Clang/GCC diagnostics**, too.
 
### Available checkers

The avalaible checkers avaible can be found with the command `list-checks`.
 
```sh
clang-tidy --list-checks
```
 
To automatically **fix** issues, you can use the `-fix` option. However, not all checkers support automatic fixes. For example for `modernize` and `performance`.
 
#### Modernize

To modernize code we have the following options:

 ```console
 clang-tidy --list-checks -checks='*' | grep "modernize"
    modernize-avoid-bind
    modernize-avoid-c-arrays
    modernize-concat-nested-namespaces
    modernize-deprecated-headers
    modernize-deprecated-ios-base-aliases
    modernize-loop-convert
    modernize-make-shared
    modernize-make-unique
    modernize-pass-by-value
    modernize-raw-string-literal
    modernize-redundant-void-arg
    modernize-replace-auto-ptr
    modernize-replace-disallow-copy-and-assign-macro
    modernize-replace-random-shuffle
    modernize-return-braced-init-list
    modernize-shrink-to-fit
    modernize-unary-static-assert
    modernize-use-auto
    modernize-use-bool-literals
    modernize-use-default-member-init
    modernize-use-emplace
    modernize-use-equals-default
    modernize-use-equals-delete
    modernize-use-nodiscard
    modernize-use-noexcept
    modernize-use-nullptr
    modernize-use-override
    modernize-use-trailing-return-type
    modernize-use-transparent-functors
    modernize-use-uncaught-exceptions
    modernize-use-using
 ```

##### **Example: Adding `override` Specifiers**

For example to add the missing overrides in the following code:
 
 ```c
 struct Base {
    virtual void reimplementMe(int a) {}
};
struct Derived : public Base  {
    virtual void reimplementMe(int a) {}
};
 ```
 
 We use the next command:
 
 ```console
 clang-tidy -checks='modernize-use-override' -fix test.cpp -- -std=c++11
 ```
 
#### Performance

For performance checks we have the following options:

 ```console
clang-tidy --list-checks -checks='*' | grep "performance"                                                                                                                                                                   
    clang-analyzer-optin.performance.GCDAntipattern
    clang-analyzer-optin.performance.Padding
    performance-faster-string-find
    performance-for-range-copy
    performance-implicit-conversion-in-loop
    performance-inefficient-algorithm
    performance-inefficient-string-concatenation
    performance-inefficient-vector-operation
    performance-move-const-arg
    performance-move-constructor-init
    performance-no-automatic-move
    performance-no-int-to-ptr
    performance-noexcept-move-constructor
    performance-trivially-destructible
    performance-type-promotion-in-math-fn
    performance-unnecessary-copy-initialization
    performance-unnecessary-value-param
 ```

Some commonly used performance checkers include:
- `performance-faster-string-find`
- `performance-inefficient-string-concatenation`
- `performance-move-const-arg`
- `performance-unnecessary-copy-initialization`

For example, to **optimize unnecessary copy initialization**, use:

```sh
clang-tidy -checks='performance-unnecessary-copy-initialization' -fix test.cpp
```

### **Applying fixes to an entire project**
The previous examples work on **single files**, but for larger projects, we need **Clang-Tidy to analyze all source files**.

#### **1. Generate `compile_commands.json`**
To allow Clang-Tidy to understand the project structure, enable `CMAKE_EXPORT_COMPILE_COMMANDS` in **CMake** by adding this line to your `configure.sh`:

```sh
-DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

This generates a file called **`compile_commands.json`**, which Clang-Tidy will use to analyze the project.

#### **2. Run Clang-Tidy on the Whole Project**
Download the **Clang-Tidy Python script** from the [LLVM project](https://github.com/llvm/llvm-project/blob/main/clang-tools-extra/clang-tidy/tool/run-clang-tidy.py).

Then, run:

```sh
run-clang-tidy.py -header-filter='.*' -checks='-*,modernize-use-override' -fix
```

This command:
✅ **Checks** all files in the project  
✅ **Applies** `modernize-use-override` fixes  
✅ **Uses** the compilation database (`compile_commands.json`)

#### **3. Run All Modernization Fixes**
To apply **all modernize transformations** at once, use the following **shell script**:

```sh
sh modernize.sh
```

This script is available at:  
👉 [Modernize Shell Script](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Resources_files/Clang-tidy%20modernize/modernize.sh)

 
