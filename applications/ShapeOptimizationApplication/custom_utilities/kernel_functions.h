// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala https://github.com/adityaghantasala
//
//   Based on the previous implementations of filter functions and damping functions.
//   This just unifies them so they can be used everywhere in the code from the same source.
//
// ==============================================================================

#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "shape_optimization_application.h"

// ==============================================================================

namespace Kratos
{

/// Short class definition.
/**
 * Base class for all the kernel functions
*/
class KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{
    typedef array_1d<double, 3> array_3d;

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(KernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    explicit KernelFunction(const double Radius)
        : mRadius(Radius)
    {}
    /// Destructor.
    virtual ~KernelFunction()
    {}
    ///@}
    ///@name Operations
    ///@{
    /**
     * @brief Creates the object of the kernel function and returns a shared pointer to it.
     * @param FunctionType Type of the kernel function which is to be created.
     * @param Radius The radius of the kernel function.
     * @return Shared pointer to the kernel function.
     */
    static KernelFunction::Pointer New(const std::string FunctionType, const double Radius)
    {
        if(FunctionType.compare("linear"))
            return Kratos::make_shared<LinearKernelFunction>(Radius);
        else if (FunctionType.compare("cosine"))
            return Kratos::make_shared<CosineKernelFunction>(Radius);
        else if (FunctionType.compare("quartic"))
            return Kratos::make_shared<QuarticKernelFunction>(Radius);
        else if (FunctionType.compare("gaussian"))
            return Kratos::make_shared<GaussianKernelFunction>(Radius);
        else if (FunctionType.compare("constant"))
            return Kratos::make_shared<ConstantKernelFunction>(Radius);
        else
            KRATOS_ERROR<<"Unknown kernel function of type : "<<FunctionType<<std::endl;
    }

    /**
     * @brief Computes the weight of point J with the function center as point I.
     * @param iCoord Coordinates of the point I.
     * @param jCoord Coordinates of the point J.
     * @return Weight of the point J
     */
    virtual double ComputeWeight(const array_3d iCoord, const array_3d jCoord) const = 0;

    /**
     * @brief Function to get the radius of the current filter.
     * @return Radius of the filter function.
     */
    virtual inline double GetRadius() const {return mRadius;}

    /**
     * @brief Function to set the radius of the current filter.
     * @param Radius of the filter function.
     */
    virtual inline void GetRadius(const double NewRadius) {mRadius = NewRadius;}
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BaseKernelFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "BaseKernelFunction";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        PrintInfo(rOStream);
    }
    ///@}
  protected:
    ///@name Operations
    ///@{
    /**
     * @brief Gets the square of distance between to points in 3D space
     * @param iCoord Coordinates of the point I.
     * @param jCoord Coordinates of the point J.
     * @return Square of the distance between point I and point J.
     */
    double GetSquaredScalarDistance(const array_3d iCoord, const array_3d jCoord) const
    {
        array_3d dist_vector = iCoord - jCoord;
        return dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2];
    }
    ///@}
    ///@name Member Variables
    ///@{
    double mRadius;
    ///@}
}; // Class KernelFunction


/**
 * Cosine kernel Function
*/
class CosineKernelFunction : public KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(CosineKernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    explicit CosineKernelFunction(const double Radius)
        : KernelFunction(Radius)
    {}
    ///@}
    ///@name Operations
    ///@{
    double ComputeWeight(const array_3d iCoord, const array_3d jCoord) const override
    {
        const double distance = sqrt(GetSquaredScalarDistance(iCoord, jCoord));
        // Compute damping factor
        const double weight = (distance < mRadius) ? 1-(0.5*(1-cos((distance/mRadius)*PI))) : 0.0;

        return weight;
    }
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CosineKernelFunction";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "CosineKernelFunction";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        PrintInfo(rOStream);
    }
    ///@}
}; // Class CosineKernelFunction


/**
 * Linear kernel Function
*/
class LinearKernelFunction : public KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(LinearKernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    explicit LinearKernelFunction(const double Radius)
        : KernelFunction(Radius)
    {}
    ///@}
    ///@name Operations
    ///@{
    double ComputeWeight(const array_3d iCoord, const array_3d jCoord) const override
    {
        const double distance = sqrt(GetSquaredScalarDistance(iCoord, jCoord));
        // Compute weight
        const double weight = (distance < mRadius) ? (1-(distance/mRadius)) : 0.0;

        return weight;
    }
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "LinearKernelFunction";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "LinearKernelFunction";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        PrintInfo(rOStream);
    }
    ///@}
}; // Class CosineKernelFunction


/**
 * Quartic kernel Function
*/
class QuarticKernelFunction : public KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(QuarticKernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    explicit QuarticKernelFunction(const double Radius)
        : KernelFunction(Radius)
    {}
    ///@}
    ///@name Operations
    ///@{
    double ComputeWeight(const array_3d iCoord, const array_3d jCoord) const override
    {
        const double distance = sqrt(GetSquaredScalarDistance(iCoord, jCoord));
        const double numerator = distance-mRadius;
        // Compute weight
        const double weight = (distance < mRadius) ? (1-pow(numerator/mRadius,4.0)) : 0.0;

        return weight;
    }
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "QuarticKernelFunction";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "QuarticKernelFunction";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        PrintInfo(rOStream);
    }
    ///@}
}; // Class QuarticKernelFunction


/**
 * Gaussian kernel Function
*/
class GaussianKernelFunction : public KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(GaussianKernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    explicit GaussianKernelFunction(const double Radius)
        : KernelFunction(Radius)
    {}
    ///@}
    ///@name Operations
    ///@{
    double ComputeWeight(const array_3d iCoord, const array_3d jCoord) const override
    {
        const double distance = sqrt(GetSquaredScalarDistance(iCoord, jCoord));
        const double ratio = distance/mRadius;
        // Compute weight
        const double weight = (distance < mRadius) ? exp(-ratio*ratio*4.5) : 0.0;

        return weight;
    }
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GaussianKernelFunction";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "GaussianKernelFunction";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        PrintInfo(rOStream);
    }
    ///@}
}; // Class GaussianKernelFunction


/**
 * Constant kernel Function
*/
class ConstantKernelFunction : public KernelFunction
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(ConstantKernelFunction);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    explicit ConstantKernelFunction(const double Radius)
        : KernelFunction(Radius)
    {}
    ///@}
    ///@name Operations
    ///@{
    double ComputeWeight(const array_3d iCoord, const array_3d jCoord) const override
    {
        const double distance = sqrt(GetSquaredScalarDistance(iCoord, jCoord));
        // Compute weight
        const double weight = (distance < mRadius) ? 1.0 : 0.0;

        return weight;
    }
    ///@}

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ConstantKernelFunction";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ConstantKernelFunction";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        PrintInfo(rOStream);
    }
    ///@}
}; // Class ConstantKernelFunction


} // namespace Kratos.

#endif // KERNEL_FUNCTIONS_H
