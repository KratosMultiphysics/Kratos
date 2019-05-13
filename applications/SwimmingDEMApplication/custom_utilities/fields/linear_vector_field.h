#if !defined(KRATOS_LINEAR_VECTOR_FIELD_H)
#define KRATOS_LINEAR_VECTOR_FIELD_H

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "vector_field.h"
#include "includes/variables.h"
#include "includes/model_part.h"

namespace Kratos
{
class LinearVectorField : public VectorField3
{
public:

KRATOS_CLASS_POINTER_DEFINITION(LinearVectorField);

/// Default constructor.

LinearVectorField(): VectorField3()
{
    mb = ZeroVector(3);
    for (auto i = 0; i < 3; ++i){
        for (auto j = 0; j < 3; ++j){
            mA(i, j) = 0.0;
        }
    }
}

LinearVectorField(Parameters rParameters)
{
    Parameters default_parameters( R"(
        {
            "A": [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "b" : [0.0, 0.0, 0.0]
        }
        )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(rParameters["A"].IsArray()
                    && (rParameters["A"].GetMatrix()).size1() == 3
                    && (rParameters["A"].GetMatrix()).size2() == 3) << "A 3 by 3 matrix should be provided in the parameters to define A." << std::endl;


    const auto A = rParameters["A"].GetMatrix();

    KRATOS_ERROR_IF_NOT(rParameters["b"].IsVector()
                    && (rParameters["b"].GetVector()).size() == 3)  << "A 3-vector should be provided in the parameters to define b." << std::endl;

    const auto b = rParameters["b"].GetVector();

    noalias(mb) = b;
    for (auto i = 0; i < 3; ++i){
        for (auto j = 0; j < 3; ++j){
            mA(i, j) = A(i, j);
        }
    }
}
/// Destructor.

virtual ~LinearVectorField(){}


//***************************************************************************************************************
//***************************************************************************************************************

void Evaluate(const double time,
              const array_1d<double, 3>& coor,
              array_1d<double, 3>& vector,
              const int i_thread = 0) override;

void CalculateGradient(const double time,
                       const array_1d<double, 3>& coor,
                       array_1d< array_1d<double, 3>, 3>& gradient,
                       const int i_thread = 0) override;

double CalculateDivergence(const double time,
                           const array_1d<double, 3>& coor,
                           const int i_thread = 0) override;

void CalculateRotational(const double time,
                         const array_1d<double, 3>& coor,
                         array_1d<double, 3>& rot,
                         const int i_thread = 0) override;

void Evaluate(const double time,
              const DenseVector<double>& coor,
              DenseVector<double>& vector,
              const int i_thread = 0) override;

void CalculateTimeDerivative(const double time,
                             const DenseVector<double>& coor,
                             DenseVector<double>& result,
                             const int i_thread = 0) override {}

double CalculateDivergence(const double time,
                           const DenseVector<double>& coor,
                           const int i_thread = 0) override;

void CalculateRotational(const double time,
                         const DenseVector<double>& coor,
                         DenseVector<double>& rot,
                         const int i_thread = 0) override;

//virtual void ImposeVelocityOnNodes(ModelPart& r_model_part, const VariableData& container_variable){}

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const override
{
    return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const override
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const override
{
}


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{


///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

private:

///@name Static Member r_variables
///@{


///@}
///@name Member r_variables
///@{

BoundedMatrix<double, 3, 3> mA;
array_1d<double, 3> mb;

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{


///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessible methods
///@{

/// Assignment operator.
LinearVectorField & operator=(LinearVectorField const& rOther);


///@}

}; // Class LinearVectorField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_LINEAR_VECTOR_FIELD_H defined
