#if !defined(KRATOS_VELOCITY_FIELD_H)
#define KRATOS_VELOCITY_FIELD_H

// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

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
#include "real_functions.h"
#include "vector_field.h"
#include "swimming_DEM_application.h"


namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) VelocityField : public VectorField<3>
{
public:

KRATOS_CLASS_POINTER_DEFINITION(VelocityField);

/// Default constructor.

VelocityField():VectorField<3>(){}
/// Destructor.

virtual ~VelocityField(){}


//***************************************************************************************************************
//***************************************************************************************************************
void Evaluate(const double time,
              const array_1d<double, 3>& coor,
              array_1d<double, 3>& vector,
              const int i_thread = 0) override;

void CalculateTimeDerivative(const double time,
                             const array_1d<double, 3>& coor,
                             array_1d<double, 3>& deriv,
                             const int i_thread = 0) override;

void CalculateGradient(const double time,
                       const array_1d<double, 3>& coor,
                       array_1d< array_1d<double, 3>, 3>& gradient,
                       const int i_thread = 0) override;

void CalculateGradient(const double time,
                       const array_1d<double, 3>& coor,
                       DenseVector< double>& gradient_x,
                       DenseVector< double>& gradient_y,
                       DenseVector< double>& gradient_z,
                       const int i_thread);

double CalculateDivergence(const double time,
                           const array_1d<double, 3>& coor,
                           const int i_thread = 0) override;

void CalculateRotational(const double time,
                         const array_1d<double, 3>& coor,
                         array_1d<double, 3>& rot,
                         const int i_thread = 0) override;

void CalculateLaplacian(const double time,
                        const array_1d<double, 3>& coor,
                        array_1d<double, 3>& lapl,
                        const int i_thread = 0) override;

virtual void CalculateMaterialAcceleration(const double time,
                                           const array_1d<double, 3>& coor,
                                           array_1d<double, 3>& accel,
                                           const int i_thread = 0);

virtual void CalculateConvectiveDerivative(const double time,
                                           const array_1d<double, 3>& coor,
                                           array_1d<double, 3>& accel,
                                           const int i_thread = 0);

virtual void CalculateAccelerationFollowingTheParticle(const double time,
                                                       const array_1d<double, 3>& coor,
                                                       array_1d<double, 3>& accel,
                                                       const array_1d<double, 3>& particle_vel,
                                                       const int i_thread);

virtual void UpdateCoordinates(const double time,
                               const array_1d<double, 3>& coor,
                               const int i_thread = 0){}

virtual void LockCoordinates(const int i_thread = 0){(void)i_thread;}

virtual void UnlockCoordinates(const int i_thread = 0){(void)i_thread;}

void Evaluate(const double time,
              const DenseVector<double>& coor,
              DenseVector<double>& result,
              const int i_thread = 0) override;

void CalculateTimeDerivative(const double time,
                             const DenseVector<double>& coor,
                             DenseVector<double>& result,
                             const int i_thread = 0) override;

double CalculateDivergence(const double time, const DenseVector<double>& coor, const int i_thread = 0) override;

void CalculateRotational(const double time,
                         const DenseVector<double>& coor,
                         DenseVector<double>& result,
                         const int i_thread = 0) override;

void CalculateLaplacian(const double time,
                        const DenseVector<double>& coor,
                        DenseVector<double>& result,
                        const int i_thread = 0) override;

virtual void CalculateMaterialAcceleration(const double time,
                                           const DenseVector<double>& coor,
                                           DenseVector<double>& result,
                                           const int i_thread = 0);

virtual void CalculateConvectiveDerivative(const double time,
                                           const DenseVector<double>& coor,
                                           DenseVector<double>& result,
                                           const int i_thread = 0);

void ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed) override;

virtual void ImposeVelocityOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& container_variable);

virtual void UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread = 0){}
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
// Values

virtual double U0(const int i_thread = 0){return 0.0;}

virtual double U1(const int i_thread = 0){return 0.0;}

virtual double U2(const int i_thread = 0){return 0.0;}

// First-order derivatives

virtual double U0DT(const int i_thread = 0){return 0.0;}
virtual double U0D0(const int i_thread = 0){return 0.0;}
virtual double U0D1(const int i_thread = 0){return 0.0;}
virtual double U0D2(const int i_thread = 0){return 0.0;}

virtual double U1DT(const int i_thread = 0){return 0.0;}
virtual double U1D0(const int i_thread = 0){return 0.0;}
virtual double U1D1(const int i_thread = 0){return 0.0;}
virtual double U1D2(const int i_thread = 0){return 0.0;}

virtual double U2DT(const int i_thread = 0){return 0.0;}
virtual double U2D0(const int i_thread = 0){return 0.0;}
virtual double U2D1(const int i_thread = 0){return 0.0;}
virtual double U2D2(const int i_thread = 0){return 0.0;}

// Second-order derivatives

virtual double U0DTDT(const int i_thread = 0){return 0.0;}
virtual double U0DTD0(const int i_thread = 0){return 0.0;}
virtual double U0DTD1(const int i_thread = 0){return 0.0;}
virtual double U0DTD2(const int i_thread = 0){return 0.0;}
virtual double U0D0D0(const int i_thread = 0){return 0.0;}
virtual double U0D0D1(const int i_thread = 0){return 0.0;}
virtual double U0D0D2(const int i_thread = 0){return 0.0;}
virtual double U0D1D1(const int i_thread = 0){return 0.0;}
virtual double U0D1D2(const int i_thread = 0){return 0.0;}
virtual double U0D2D2(const int i_thread = 0){return 0.0;}

virtual double U1DTDT(const int i_thread = 0){return 0.0;}
virtual double U1DTD0(const int i_thread = 0){return 0.0;}
virtual double U1DTD1(const int i_thread = 0){return 0.0;}
virtual double U1DTD2(const int i_thread = 0){return 0.0;}
virtual double U1D0D0(const int i_thread = 0){return 0.0;}
virtual double U1D0D1(const int i_thread = 0){return 0.0;}
virtual double U1D0D2(const int i_thread = 0){return 0.0;}
virtual double U1D1D1(const int i_thread = 0){return 0.0;}
virtual double U1D1D2(const int i_thread = 0){return 0.0;}
virtual double U1D2D2(const int i_thread = 0){return 0.0;}

virtual double U2DTDT(const int i_thread = 0){return 0.0;}
virtual double U2DTD0(const int i_thread = 0){return 0.0;}
virtual double U2DTD1(const int i_thread = 0){return 0.0;}
virtual double U2DTD2(const int i_thread = 0){return 0.0;}
virtual double U2D0D0(const int i_thread = 0){return 0.0;}
virtual double U2D0D1(const int i_thread = 0){return 0.0;}
virtual double U2D0D2(const int i_thread = 0){return 0.0;}
virtual double U2D1D1(const int i_thread = 0){return 0.0;}
virtual double U2D1D2(const int i_thread = 0){return 0.0;}
virtual double U2D2D2(const int i_thread = 0){return 0.0;}

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

bool VariableIsInList(const VariablesList var_list, const VariableData& var);

/// Assignment operator.
VelocityField & operator=(VelocityField const& rOther);


///@}

}; // Class VelocityField

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_VELOCITY_FIELD_H defined
