#ifndef KRATOS_TIME_DEPENDANT_1D_POROSITY_FIELD_H
#define KRATOS_TIME_DEPENDANT_1D_POROSITY_FIELD_H

//#include "real_functions.h"
#include "real_field.h"
#include "real_functions.h"
#include "vector_field.h"

namespace Kratos
{

// This class implements the necessary fields for the imposition of a 1D manufactured solution

class TimeDependantPorosityField: public RealField
    {
     friend class TimeDependantForceField;
     public:

     KRATOS_CLASS_POINTER_DEFINITION(TimeDependantPorosityField);

      /// Default constructor.

     TimeDependantPorosityField(const double& max_time): mC(2 * max_time){}

      /// Destructor.

      virtual ~TimeDependantPorosityField(){}

      //***************************************************************************************************************
      //***************************************************************************************************************

      double Evaluate(const double time, const array_1d<double, 3>& coor)
      {
        return ((coor[1] - 2) / (2 * time - mC));
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor)
      {
        return (2 * (2 - coor[1]) / ((2 * time - 4) * (2 * time - mC)));
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient)
      {
        gradient[0] = 0.0;
        gradient[1] = 1.0 / (2 * time - mC);
        gradient[2] = 0.0;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      void CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& laplacian)
      {
        laplacian[0] = 0.0;
        laplacian[1] = 0.0;
        laplacian[2] = 0.0;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a stemplate<class T, std::size_t dim> tring.

        virtual std::string Info() const
        {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
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
      double mC;
     // RealFunction& mF;
     // RealFunction& mG;
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

        /// Assignment operator.
        TimeDependantPorosityField & operator=(TimeDependantPorosityField const& rOther);


        ///@}

    }; // Class TimeDependantPorosityField

class TimeDependantForceField: public VectorField<3>
{

public:

TimeDependantForceField(const double max_time): mAlpha(max_time){}


virtual ~TimeDependantForceField(){}

//***************************************************************************************************************
//***************************************************************************************************************

void Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector)
{
    array_1d<double, 3> porosity_grad;
    double porosity = mAlpha.Evaluate(time, coor);

    mAlpha.CalculateGradient(time, coor, porosity_grad);
    vector[0] = 0.0;
    vector[1] = - porosity * porosity_grad[1];
    vector[2] = 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& deriv){}
void CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d< array_1d<double, 3>, 3>& gradient){}
void CalculateDivergence(const double time, const array_1d<double, 3>& coor, double& div){}
void CalculateRotational(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& rot){}

//***************************************************************************************************************
//***************************************************************************************************************


TimeDependantPorosityField GetPorosityField(){return(mAlpha);}

virtual std::string Info() const
{
  return "";
}


virtual void PrintInfo(std::ostream& rOStream) const{}


virtual void PrintData(std::ostream& rOStream) const{}


protected:

private:

TimeDependantPorosityField mAlpha;

TimeDependantForceField & operator=(TimeDependantForceField const& rOther);

}; // Class TimeDependantForceField
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_TIME_DEPENDANT_1D_POROSITY_FIELD_H
