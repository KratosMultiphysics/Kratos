//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: c.karacaova
//   Date:                $Date: 2012-08-11  $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_SPHPARTICLE_H_INCLUDED)
#define  KRATOS_SPHPARTICLE_H_INCLUDED 

// System includes 
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>


// Project includes
#include "meshless_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

#include "includes/model_part.h"
#include "includes/node.h"
#include "geometries/point_3d.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/neighbours_calculator_SPH.h"
#include "includes/serializer.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"

#include "includes/kratos_application.h"  /* ?? */

#include "meshless_application_variables.h"





// External includes 
#include "boost/smart_ptr.hpp"


// Project includes

#define PI 3.14159265

namespace Kratos
{

class KernelGaussian
{
public:
    static double ComputeKernel(double r, double h){

        double q;
        q=r/h;
        double Kernel;
        Kernel = ( 1.0/ (PI*h*h) ) * exp(-1.0*q*q);
        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);
        double q;
        q=r/h;

        double KernelValue;

        KernelValue=( 1.0/ (PI*h*h) ) * exp(-1.0*q*q) * (-2.0 * q * (1.0/h) );


        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }

};

class KernelPoly6
{
public:
    //static const unsigned int DomainSize = 2;

    static double ComputeKernel(double r, double h){
        double Kernel=0.0;
        if (r<0.0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=h) {Kernel=(4.0/(PI*pow(h,8))) * (pow((h*h-r*r),3)); }
        else if (r>h) {Kernel=0.0;}
        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);

        double KernelValue=0.0;
        if (r<0.0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=h) {KernelValue=( (4.0/(PI*pow(h,8))) * (-6.0 * r * pow((h*h-r*r),2)) ); }
        else if (r>h) {KernelValue=0.0;}

        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }

};

class KernelQuadratic
{
public:
    static double ComputeKernel(double r, double h){

        double q;
        q=r/h;

        double Kernel=0.0;

        if (q<0.0) {std::cout<<"WARNING: q is negative";}
        else if (q>=0.0 && q<=2.0) {Kernel=(2.0/(PI*h*h)) * ((3.0/16.0)*q*q - 0.75*q + 0.75);}
        else if (q>2.0) {Kernel=0.0;}

        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);

        double KernelValue=0.0;
        double q;
        q=r/h;

        if (q<0.0) {std::cout<<"WARNING: q is negative";}
        else if (q>=0.0 && q<=2.0) {KernelValue=(2.0/(PI*h*h)) * ((3.0/8.0)*q*(1.0/h) - 0.75*(1.0/h));}
        else if (q>2.0) {KernelValue=0.0;}

        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }

};

class KernelQuintic
{
public:
    static double ComputeKernel(double r, double h){

        double q;
        q=r/h;

        //Calculate the Value
        double Kernel=0.0;



        if (q<0.0) {std::cout<<"WARNING: q is negative";}
        else if (q>=0.0 && q<1.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5) - 6.0*pow(2.0-q,5) + 15.0*pow(1.0-q,5)); }
        else if (q>=1.0 && q<2.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5) - 6.0*pow(2.0-q,5)); }
        else if (q>=2.0 && q<3.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5)); }
        else if (q>=3.0) {Kernel=0.0;}

        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);

        double KernelValue=0.0;
        double q;
        q=r/h;

        if (q<0.0) {std::cout<<"WARNING: q is negative (GRADIENT) ";}
        else if (q>=0.0 && q<1.0) {KernelValue=(7.0/(478.0*PI*h*h)) * ( ( 5.0 * pow(3.0-q,4) * (-1.0/h) ) - ( 30.0 * pow(2.0-q,4) * (-1.0/h) ) + ( 75.0 * pow(1.0-q,4) * (-1.0/h) ) ); }
        else if (q>=1.0 && q<2.0) {KernelValue=(7.0/(478.0*PI*h*h)) * ( ( 5.0 * pow(3.0-q,4) * (-1.0/h) ) - ( 30.0 * pow(2.0-q,4) * (-1.0/h) ) ) ; }
        else if (q>=2.0 && q<3.0) {KernelValue=(7.0/(478.0*PI*h*h)) * ( ( 5.0 * pow(3.0-q,4) * (-1.0/h) ) ) ;}
        else if (q>=3.0) {KernelValue=0;}

        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }

};

class KernelC2
{
public:
    static double ComputeKernel(double r, double h){
        double Kernel=0.0;
        double q;
        q= r/h;
        if (q<0.0) {std::cout<<"WARNING: q is negative";}
        else if (q>=0.0 && q<1.0) {Kernel=(10.0/(7.0*PI*h*h)) * (1.0-1.5*pow(q,2)+0.75*pow(q,3)); }
        else if (q>=1.0 && q<2.0) {Kernel=(10.0/(7.0*PI*h*h)) * (0.25 * pow(2.0-q,3));}
        else if (q>=2.0) {Kernel=0.0;}


        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);

        double q;
        q=r/h;
        double KernelValue=0.0;

        if (q<0.0) {std::cout<<"WARNING: q is negative (GRADIENT) ";}
        else if (q>=0.0 && q<1.0) {KernelValue=(10.0/(7.0*PI*h*h)) * (-3.0 * q * (1.0/h) + 2.25 * q * q * (1.0/h)); }
        else if (q>=1.0 && q<2.0) {KernelValue=(10.0/(7.0*PI*h*h)) * (-0.75 * (1.0/h) * pow(2.0-q,2) );}
        else if (q>=2.0) {KernelValue=0;}

        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }

};

class KernelSpiky
{
public:
    static double ComputeKernel(double r, double h){
        double Kernel=0.0;


        if (r<0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=h) {Kernel=(10.0/(PI*pow(h,5))) * (pow((h-r),3)); }
        else if (r>h) {Kernel=0;}


        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);


        double KernelValue=0.0;

        if (r<0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0 && r<=h) {KernelValue=(10.0/(PI*pow(h,5))) * (-3* (h-r) * (h-r)) ; }
        else if (r>h) {KernelValue=0;}


        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }
};

class KernelSpikyAlt
{
public:
    static double ComputeKernel(double r, double h){
        double Kernel=0.0;

        if (r<0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=(2*h)) {Kernel=(5.0/(PI*pow((4.0*h),2))) * (pow((2.0-(r/h)),3)); }
        else if (r>(2*h)) {Kernel=0;}


        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);

        double KernelValue=0.0;

        if (r<0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=(2*h)) {KernelValue=(5.0/(PI*pow((4.0*h),2))) * (3 * (-1/h) * pow((2.0-(r/h)),2)); }
        else if (r>(2*h)) {KernelValue=0;}


        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }


};

class KernelMullerViscous
{
public:
    static double ComputeKernel(double r, double h){
        double Kernel=0.0;

        if (r<0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=(h)) {Kernel=(10.0/(3.0 * PI * h * h )) * ( -0.5 * (pow(r/h,3) + pow(r/h,2) + h/(2*r) - 1  ) ); }
        else if (r>(h)) {Kernel=0;}


        return Kernel;
    }

    static array_1d<double,3> ComputeKernelDerivative(array_1d<double,3> rvec, double h){
        double r = norm_2(rvec);

        double KernelValue=0.0;

        if (r<0) {std::cout<<"WARNING: distance is negative";}
        else if (r>=0.0 && r<=h) {KernelValue=(10.0/(3.0 * PI * h * h )) * ( -1.5 * ((r*r)/(h*h*h)) + 2 * (r/(h*h)) + -0.5 * (h/(r*r)) ); }
        else if (r>h) {KernelValue=0;}


        array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        return w;
    }


};

template< class K_dens,class K_pres,class K_visc >
class SPHparticle
        : public Element
{
public:

    /// Counted pointer of SPHparticle
    KRATOS_CLASS_POINTER_DEFINITION(SPHparticle);

    typedef SPHparticle                     ParticleType;
    typedef Neighbours_Calculator_SPH<ParticleType> NeighboursCalculatorType;

    typedef WeakPointerVector<Element> ParticleWeakVectorType;
    typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
    typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

    typedef Node < 3 > PointType;
    typedef PointerVector<PointType > PointVector;
    typedef PointVector::iterator PointIterator;

    typedef ModelPart::NodesContainerType::iterator NodeIterator;




    /// Constructors.
    //************************************************************************************
    //************************************************************************************
    SPHparticle(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    SPHparticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {

    }
    /// Destructor.
    virtual ~ SPHparticle()
    {}



    /// Create a new Element of this type

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new SPHparticle<K_dens, K_pres, K_visc>(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }



    /// Element Functions
    //************************************************************************************
    //************************************************************************************
    void Initialize()
    {
        KRATOS_TRY;


//        Collacative Density approach

//        ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
//        double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);

//        //compute the initial density
//        double& density = GetGeometry()(0)->GetValue(DENSITY);
//        density = 0;
//        array_1d<double,3> rvec;
//        for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
//        {
//            double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;
//            noalias(rvec) = (neighbour_iterator->GetGeometry()[0].Coordinates());
//            noalias(rvec) -= this->GetGeometry()(0)->Coordinates();
//            double r = norm_2(rvec);

//            //Calculate the Value
//            double Kernel;
//            Kernel = K_dens::ComputeKernel(r,h);

//            density += mass * Kernel;
//        }


//        this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) = density;


        //Continuum density approach

                this->GetGeometry()(0)->GetValue(DENSITY) = this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;



        if ( this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE)==0.0 )
        {

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);
            const array_1d<double,3>& vel1 = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
            const double density1=this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
            const double p1=this->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) ;
            const double v1=this->GetGeometry()(0)->FastGetSolutionStepValue(VISCOSITY) ;


            array_1d<double,3> aux_sum = ZeroVector(3);
            array_1d<double,3> aux_sum_visc = ZeroVector(3);
            array_1d<double,3> rvec;
            array_1d<double,3> w1;
            array_1d<double,3> w2;

            const double parameter = 0.01 * h * h;



            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                const double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;
                const double p2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) ;
                const double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
                const double v2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(VISCOSITY) ;
                const array_1d<double,3>& vel2 = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

                array_1d<double,3> u = vel1;
                noalias(u) -= vel2 ;


                noalias(rvec) = this->GetGeometry()(0)->Coordinates();
                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());
                const double r = norm_2(rvec);



                w1 = K_pres::ComputeKernelDerivative(rvec,h);


                w2 = K_visc::ComputeKernelDerivative(rvec,h);



                //Artificial visc

                double art_visc = 0;


                double vDOTr = u[0]*rvec[0] + u[1]*rvec[1] + u[2]*rvec[2] ;

                if (vDOTr < 0.0){
                    double phi = h * vDOTr / (r*r + parameter);
                    double dens_aver = (density2 + density1)/2.0;
                    art_visc = 1.0 * phi * phi / dens_aver;
                }

                //PRESSURE
                const double same_valueP = mass * ( p1/(density1*density1) + p2/(density2*density2) + art_visc);

                noalias(aux_sum ) += same_valueP*w1;





                //VISCOSITY
                //                double same_valueV = mass * (8 * ((v1+v2)/(density1+density2)) * (( u[0]*rvec[0] +u[1]*rvec[1]+u[2]*rvec[2]) / (r*r + parameter)) );
                //                noalias(aux_sum_visc ) += same_valueV*w2;


                const double same_valueV = mass * ( ( (v1+v2) / (density1*density2) ) * ( ( w2[0]*rvec[0] +w2[1]*rvec[1]+w2[2]*rvec[2]) / (r*r + parameter) ) );
                noalias(aux_sum_visc ) += same_valueV*u;


            }





            //// ASSING THE PRESSURE ACCELERATION
            array_1d<double,3>& pres_acc = this->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_ACC);
            noalias(pres_acc) = -density1*aux_sum;

            //// ASSING THE VISCOUS ACCELERATION
            array_1d<double,3>& visc_acc = this->GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_ACC);
            noalias(visc_acc) = density1*aux_sum_visc;


            // ADD EVERYTHING UP AND ASSIGN TO RHS VARIABLE
            const array_1d<double,3>& bodyf_acc = this->GetGeometry()[0].FastGetSolutionStepValue(BODYFORCE_ACC);
            array_1d<double,3>& rhs = this->GetGeometry()[0].FastGetSolutionStepValue(RHS);
            noalias(rhs) = density1*bodyf_acc;
            noalias(rhs) += visc_acc;
            noalias(rhs) += pres_acc;
            rhs /= density1;







        }

        ApplyBoundaryForces();

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void NormalizeRightHandSide()
    {
        KRATOS_TRY;

        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE)==0)
        {

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);

            const double density1=this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;



            Matrix normalizer_RHS=ZeroMatrix(3,3);

            array_1d<double,3> rvec;




            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                const double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;
                const double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;



                noalias(rvec) = (neighbour_iterator->GetGeometry()[0].Coordinates());
                noalias(rvec) -= this->GetGeometry()(0)->Coordinates();

                array_1d<double,3> w;
                w = K_dens::ComputeKernelDerivative(rvec,h);

                const double massDBdensity2 = mass/density2;
                normalizer_RHS(0,0) += massDBdensity2 * rvec[0]*w[0];
                normalizer_RHS(0,1) += massDBdensity2 * rvec[0]*w[1];
                normalizer_RHS(0,2) += massDBdensity2 * rvec[0]*w[2];
                normalizer_RHS(1,0) += massDBdensity2 * rvec[1]*w[0];
                normalizer_RHS(1,1) += massDBdensity2 * rvec[1]*w[1];
                normalizer_RHS(1,2) += massDBdensity2 * rvec[1]*w[2];
                normalizer_RHS(2,0) += massDBdensity2 * rvec[2]*w[0];
                normalizer_RHS(2,1) += massDBdensity2 * rvec[2]*w[1];
                normalizer_RHS(2,2) += massDBdensity2 * rvec[2]*w[2];

            }

            normalizer_RHS = -normalizer_RHS;

            // Chapuza
            Matrix normalizer_RHS_2D = ZeroMatrix(2,2);

            normalizer_RHS_2D(0,0) = normalizer_RHS(0,0);
            normalizer_RHS_2D(0,1) = normalizer_RHS(0,1);
            normalizer_RHS_2D(1,0) = normalizer_RHS(1,0);
            normalizer_RHS_2D(1,1) = normalizer_RHS(1,1);



            const double tol = 0.001;

            double determinant_2D = MathUtils<double>::Det2(normalizer_RHS_2D);

            Matrix Inverted_normalizer_RHS_2D = ZeroMatrix(2,2);

            if ( fabs(determinant_2D) < tol){
                Inverted_normalizer_RHS_2D(0,0) = 1.0;
                Inverted_normalizer_RHS_2D(1,1) = 1.0;
            }
            else{

                MathUtils<double>::InvertMatrix2(normalizer_RHS_2D,Inverted_normalizer_RHS_2D,determinant_2D);
            }


            //            NORMALIZE THE RHS


            const array_1d<double,3>& bodyf_acc = this->GetGeometry()[0].FastGetSolutionStepValue(BODYFORCE_ACC);


            array_1d<double,3>& rhs      = this->GetGeometry()[0].FastGetSolutionStepValue(RHS);


            noalias(rhs) -= bodyf_acc;

            rhs *= density1;


            array_1d<double,2> rhs_2D;

            rhs_2D[0] = rhs[0];   rhs_2D[1] = rhs[1];

            Vector rhs_2D_normalized = prod(Inverted_normalizer_RHS_2D,rhs_2D);

            rhs[0] = rhs_2D_normalized[0];   rhs[1] = rhs_2D_normalized[1];


            rhs /= density1;

            noalias(rhs) += bodyf_acc;





            // How it should be, but I can"t because determinant is 0
            //            if (this->Id()==2162){KRATOS_WATCH(normalizer_RHS);}

            //            double determinant = MathUtils<double>::Det3(normalizer_RHS);

            //            if (this->Id()==2162){KRATOS_WATCH(determinant);}
            //            Matrix Inverted_normalizer_RHS = ZeroMatrix(3,3);

            //            MathUtils<double>::InvertMatrix3(normalizer_RHS,Inverted_normalizer_RHS,determinant);

            //            if (this->Id()==2162){KRATOS_WATCH(Inverted_normalizer_RHS);}

            //            //NORMALIZE THE RHS


            //            array_1d<double,3>& bodyf_acc = this->GetGeometry()[0].FastGetSolutionStepValue(BODYFORCE_ACC);

            //            array_1d<double,3>& rhs      = this->GetGeometry()[0].FastGetSolutionStepValue(RHS);

//            noalias(rhs) -= bodyf_acc;

//            rhs *= density1;


            //            rhs = prod(Inverted_normalizer_RHS,rhs); //*******/*///*/*/*

//            rhs /= density1;

//            noalias(rhs) += bodyf_acc;




        }



        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void CopyBoundaryPressures()
    {
        KRATOS_TRY;

        if ( (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) && (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_BOUNDARY) == 0.0) ){

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            double pres_sum=0;
            double counter=0;

            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {

                if ( neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0 ){


                    double neighbour_pressure=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) ;

                    pres_sum += neighbour_pressure;
                    counter+=1;
                }
            }

            if (counter != 0){this->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE)=pres_sum/counter;}

        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void ApplyXSPH()
    {
        KRATOS_TRY;

        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0){

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
            const array_1d<double,3>& vel1 = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);


            const double density1=this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;

            const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;

            array_1d<double,3> velocity_correction = ZeroVector(3);


            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0){
                    const double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

                    const double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;

                    const array_1d<double,3> vel_diff = vel1 - neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

                    array_1d<double,3> rvec = (neighbour_iterator->GetGeometry()[0].Coordinates());
                    noalias(rvec) -= this->GetGeometry()(0)->Coordinates();
                    const double r = norm_2(rvec);

                    //Calculate the Value
                    double Kernel;
                    Kernel = K_dens::ComputeKernel(r,h);

                    const double dens_aver = (density1 + density2)/2.0;

                    velocity_correction += (Kernel*mass/dens_aver) *  vel_diff;


                }


            }


            velocity_correction = 0.5 * velocity_correction ;
            array_1d<double,3>& temp_vel1 = this->GetGeometry()(0)->FastGetSolutionStepValue(XSPH_VELOCITY);
            temp_vel1 = vel1 - velocity_correction;
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    inline void ApplyBoundaryForces()
    {
        KRATOS_TRY;

        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0){

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;

            double D= 5.0* 0.3 * 9.81;

            double spacing = h;

            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0){

                    array_1d<double,3> rvec = this->GetGeometry()(0)->Coordinates();
                    noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());
                    const double r = norm_2(rvec);
                    //Calculate the Value
                    const double ratio = r/spacing;

                    if (ratio <= 0.9){
                        array_1d<double,3>& bound_acc = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(BOUNDARY_ACC);
                        const double SameValue = D* (pow(ratio,12) - pow(ratio,4) ) ;
                        bound_acc[0] += SameValue * ( rvec[0] /(r*r) ) ;
                        bound_acc[1] += SameValue * ( rvec[1] /(r*r) ) ;
                        bound_acc[2] += SameValue * ( rvec[2] /(r*r) ) ;

                    }

                }

            }
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void CatchFreeSurface()
    {
        KRATOS_TRY;
        this->GetGeometry()(0)->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0.0;
        ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);


        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0){

            array_1d<double,3> rvec ;


            double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;
            double aux=0;
            double lower_catcher =0;


            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

                double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;

                noalias(rvec) = this->GetGeometry()(0)->Coordinates();
                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());
                array_1d<double,3> w;
                w = K_dens::ComputeKernelDerivative(rvec,h);
                aux += (mass/density2) * (rvec[0]*w[0] + rvec[1]*w[1] + rvec[2]*w[2]);
                if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) ==1.0 ){
                    lower_catcher=1;
                }
            }

            aux = -aux;
            //KRATOS_WATCH(aux);
            if (aux<1.6){this->GetGeometry()(0)->FastGetSolutionStepValue(IS_FREE_SURFACE) = 1.0;}
            if (lower_catcher==1){this->GetGeometry()(0)->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0.0;}

        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void CalculateIntermediateRHS()
    {
        KRATOS_TRY;

        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0){

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);
            const array_1d<double,3>& vel1 = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);

            double density1=this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
            double v1=this->GetGeometry()(0)->FastGetSolutionStepValue(VISCOSITY) ;

            array_1d<double,3> aux_sum_visc = ZeroVector(3);
            array_1d<double,3> rvec;
            array_1d<double,3> w2;

            const double parameter = 0.01 * h * h;

            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;
                double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
                double v2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(VISCOSITY) ;
                const array_1d<double,3>& vel2 = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

                array_1d<double,3> u = vel1;
                noalias(u) -= vel2 ;


                //                noalias(rvec) = this->GetGeometry()(0)->Coordinates();
                //                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());

                noalias(rvec) = this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_POS);
                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(TEMP_POS));
                double r = norm_2(rvec);


                w2 = K_visc::ComputeKernelDerivative(rvec,h);

                //VISCOSITY


                double same_valueV = 4*mass * ((v1+v2) / ((density1+density2)*(density1+density2)) * (( w2[0]*rvec[0] +w2[1]*rvec[1]+w2[2]*rvec[2]) / (r*r + parameter)) );
                noalias(aux_sum_visc ) += same_valueV*u;


            }


            //// ASSING THE VISCOUS ACCELERATION
            array_1d<double,3>& visc_acc = this->GetGeometry()[0].FastGetSolutionStepValue(VISCOUS_ACC);

            // ADD EVERYTHING UP AND ASSIGN TO RHS VARIABLE
            array_1d<double,3>& bodyf_acc = this->GetGeometry()[0].FastGetSolutionStepValue(BODYFORCE_ACC);
            array_1d<double,3>& temp_rhs = this->GetGeometry()[0].FastGetSolutionStepValue(TEMP_RHS);
            noalias(temp_rhs) = bodyf_acc;
            noalias(temp_rhs) += visc_acc;
            //            temp_rhs /= density1;

        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        double time_step = rCurrentProcessInfo[DELTA_TIME];
        const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;
        double& density1=this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
        //        double ref_density1=this->GetGeometry()(0)->GetValue(DENSITY) ;
        const double parameter = 0.01 * h * h;

        ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);




        //        KRATOS_WATCH(this->Id());

        //        KRATOS_WATCH(neighbours.size());



        //The right hand side vector

        //const array_1d<double,3>& temp_vel1= this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_VEL);

        //        density1 = 0;

        //        array_1d<double,3> rvec;

        //        for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
        //        {
        //            double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

        //            noalias(rvec) = this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_POS);
        //            noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(TEMP_POS));
        //            double r = norm_2(rvec);

        //            //Calculate the Value
        //            double Kernel = K_dens::ComputeKernel(r,h);


        //            density1 += mass * Kernel;
        //        }


        //        rRightHandSideVector(0)= (ref_density1-density1) / (ref_density1 * time_step * time_step);

        //The right hand side vector

        array_1d<double,3>& temp_vel1= this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_VEL);


        double div_of_velocity = 0;


        array_1d<double,3> rvec;

        for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
        {
            if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 ){

                double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;
                double& density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;


                array_1d<double,3>& temp_vel2 = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(TEMP_VEL);

                temp_vel1 = (1/(density1 * density1) ) * temp_vel1;
                temp_vel2 = (1/(density2 * density2) ) * temp_vel2;

                array_1d<double,3> u;
                u = temp_vel1 + temp_vel2 ;



                //                noalias(rvec) = this->GetGeometry()(0)->Coordinates();
                //                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());

                noalias(rvec) = this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_POS);
                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(TEMP_POS));

                array_1d<double,3> w;
                w = K_dens::ComputeKernelDerivative(rvec,h);

                div_of_velocity += mass * (u[0]*w[0] + u[1]*w[1] + u[2]*w[2]);
            }
        }

        div_of_velocity = density1 * div_of_velocity;

        rRightHandSideVector(0)= div_of_velocity / time_step;

        //        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) ==1.0  ){
        //            KRATOS_WATCH(rRightHandSideVector);
        //        }





        double aux=0;
        for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
        {
            if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(OUT_OF_SYSTEM) == 0.0 ){
                aux++;
            }
        }




        double entry=0;
        double diagonal_entry=0;
        int diagonal_entry_location=0;
        int counter=0;

        Vector pressures(aux);

        rLeftHandSideMatrix = ZeroMatrix(2,aux);

        //        KRATOS_WATCH(aux);


        for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
        {
            if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(OUT_OF_SYSTEM) == 0.0 ){


                if ( this->Id() != neighbour_iterator->Id() ){

                    double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

                    double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;

                    noalias(rvec) = this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_POS);
                    noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(TEMP_POS));

                    double r = norm_2(rvec);

                    array_1d<double,3> w;
                    w = K_pres::ComputeKernelDerivative(rvec,h);

                    double density_sum = (density1+density2);
                    entry = ( ( mass / density2 ) * ( 4/(density_sum)) ) * ( (rvec[0]*w[0] + rvec[1]*w[1] + rvec[2]*w[2]) / (r*r + parameter) );

                    diagonal_entry += entry;

                    rLeftHandSideMatrix(0,counter) =entry;
                    rLeftHandSideMatrix(1,counter) =neighbour_iterator->Id() - 1;

                    pressures(counter) = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE);
                    counter++;
                }
                else{
                    diagonal_entry_location=counter;
                    pressures(counter) = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE);
                    counter++;

                }
            }
        }

        //        KRATOS_WATCH(counter);
        //        KRATOS_WATCH(pressures);




        rLeftHandSideMatrix(0,diagonal_entry_location) = -diagonal_entry;

        //        KRATOS_WATCH(diagonal_entry_location);


        rLeftHandSideMatrix(1,diagonal_entry_location) = this->Id() - 1;




        Matrix rLeftHandSideMatrix_FirstLine( 1, rLeftHandSideMatrix.size2() );

        for (unsigned int i =0;i<rLeftHandSideMatrix.size2();i++){
            rLeftHandSideMatrix_FirstLine(0,i) = rLeftHandSideMatrix(0,i);
        }


        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix_FirstLine,pressures);

        //        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) ==1.0  ){
        //            KRATOS_WATCH(rRightHandSideVector);
        //        }

        //        KRATOS_WATCH(rRightHandSideVector);






        KRATOS_CATCH("");

    }

    //************************************************************************************
    //************************************************************************************
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY;

        if (this->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0){

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);

            double density1=this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
            double p1=this->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) ;

            array_1d<double,3> aux_sum = ZeroVector(3);
            array_1d<double,3> rvec;
            array_1d<double,3> w1;
            if (this->Id() == 1796){
                KRATOS_WATCH(p1);
            }

            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;
                double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
                double p2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) ;


                //                noalias(rvec) = this->GetGeometry()(0)->Coordinates();
                //                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());

                noalias(rvec) = this->GetGeometry()(0)->FastGetSolutionStepValue(TEMP_POS);
                noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(TEMP_POS));

                w1 = K_visc::ComputeKernelDerivative(rvec,h);

                //PRESSURE
                double same_valueP = mass * ( p1/(density1*density1) + p2/(density2*density2));

                noalias(aux_sum ) += same_valueP*w1;




            }


            //// ASSING THE PRESSURE ACCELERATION
            array_1d<double,3>& pres_acc = this->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_ACC);
            noalias(pres_acc) =-aux_sum;

            // USE TEMP_RHS to reach final RHS
            array_1d<double,3>& temp_rhs = this->GetGeometry()[0].FastGetSolutionStepValue(TEMP_RHS);
            array_1d<double,3>& rhs = this->GetGeometry()[0].FastGetSolutionStepValue(RHS);

            noalias(rhs) = temp_rhs + pres_acc;





        }




        KRATOS_CATCH("");
    }


    //************************************************************************************
    //************************************************************************************

    void Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY;

        double soundspeed = rCurrentProcessInfo[SOUND_VELOCITY];

        if ( rVariable == DUMMY_NORMALIZE_RHS)
        {
            this->NormalizeRightHandSide();
        }

        else if ( rVariable == DUMMY_CATCH_FREESURFACE)
        {
            this->CatchFreeSurface();
        }

        else if ( rVariable == DUMMY_INTERMEDIATE_RHS)
        {
            this->CalculateIntermediateRHS();
        }

        else if ( rVariable == DUMMY_APPLY_XSPH)
        {
            this->ApplyXSPH();
        }

        else if ( rVariable == DUMMY_BOUNDARY_PRESSURES)
        {
            this->CopyBoundaryPressures();
        }

        else if ( rVariable == DENSITY )
        {


//            Collacative Method

//            const  double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;

//            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

//            double& density = this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);
//            density = 0.0;

//            array_1d<double,3> rvec;

//            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
//            {
//                const double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

//                noalias(rvec) = (neighbour_iterator->GetGeometry()[0].Coordinates());
//                noalias(rvec) -= this->GetGeometry()(0)->Coordinates();
//                const double r = norm_2(rvec);

//                //Calculate the Value
//                double Kernel = K_dens::ComputeKernel(r,h);


//                density += mass * Kernel;
//            }



//            const double ref_dens=this->GetGeometry()(0)->GetValue(DENSITY);

//            this->GetGeometry()(0)->FastGetSolutionStepValue(DENS_VARIATION) = ((density-ref_dens) / ref_dens);




            //Continuum Method

            double time_step = rCurrentProcessInfo[DELTA_TIME];


            double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;

            double& density1 = this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);

            ParticleWeakVectorType& neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

            const array_1d<double,3>& vel1 = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);


            double div_of_velocity = 0;

            array_1d<double,3> rvec;

            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                if (neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 ){

                    const double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

                    const array_1d<double,3>& vel2 = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

                    array_1d<double,3> u = vel1;
                    noalias(u) -= vel2 ;

                    noalias(rvec) = this->GetGeometry()(0)->Coordinates();
                    noalias(rvec) -= (neighbour_iterator->GetGeometry()[0].Coordinates());

                    array_1d<double,3> w;
                    w = K_dens::ComputeKernelDerivative(rvec,h);

                    div_of_velocity += mass * (u[0]*w[0] + u[1]*w[1] + u[2]*w[2]);
                }
            }

            div_of_velocity = -div_of_velocity;

            //Update the density

            density1 -= div_of_velocity * time_step;

            const double ref_dens=this->GetGeometry()(0)->GetValue(DENSITY);

            this->GetGeometry()(0)->FastGetSolutionStepValue(DENS_VARIATION) = ((density1-ref_dens) / ref_dens);
            if (this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) < 0.0){
                KRATOS_WATCH(div_of_velocity);
                KRATOS_WATCH(this->GetGeometry()(0)->Id());
                KRATOS_WATCH(this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY));
                KRATOS_WATCH(this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY,1));
            }


        }

        else if ( rVariable == PRESSURE )
        {
            //State Equation 1
            //            this->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) =
            //                    soundspeed * soundspeed * (this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) - this->GetGeometry()(0)->GetValue(DENSITY));

            // ANTONIA STATE EQUATION

            //                        this->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) +=
            //                                soundspeed * soundspeed * (this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) - this->GetGeometry()(0)->GetValue(DENSITY));

            //State Equation 2
            const double density = this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);
            const double ref_dens=this->GetGeometry()(0)->GetValue(DENSITY);
            const double ratio = density/ref_dens;
            this->GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE) =
                    (ref_dens*soundspeed*soundspeed/7.0) * (pow(ratio,7) - 1);
        }

        else if ( rVariable == DELTA_TIME )
        {

            const double radius = this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);
            const double viscosity = this->GetGeometry()(0)->FastGetSolutionStepValue(VISCOSITY);
            const double density = this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);
            const double ref_density = this->GetGeometry()(0)->GetValue(DENSITY);

            const array_1d<double, 3 > particle_acc= this->GetGeometry()[0].FastGetSolutionStepValue(RHS);
            const double particle_accNorm=norm_2(particle_acc);
            const array_1d<double, 3 > velocity= this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const double velocity_Norm=norm_2(velocity);



            //The CFL condition

            const double CFL=0.25*radius/( soundspeed + velocity_Norm );

            // The MFC ( Mass force condition )


            double MFC = CFL;
            if (particle_accNorm != 0.0){

                MFC= 0.25 * sqrt(radius/fabs(particle_accNorm));
            }


            // The VFC ( Viscous force condition )

            const double VFC=0.125 * (radius*radius/ (viscosity/ref_density) );

            //Put the smallest as Output
            Output = fmin(CFL,fmin(MFC,VFC));

            if (Output < 0.0){
                KRATOS_WATCH(CFL);
                KRATOS_WATCH(density);
                KRATOS_WATCH(this->Id());
                KRATOS_WATCH(this->GetGeometry()(0)->Id());

            }

        }

        else if ( rVariable == DELTA_TIME_ISPH )
        {

            double radius = this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);
            double viscosity = this->GetGeometry()(0)->FastGetSolutionStepValue(VISCOSITY);
            double density = this->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);

            array_1d<double, 3 > velocity= this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            double velocity_norm = norm_2(velocity);
            double CFL=1000;
            //The CFL condition
            if (velocity_norm != 0.0){

                CFL=0.25*radius/velocity_norm;
            }
            // The MFC ( Mass force condition )

            array_1d<double, 3 > particle_acc= this->GetGeometry()[0].FastGetSolutionStepValue(RHS);
            double particle_accNorm=norm_2(particle_acc);
            double MFC = CFL;
            if (particle_accNorm != 0.0){

                MFC= 0.25 * sqrt(radius/fabs(particle_accNorm));
            }


            // The VFC ( Viscous force condition )

            double VFC=0.125 * (radius*radius/ (viscosity/density) );

            //Put the smallest as Output
            Output = fmin(CFL,fmin(MFC,VFC));

        }

        KRATOS_CATCH("");

    }


protected:

private:
    friend class Serializer;

    SPHparticle() : Element()
    {
    }


}; // Class SPHparticle
}  // namespace Kratos.

#endif // KRATOS_SPHPARTICLE_H_INCLUDED  defined
