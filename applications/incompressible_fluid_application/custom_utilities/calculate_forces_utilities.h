/*
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_CALCULATE_FORCE_UTILITIES_INCLUDED )
#define  KRATOS_CALCULATE_FORCE_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
class CalculateForcesUtils
{
public:

    //**********************************************************************************************
    //**********************************************************************************************
    //calculates the mesh motion spatial acceleration, taking in account both the spatial and the time dependent
    //componenents acc = (v-vm)*dv/dx
    void CalculateForces3D(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
        array_1d<double,6> sigma;
        array_1d<double,4> N;
        array_1d<double,3> aux;
        array_1d<double,3> aux_stab_term;
        array_1d<double,3> aux_conv_proj;
        array_1d<double,3> vconv;
        const array_1d<double,3> zero3 = ZeroVector(3);
        boost::numeric::ublas::bounded_matrix<double,3,3> dv_dx;
        array_1d<double,4> ms_u_DN;
        const int nnodes = 4;
        int dim = 3;

        const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
        const Vector& BDFcoeffs = CurrentProcessInfo[BDF_COEFFICIENTS];

        //reset the acceleration on the nodes
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            noalias(in->FastGetSolutionStepValue(FORCE)) = zero3;
        }

        //compute the projection (first step
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node<3> >& geom = i->GetGeometry();
            double Volume;

            GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

            //calculation of the average density and viscosity
            double density = 0.00;
            double viscosity = 0.00;
            for(int I = 0; I<nnodes; I++)
            {
                density += geom[I].FastGetSolutionStepValue(DENSITY);
                viscosity += geom[I].FastGetSolutionStepValue(VISCOSITY);
            }
            density *= 0.25;
            viscosity *= 0.25;
            double coeff = density * viscosity;

            //constructing the velocity gradient
            array_1d<double,3>& v = geom[0].FastGetSolutionStepValue(VELOCITY);

            noalias(vconv) = v;
            noalias(vconv) -= geom[0].FastGetSolutionStepValue(MESH_VELOCITY);

            for(int i = 0; i<dim; i++)
                for(int j = 0; j<dim; j++)
                    dv_dx(i,j) = DN_DX(0,j) * v[i];
            for(int I = 1; I<nnodes; I++)
            {
                array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY);

                noalias(vconv) += v;
                noalias(vconv) -= geom[0].FastGetSolutionStepValue(MESH_VELOCITY);

                for(int i = 0; i<dim; i++)
                    for(int j = 0; j<dim; j++)
                        dv_dx(i,j) += DN_DX(I,j) * v[i];
            }
            vconv *= 0.25;

            //calculating avg pressure
            double pavg = geom[0].FastGetSolutionStepValue(PRESSURE);
            for(int I = 1; I<nnodes; I++)
            {
                pavg += geom[I].FastGetSolutionStepValue(PRESSURE);
            }
            pavg *= 0.25;

            //calculating the stress
            sigma[0] = -dv_dx(0,0) * coeff *2.0 + pavg;
            sigma[1] = -dv_dx(1,1) * coeff *2.0 + pavg;
            sigma[2] = -dv_dx(2,2) * coeff *2.0 + pavg;
            sigma[3] = -(dv_dx(0,1) + dv_dx(1,0) )* coeff;
            sigma[4] = -(dv_dx(1,2) + dv_dx(2,1) )* coeff;
            sigma[5] = -(dv_dx(2,0) + dv_dx(0,2) )* coeff;

            //calculating the stabilization parameter
            double c1 = 4.00;
            double c2 = 2.00;
            double h = pow(6.00*Volume,0.3333333);
            double norm_u = vconv[0]*vconv[0] + vconv[1]*vconv[1] + vconv[2]*vconv[2];
            norm_u = sqrt(norm_u);
            double tau = 1.00 / ( c1*viscosity/(h*h) + c2*norm_u/h );

            noalias(ms_u_DN) = prod(DN_DX , vconv);

            noalias(aux_stab_term) = ms_u_DN[0]*geom[0].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[1]*geom[1].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[2]*geom[2].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[3]*geom[3].FastGetSolutionStepValue(VELOCITY);

            noalias(aux_conv_proj) = geom[0].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[1].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[2].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[3].FastGetSolutionStepValue(CONV_PROJ);
            aux_conv_proj *= 0.25;

            //adding the force to the nodes
            for(int I = 0; I<nnodes; I++)
            {
                array_1d<double,3>& force = geom[I].FastGetSolutionStepValue(FORCE);

                //adding the viscous and pressure contribution
                force[0] += Volume * ( DN_DX(I,0)*( sigma[0] )
                                       + DN_DX(I,1)*sigma[3] + DN_DX(I,2)*sigma[4]	);

                force[1] += Volume * ( DN_DX(I,1)*( sigma[1] )
                                       + DN_DX(I,0)*sigma[3] + DN_DX(I,2)*sigma[5]	);

                force[2] += Volume * ( DN_DX(I,2)*( sigma[2] )
                                       + DN_DX(I,0)*sigma[4] + DN_DX(I,1)*sigma[5]	);

                //calculating acceleration (direct)
                if(BDFcoeffs.size() != 0)
                {
                    noalias(aux) =  zero3;
                    for(unsigned int step = 0; step<BDFcoeffs.size(); step++)
                    {
                        const array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY,step);
                        noalias(aux) += BDFcoeffs[step]*v;
                    }
                    //acceleration (spatial part)
                    noalias(aux) += aux_stab_term;
                }
                else
                    noalias(aux) = zero3;

                noalias(force) -= (0.25 * Volume * density) * aux;

            }

            //adding stabilization terms
            noalias(aux) = aux_conv_proj - aux_stab_term;
            aux *= tau;

            for(int I = 0; I<nnodes; I++)
            {
                array_1d<double,3>& force = geom[I].FastGetSolutionStepValue(FORCE);
                noalias(force) += 	(ms_u_DN[I] * Volume)* ( aux);
            }
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    //calculates the mesh motion spatial acceleration, taking in account both the spatial and the time dependent
    //componenents acc = (v-vm)*dv/dx
    void CalculateForces2D(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
        array_1d<double,3> sigma;
        array_1d<double,3> N;
        array_1d<double,3> aux;
        array_1d<double,3> aux_stab_term;
        array_1d<double,3> aux_conv_proj;
        array_1d<double,2> vconv;
        const array_1d<double,3> zero3 = ZeroVector(3);
        boost::numeric::ublas::bounded_matrix<double,2,2> dv_dx;
        array_1d<double,3> ms_u_DN;
        const int nnodes = 3;
        int dim = 2;

        const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
        const Vector& BDFcoeffs = CurrentProcessInfo[BDF_COEFFICIENTS];

        //reset the acceleration on the nodes
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            noalias(in->FastGetSolutionStepValue(FORCE)) = zero3;
        }

        //compute the projection (first step
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node<3> >& geom = i->GetGeometry();
            double Volume;

            GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

            //calculation of the average density and viscosity
            double density = 0.00;
            double viscosity = 0.00;
            for(int I = 0; I<nnodes; I++)
            {
                density += geom[I].FastGetSolutionStepValue(DENSITY);
                viscosity += geom[I].FastGetSolutionStepValue(VISCOSITY);
            }
            density *= 0.3333333333333333333333;
            viscosity *= 0.3333333333333333333333;
            double coeff = density * viscosity;

            //constructing the velocity gradient
            array_1d<double,3>& v = geom[0].FastGetSolutionStepValue(VELOCITY);
            array_1d<double,3>& mesh_vel = geom[0].FastGetSolutionStepValue(MESH_VELOCITY);

            vconv[0] = v[0]-mesh_vel[0];
            vconv[1] = v[1]-mesh_vel[1];

            for(int i = 0; i<dim; i++)
                for(int j = 0; j<dim; j++)
                    dv_dx(i,j) = DN_DX(0,j) * v[i];
            for(int I = 1; I<nnodes; I++)
            {
                array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY);
                array_1d<double,3>& mesh_vel = geom[I].FastGetSolutionStepValue(MESH_VELOCITY);

                vconv[0] += v[0]-mesh_vel[0];
                vconv[1] += v[1]-mesh_vel[1];

                for(int i = 0; i<dim; i++)
                    for(int j = 0; j<dim; j++)
                        dv_dx(i,j) += DN_DX(I,j) * v[i];
            }
            vconv *= 0.3333333333333333333333;

            //calculating avg pressure
            double pavg = geom[0].FastGetSolutionStepValue(PRESSURE);
            for(int I = 1; I<nnodes; I++)
            {
                pavg += geom[I].FastGetSolutionStepValue(PRESSURE);
            }
            pavg *= 0.3333333333333333333333;

            //calculating the stress
            sigma[0] = -dv_dx(0,0) * coeff *2.0 + pavg;
            sigma[1] = -dv_dx(1,1) * coeff *2.0 + pavg;
            sigma[2] = -(dv_dx(0,1) + dv_dx(1,0) )* coeff;

            //calculating the stabilization parameter
            double c1 = 4.00;
            double c2 = 2.00;
            double h = sqrt(2.00*Volume);
            double norm_u = vconv[0]*vconv[0] + vconv[1]*vconv[1] + vconv[2]*vconv[2];
            norm_u = sqrt(norm_u);
            double tau = 1.00 / ( c1*viscosity/(h*h) + c2*norm_u/h );

            noalias(ms_u_DN) = prod(DN_DX , vconv);

            noalias(aux_stab_term) = ms_u_DN[0]*geom[0].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[1]*geom[1].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[2]*geom[2].FastGetSolutionStepValue(VELOCITY);

            noalias(aux_conv_proj) = geom[0].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[1].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[2].FastGetSolutionStepValue(CONV_PROJ);
            aux_conv_proj *= 0.3333333333333333333333;

            //adding the force to the nodes
            for(int I = 0; I<nnodes; I++)
            {
                array_1d<double,3>& force = geom[I].FastGetSolutionStepValue(FORCE);

                //adding the viscous and pressure contribution
                force[0] += Volume * ( DN_DX(I,0)*( sigma[0] )
                                       + DN_DX(I,1)*sigma[2]	);

                force[1] += Volume * ( DN_DX(I,1)*( sigma[1] )
                                       + DN_DX(I,0)*sigma[2]	);


                //calculating acceleration (direct)
                if(BDFcoeffs.size() != 0)
                {
                    noalias(aux) =  zero3;
                    for(unsigned int step = 0; step<BDFcoeffs.size(); step++)
                    {
                        const array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY,step);
                        noalias(aux) += BDFcoeffs[step]*v;
                    }
                    //acceleration (spatial part)
                    noalias(aux) += aux_stab_term;
                }
                else
                    noalias(aux) = zero3;

                noalias(force) -= (0.3333333333333333333333 * Volume * density) * aux;

            }

            //adding stabilization terms
            noalias(aux) = aux_conv_proj - aux_stab_term;
            aux *= tau;

            for(int I = 0; I<nnodes; I++)
            {
                array_1d<double,3>& force = geom[I].FastGetSolutionStepValue(FORCE);
                noalias(force) += 	(ms_u_DN[I] * Volume)* ( aux);
            }
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    //calculates the mesh motion spatial acceleration, taking in account both the spatial and the time dependent
    //componenents acc = (v-vm)*dv/dx
    void CalculatePressureForces2D(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
        array_1d<double,3> sigma;
        array_1d<double,3> N;
        array_1d<double,3> aux;
        array_1d<double,3> aux_stab_term;
        array_1d<double,3> aux_conv_proj;
        array_1d<double,2> vconv;
        const array_1d<double,3> zero3 = ZeroVector(3);
        boost::numeric::ublas::bounded_matrix<double,2,2> dv_dx;
        array_1d<double,3> ms_u_DN;
        const int nnodes = 3;
        int dim = 2;

        const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
        const Vector& BDFcoeffs = CurrentProcessInfo[BDF_COEFFICIENTS];

        //reset the acceleration on the nodes
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            noalias(in->FastGetSolutionStepValue(FORCE)) = zero3;
        }

        //compute the projection (first step
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node<3> >& geom = i->GetGeometry();
            double Volume;

            GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

            //calculation of the average density and viscosity
            double density = 0.00;
            double viscosity = 0.00;
            for(int I = 0; I<nnodes; I++)
            {
                density += geom[I].FastGetSolutionStepValue(DENSITY);
                viscosity += geom[I].FastGetSolutionStepValue(VISCOSITY);
            }
            density *= 0.3333333333333333333333;
            viscosity *= 0.3333333333333333333333;
            //double coeff = 0.0;

            //constructing the velocity gradient
            array_1d<double,3>& v = geom[0].FastGetSolutionStepValue(VELOCITY);
            array_1d<double,3>& mesh_vel = geom[0].FastGetSolutionStepValue(MESH_VELOCITY);

            vconv[0] = v[0]-mesh_vel[0];
            vconv[1] = v[1]-mesh_vel[1];

            for(int i = 0; i<dim; i++)
                for(int j = 0; j<dim; j++)
                    dv_dx(i,j) = DN_DX(0,j) * v[i];
            for(int I = 1; I<nnodes; I++)
            {
                array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY);
                array_1d<double,3>& mesh_vel = geom[I].FastGetSolutionStepValue(MESH_VELOCITY);

                vconv[0] += v[0]-mesh_vel[0];
                vconv[1] += v[1]-mesh_vel[1];

                for(int i = 0; i<dim; i++)
                    for(int j = 0; j<dim; j++)
                        dv_dx(i,j) += DN_DX(I,j) * v[i];
            }
            vconv *= 0.3333333333333333333333;

            //calculating avg pressure
            double pavg = geom[0].FastGetSolutionStepValue(PRESSURE);
            for(int I = 1; I<nnodes; I++)
            {
                pavg += geom[I].FastGetSolutionStepValue(PRESSURE);
            }
            pavg *= 0.3333333333333333333333;

            //calculating the stress
            sigma[0] =  pavg;
            sigma[1] =  pavg;
            sigma[2] = 0.0;

            //calculating the stabilization parameter
            double c1 = 4.00;
            double c2 = 2.00;
            double h = sqrt(2.00*Volume);
            double norm_u = vconv[0]*vconv[0] + vconv[1]*vconv[1] + vconv[2]*vconv[2];
            norm_u = sqrt(norm_u);
            double tau = 1.00 / ( c1*viscosity/(h*h) + c2*norm_u/h );

            noalias(ms_u_DN) = prod(DN_DX , vconv);

            noalias(aux_stab_term) = ms_u_DN[0]*geom[0].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[1]*geom[1].FastGetSolutionStepValue(VELOCITY);
            noalias(aux_stab_term) += ms_u_DN[2]*geom[2].FastGetSolutionStepValue(VELOCITY);

            noalias(aux_conv_proj) = geom[0].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[1].FastGetSolutionStepValue(CONV_PROJ);
            noalias(aux_conv_proj) += geom[2].FastGetSolutionStepValue(CONV_PROJ);
            aux_conv_proj *= 0.3333333333333333333333;

            //adding the force to the nodes
            for(int I = 0; I<nnodes; I++)
            {
                array_1d<double,3>& force = geom[I].FastGetSolutionStepValue(FORCE);

                //adding the viscous and pressure contribution
                force[0] += Volume * ( DN_DX(I,0)*( sigma[0] )
                                       + DN_DX(I,1)*sigma[2]	);

                force[1] += Volume * ( DN_DX(I,1)*( sigma[1] )
                                       + DN_DX(I,0)*sigma[2]	);


                //calculating acceleration (direct)
                if(BDFcoeffs.size() != 0)
                {
                    noalias(aux) =  zero3;
                    for(unsigned int step = 0; step<BDFcoeffs.size(); step++)
                    {
                        const array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY,step);
                        noalias(aux) += BDFcoeffs[step]*v;
                    }
                    //acceleration (spatial part)
                    noalias(aux) += aux_stab_term;
                }
                else
                    noalias(aux) = zero3;

                noalias(force) -= (0.3333333333333333333333 * Volume * density) * aux;

            }

            //adding stabilization terms
            noalias(aux) = aux_conv_proj - aux_stab_term;
            aux *= tau;

            for(int I = 0; I<nnodes; I++)
            {
                array_1d<double,3>& force = geom[I].FastGetSolutionStepValue(FORCE);
                noalias(force) += 	(ms_u_DN[I] * Volume)* ( aux);
            }
        }

        KRATOS_CATCH("")
    }

private:


};

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_FORCE_UTILITIES_INCLUDED  defined


