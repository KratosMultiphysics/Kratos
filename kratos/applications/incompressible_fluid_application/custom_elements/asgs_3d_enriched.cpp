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
//   Last modified by:    $Author: kazem, pavel $
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//

//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/asgs_3d_enriched.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h"
#include "utilities/split_tetrahedra.h"
#include "utilities/enrichment_utilities.h"

namespace Kratos
{



//************************************************************************************
//************************************************************************************

ASGS3D_ENR::ASGS3D_ENR(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ASGS3D_ENR::ASGS3D_ENR(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{

}

Element::Pointer ASGS3D_ENR::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new ASGS3D_ENR(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ASGS3D_ENR::~ASGS3D_ENR()
{
}

//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int nodes_number = 4;
    int dim = 3;
    unsigned int matsize = nodes_number * (dim + 1);

    if (rLeftHandSideMatrix.size1() != matsize)
        rLeftHandSideMatrix.resize(matsize, matsize,false); //false says not to preserve existing storage!!

    if (rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize,false); //false says not to preserve existing storage!!


    noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
    noalias(rRightHandSideVector) = ZeroVector(matsize);

    ///////////////////////////////////////////////////////////////PAVEL//////////////////////////////////////////////////////
    Matrix Points(4,3);
    double Volume;								   //
    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;		   //
    array_1d<double,4> N;							   //
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);	   //
    array_1d<double,4> Distances;

    Vector PartitionVolumes(6);
    Vector PartitionSigns(6);
    Matrix Partition_Shape_Functions(6,4);
    //one shape function per Gauss points
    Matrix Enriched_Shape_Functions(6,1);
    //vector of vectors (6 rows of matrices 1*3) (vector of 3 = matrix 1*3)
    std::vector<Matrix> Partition_Gradients(6, Matrix(1,3));
    /////////////////////////////////////////////////////////////////////////////
    for (int i=0; i<nodes_number; i++)
    {
        Distances[i]=GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        // and now we store the verices of the element
        Points(i,0)=GetGeometry()[i].X();
        Points(i,1)=GetGeometry()[i].Y();
        Points(i,2)=GetGeometry()[i].Z();
    }

    array_1d<double,6> edge_areas;
    unsigned int n_subdivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions<MatrixType,VectorType,boost::numeric::ublas::bounded_matrix<double,4,3> >(Points, DN_DX, Distances, PartitionVolumes, 																	
                                Partition_Shape_Functions,PartitionSigns, Partition_Gradients, Enriched_Shape_Functions,edge_areas );
    array_1d<double,4> N_at_igauss = ZeroVector(4);

    //this loops is in progress
    for(unsigned int igauss=0; igauss<n_subdivisions; igauss++)
    {
        double density;

        N_at_igauss[0] = Partition_Shape_Functions(igauss,0);
        N_at_igauss[1] = Partition_Shape_Functions(igauss,1);
        N_at_igauss[2] = Partition_Shape_Functions(igauss,2);
        N_at_igauss[3] = Partition_Shape_Functions(igauss,3);

        Volume=PartitionVolumes[igauss];

        EvaluateAtGaussPoint(density, DENSITY, N_at_igauss);
        AddBodyForceAndMomentum(rRightHandSideVector, N_at_igauss, Volume, density);
    }


    KRATOS_CATCH("")
}
//***********************************************************************************++
//**************************************************************************************++

void ASGS3D_ENR::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MatrixType temp = Matrix();
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
//checked
void ASGS3D_ENR::CalculateMassContribution(MatrixType& K, const array_1d<double, 4 > & N_at_igauss, const double volume, const double density)
{
    KRATOS_TRY
    boost::numeric::ublas::bounded_matrix<double, 4, 4 > ConsistentMass = ZeroMatrix(4, 4);
    boost::numeric::ublas::bounded_matrix<double, 16, 16 > ConsistentMassFull = ZeroMatrix(16, 16);

    ConsistentMass=volume*density*outer_prod(N_at_igauss, trans(N_at_igauss));

    //KRATOS_WATCH(ConsistentMass);
    //double lump_mass_fac = volume * 0.25;

    //int nodes_number = 4;
    //int dof = 3;

    for (int i = 0; i < 4; i++)
    {
        int row = 4* i;
        for (int j = 0; j < 4; j++)
	    {
	    int column=4*j;
            ConsistentMassFull(row,   column) += ConsistentMass(i,j);   ConsistentMassFull(row,   column+1) += ConsistentMass(i,j);   ConsistentMassFull(row,   column+2) += ConsistentMass(i,j);
            ConsistentMassFull(row+1, column) += ConsistentMass(i,j);   ConsistentMassFull(row+1, column+1) += ConsistentMass(i,j);   ConsistentMassFull(row+1, column+2) += ConsistentMass(i,j);
            ConsistentMassFull(row+2, column) += ConsistentMass(i,j);   ConsistentMassFull(row+1, column+1) += ConsistentMass(i,j);   ConsistentMassFull(row+2, column+2) += ConsistentMass(i,j);

	    }
    }
    //KRATOS_WATCH(ConsistentMassFull)
    
    double lumped_term=0.0;
    for (int i = 0; i < 16; i++)
    {
       lumped_term=0.0;
       for (int j=0;j<16;j++)
  	  {
          lumped_term+=ConsistentMassFull(i,j);		
	  }
       for (int j=0;j<16;j++)
	  {
	     if(i==j)
	          K(i,j)+=lumped_term;		
	  }
    }
    
    /*
    for (int nd = 0; nd < nodes_number; nd++)
    {
        int row = nd * (dof + 1);
        for (int jj = 0; jj < dof; jj++)
            K(row + jj, row + jj) += density / 1.0 * lump_mass_fac;
    }
    */
    KRATOS_CATCH("")

}


//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int nodes_number = GetGeometry().size();
    unsigned int MatSize = (dimension + 1) * nodes_number;
    if (rMassMatrix.size1() != MatSize)
        rMassMatrix.resize(MatSize, MatSize, false);

    rMassMatrix = ZeroMatrix(MatSize, MatSize);
    double delta_t = rCurrentProcessInfo[DELTA_TIME];
    ///////////////////////////////////////////////////////////////PAVEL//////////////////////////////////////////////////////
    Matrix Points(4,3);
    double Volume;								   //
    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;		   //
    array_1d<double,4> N;							   //
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);	   //
    array_1d<double,4> Distances;

    Vector PartitionVolumes(6);
    Vector PartitionSigns(6);
    Matrix Partition_Shape_Functions(6,4);
    //one shape function per Gauss points
    Matrix Enriched_Shape_Functions(6,1);
    //vector of vectors (6 rows of matrices 1*3) (vector of 3 = matrix 1*3)
    std::vector<Matrix> Partition_Gradients(6, Matrix(1,3));
    /////////////////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i<nodes_number; i++)
    {
        Distances[i]=GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        // and now we store the verices of the element
        Points(i,0)=GetGeometry()[i].X();
        Points(i,1)=GetGeometry()[i].Y();
        Points(i,2)=GetGeometry()[i].Z();
    }

    array_1d<double,6> edge_areas;
    unsigned int n_subdivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions<MatrixType,VectorType,boost::numeric::ublas::bounded_matrix<double,4,3> >(Points, DN_DX, Distances, PartitionVolumes, 																	
                                Partition_Shape_Functions,PartitionSigns, Partition_Gradients, Enriched_Shape_Functions,edge_areas );
    array_1d<double,4> N_at_igauss = ZeroVector(4);


    boost::numeric::ublas::bounded_matrix<double, 16, 16 > ConsistentMass = ZeroMatrix(16, 16);
    //this loops is in progress
    for(unsigned int igauss=0; igauss<n_subdivisions; igauss++)
    {
        double density;
        double viscosity;
        double tauone;


        N_at_igauss[0] = Partition_Shape_Functions(igauss,0);
        N_at_igauss[1] = Partition_Shape_Functions(igauss,1);
        N_at_igauss[2] = Partition_Shape_Functions(igauss,2);
        N_at_igauss[3] = Partition_Shape_Functions(igauss,3);

        Volume=PartitionVolumes[igauss];

        EvaluateAtGaussPoint(density, DENSITY, N_at_igauss);
        EvaluateAtGaussPoint(viscosity, VISCOSITY, N_at_igauss);
        CalculateTau(N_at_igauss,tauone, delta_t, Volume, density, viscosity);

        //the next three functions are correct
        CalculateMassContribution(rMassMatrix, N_at_igauss, Volume, density);
        //add stablilization terms due to advective term (a)grad(V) * ro*Acce
        CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N_at_igauss, tauone, Volume, density);
        //add stablilization terms due to grad term grad(q) * ro*Acce
        CalculateGradMassStblTerms(rMassMatrix, DN_DX, N_at_igauss, tauone, Volume, density);
    }

    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    int nodes_number = 4;
    int dim = 3;
    unsigned int matsize = nodes_number * (dim + 1);

    if (rDampingMatrix.size1() != matsize)
        rDampingMatrix.resize(matsize, matsize, false); //false says not to preserve existing storage!!


    noalias(rDampingMatrix) = ZeroMatrix(matsize, matsize);

    double delta_t = rCurrentProcessInfo[DELTA_TIME];


    ///////////////////////////////////////////////////////////////PAVEL//////////////////////////////////////////////////////
    Matrix Points(4,3);
    double Volume;								   //
    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;		   //
    array_1d<double,4> N;							   //
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);	   //
    array_1d<double,4> Distances;

    Vector PartitionVolumes(6);
    Vector PartitionSigns(6);
    Matrix Partition_Shape_Functions(6,4);
    //one shape function per Gauss points
    Matrix Enriched_Shape_Functions(6,1);
    //vector of vectors (6 rows of matrices 1*3) (vector of 3 = matrix 1*3)
    std::vector<Matrix> Partition_Gradients(6, Matrix(1,3));
    /////////////////////////////////////////////////////////////////////////////
    for (int i=0; i<nodes_number; i++)
    {
        Distances[i]=GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        // and now we store the verices of the element
        Points(i,0)=GetGeometry()[i].X();
        Points(i,1)=GetGeometry()[i].Y();
        Points(i,2)=GetGeometry()[i].Z();
    }

    array_1d<double,6> edge_areas;
    unsigned int n_subdivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions<MatrixType,VectorType,boost::numeric::ublas::bounded_matrix<double,4,3> >(Points, DN_DX, Distances, PartitionVolumes, 																	
                                    Partition_Shape_Functions,PartitionSigns, Partition_Gradients, Enriched_Shape_Functions,edge_areas );
    array_1d<double,4> N_at_igauss = ZeroVector(4);

    //this loops is in progress
    for(unsigned int igauss=0; igauss<n_subdivisions; igauss++)
    {
        double density;
        double viscosity;
        double tauone;

        N_at_igauss[0] = Partition_Shape_Functions(igauss,0);
        N_at_igauss[1] = Partition_Shape_Functions(igauss,1);
        N_at_igauss[2] = Partition_Shape_Functions(igauss,2);
        N_at_igauss[3] = Partition_Shape_Functions(igauss,3);
        const double partition_volume=PartitionVolumes[igauss];

        EvaluateAtGaussPoint(density, DENSITY, N_at_igauss);
        EvaluateAtGaussPoint(viscosity, VISCOSITY, N_at_igauss);
        CalculateTau(N_at_igauss,tauone, delta_t, partition_volume, density, viscosity);
        //viscous term
        CalculateViscousTerm(rDampingMatrix, DN_DX, partition_volume, viscosity, density);
        //Advective term
        CalculateAdvectiveTerm(rDampingMatrix, DN_DX, N_at_igauss, partition_volume, density);
        //calculate pressure term
        CalculatePressureTerm(rDampingMatrix, DN_DX, N_at_igauss, partition_volume, density);

        //stabilization terms PAVEL: REMOVED all the terms multiplied by tautwo
        ////CalculateDivStblTerm(rDampingMatrix, DN_DX, tautwo, Volume);
        CalculateAdvStblAllTerms(rDampingMatrix, rRightHandSideVector, DN_DX, N_at_igauss, tauone, partition_volume, density);
        CalculateGradStblAllTerms(rDampingMatrix, rRightHandSideVector, DN_DX, N_at_igauss, tauone, partition_volume, density);

    }

    //Calculate the term corresponding to the enrichment of pressure in the split elements
    //and add the term corresponding to the modification of the RHS due to enrichment
    //KRATOS_WATCH("Computing enrichment terms")

    CalculateEnrichmentTerms(rDampingMatrix, rRightHandSideVector, delta_t);
    //KRATOS_WATCH("Finished computing enrichment terms")

    CalculateResidual(rDampingMatrix, rRightHandSideVector);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateViscousTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const double volume, const double nu, const double rho)
{

    KRATOS_TRY
    double mu = nu*rho;

    int nodes_number = 4;
    int dof = 3;

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1);
            K(row, column) += mu * 1.0 * volume * (DN_DX(ii, 0) * DN_DX(jj, 0) + DN_DX(ii, 1) * DN_DX(jj, 1) + DN_DX(ii, 2) * DN_DX(jj, 2));
            K(row + 1, column + 1) += mu * 1.0 * volume * (DN_DX(ii, 0) * DN_DX(jj, 0) + DN_DX(ii, 1) * DN_DX(jj, 1) + DN_DX(ii, 2) * DN_DX(jj, 2));
            K(row + 2, column + 2) += mu * 1.0 * volume * (DN_DX(ii, 0) * DN_DX(jj, 0) + DN_DX(ii, 1) * DN_DX(jj, 1) + DN_DX(ii, 2) * DN_DX(jj, 2));
        }
    }


    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateAdvectiveTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const array_1d<double, 4 > & N, const double volume, const double density)
{
    KRATOS_TRY
    //calculate mean advective velocity and taus
    const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);

    /*
    array_1d<double,3> ms_adv_vel;
        ms_adv_vel[0] = (adv_vel0[0] - mesh_vel0[0]) + (adv_vel1[0] - mesh_vel1[0]) + (adv_vel2[0] - mesh_vel2[0]) + (adv_vel3[0] - mesh_vel3[0]);
        ms_adv_vel[1] = (adv_vel0[1] - mesh_vel0[1]) + (adv_vel1[1] - mesh_vel1[1]) + (adv_vel2[1] - mesh_vel2[1]) + (adv_vel3[1] - mesh_vel3[1]);
        ms_adv_vel[2] = (adv_vel0[2] - mesh_vel0[2]) + (adv_vel1[2] - mesh_vel1[2]) + (adv_vel2[2] - mesh_vel2[2]) + (adv_vel3[2] - mesh_vel3[2]);
    ms_adv_vel *= 0.25;
    */
    array_1d<double,3> ms_adv_vel;
    ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]) + N[3]*(adv_vel3[0] - mesh_vel3[0]);
    ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]) + N[3]*(adv_vel3[1] - mesh_vel3[1]);
    ms_adv_vel[2] = N[0]*(adv_vel0[2] - mesh_vel0[2]) + N[1]*(adv_vel1[2] - mesh_vel1[2]) + N[2]*(adv_vel2[2] - mesh_vel2[2]) + N[3]*(adv_vel3[2] - mesh_vel3[2]);

    //calculate convective term
    int nodes_number = 4;
    int dof = 3;
    int matsize = dof*nodes_number;

    boost::numeric::ublas::bounded_matrix<double, 3, 12 > conv_opr = ZeroMatrix(dof, matsize);
    boost::numeric::ublas::bounded_matrix<double, 12, 3 > shape_func = ZeroMatrix(matsize, dof);


    for (int ii = 0; ii < nodes_number; ii++)
    {
        int column = ii*dof;
        conv_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1] + DN_DX(ii, 2) * ms_adv_vel[2];
        conv_opr(1, column + 1) = conv_opr(0, column);
        conv_opr(2, column + 2) = conv_opr(0, column);

        shape_func(column, 0) = N[ii];
        shape_func(column + 1, 1) = shape_func(column, 0);
        shape_func(column + 2, 2) = shape_func(column, 0);
    }
    boost::numeric::ublas::bounded_matrix<double, 12, 12 > temp_convterm = ZeroMatrix(matsize, matsize);
    temp_convterm = prod(shape_func, conv_opr);

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        int loc_row = ii*dof;
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1);
            int loc_column = jj*dof;

            K(row, column) += volume * density * temp_convterm(loc_row, loc_column);
            K(row + 1, column + 1) += volume * density * temp_convterm(loc_row + 1, loc_column + 1);
            K(row + 2, column + 2) += volume * density * temp_convterm(loc_row + 2, loc_column + 2);
        }
    }

    KRATOS_CATCH("")

}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const array_1d<double, 4 > & N, const double volume, const double density)
{
    KRATOS_TRY
    int nodes_number = 4;
    int dof = 3;

    //Pavel: is equivalent to not multiplying the divergence op. by rho.. so I also use tau=tau/rho for the continuity equation???

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;

            K(row, column) += -1 * volume * N(jj) * DN_DX(ii, 0);
            K(column, row) += 1 * volume * N(jj) * DN_DX(ii, 0);

            K(row + 1, column) += -1 * volume * N(jj) * DN_DX(ii, 1);
            K(column, row + 1) += 1 * volume * N(jj) * DN_DX(ii, 1);

            K(row + 2, column) += -1 * volume * N(jj) * DN_DX(ii, 2);
            K(column, row + 2) += 1 * volume * N(jj) * DN_DX(ii, 2);
        }
    }

    /*
    for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                K(row, column) += -1 * volume * N(jj) * DN_DX(ii, 0);
                K(column, row) += 1 * volume * density * N(jj) * DN_DX(ii, 0);

                K(row + 1, column) += -1 * volume * N(jj) * DN_DX(ii, 1);
                K(column, row + 1) += 1 * volume * density * N(jj) * DN_DX(ii, 1);

                K(row + 2, column) += -1 * volume * N(jj) * DN_DX(ii, 2);
                K(column, row + 2) += 1 * volume * density * N(jj) * DN_DX(ii, 2);
            }
        }
    */

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateAdvStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const array_1d<double, 4 > & N, const double tauone, const double volume, const double density)
{
    KRATOS_TRY
    const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);

    array_1d<double,3> ms_adv_vel;
    ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]) + N[3]*(adv_vel3[0] - mesh_vel3[0]);
    ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]) + N[3]*(adv_vel3[1] - mesh_vel3[1]);
    ms_adv_vel[2] = N[0]*(adv_vel0[2] - mesh_vel0[2]) + N[1]*(adv_vel1[2] - mesh_vel1[2]) + N[2]*(adv_vel2[2] - mesh_vel2[2]) + N[3]*(adv_vel3[2] - mesh_vel3[2]);

    //calculate convective term
    int nodes_number = 4;
    int dof = 3;
    int matsize = dof*nodes_number;

    boost::numeric::ublas::bounded_matrix<double, 3, 12 > conv_opr = ZeroMatrix(dof, matsize);
    boost::numeric::ublas::bounded_matrix<double, 12, 3 > shape_func = ZeroMatrix(matsize, dof);

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int column = ii*dof;
        conv_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1] + DN_DX(ii, 2) * ms_adv_vel[2];
        conv_opr(1, column + 1) = conv_opr(0, column);
        conv_opr(2, column + 2) = conv_opr(0, column);

        shape_func(column, 0) = N[ii];
        shape_func(column + 1, 1) = shape_func(column, 0);
        shape_func(column + 2, 2) = shape_func(column, 0);
    }

    //build (a.grad V)(ro*a.grad U) stabilization term & assemble
    boost::numeric::ublas::bounded_matrix<double, 12, 12 > adv_stblterm = ZeroMatrix(matsize, matsize);
    adv_stblterm = tauone * prod(trans(conv_opr), conv_opr);
    //Pavel: due to using tau divided by rho, must multiply it here
    adv_stblterm*=density;

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        int loc_row = ii*dof;
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1);
            int loc_column = jj*dof;

            K(row, column) += 1.0 * volume * density * adv_stblterm(loc_row, loc_column);
            K(row, column + 1) += 1.0 * volume * density * adv_stblterm(loc_row, loc_column + 1);
            K(row, column + 2) += 1.0 * volume * density * adv_stblterm(loc_row, loc_column + 2);

            K(row + 1, column) += 1.0 * volume * density * adv_stblterm(loc_row + 1, loc_column);
            K(row + 1, column + 1) += 1.0 * volume * density * adv_stblterm(loc_row + 1, loc_column + 1);
            K(row + 1, column + 2) += 1.0 * volume * density * adv_stblterm(loc_row + 1, loc_column + 2);

            K(row + 2, column) += 1.0 * volume * density * adv_stblterm(loc_row + 2, loc_column);
            K(row + 2, column + 1) += 1.0 * volume * density * adv_stblterm(loc_row + 2, loc_column + 1);
            K(row + 2, column + 2) += 1.0 * volume * density * adv_stblterm(loc_row + 2, loc_column + 2);
        }
    }

    //build 1*tau1*(a.grad V)(grad P) & 1*tau1*(grad q)(ro*a.grad U) stabilization terms & assemble
    boost::numeric::ublas::bounded_matrix<double, 12, 4 > grad_stblterm = ZeroMatrix(matsize, nodes_number);
    grad_stblterm = tauone * prod(trans(conv_opr), trans(DN_DX));
    //Pavel: due to using tau divided by rho, must multiply it here - only the first one:tau1*(a.grad V)(grad P)
    grad_stblterm*=density;

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        int loc_row = ii*dof;
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;

            K(row, column) += 1.0 * volume * 1.0 * grad_stblterm(loc_row, jj);
            K(row + 1, column) += 1.0 * volume * 1.0 * grad_stblterm(loc_row + 1, jj);
            K(row + 2, column) += 1.0 * volume * 1.0 * grad_stblterm(loc_row + 2, jj);

            //Pavel: but the terms that go to the second equation - no! (so I divided them by density again)
            K(column, row) += 1.0 * volume * grad_stblterm(loc_row, jj);
            K(column, row + 1) += 1.0 * volume * grad_stblterm(loc_row + 1, jj);
            K(column, row + 2) += 1.0 * volume * grad_stblterm(loc_row + 2, jj);

            /*
                    K(column, row) += 1.0 * volume * density * grad_stblterm(loc_row, jj);
                    K(column, row + 1) += 1.0 * volume * density * grad_stblterm(loc_row + 1, jj);
                    K(column, row + 2) += 1.0 * volume * density * grad_stblterm(loc_row + 2, jj);
            */
        }
    }

    //build (1.0*a.grad V) (Fbody) stabilization term & assemble
    array_1d<double, 3 > bdf = ZeroVector(3);
    const array_1d<double, 3 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);


    bdf[0] = N[0]*(density * bdf0[0]) + N[1]*(density * bdf1[0]) + N[2]*(density * bdf2[0]) + N[3]*(density * bdf3[0]);
    bdf[1] = N[0]*(density * bdf0[1]) + N[1]*(density * bdf1[1]) + N[2]*(density * bdf2[1]) + N[3]*(density * bdf3[1]);
    bdf[2] = N[0]*(density * bdf0[2]) + N[1]*(density * bdf1[2]) + N[2]*(density * bdf2[2]) + N[3]*(density * bdf3[2]);


    array_1d<double, 12 > fbd_stblterm = ZeroVector(matsize);
    fbd_stblterm = tauone * 1.0 * prod(trans(conv_opr), bdf);
    //Pavel: due to using tau divided by rho, it must be  multiplied by rho here
    fbd_stblterm*=density;

    for (int ii = 0; ii < nodes_number; ++ii)
    {
        int index = ii * (dof + 1);
        int loc_index = ii*dof;
        F[index] +=  volume * fbd_stblterm[loc_index];
        F[index + 1] +=  volume * fbd_stblterm[loc_index + 1];
        F[index + 2] +=  volume * fbd_stblterm[loc_index + 2];
    }
    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateAdvMassStblTerms(MatrixType& M, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const array_1d<double, 4 > & N, const double tauone, const double volume, const double density)
{
    KRATOS_TRY

    const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);

    array_1d<double,3> ms_adv_vel;
    ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]) + N[3]*(adv_vel3[0] - mesh_vel3[0]);
    ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]) + N[3]*(adv_vel3[1] - mesh_vel3[1]);
    ms_adv_vel[2] = N[0]*(adv_vel0[2] - mesh_vel0[2]) + N[1]*(adv_vel1[2] - mesh_vel1[2]) + N[2]*(adv_vel2[2] - mesh_vel2[2]) + N[3]*(adv_vel3[2] - mesh_vel3[2]);

    //calculate convective term
    int nodes_number = 4;
    int dof = 3;
    int matsize = dof*nodes_number;

    boost::numeric::ublas::bounded_matrix<double, 3, 12 > conv_opr = ZeroMatrix(dof, matsize);
    boost::numeric::ublas::bounded_matrix<double, 12, 3 > shape_func = ZeroMatrix(matsize, dof);

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int column = ii*dof;
        conv_opr(0, column) = DN_DX(ii, 0) * ms_adv_vel[0] + DN_DX(ii, 1) * ms_adv_vel[1] + DN_DX(ii, 2) * ms_adv_vel[2];
        conv_opr(1, column + 1) = conv_opr(0, column);
        conv_opr(2, column + 2) = conv_opr(0, column);

        shape_func(column, 0) = N[ii];
        shape_func(column + 1, 1) = shape_func(column, 0);
        shape_func(column + 2, 2) = shape_func(column, 0);
    }

    //tau1*ro*Nacc.(1.0*a.grad V)
    boost::numeric::ublas::bounded_matrix<double, 12, 12 > temp_convterm = ZeroMatrix(matsize, matsize);
    temp_convterm = prod(shape_func, conv_opr);

    double fac = tauone*density;
    //PAVEL: as I am using tau divided by rho, I must multiply it here by rho	(since this is the momentum eq. stabilization)
    fac *= density;

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        int loc_row = ii*dof;
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1);
            int loc_column = jj*dof;

            M(row, column) += volume * fac * temp_convterm(loc_row, loc_column);
            M(row + 1, column + 1) += volume * fac * temp_convterm(loc_row + 1, loc_column + 1);
            M(row + 2, column + 2) += volume * fac * temp_convterm(loc_row + 2, loc_column + 2);
        }
    }


    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateGradStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const array_1d<double,4>& N, const double tauone, const double volume, const double density)
{
    KRATOS_TRY

    int nodes_number = 4;
    int dof = 3;


    //build 1*(grad q . grad p) stabilization term & assemble
    boost::numeric::ublas::bounded_matrix<double, 4, 4 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
    gard_opr = 1.0 * tauone * prod(DN_DX, trans(DN_DX));

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1) + dof;

        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;

            K(row, column) += volume * 1.0 * gard_opr(ii, jj);

        }
    }

    //build 1*(grad q) (Fbody ) stabilization term & assemble
    array_1d<double, 3 > bdf = ZeroVector(3);
    const array_1d<double, 3 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);


    bdf[0] = N[0]*(density * bdf0[0]) + N[1]*(density * bdf1[0]) + N[2]*(density * bdf2[0]) + N[3]*(density * bdf3[0]);
    bdf[1] = N[0]*(density * bdf0[1]) + N[1]*(density * bdf1[1]) + N[2]*(density * bdf2[1]) + N[3]*(density * bdf3[1]);
    bdf[2] = N[0]*(density * bdf0[2]) + N[1]*(density * bdf1[2]) + N[2]*(density * bdf2[2]) + N[3]*(density * bdf3[2]);

    array_1d<double, 4 > fbd_stblterm = ZeroVector(nodes_number);
    fbd_stblterm = tauone * prod(DN_DX, bdf);


    for (int ii = 0; ii < nodes_number; ++ii)
    {
        int index = ii * (dof + 1) + dof;
        F[index] += 1.0 * volume * fbd_stblterm[ii];
    }

    KRATOS_CATCH("")

}
//************************************************************************************
//************************************************************************************
//checked
void ASGS3D_ENR::CalculateGradMassStblTerms(MatrixType& M, const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, const array_1d<double,4>& N,  const double tauone, const double volume, const double density)
{
    KRATOS_TRY

    int nodes_number = 4;
    int dof = 3;

    //build 1*tau1*ro Nacc grad q)
    double fac = tauone*density;

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;


            M(column, row) += volume * fac * N[ii] * DN_DX(jj, 0);

            M(column, row + 1) += volume * fac * N[ii] * DN_DX(jj, 1);

            M(column, row + 2) += volume * fac * N[ii] * DN_DX(jj, 2);
        }
    }

    KRATOS_CATCH("")

}
//************************************************************************************
//************************************************************************************
//checked
void ASGS3D_ENR::AddBodyForceAndMomentum(VectorType& F,const array_1d<double, 4 > & N, const double volume, const double density)
{
    KRATOS_TRY
    int nodes_number = 4;
    int dof = 3;

    //body  & momentum term force
    for (int ii = 0; ii < nodes_number; ii++)
    {
        int index = ii * (dof + 1);
        const array_1d<double, 3 > bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);


        F[index] += volume * N[ii] * density * bdf[0];
        F[index + 1] += volume * N[ii] * density * bdf[1];
        F[index + 2] += volume * N[ii] * density * bdf[2];

    }


    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::CalculateResidual(const MatrixType& K, VectorType& F)
{
    KRATOS_TRY

    int nodes_number = 4;
    int dof = 3;


    array_1d<double, 16 > UP = ZeroVector(16);
    for (int ii = 0; ii < nodes_number; ++ii)
    {
        int index = ii * (dof + 1);
        UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[0];
        UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[1];
        UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[2];
        UP[index + 3] = GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE, 0);
    }

    F -= prod(K, UP);

    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 3;
    unsigned int node_size = dim + 1;


    if (rResult.size() != number_of_nodes * node_size)
        rResult.resize(number_of_nodes * node_size, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rResult[i * node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i * node_size + 1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i * node_size + 2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
        rResult[i * node_size + 3] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 3;
    unsigned int node_size = dim + 1;


    if (ElementalDofList.size() != number_of_nodes * node_size)
        ElementalDofList.resize(number_of_nodes * node_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        ElementalDofList[i * node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
        ElementalDofList[i * node_size + 1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
        ElementalDofList[i * node_size + 2] = GetGeometry()[i].pGetDof(VELOCITY_Z);
        ElementalDofList[i * node_size + 3] = GetGeometry()[i].pGetDof(PRESSURE);
    }
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & Output,
                           const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_WATCH("Calculate is an Empty function for this element")

}

//************************************************************************************
//************************************************************************************
void ASGS3D_ENR::CalculateTau(const array_1d<double,4>& N, double& tauone, const double time, const double volume, const double density, const double viscosity)
{
    KRATOS_TRY
    //calculate mean advective velocity and taus

    const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);

    array_1d<double,3> ms_adv_vel;
    ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]) + N[3]*(adv_vel3[0] - mesh_vel3[0]);
    ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]) + N[3]*(adv_vel3[1] - mesh_vel3[1]);
    ms_adv_vel[2] = N[0]*(adv_vel0[2] - mesh_vel0[2]) + N[1]*(adv_vel1[2] - mesh_vel1[2]) + N[2]*(adv_vel2[2] - mesh_vel2[2]) + N[3]*(adv_vel3[2] - mesh_vel3[2]);

    //ms_adv_vel[0] = 0.0;
    //ms_adv_vel[1] = 0.0;


    double advvel_norm = ms_adv_vel[0] * ms_adv_vel[0] + ms_adv_vel[1] * ms_adv_vel[1]+ ms_adv_vel[2] * ms_adv_vel[2];
    advvel_norm = sqrt(advvel_norm);

    double ele_length = pow(12.0*volume,0.333333333333333333333);
    ele_length = 0.666666667 * ele_length * 1.732;//2.0/3.0 * ele_length * sqrt(3.00)

    //double mu=density*viscosity;

    //const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
    //tauone = 1.0 / (1 / time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm / ele_length);
    //tauone = 1.0 / ( (dyn_st_beta)/time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm / ele_length);
    tauone=time/density;


    //tauone*=0.1;
    //tauone = time;
    KRATOS_CATCH("")


}


//*************************************************************************************
//*************************************************************************************

void ASGS3D_ENR::GetFirstDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dim + 1);
    if (values.size() != MatSize) values.resize(MatSize, false);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int index = i * (dim + 1);
        values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X, Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y, Step);
        values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z, Step);
        values[index + 3] = GetGeometry()[i].GetSolutionStepValue(PRESSURE, Step);

    }
}
//************************************************************************************
//************************************************************************************

void ASGS3D_ENR::GetSecondDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dim + 1);
    if (values.size() != MatSize) values.resize(MatSize, false);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int index = i * (dim + 1);
        values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X, Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y, Step);
        values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z, Step);
        values[index + 3] = 0.0;
    }
}
//************************************************************************************
//************************************************************************************
void ASGS3D_ENR::EvaluateAtGaussPoint(double& rResult, const Variable< double >& rVariable, const array_1d< double, 4 >& rShapeFunc)
{
    //compute sign of distance on gauss point
    double dist = 0.0;
    for (unsigned int i = 0; i < 4; i++)
        dist += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);

    double navg = 0.0;
    double value = 0.0;
    for (unsigned int i = 0; i < 4; i++)
    {
        if ( (dist * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)) > 0.0)
        {
            navg += 1.0;
            value += this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);

        }
    }
    if(navg != 0)
        value /= navg;
    else
    {
        rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
        for (unsigned int iNode = 1; iNode < 4; ++iNode)
            rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
    }
    rResult = value;
}



//*************************************************************************************
//*************************************************************************************
void ASGS3D_ENR::CalculateEnrichmentTerms(MatrixType& DampingMatrix, VectorType& rRightHandSideVector, const double delta_t)
{

    KRATOS_TRY

    const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& ff3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& v0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& v1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& v2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& v3_old = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY,1);

    const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    
    unsigned int nodes_number = 4;
    int dof = 3;



    //next we need to add the enrichment terms accounting for pressure disconuity in the "cut" elements, if such exist.
    //The element is "cut" if it contains both positive and negative vales of the distance function
    double distance=GetGeometry()[0].FastGetSolutionStepValue(DISTANCE); //Make sure that the case of the "cut" passing through the nodes is also accounted for correctly.. TO DO!!!
    //double temp=0.0;
    bool cut_element=false;

    ///////////////////////////////////////////////////////////
    for (unsigned int i =1; i<nodes_number; i++)
    {
        if (distance*GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0)
        {
            cut_element=true;
            break;
        }
    }
    if (cut_element==true)
    {
        /////////////////////////////////////////////////////////////////////////////
        //if the element is cut, we get the information of the subdivision	   //

        //boost::numeric::ublas::bounded_matrix<double, 4, 3 > Points;		   //
        Matrix Points(4,3);
        double Volume;								   //
        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;		   //
        //Matrix DN_DX (4,3);
        array_1d<double,4> N;							   //
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);	   //
        array_1d<double,4> Distances;

        //**************************************************************************

        Vector PartitionVolumes(6);
        Vector PartitionSigns(6);
        Matrix Partition_Shape_Functions(6,4);
        //one shape function per Gauss points
        Matrix Enriched_Shape_Functions(6,1);
        //vector of vectors (6 rows of matrices 1*3) (vector of 3 = matrix 1*3)
        std::vector<Matrix> Partition_Gradients(6, Matrix(1,3));
        /////////////////////////////////////////////////////////////////////////////
        for (unsigned int i=0; i<nodes_number; i++)
        {
            Distances[i]=GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            // and now we store the verices of the element
            Points(i,0)=GetGeometry()[i].X();
            Points(i,1)=GetGeometry()[i].Y();
            Points(i,2)=GetGeometry()[i].Z();
        }

        array_1d<double,6> edge_areas;
        unsigned int n_subdivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions<MatrixType,VectorType,boost::numeric::ublas::bounded_matrix<double,4,3> >(Points, DN_DX, Distances, PartitionVolumes, 																	
                                    Partition_Shape_Functions,PartitionSigns, Partition_Gradients, Enriched_Shape_Functions,edge_areas );
        double Lap = 0.0;	//laplacian is a 1x1 matrix since we have one enrichment shape fct per element
        boost::numeric::ublas::bounded_matrix<double,3,3> Grad_u; //gradient of velocity
        boost::numeric::ublas::bounded_matrix<double, 1, 4 >Lap_star = ZeroMatrix(1,4);  //gradN* gradN
        boost::numeric::ublas::bounded_matrix<double, 4, 1 >Lap_starT = ZeroMatrix(4,1);  //gradN gradN*
        boost::numeric::ublas::bounded_matrix<double, 4, 4 >Lap_starT_Lap_star = ZeroMatrix(4,4);
        boost::numeric::ublas::bounded_matrix<double, 4, 1 >Lap_starT_f_star = ZeroMatrix(4,1);

        boost::numeric::ublas::bounded_matrix<double, 12, 4 > Gstar_Lap_star = ZeroMatrix(12, 4);
        //the one below is not necessary- it is a transpose of the one above
        boost::numeric::ublas::bounded_matrix<double, 4, 12 > Lap_starT_Dstar = ZeroMatrix(4, 12);

        boost::numeric::ublas::bounded_matrix<double, 1, 12 > temp = ZeroMatrix(1, 12);
        boost::numeric::ublas::bounded_matrix<double, 1, 12 > Dstar = ZeroMatrix(1, 12);
        boost::numeric::ublas::bounded_matrix<double, 12, 1 > Gstar = ZeroMatrix(12, 1);
        boost::numeric::ublas::bounded_matrix<double, 12, 12 > Gstar_Dstar = ZeroMatrix(12, 12);
        boost::numeric::ublas::bounded_matrix<double, 12, 1 > Gstar_fstar = ZeroMatrix(12, 1);
        boost::numeric::ublas::bounded_matrix<double, 1, 3 > Grad_at_igauss = ZeroMatrix(1, 3);

        array_1d<double,4> N_at_igauss = ZeroVector(4);
        
        //the gradient is constant in the element
        for (unsigned int ii = 0; ii < nodes_number; ii++)
        {
            int index = dof*ii;
            temp(0, index) = DN_DX(ii, 0);
            temp(0, index + 1) = DN_DX(ii, 1);
            temp(0, index + 2) = DN_DX(ii, 2);
        }
        //KRATOS_WATCH(temp)
        //KRATOS_WATCH(tau)
        double f_star=0.0;
        double tau=0.0;
        //FOR CONVECTIVE CONTRIBUTION TO THE RHS WE WILL NEED THE GRADIENT OF VELOCITY!
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //KRATOS_WATCH("___________________________________________________________________________________________________________________")
        //Now we shall loop over the Gauss points of the sub-elements

        for(unsigned int igauss=0; igauss<n_subdivisions; igauss++)
        {
            double density;
            double viscosity;

            Grad_at_igauss = Partition_Gradients[igauss];
            N_at_igauss[0] = Partition_Shape_Functions(igauss,0);
            N_at_igauss[1] = Partition_Shape_Functions(igauss,1);
            N_at_igauss[2] = Partition_Shape_Functions(igauss,2);
            N_at_igauss[3] = Partition_Shape_Functions(igauss,3);

            EvaluateAtGaussPoint(density, DENSITY, N_at_igauss);
            EvaluateAtGaussPoint(viscosity, VISCOSITY, N_at_igauss);
            CalculateTau(N_at_igauss, tau , delta_t, Volume, density, viscosity);

            //Grad_at_igauss:(1, 3); DN_DX: (4,3)
            Lap_star += tau*PartitionVolumes[igauss] * prod(Grad_at_igauss, trans(DN_DX) );

            //scalar product of Grad_at_igauss*Grad_at_igauss
            Lap += tau*PartitionVolumes[igauss]*(Grad_at_igauss(0,0)*Grad_at_igauss(0,0)+
                                                 Grad_at_igauss(0,1)*Grad_at_igauss(0,1)+
                                                 Grad_at_igauss(0,2)*Grad_at_igauss(0,2));
            //KRATOS_WATCH(Lap)

            //and the enriched gradient operator: D*=N* DN_DX... since N* is a scalar (we just have one enrichment function), its just a product
            //now Enriched_Shape_Functions is a matrix (6,1)
            //NOTE THAT DIV must be multilpied by density, while Grad - NO! Put the densities of the Gauss points later! TO DO!

            //COMPUTING CONVECTION OPERATOR///////////////////////vk*dNi/dxk
            array_1d<double,3> adv_vel = ZeroVector(3);
            adv_vel=N_at_igauss[0]*fv0 + N_at_igauss[1]*fv1 + N_at_igauss[2]*fv2 + N_at_igauss[3]*fv3;
            //KRATOS_WATCH(adv_vel)
            Matrix conv_opr=ZeroMatrix(3,12);
            Matrix shape_func=ZeroMatrix(12,3);
            Matrix conv_contrib=ZeroMatrix(1,12);
            //inertia contribution:
            Matrix acc_contrib=ZeroMatrix(1,12);            
            int dof=3;
            for (unsigned int ii = 0; ii < nodes_number; ii++)
            {
                int column = ii*dof;
                conv_opr(0, column) = DN_DX(ii, 0) * adv_vel[0] + DN_DX(ii, 1) * adv_vel[1] + DN_DX(ii, 2) * adv_vel[2];
                conv_opr(1, column + 1) = conv_opr(0, column);
                conv_opr(2, column + 2) = conv_opr(0, column);

                shape_func(column, 0) = N_at_igauss[ii];
                shape_func(column + 1, 1) = shape_func(column, 0);
                shape_func(column + 2, 2) = shape_func(column, 0);
            }
            //tau*grad q*(conv_op)
            conv_contrib=tau*density*PartitionVolumes[igauss]*prod(Grad_at_igauss, conv_opr);
            acc_contrib=tau*density*PartitionVolumes[igauss]*prod(Grad_at_igauss, trans(shape_func))/delta_t;

            //Dstar+=density*PartitionVolumes[igauss]*Enriched_Shape_Functions(igauss,0)*temp;
            Dstar+=PartitionVolumes[igauss]*Enriched_Shape_Functions(igauss,0)*temp;
            Dstar+=conv_contrib;
            Dstar+=acc_contrib;

            //gradient of pressure integrated by part: -(nabla v)*pstar??
            Gstar-=PartitionVolumes[igauss]*Enriched_Shape_Functions(igauss,0)*trans(temp);
            //v gradp* contribution to the momentum equation
            Gstar+=trans(conv_contrib);
            ///////////////////////////NOW THE RHS TERMS!////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////
            //Vector body_force_and_old_acc = ZeroVector(3);
            Vector body_force = ZeroVector(3);
            Vector old_acc = ZeroVector(3);
            //body force
            //body_force_and_old_acc=N_at_igauss[0]*(ff0+v0_old/delta_t) + N_at_igauss[1]*(ff1+v1_old/delta_t) + N_at_igauss[2]*(ff2+v2_old/delta_t) + N_at_igauss[3]*(ff3+v3_old/delta_t);
            body_force=N_at_igauss[0]*ff0 + N_at_igauss[1]*ff1 + N_at_igauss[2]*ff2 + N_at_igauss[3]*ff3;
            old_acc=(N_at_igauss[0]*v0_old + N_at_igauss[1]*v1_old + N_at_igauss[2]*v2_old + N_at_igauss[3]*v3_old)/delta_t;
            //ms_aux1 - is the product of: (nabla q, f)
            //f*=tau*nabla q* (-nu lap u - grad p - f - du/dt-v*grad v) ... first product is zero for the linear elements
            //check if it should be multiplied by density
            double f_term=Grad_at_igauss(0,0)*(body_force[0]+old_acc[0]) + Grad_at_igauss(0,1)*(body_force[1]+old_acc[1]) + Grad_at_igauss(0,2)*(body_force[2]+old_acc[2]);

            f_star+=tau*density*PartitionVolumes[igauss]*f_term;
        }


        //now the auxilliary matrix, the product of enriched grad and div (G*D*), we shall need it afterwards to modify the first equation
        Gstar_Dstar=prod(Gstar, Dstar)/Lap;
        Gstar_Lap_star= prod(Gstar, Lap_star)/Lap;
        //gradN gradN*D* to be added to the momentum equation

        Lap_starT=trans(Lap_star);
        Lap_starT_Dstar=prod (Lap_starT, Dstar)/Lap;
        Lap_starT_Lap_star=prod(Lap_starT, Lap_star)/Lap;
        Gstar_fstar=(f_star*Gstar)/Lap;

        //this is the term that corresponds to: gradN gradN*F*
        Lap_starT_f_star=f_star*Lap_starT;
        Lap_starT_f_star/=Lap;

        //Now we modify the left hand*-side matrix K and the right-hand-side vector
        //now we add the computed values to the "Damp Matrix" (which is the LHS matrix (without mass term)
        //as the "Damp matrix" is monolithic, it is 16x16. Following the structure of: v1_x v1_y v1_z p1. v2_x v2_y etc..
        //so every fourth term should not be touched (i.e. it corresponds to the continuity equation contribution) (Gstar_Dstar is 12x12)
        for (int i=0; i<4; i++)
        {
            int row=i*4;
            int row_small=i*3;

            rRightHandSideVector[row]-=Gstar_fstar(row_small, 0);
            rRightHandSideVector[row+1]-=Gstar_fstar(row_small+1,0);
            rRightHandSideVector[row+2]-=Gstar_fstar(row_small+2,0);
            rRightHandSideVector[row+3]-=Lap_starT_f_star(i, 0);//this term was checked and seems to be OK


            for (int j=0; j<4; j++)
            {
                int column=j*4;
                int column_small=j*3;

                DampingMatrix(row,column)-=Gstar_Dstar(row_small,column_small);
                DampingMatrix(row,column+1)-=Gstar_Dstar(row_small,column_small+1);
                DampingMatrix(row,column+2)-=Gstar_Dstar(row_small,column_small+2);
                /////////////////////////////////////////////////////////////////
                DampingMatrix(row,column+3)-=Gstar_Lap_star(row_small, j);	//this term was checked and seems to be OK

                DampingMatrix(row+1,column)-=Gstar_Dstar(row_small+1,column_small);
                DampingMatrix(row+1,column+1)-=Gstar_Dstar(row_small+1,column_small+1);
                DampingMatrix(row+1,column+2)-=Gstar_Dstar(row_small+1,column_small+2);
                /////////////////////////////////////////////////////////////////
                DampingMatrix(row+1,column+3)-=Gstar_Lap_star(row_small+1, j);//this term was checked and seems to be OK

                DampingMatrix(row+2,column)-=Gstar_Dstar(row_small+2,column_small);
                DampingMatrix(row+2,column+1)-=Gstar_Dstar(row_small+2,column_small+1);
                DampingMatrix(row+2,column+2)-=Gstar_Dstar(row_small+2,column_small+2);
                /////////////////////////////////////////////////////////////////
                DampingMatrix(row+2,column+3)-=Gstar_Lap_star(row_small+2, j);//this term was checked and seems to be OK

                DampingMatrix(row+3,column)-=Lap_starT_Dstar(i, column_small);//this term was ALSO checked and ALSO seems to be OK
                DampingMatrix(row+3,column+1)-=Lap_starT_Dstar(i, column_small+1);//this term was ALSO checked and ALSO seems to be OK
                DampingMatrix(row+3,column+2)-=Lap_starT_Dstar(i, column_small+2);//this term was ALSO checked and ALSO seems to be OK
                ////////////////////////////////////////////////////////////////
                DampingMatrix(row+3,column+3)-=Lap_starT_Lap_star(i, j);//this term was checked and ALSO seems to be OK
            }


        }

        //here we end the part of the code that is executed only if the element is a cut element
        //end of "if cut_element ==true
    }
    //else if (cut_element==false)
    //KRATOS_WATCH("Uncut element - not doing anything")

    KRATOS_CATCH("")

}
} // Namespace Kratos


