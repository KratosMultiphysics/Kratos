#include "fractional_step_discontinuous.h"
#include "utilities/geometry_utilities.h"
#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

/**The FractionalStepDiscontinuous element is designed to allow including embedded objects in the computations.
     * The element employes discontinuous shape functions as described by Ausas [.....].
     * This allows considering an independent flow field on the two sides of the element. The practical problem is however that
     * by construction the normal gradient (with respect to the embedded surface) of the finite element variable is ZERO.
     *
     * The position and direction of the discontinuity is obtained by reading the value
     *
     *      const Vector& distances = this->GetValue(ELEMENTAL_DISTANCES);
     *
     * which defines a level set function at element level. Note that
     * this function is allowed to be discontinuous from one element to the next.
     *
     * cut elements are identified by the boolean variable "SPLIT_ELEMENT" as
     *
     *      bool split_element = this->GetValue(SPLIT_ELEMENT);
     *
     * which has to be set to "true" in order for the embedded obstacle to be considered.
     * ELEMENTAL_DISTANCES will NOT be read unless SPLIT_ELEMENT is set to "true"
     *
     * normal velocity at the cut surface is imposed weakly while solving the pressure system.
     * tangential velocity is also imposed weakly in solving the momentum equation. Please note that
     * since the normal gradient of velocity is zero, the imposition of the tangential component of the velocity can
     * be understood as a numerical wall law.
     * In both cases the value of the velocity to be imposed is given in the elemental variable EMBEDDED_VELOCITY
     * and can be fetched as
     *
     *         const array_1d<double,3>& vel = this->GetValue(EMBEDDED_VELOCITY);
     *
     * Note that the element assumes ASSIGNED the variables ELEMENTAL_DISTANCES, SPLIT_ELEMENT, EMBEDDED_VELOCITY
     * hence this assignment has to be done by the user (automatic tools are available) prior to the computational loop
     * Note also that if the element is to be used only for a one-sided flow, the flow variables
     * shall be fixed ( by the user, not by this element) in the area to be deactivated.
     *
     * As a final comment, the element allows as a feature to ease FSI coupling with lightweight structure:
     * if an estimate of the structural mass PER UNIT AREA is available, this data will be used in the
     * computation of the pressure. This follows the ideas in [.....] and allows avoiding added-mass effects.
     * The estimate of the structural mass can be given as:
     *
     * another feature of the element is the capability of adding a mass-like matrix epsilon*1/nu*Mass to the pressure equation
     * this is useful for cases in which there may be domains without any Neumann condition, where the pressure would be otherwise undefined.
     * the value of the coefficient epsilon (which should be a very small number in the range 1e-6-1e-15) is controlled by
     * epsilon = rCurrentProcessInfo[PRESSURE_MASSMATRIX_COEFFICIENT]
     *
     *
     */

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::AddMomentumSystemTerms(Matrix& rLHSMatrix,
        Vector& rRHSVector,
        const double Density,
        const Vector& rConvOperator,
        const array_1d<double, 3 > & rBodyForce,
        const double OldPressure,
        const double TauOne,
        const double TauTwo,
        const array_1d<double, 3 > & rMomentumProjection,
        const double MassProjection,
        const ShapeFunctionsType& rN,
        const ShapeFunctionDerivativesType& rDN_DX,
        const double Weight)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;

    array_1d<double, TDim + 1 > OldPressures;
    for (SizeType i = 0; i < NumNodes; ++i)
        OldPressures[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

    array_1d<double, TDim> grad_p = prod(trans(rDN_DX), OldPressures);

    for (SizeType i = 0; i < NumNodes; ++i)
    {

        // Build RHS
        for (SizeType d = 0; d < TDim; ++d)
        {
            // Body force
            double RHSi = Density * rN[i] * rBodyForce[d];
            // Pressure gradient (integrated by parts)
            //RHSi += rDN_DX(i,d) * OldPressure;
            RHSi -= rN[i] * grad_p[d];

            // Momentum Stabilization
            RHSi += Density * rConvOperator[i] * TauOne * (/*rBodyForce[d] + rOldPressureGradient[d]*/ -rMomentumProjection[d]);
            // Mass Stabilization
            RHSi -= rDN_DX(i, d) * TauTwo * MassProjection;

            rRHSVector[FirstRow + d] += Weight * RHSi;
        }

        // Build LHS
        for (SizeType j = 0; j < NumNodes; ++j)
        {
            // Convective term
            double Kij = Density * rN[i] * rConvOperator[j];

            // Streamline stabilization
            Kij += Density * rConvOperator[i] * TauOne * Density * rConvOperator[j];

            Kij *= Weight;

            for (SizeType d = 0; d < TDim; ++d)
                rLHSMatrix(FirstRow + d, FirstCol + d) += Kij;

            // Mass-GLS (TauTwo) stabiliziation term
            for (SizeType m = 0; m < TDim; ++m)
                for (SizeType n = 0; n < TDim; ++n)
                    rLHSMatrix(FirstRow + m, FirstCol + n) += Weight * rDN_DX(i, m) * TauTwo * rDN_DX(j, n);

            FirstCol += TDim;
        }
        FirstRow += TDim;
        FirstCol = 0;
    }
}

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    Matrix NContainer;
    ShapeFunctionDerivativesArrayType DN_DX;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const SizeType NumGauss = NContainer.size1();

    // Stabilization parameters
    double ElemSize = this->ElementSize();
    double TauOne;
    double TauTwo;

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& N = row(NContainer, g);
        const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

        // Evaluate required variables at the integration point
        double Density;
        double pgauss;
        //        array_1d<double,3> Velocity = ZeroVector(3);
        //        array_1d<double,3> MeshVelocity = ZeroVector(3);
        array_1d<double, 3 > BodyForce = ZeroVector(3);
        array_1d<double, 3 > MomentumProjection = ZeroVector(3);

        this->EvaluateInPoint(Density, DENSITY, N);
        this->EvaluateInPoint(pgauss, PRESSURE, N);
        //        this->EvaluateInPoint(Velocity,VELOCITY,N);
        //        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
        this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
        //        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MomentumProjection, PRESS_PROJ, N);

        //        // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
        //        double OldPressure;
        //        this->EvaluateInPoint(OldPressure,PRESSURE,N,0);

        array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
        this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

        //        // For ALE: convective velocity
        //        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

        // Stabilization parameters
        array_1d<double, 3 > ConvVel = ZeroVector(3);
        this->EvaluateConvVelocity(ConvVel, N);
        double Viscosity = this->EffectiveViscosity(Density,N, rDN_DX, ElemSize, rCurrentProcessInfo);
        this->CalculateTau(TauOne, TauTwo, ElemSize, ConvVel, Density, Viscosity, rCurrentProcessInfo);
        //TauOne = 0.0;

        //        // Evaluate convection operator Velocity * Grad(N)
        //        Vector UGradN(NumNodes);
        //        this->EvaluateConvection(UGradN,ConvVel,mDN_DX);

        //double DivU;
        //this->EvaluateDivergenceInPoint(DivU,VELOCITY,rDN_DX);
        array_1d<double, 3 > ugauss;
        this->EvaluateInPoint(ugauss, VELOCITY, N);
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            double aux = 0.0;
            for (SizeType d = 0; d < TDim; ++d)
                aux += ugauss[d] * rDN_DX(i, d);

            rRightHandSideVector[i] += GaussWeight * aux;
        }

        // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
        const double LaplacianCoeff = 1.0 / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);

        // Add convection, stabilization and RHS contributions to the local system equation
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            // LHS contribution
            for (SizeType j = 0; j < NumNodes; ++j)
            {
                double Lij = 0.0;
                for (SizeType d = 0; d < TDim; ++d)
                    Lij += rDN_DX(i, d) * rDN_DX(j, d);
                Lij *= (LaplacianCoeff + TauOne);

                rLeftHandSideMatrix(i, j) += GaussWeight * Lij;
            }

            // RHS contribution

            // Velocity divergence
            //double RHSi = - N[i] * DivU;
            double RHSi = 0.0;
            for (SizeType d = 0; d < TDim; ++d)
            {
                //                double Conv = UGradN[0] * rGeom[0].FastGetSolutionStepValue(VELOCITY)[d];
                //                for (SizeType j = 1; j < NumNodes; ++j) Conv += UGradN[i] * rGeom[i].FastGetSolutionStepValue(VELOCITY)[d];
                // Momentum stabilization
                RHSi += rDN_DX(i, d) * TauOne * (Density * (BodyForce[d]/* - Conv*/) - OldPressureGradient[d] - MomentumProjection[d]);
            }

            rRightHandSideVector[i] += GaussWeight * RHSi;
        }

        //adding a penalization term to ask for the pressure to be zero in average
        double epsilon = 1.0e-15/Viscosity;
        if (rCurrentProcessInfo.Has(PRESSURE_MASSMATRIX_COEFFICIENT))
            epsilon = rCurrentProcessInfo[PRESSURE_MASSMATRIX_COEFFICIENT]/Viscosity;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            for (SizeType j = 0; j < NumNodes; ++j)
            {
                rLeftHandSideMatrix(i,i) += GaussWeight*epsilon*N[i]*N[j];
            }
            rRightHandSideVector[i] -= GaussWeight*epsilon*N[i]*pgauss;
        }
    }

    //
    bool split_element = this->GetValue(SPLIT_ELEMENT);
    if (split_element == true)
    {
        const array_1d<double, 3 > & vel = this->GetValue(EMBEDDED_VELOCITY);
        const Vector& distances = this->GetValue(ELEMENTAL_DISTANCES);

        double Volume_tot;
        BoundedMatrix<double, 4, 3 > DN_DXcontinuous;
        array_1d<double, 4 > Ncontinuous;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DXcontinuous, Ncontinuous, Volume_tot);

        array_1d<double, TDim> grad_d = prod(trans(DN_DXcontinuous), distances);
        grad_d /= norm_2(grad_d);

        double vn = grad_d[0] * vel[0];
        for (unsigned int i = 1; i < TDim; i++) vn += grad_d[i] * vel[i];

        //KRATOS_WATCH(medge_areas);

        //loop on edges
        unsigned int edge_counter = 0.0;
        for (unsigned int i = 0; i < TDim; i++)
        {
            for (unsigned int j = i + 1; j < TDim + 1; j++)
            {
                const double aux = medge_areas[edge_counter] * vn;
                if (distances[i] * distances[j] < 0.0) //cut edge
                {
                    if (distances[i] > 0)
                    {
                        rRightHandSideVector[i] += aux;
                        rRightHandSideVector[j] -= aux;
                    }
                    else
                    {
                        rRightHandSideVector[i] -= aux;
                        rRightHandSideVector[j] += aux;
                    }
                }
                edge_counter++;
            }
        }

        //FSI STRUCTURAL MASS TERM
        if (rCurrentProcessInfo.Has(DENSITY))
        {
            const double structural_density = rCurrentProcessInfo[DENSITY];
            edge_counter = 0.0;
            for (unsigned int i = 0; i < TDim; i++)
            {
                for (unsigned int j = i + 1; j < TDim + 1; j++)
                {
                    if (distances[i] * distances[j] < 0.0) //cut edge
                    {
                        //double li = fabs(distances[i]) + fabs(distances[j]);
                        double Nj = 1.0; //distances[i]/li;
                        double Ni = 1.0; // - Nj;

                        const double aux = 1.0 / (structural_density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]) * medge_areas[edge_counter];
                        rLeftHandSideMatrix(i, i) += Ni * aux;
                        rLeftHandSideMatrix(j, j) += Nj * aux;
                    }
                    edge_counter++;
                }
            }
        }

    }
}

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::Calculate(const Variable<array_1d<double, 3 > > &rVariable,
        array_1d<double, 3 > &rOutput,
        const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == ADVPROJ)
    {
        double Tmp = 0.0;
        this->Calculate(DIVPROJ, Tmp, rCurrentProcessInfo);
    }
    else if (rVariable == CONV_PROJ)
    {
        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        const unsigned int LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();

        VectorType ConvTerm = ZeroVector(LocalSize);
        VectorType PresTerm = ZeroVector(LocalSize);
        VectorType DivTerm = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer, g);
            const double GaussWeight = GaussWeights[g];

            for (unsigned int i = 0; i < NumNodes; i++)
                NodalArea[i] += N[i] * GaussWeight;

            this->CalculateProjectionRHS(ConvTerm, PresTerm, DivTerm, N, DN_DX[g], GaussWeight);
        }

        // Carefully write results to nodal variables, to avoid parallelism problems
        unsigned int RowIndex = 0;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double, 3 > & rConvVal = this->GetGeometry()[i].FastGetSolutionStepValue(CONV_PROJ);
            array_1d<double, 3 > & rPresVal = this->GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rConvVal[d] += ConvTerm[RowIndex];
                rPresVal[d] += PresTerm[RowIndex];
                ++RowIndex;
            }
            this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += DivTerm[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        }
    }
    else if (rVariable == VELOCITY)
    {
        GeometryType& rGeom = this->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const SizeType LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();

        VectorType NodalVelCorrection = ZeroVector(LocalSize);



        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; ++g)
        {
            const ShapeFunctionsType& N = row(NContainer, g);
            const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

            double Density;

            this->EvaluateInPoint(Density, DENSITY, N);

            const double Coeff = GaussWeights[g] / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);

            array_1d<double, TDim> DeltaPressureGradient = ZeroVector(TDim);
            this->EvaluateGradientInPoint(DeltaPressureGradient, PRESSURE_OLD_IT, rDN_DX);

            // Calculate contribution to the gradient term (RHS)
            //             double DeltaPressure;
            //             this->EvaluateInPoint(DeltaPressure,PRESSURE_OLD_IT,N);

            SizeType RowIndex = 0;

            for (SizeType i = 0; i < NumNodes; ++i)
            {
                for (SizeType d = 0; d < TDim; ++d)
                {
                    //                     NodalVelCorrection[RowIndex++] += Coeff * rDN_DX(i,d) * DeltaPressure;
                    NodalVelCorrection[RowIndex++] -= Coeff * N(i) * DeltaPressureGradient[d];
                }
            }
        }

        SizeType Index = 0;

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double, 3 > & rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < TDim; ++d)
            {
                rTemp[d] += NodalVelCorrection[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }
}

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::Calculate(const Variable<double> &rVariable,
        double &rOutput,
        const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DIVPROJ)
    {
        const GeometryType& rGeom = this->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const unsigned int LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();

        VectorType MomentumRHS = ZeroVector(LocalSize);
        VectorType MassRHS = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer, g);
            const double GaussWeight = GaussWeights[g];

            for (unsigned int i = 0; i < NumNodes; i++)
                NodalArea[i] += N[i] * GaussWeight;

            this->CalculateProjectionRHS(MomentumRHS, MassRHS, N, DN_DX[g], GaussWeight);
        }

        // Carefully write results to nodal variables, to avoid parallelism problems
        unsigned int RowIndex = 0;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double, 3 > & rMomValue = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
            for (unsigned int d = 0; d < TDim; ++d)
                rMomValue[d] += MomentumRHS[RowIndex++];
            this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        }
    }
}

/*
 * protected FractionalStepDiscontinuous<TDim> functions
 */

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
        Matrix& rNContainer,
        Vector& rGaussWeights)
{
    bool split_element = this->GetValue(SPLIT_ELEMENT);
    /*	KRATOS_WATCH(this->Id())
        KRATOS_WATCH(this->GetValue(SPLIT_ELEMENT))*/
    if (split_element == true)
    {

        array_1d<double, (TDim + 1) > distances = this->GetValue(ELEMENTAL_DISTANCES);
        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < TDim + 1; i++)
            if (distances[i] >= 0) npos++;
            else nneg++;

        if (nneg != 0 && npos != 0)
        {
            BoundedMatrix<double, 3 * (TDim - 1), (TDim + 1) > Nenriched;
            array_1d<double, (3 * (TDim - 1)) > volumes;
            BoundedMatrix<double, (TDim + 1), TDim > coords;
            BoundedMatrix<double, 3 * (TDim - 1), (TDim + 1) > Ngauss;
            array_1d<double, (3 * (TDim - 1)) > signs;
            std::vector< Matrix > gauss_gradients(3 * (TDim - 1));

            BoundedMatrix<double, (TDim + 1), TDim > DN_DX;
            array_1d<double, (TDim + 1) > N;



            for (unsigned int i = 0; i < TDim + 1; i++)
            {
                const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
                volumes[i] = 0.0;
                for (unsigned int j = 0; j < TDim; j++) coords(i, j) = xyz[j];
            }

            for (unsigned int i = 0; i < 3 * (TDim - 1); i++) gauss_gradients[i].resize(TDim + 1, TDim, false); //2 values of the 2 shape functions, and derivates in (xy) direction).
            unsigned int ndivisions = DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched, medge_areas);

            if (rGaussWeights.size() != ndivisions)
                rGaussWeights.resize(ndivisions, false);
            if (rNContainer.size1() != ndivisions || rNContainer.size2() != TDim + 1)
                rNContainer.resize(ndivisions, TDim + 1, false);
            rDN_DX.resize(ndivisions);

            //KRATOS_WATCH(this->Id())
            for (unsigned int g = 0; g < ndivisions; g++)
            {
                for (unsigned int d = 0; d < TDim + 1; d++)
                    rNContainer(g, d) = Nenriched(g, d);
                rGaussWeights[g] = volumes[g];
                rDN_DX[g] = gauss_gradients[g];
                /*			KRATOS_WATCH(g)
                                        KRATOS_WATCH(rDN_DX[g])
                                        KRATOS_WATCH(rGaussWeights[g])*/
            }
            //			KRATOS_WATCH(rNContainer)


        }
        else
        {
            const GeometryType& rGeom = this->GetGeometry();
            Vector DetJ;
            rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_2);
            rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

            rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2), false);

            for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
                rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
        }
    }
    else
    {
        const GeometryType& rGeom = this->GetGeometry();
        Vector DetJ;
        rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_2);
        rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

        rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2), false);

        for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
            rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }
    /*
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const unsigned int NumGauss = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);

    // Initialize arrays to proper size
    rDN_DX.resize(NumGauss);
    rDetJ.resize(NumGauss);

    const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( GeometryData::GI_GAUSS_2 );

    // Temporary container for inverse of J
    Matrix InvJ;

    GeometryType::JacobiansType J;
    rGeom.Jacobian( J, GeometryData::GI_GAUSS_2 );

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        // calculate inverse of the jacobian and its determinant
        MathUtils<double>::InvertMatrix( J[g], InvJ, rDetJ[g] );

        // calculate the shape function derivatives in global coordinates
        rDN_DX[g].resize(NumNodes,TDim);
        noalias( rDN_DX[g] ) = prod( DN_De[g], InvJ );
    }*/
}

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    FractionalStep<TDim>::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    //add a numeric wall law to take into account the tangential velocity at the boundaries
    bool split_element = this->GetValue(SPLIT_ELEMENT);
    if (split_element == true)
    {


        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
        {

            double Volume_tot;
            BoundedMatrix<double, 4, 3 > DN_DXcontinuous;
            array_1d<double, 4 > Ncontinuous;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DXcontinuous, Ncontinuous, Volume_tot);

            const array_1d<double, 3 > & embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
            //        KRATOS_WATCH(embedded_vel)

            const Vector& distances = this->GetValue(ELEMENTAL_DISTANCES);

            double Density;
            this->EvaluateInPoint(Density, DENSITY, Ncontinuous);
            double ElemSize = this->ElementSize();
            const double Viscosity = this->EffectiveViscosity(Density, Ncontinuous, DN_DXcontinuous, ElemSize, rCurrentProcessInfo);


            const double min_dist_in_viscous_term = 0.001 * ElemSize; //the 0.001 it is arbitrary. This terms is only to avoid divisions by zero



            array_1d<double, TDim> normal = prod(trans(DN_DXcontinuous), distances);
            normal /= norm_2(normal);
            //KRATOS_WATCH(normal)
            //compute the block diagonal parallel projection
            //defined as the operator which extracts the part of the velocity
            //tangent to the embedded wall
            #ifdef KRATOS_USE_AMATRIX
            BoundedMatrix<double, TDim, TDim> block = IdentityMatrix(TDim);
            #else
            BoundedMatrix<double, TDim, TDim> block = IdentityMatrix(TDim, TDim);
            #endif
            BoundedMatrix<double, TDim, TDim> nn_matrix = outer_prod(normal, normal);
            noalias(block) -= nn_matrix;
            //KRATOS_WATCH(block)

            array_1d<double, 3 > tangent_vel_dirichlet = prod(block, embedded_vel);
            //          KRATOS_WATCH(tangent_vel_dirichlet)

            //add tangential component
            //loop on cut edges
            unsigned int edge_counter;

            if (this->IsDefined(SLIP) == true && this->IsNot(SLIP))
            {
                edge_counter = 0;
                for (unsigned int i = 0; i < TDim; i++)
                {
                    unsigned int base_i = i * (TDim);

                    const array_1d<double, 3 > tangent_vi = prod(block, this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY));

                    const double di = std::max(min_dist_in_viscous_term, fabs(distances[i])); //avoid divisions by zero

                    const double coeff_i = Viscosity / di;

                    for (unsigned int j = i + 1; j < TDim + 1; j++)
                    {
                        array_1d<double, 3 > tangent_vj = prod(block, this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY));

                        const double dj = std::max(min_dist_in_viscous_term, fabs(distances[j])); //avoid divisions by zero
                        const double coeff_j = Viscosity / dj;

                        unsigned int base_j = j * (TDim);

                        if (distances[i] * distances[j] < 0.0) //cut edge
                        {
                            for (unsigned int kkk = 0; kkk < TDim; kkk++)
                            {
                                rRightHandSideVector(base_i + kkk) += coeff_i * medge_areas[edge_counter] * (tangent_vel_dirichlet[kkk] - tangent_vi[kkk]);
                                rRightHandSideVector(base_j + kkk) += coeff_j * medge_areas[edge_counter] * (tangent_vel_dirichlet[kkk] - tangent_vj[kkk]);

                                for (unsigned int lll = 0; lll < TDim; lll++)
                                {
                                    //node I contribution
                                    rLeftHandSideMatrix(base_i + kkk, base_i + lll) += coeff_i * medge_areas[edge_counter] * block(kkk, lll);

                                    //node J contribution
                                    rLeftHandSideMatrix(base_j + kkk, base_j + lll) += coeff_j * medge_areas[edge_counter] * block(kkk, lll);

                                }
                            }
                        }
                        edge_counter++;
                    }
                }
            }

            //compute penalty_coeff
            const double embedded_vnorm = norm_2(embedded_vel);
            double penalty_coeff = Density*(Viscosity + embedded_vnorm * ElemSize ) / ElemSize;
//                 penalty_coeff = 0.0
//                 KRATOS_WATCH(penalty_coeff);

            //add normal component
            edge_counter = 0.0;
            for (unsigned int i = 0; i < TDim; i++)
            {
                unsigned int base_i = i * (TDim);
                const array_1d<double,3>& v_i = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                array_1d<double,TDim> v_i_diff;
                for (unsigned int kkk = 0; kkk < TDim; kkk++)  v_i_diff[kkk] = v_i[kkk] - embedded_vel[kkk];

                array_1d<double,TDim> rhs_penalty_I = prod(nn_matrix,v_i_diff);

                for (unsigned int j = i + 1; j < TDim + 1; j++)
                {
                    if (distances[i] * distances[j] < 0.0) //cut edge
                    {
                        unsigned int base_j = j * (TDim);

                        const array_1d<double,3>& v_j = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                        array_1d<double,TDim> v_j_diff;
                        for (unsigned int kkk = 0; kkk < TDim; kkk++)  v_j_diff[kkk] = v_j[kkk] - embedded_vel[kkk];


                        array_1d<double,TDim> rhs_penalty_J = prod(nn_matrix,v_j_diff);

                        for (unsigned int kkk = 0; kkk < TDim; kkk++)
                        {
                            rRightHandSideVector(base_i + kkk) -= penalty_coeff * medge_areas[edge_counter] * (rhs_penalty_I[kkk]); //contribution to I
                            rRightHandSideVector(base_j + kkk) -= penalty_coeff * medge_areas[edge_counter] * (rhs_penalty_J[kkk]); //contribution to J

                            for (unsigned int lll = 0; lll < TDim; lll++)
                            {
                                //node I contribution
                                rLeftHandSideMatrix(base_i + kkk, base_i + lll) += penalty_coeff * medge_areas[edge_counter] * nn_matrix(kkk, lll);

                                //node J contribution
                                rLeftHandSideMatrix(base_j + kkk, base_j + lll) += penalty_coeff * medge_areas[edge_counter] * nn_matrix(kkk, lll);

                            }
                        }
                    }
                    edge_counter++;
                }

            }
        }

    }

    KRATOS_CATCH("");
}


///*
// * FractionalStepDiscontinuous<TDim> static members
// */

//template< unsigned int TDim >
//static const FractionalStepDiscontinuous<TDim>::ShapeFunctionsType FractionalStepDiscontinuous<TDim>::msNg = FractionalStepDiscontinuous<TDim>::InitializeShapeFunctions();

/*
 * Template class definition (this should allow us to compile the desired template instantiations)
 */

template class FractionalStepDiscontinuous < 2 >;
template class FractionalStepDiscontinuous < 3 >;

}
