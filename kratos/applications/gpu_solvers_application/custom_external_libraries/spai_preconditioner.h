//
// SPAI preconditioner on GPU
//
// Computes the values of the SPAI preconditioner on GPU
// Written by: Borja Servan
// Optimized and styled by: Farshid Mossaiby

__global__ void GPU_SPAIPreconditioner_CSR_Kernel(const size_t Rows, const size_t *A_Columns, const size_t *A_RowIndices, const double *A_Values, double *M_Values)
{
	const size_t MaxNZ = 40;
	const size_t MaxIterates = 300;
	
    double Alpha, Beta, BNorm,
		Rho, RhoP,
		P[MaxNZ], P2[MaxNZ], 
		Pp[MaxNZ], Pp2[MaxNZ], 
		Q[MaxNZ], Q2[MaxNZ], 
		R[MaxNZ], R2[MaxNZ], 
		B[MaxNZ], X[MaxNZ],
		Matrix[MaxNZ * MaxNZ];

	size_t Idx, i_begin, i_end, j_begin, j_end, n;

	Idx = GlobalIdx(); 

	if (Idx < Rows)
	{
		i_begin = A_RowIndices[Idx];
		i_end = A_RowIndices[Idx + 1];
		
		n = i_end - i_begin;

		for (size_t i1 = i_begin; i1 < i_end; i1++)
		{
			size_t j = A_Columns[i1];
			
			j_begin = A_RowIndices[j];
			j_end = A_RowIndices[j + 1];
			
			for (size_t i2 = i_begin; i2 < i_end; i2++)
			{
				Matrix[(i1 - i_begin) * n + (i2 - i_begin)] = 0.00;
				
				for (size_t i3 = j_begin; i3 < j_end; i3++)
					if (A_Columns[i3] == A_Columns[i2])
					{
						Matrix[(i1 - i_begin) * n + (i2 - i_begin)] = A_Values[i3];
						break;
					}
			}
			
			B[i1 - i_begin] = j == Idx ? 1.00 : 0.00;
		}

		BNorm = 0.00;
		
		for (size_t i1 = 0; i1 < n; i1++)
		{
			BNorm += B[i1] * B[i1];
			X[i1] = B[i1];
			
			double T = 0.00;

			for (size_t i2 = 0; i2 < n; i2++)
				T += Matrix[i1 * n + i2] * X[i2];
			
			T -= B[i1];
			R2[i1] = R[i1] = T;
		}

		size_t k1 = 0;
		
		while (k1 < MaxIterates)
		{
			Rho = 0.00;
			
			for (size_t j = 0; j < n; j++)
			{
				Rho += R[j] * R2[j];
				P[j] = R[j];
				P2[j] = R2[j];
			}

			// TODO: Should this be exact?
			if (Rho == 0.00)
                break;
			else
			{			    if (k1 > 0)
			    {
					Beta = Rho / RhoP;
				    
				    for (size_t j = 0; j < n; j++)
				    {
					    P[j] += Beta * Pp[j];
					    P2[j] += Beta * Pp2[j];
				    }
			    }
    	                
			    for (size_t i2 = 0; i2 < n; i2++)
			    {
					double T = 0.00, T2 = 0.00;

					for (size_t i3 = 0; i3 < n; i3++)
					{
						T += Matrix[i2 * n + i3] * P[i3];
					    T2 += Matrix[i3 * n + i2] * P2[i3];
					}
						
					Q[i2] = T;
					Q2[i2] = T2;
			    }

			    Alpha = 0.00;
			    
			    for (size_t j = 0; j < n; j++)
				    Alpha += P2[j] * Q[j];

			    Alpha = Rho / Alpha;

			    for (size_t j = 0; j < n; j++)
			    {
				    X[j] -= Alpha * P[j];
				    R[j] -= Alpha * Q[j];
				    R2[j] -= Alpha * Q2[j];
				}

			    RhoP = 0.00;
			    
			    for (size_t j = 0; j < n; j++)
					RhoP += R[j] * R[j];
			    
				RhoP /= BNorm;

				if (RhoP < 1.0e-10)
				    break;
				
				k1++;
				
				RhoP = Rho;

				for (size_t j = 0; j < n; j++)
				{
					Pp[j] = P[j];
					Pp2[j] = P2[j];
				}
		    }
        }

		for (size_t i1 = 0; i1 < n; i1++)
			M_Values[i_begin + i1] = X[i1];
    }
}
