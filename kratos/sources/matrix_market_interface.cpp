#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

#include "includes/matrix_market_interface.h"

namespace Kratos {

// Type checks
template<>
constexpr bool IsCorrectType<double>(MM_typecode& mm_code)
{
    return mm_is_real(mm_code);
}

template<>
constexpr bool IsCorrectType<std::complex<double>>(MM_typecode& mm_code)
{
    return mm_is_complex(mm_code);
}

// Matrix I/O routines
bool ReadMatrixMarketMatrixEntry(FILE *f, int& I, int& J, double& V)
{
    return fscanf(f, "%d %d %lg", &I, &J, &V) == 3;
}

bool ReadMatrixMarketMatrixEntry(FILE *f, int& I, int& J, std::complex<double>& V)
{
    double real;
    double imag;

    const int i = fscanf(f, "%d %d %lg %lg", &I, &J, &real, &imag);
    V = std::complex<double>(real, imag);
    return i == 4;
}

template <typename CompressedMatrixType> 
bool ReadMatrixMarketMatrix(const char *FileName, CompressedMatrixType &M)
{
    typedef typename CompressedMatrixType::value_type ValueType;

    // Open MM file for reading
    FILE *f = fopen(FileName, "r");

    if (f == NULL)
    {
        printf("ReadMatrixMarketMatrix(): unable to open %s.\n", FileName);
        return false;
    }

    // Process MM file header
    MM_typecode mm_code;

    if (mm_read_banner(f, &mm_code) != 0)
    {
        printf("ReadMatrixMarketMatrix(): unable to read MatrixMarket banner.\n");
        fclose(f);
        return false;
    }

    if (!mm_is_valid(mm_code))
    {
        printf("ReadMatrixMarketMatrix(): invalid MatrixMarket banner.\n");
        fclose(f);
        return false;
    }

    // Check for supported types of MM file
    if (!(mm_is_coordinate(mm_code) && mm_is_sparse(mm_code)))
    {
        printf("ReadMatrixMarketMatrix(): unsupported MatrixMarket type, \"%s\".\n",  mm_typecode_to_str(mm_code));
        fclose(f);
        return false;
    }

    // Read MM dimensions and NNZ
    int size1, size2, nnz;

    if (mm_read_mtx_crd_size(f, &size1, &size2, &nnz) != 0)
    {
        printf("ReadMatrixMarketMatrix(): cannot read dimensions and NNZ.\n");
        fclose(f);
        return false;
    }

    // Allocate temporary arrays
    int *I = new int[nnz];
    int *J = new int[nnz];
    ValueType *V = new ValueType[nnz];

    // Check if matrix type matches MM file
    if (!IsCorrectType<ValueType>(mm_code))
    {
        printf("ReadMatrixMarketMatrix(): MatrixMarket type, \"%s\" does not match provided matrix type.\n",  mm_typecode_to_str(mm_code));
        fclose(f);
        return false;
    }

    // Read MM file
    // Pattern file, only non-zero structure
    if (mm_is_pattern(mm_code))
        for (int i = 0; i < nnz; i++)
        {
            if (fscanf(f, "%d %d", &I[i], &J[i]) != 2)
            {
                printf("ReadMatrixMarketMatrix(): invalid data.\n");
                fclose(f);

                delete[] I;
                delete[] J;
                delete[] V;

                return false;
            }

            // Adjust to 0-based
            I[i]--;
            J[i]--;

            // Set all values to 1.00
            V[i] = 1.00;
        }
    else
        for (int i = 0; i < nnz; i++)
        {
            if (! ReadMatrixMarketMatrixEntry(f, I[i], J[i], V[i]))
            {
                printf("ReadMatrixMarketMatrix(): invalid data.\n");
                fclose(f);

                delete[] I;
                delete[] J;
                delete[] V;

                return false;
            }

            // Adjust to 0-based
            I[i]--;
            J[i]--;
        }

    fclose(f);

    // Second stage
    int *nz = new int[size1];

    for (int i = 0; i < size1; i++)
        nz[i] = 0;

    // Count non-zeros on each line
    if (mm_is_symmetric(mm_code))
        for (int i = 0; i < nnz; i++)
        {
            if (I[i] == J[i])
                nz[I[i]]++;
            else
            {
                nz[I[i]]++;
                nz[J[i]]++;
            }
        }
    else
        for (int i = 0; i < nnz; i++)
            nz[I[i]]++;

    // Find out total number of non-zeros
    int nnz2;

    if (mm_is_symmetric(mm_code))
    {
        int diagonals = 0;

        for (int i = 0; i < nnz; i++)
            if (I[i] == J[i])
                diagonals++;

        nnz2 = diagonals + 2 * (nnz - diagonals);
    }
    else
        nnz2 = nnz;

    // Fill in an almost-CSR data structure
    int *filled = new int[size1];
    int *indices = new int[size1];
    int *columns = new int[nnz2];
    ValueType *values = new ValueType[nnz2];

    indices[0] = 0;
    for (int i = 1; i < size1; i++)
        indices[i] = indices[i - 1] + nz[i - 1];

    for (int i = 0; i < size1; i++)
        filled[i] = 0;

    if (mm_is_symmetric(mm_code))
        for (int i = 0; i < nnz; i++)
            if (I[i] == J[i])
            {
                int index;

                index = indices[I[i]] + filled[I[i]];
                columns[index] = J[i];
                values[index] = V[i];
                filled[I[i]]++;
            }
            else
            {
                int index;

                index = indices[I[i]] + filled[I[i]];
                columns[index] = J[i];
                values[index] = V[i];
                filled[I[i]]++;

                index = indices[J[i]] + filled[J[i]];
                columns[index] = I[i];
                values[index] = V[i];
                filled[J[i]]++;
            }
    else
        for (int i = 0; i < nnz; i++)
        {
            int index;

            index = indices[I[i]] + filled[I[i]];
            columns[index] = J[i];
            values[index] = V[i];
            filled[I[i]]++;
        }

    // Create the matrix
    CompressedMatrixType *m = new CompressedMatrixType(size1, size2, nnz2);

    int k = 0;

    for (int i = 0; i < size1; i++)
        for (int j = 0; j < nz[i]; j++)
            (*m)(i, columns[indices[i] + j]) = values[k++];

    M.resize(m->size1(), m->size2(), false);

    M = *m;

    delete[] I;
    delete[] J;
    delete[] V;

    delete[] filled;
    delete[] indices;
    delete[] columns;
    delete[] values;
    delete[] nz;

    delete m;

    return true;
}

void SetMatrixMarketValueTypeCode(MM_typecode& mm_code, const double& value)
{
    mm_set_real(&mm_code);
}

void SetMatrixMarketValueTypeCode(MM_typecode& mm_code, const std::complex<double>& value)
{
    mm_set_complex(&mm_code);
}

int WriteMatrixMarketMatrixEntry(FILE *f, int I, int J, const double& entry)
{
    return fprintf(f, "%d %d %.12e\n", I, J, entry);
}

int WriteMatrixMarketMatrixEntry(FILE *f, int I, int J, const std::complex<double>& entry)
{
    return fprintf(f, "%d %d %.12e %.12e\n", I, J, std::real(entry), std::imag(entry));
}

template <typename CompressedMatrixType> 
bool WriteMatrixMarketMatrix(const char *FileName, CompressedMatrixType &M, bool Symmetric)
{
    // Open MM file for writing
    FILE *f = fopen(FileName, "w");

    if (f == NULL)
    {
        printf("WriteMatrixMarketMatrix(): unable to open %s.\n", FileName);
        return false;
    }

    // Write MM file header
    MM_typecode mm_code;

    mm_initialize_typecode(&mm_code);

    mm_set_matrix(&mm_code);
    mm_set_coordinate(&mm_code);
    SetMatrixMarketValueTypeCode(mm_code, *M.begin1().begin());

    if (Symmetric)
        mm_set_symmetric(&mm_code);
    else
        mm_set_general(&mm_code);

    mm_write_banner(f, mm_code);

    // Find out the actual number of non-zeros in case of a symmetric matrix
    int nnz;

    if (Symmetric)
    {
        nnz = 0;

        typename CompressedMatrixType::iterator1 a_iterator = M.begin1();

        for (unsigned int i = 0; i < M.size1(); i++)
        {
    #ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (typename CompressedMatrixType::iterator2 row_iterator = a_iterator.begin(); row_iterator != a_iterator.end(); ++row_iterator)
    #else
            for (typename CompressedMatrixType::iterator2 row_iterator = begin(a_iterator, iterator1_tag()); row_iterator != end(a_iterator, iterator1_tag()); ++row_iterator)
    #endif
            {
                if (a_iterator.index1() >= row_iterator.index2())
                    nnz++;
            }

            a_iterator++;
        }
    }
    else
        nnz = M.nnz();

    // Write MM file sizes
    mm_write_mtx_crd_size(f, M.size1(), M.size2(), nnz);

    if (Symmetric)
    {
        typename CompressedMatrixType::iterator1 a_iterator = M.begin1();

        for (unsigned int i = 0; i < M.size1(); i++)
        {
    #ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (typename CompressedMatrixType::iterator2 row_iterator = a_iterator.begin(); row_iterator != a_iterator.end(); ++row_iterator)
    #else
            for (typename CompressedMatrixType::iterator2 row_iterator = begin(a_iterator, iterator1_tag()); row_iterator != end(a_iterator, iterator1_tag()); ++row_iterator)
    #endif
            {
                int I = a_iterator.index1(), J = row_iterator.index2();

                if (I >= J)
                    if (WriteMatrixMarketMatrixEntry(f, I+1, J+1, *row_iterator) < 0)
                    {
                        printf("WriteMatrixMarketMatrix(): unable to write data.\n");
                        fclose(f);
                        return false;
                    }
            }

            a_iterator++;
        }
    }
    else
    {
        typename CompressedMatrixType::iterator1 a_iterator = M.begin1();

        for (unsigned int i = 0; i < M.size1(); i++)
        {
    #ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (typename CompressedMatrixType::iterator2 row_iterator = a_iterator.begin(); row_iterator != a_iterator.end(); ++row_iterator)
    #else
            for (typename CompressedMatrixType::iterator2 row_iterator = begin(a_iterator, iterator1_tag()); row_iterator != end(a_iterator, iterator1_tag()); ++row_iterator)
    #endif
            {
                int I = a_iterator.index1(), J = row_iterator.index2();

                if (WriteMatrixMarketMatrixEntry(f, I+1, J+1, *row_iterator) < 0)
                {
                    printf("WriteMatrixMarketMatrix(): unable to write data.\n");
                    fclose(f);
                    return false;
                }
            }

            a_iterator++;
        }
    }

    fclose(f);

    return true;
}

// Vector I/O routines

bool ReadMatrixMarketVectorEntry(FILE *f, double& entry)
{
    return fscanf(f, "%lg", &entry) == 1;
}

bool ReadMatrixMarketVectorEntry(FILE *f, std::complex<double>& entry)
{
    double real;
    double imag;
    const int i = fscanf(f, "%lg %lg", &real, &imag);
    entry = std::complex<double>(real, imag);
    return i == 2;
}

template <typename VectorType> bool ReadMatrixMarketVector(const char *FileName, VectorType &V)
{
    typedef typename VectorType::value_type ValueType;

    // Open MM file for reading
    FILE *f = fopen(FileName, "r");

    if (f == NULL)
    {
        printf("ReadMatrixMarketVector(): unable to open %s.\n", FileName);
        return false;
    }

    // Process MM file header
    MM_typecode mm_code;

    if (mm_read_banner(f, &mm_code) != 0)
    {
        printf("ReadMatrixMarketVector(): unable to read MatrixMarket banner.\n");
        fclose(f);
        return false;
    }

    if (!mm_is_valid(mm_code))
    {
        printf("ReadMatrixMarketVector(): invalid MatrixMarket banner.\n");
        fclose(f);
        return false;
    }

    // Check for supported types of MM file
    if (!((!mm_is_pattern(mm_code)) && mm_is_array(mm_code)))
    {
        printf("ReadMatrixMarketVector(): unsupported MatrixMarket type, \"%s\".\n",  mm_typecode_to_str(mm_code));
        fclose(f);
        return false;
    }

    // Read MM dimensions
    int size1, size2;

    if (mm_read_mtx_array_size(f, &size1, &size2) != 0)
    {
        printf("ReadMatrixMarketVector(): cannot read dimensions.\n");
        fclose(f);
        return false;
    }

    // Check MM dimensions
    if (size2 != 1)
    {
        printf("ReadMatrixMarketVector(): not a N x 1 array.\n");
        fclose(f);
        return false;
    }

    VectorType *v = new VectorType(size1);
    ValueType T;

    // Check if vector type matches MM file
    if (!IsCorrectType<ValueType>(mm_code))
    {
        printf("ReadMatrixMarketVector(): MatrixMarket type, \"%s\" does not match provided vector type.\n",  mm_typecode_to_str(mm_code));
        fclose(f);
        return false;
    }

    // Read MM file

    for (int i = 0; i < size1; i++)
    {
        if (! ReadMatrixMarketVectorEntry(f, T))
        {
            printf("ReadMatrixMarketVector(): invalid data.\n");
            fclose(f);

            return false;
        }

        (*v)(i) = T;
    }

    fclose(f);

    V = *v;

    delete v;

    return true;
}

int WriteMatrixMarketVectorEntry(FILE *f, const double& entry)
{
    return fprintf(f, "%e\n", entry);
}

int WriteMatrixMarketVectorEntry(FILE *f, const std::complex<double>& entry)
{
    return fprintf(f, "%e %e\n", std::real(entry), std::imag(entry));
}

template <typename VectorType> 
bool WriteMatrixMarketVector(const char *FileName, VectorType &V)
{
    // Open MM file for writing
    FILE *f = fopen(FileName, "w");

    if (f == NULL)
    {
        printf("WriteMatrixMarketVector(): unable to open %s.\n", FileName);
        return false;
    }

    // Write MM file header
    MM_typecode mm_code;

    mm_initialize_typecode(&mm_code);

    mm_set_matrix(&mm_code);
    mm_set_array(&mm_code);
    SetMatrixMarketValueTypeCode(mm_code, V(0));

    mm_write_banner(f, mm_code);

    // Write MM file sizes
    mm_write_mtx_array_size(f, V.size(), 1);

    for (unsigned int i = 0; i < V.size(); i++)
        if (WriteMatrixMarketVectorEntry(f, V(i)) < 0)
        {
            printf("WriteMatrixMarketVector(): unable to write data.\n");
            fclose(f);
            return false;
        }

    fclose(f);

    return true;
}

template KRATOS_API(KRATOS_CORE) bool ReadMatrixMarketMatrix<Kratos::CompressedMatrix>(const char *FileName, Kratos::CompressedMatrix &M);
template KRATOS_API(KRATOS_CORE) bool ReadMatrixMarketMatrix<Kratos::ComplexCompressedMatrix>(const char *FileName, Kratos::ComplexCompressedMatrix &M);

template KRATOS_API(KRATOS_CORE) bool WriteMatrixMarketMatrix<Kratos::CompressedMatrix>(const char *FileName, Kratos::CompressedMatrix &M, bool Symmetric);
template KRATOS_API(KRATOS_CORE) bool WriteMatrixMarketMatrix<Kratos::ComplexCompressedMatrix>(const char *FileName, Kratos::ComplexCompressedMatrix &M, bool Symmetric);

template KRATOS_API(KRATOS_CORE) bool ReadMatrixMarketVector<Kratos::Vector>(const char *FileName, Kratos::Vector &V);
template KRATOS_API(KRATOS_CORE) bool ReadMatrixMarketVector<Kratos::ComplexVector>(const char *FileName, Kratos::ComplexVector &V);

template KRATOS_API(KRATOS_CORE) bool WriteMatrixMarketVector<Kratos::Vector>(const char *FileName, Kratos::Vector &V);
template KRATOS_API(KRATOS_CORE) bool WriteMatrixMarketVector<Kratos::ComplexVector>(const char *FileName, Kratos::ComplexVector &V);
}