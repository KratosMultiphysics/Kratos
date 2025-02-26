import numpy as np
import os

current_path = 'applications/RomApplication/tests/save_rom_coefficients_process_test_files/save_rom_coefficients_process_test_line_search/'
q = np.load(current_path+'ObtainedOutputSaveRomCoefficientsProcessLineSearch.npy')
# q = 6.08734*1.2
print(q)

path_compressible_potential = 'applications/RomApplication/tests/compressible_potential_test_files/'
basis = np.load(path_compressible_potential+'rom_data/RightBasisMatrix.npy')
print(sum(basis**2))
map = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 24, 27, 29, 31, 33, 35, 37, 36, 39, 41, 40, 43]
basis_corrected = basis[map]

# print(basis_corrected)

# print(q*basis_corrected)


print(sum(basis_corrected**2))