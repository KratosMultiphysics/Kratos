cd ..\applications\DelaunayMeshingApplication
mkdir external_modules
cd external_modules
::curl -L https://github.com/PFEM/tetgen/archive/0f0c1c95388d9a515e2383533c6b4f18524d9e76.zip -o tetgen.zip
powershell.exe -Command [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; wget https://github.com/PFEM/tetgen/archive/0f0c1c95388d9a515e2383533c6b4f18524d9e76.zip -OutFile tetgen.zip
%SEVEN_ZIP% e tetgen.zip -o"tetgen"
cd ..\..\..\cmake_build

%CMAKE% .. -G "Visual Studio 15 2017 Win64"                                       ^
-DCMAKE_CXX_FLAGS="/D EXCLUDE_EMBEDDED_PYTHON_DEBUG /Z7" 	                        ^
-DCMAKE_C_FLAGS="/D EXCLUDE_EMBEDDED_PYTHON_DEBUG /Z7" 		                        ^
-DBOOST_ROOT=%BOOST%                                                            ^
-DPYTHON_EXECUTABLE=%PYTHON%                          		 	                    ^
-DCMAKE_BUILD_TYPE="Custom"  							                                      ^
-DDEM_APPLICATION=ON                                                            ^
-DEXTERNAL_SOLVERS_APPLICATION=OFF                                              ^
-DFLUID_DYNAMICS_APPLICATION=ON                                                 ^
-DSTRUCTURAL_MECHANICS_APPLICATION=ON                                           ^
-DCONTACT_STRUCTURAL_MECHANICS_APPLICATION=ON                                   ^
-DSWIMMING_DEM_APPLICATION=ON                                                   ^
-DMESH_MOVING_APPLICATION=ON                                                    ^
-DSOLID_MECHANICS_APPLICATION=ON                                                ^
-DCONSTITUTIVE_MODELS_APPLICATION=ON                                            ^
-DDELAUNAY_MESHING_APPLICATION=ON                                               ^
-DCONTACT_MECHANICS_APPLICATION=ON                                              ^
-DPFEM_APPLICATION=ON                                                           ^
-DPFEM_SOLID_MECHANICS_APPLICATION=ON                                           ^
-DPFEM_FLUID_DYNAMICS_APPLICATION=ON                                            ^
-DMETIS_APPLICATION=OFF                                                         ^
-DPARMETIS_ROOT_DIR="UNSET"                                                     ^
-DTRILINOS_APPLICATION=OFF                                                      ^
-DTRILINOS_ROOT="UNSET"                                                         ^
-DINSTALL_EMBEDDED_PYTHON=ON                                                    ^
-DINCLUDE_FEAST=OFF                                                             
