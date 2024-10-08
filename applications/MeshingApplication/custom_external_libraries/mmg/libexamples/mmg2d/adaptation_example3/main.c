#include <stdio.h>
#include <stdlib.h>
#include <mmg/mmg2d/libmmg2d.h>
#include <math.h>

int main() {
  MMG5_pMesh mesh2 = NULL;
  MMG5_pSol met2 = NULL;
  int ier = 0;
  MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh2, MMG5_ARG_ppMet, &met2, MMG5_ARG_end);
  MMG2D_Set_iparameter(mesh2, met2, MMG2D_IPARAM_verbose, 10);

  // 4 vertices, no tris, no quads, no edges so far
  MMG2D_Set_meshSize(mesh2, 4, 0, 0, 0);

  // a rectangle [0, 2 * M_PI] x [-M_PI / 2, M_PI / 2] for spherical coordinates in usual lon/lat way
  double v[8] = {0, -M_PI/2, 2 * M_PI, -M_PI/2, 2 * M_PI, M_PI/2, 0, M_PI/2};
  MMG2D_Set_vertices(mesh2, v, NULL);

  MMG2D_Set_dparameter(mesh2, met2, MMG2D_DPARAM_hsiz, 0.01);

  // generate a regular fine mesh of the square in meshing mode
  ier = MMG2D_mmg2dmesh(mesh2, met2);
  if (ier) {
    fprintf(stderr, "error %i during meshing.\n", ier);
    exit(1);
  }

  // save the "computational geometry" mesh and the metric set by MMG2D_mmg2dmesh
  // For Medit users
  MMG2D_saveMesh(mesh2, "cg.mesh");
  MMG2D_saveSol(mesh2, met2, "cg.sol");
  // For Gmsh users
  MMG2D_saveMshMesh(mesh2, NULL, "cg.msh");

  MMG2D_Free_solutions(mesh2,met2);

  // remesh with anisotropic metric
  int np, nt, nquad, na;
  MMG2D_Get_meshSize(mesh2, &np, &nt, &nquad, &na);
  double* verts = malloc(2 * np * sizeof(double));
  MMG2D_Get_vertices(mesh2, verts, NULL, NULL, NULL);
  MMG2D_Set_solSize(mesh2, met2, MMG5_Vertex, np, MMG5_Tensor);
  for (int i = 0; i < np; i++) {

    // latitude
    double y = verts[2 * i + 1];
    // metric on the sphere, see standard textbooks on elementary differential geometry
    // Gaussian fundamental quantities
    double E = cos(y) * cos(y);
    double F = 0.0;
    double G = 1.0;

    // The symetric metric tensor (E F \\ F G) may be diagonalized in the basis of its
    // eigenvectors. The eigenvalue \lambda_i associated to the eigenvector e_i prescribes a
    // length l_i = 1/\sqrt(\lambda_i) in the direction e_i. This formula allows to compute
    // a size factor to have sufficiently small edges on the surface of the sphere of radius 1.
    //
    // For a detailed explanation, see:
    // https://forum.mmgtools.org/t/how-to-scale-the-metric-for-anisotropic-meshing-in-mmg2d-example-provided/441/2?u=algiane
    double factor = 100.0;

    // we add a small amount to the 1-1-entry to make the metric non-singular everywhere
    const double eps = 1.e-1;
    MMG2D_Set_tensorSol(met2, factor * E + eps, factor * F, factor * G, i+1);
  }

  // disable hsiz because it is incompatible with an anisotropic metric
  MMG2D_Set_dparameter(mesh2, met2, MMG2D_DPARAM_hsiz, -1);

  // set gradation to a not too restrictive value
  MMG2D_Set_dparameter(mesh2, met2, MMG2D_DPARAM_hgrad, 1.5);

  // save the "computational geometry" mesh and the setted metric
  // For Medit users
  MMG2D_saveMesh(mesh2, "cg-with-met.mesh");
  MMG2D_saveSol(mesh2, met2, "cg-with-met.sol");
  // For Gmsh users
  MMG2D_saveMshMesh(mesh2, met2, "cg-with-met.msh");

  // do it, remesh!
  ier = MMG2D_mmg2dlib(mesh2, met2);

  if (ier != 0) {
    fprintf(stderr, "error %i during remeshing.\n", ier);
    exit(1);
  }

  // Save at Gmsh file format
  MMG2D_saveMshMesh(mesh2, NULL, "out2.msh");

  // Save at Medit one
  MMG2D_saveMesh(mesh2, "out2.mesh");
  MMG2D_saveSol(mesh2,met2, "out2.sol");

  // map to 3d and save as Medit file
  MMG2D_Get_meshSize(mesh2, &np, &nt, &nquad, &na);
  verts = realloc(verts, 2 * np * sizeof(double));
  int* tris = malloc(3 * nt * sizeof(int));
  MMG2D_Get_vertices(mesh2, verts, NULL, NULL, NULL);
  MMG2D_Get_triangles(mesh2, tris, NULL, NULL);

  printf("\n\n -- Save the final mesh 'sphere-end.mesh' at Medit file format.\n");

  FILE *fp = fopen("sphere-end.mesh", "w+");

  fprintf(fp, "MeshVersionFormatted 2\n\n Dimension 3\n\n Vertices %d\n\n", np);
  for (int i = 0; i < np; i++) {
    double x = verts[2 * i];
    double y = verts[2 * i + 1];
    /* Parametrization of a sphere of radius 1 from the longitude (x) and
     * latitude (y) */
    double newx = cos(x) * cos(y);
    double newy = sin(x) * cos(y);
    double newz = sin(y);
    fprintf(fp, "%f %f %f 0\n", newx, newy, newz);
  }
  fprintf(fp, "Triangles %d\n", nt);

  for (int i = 0; i < nt; i++)
    fprintf(fp, "%i %i %i 0\n", tris[3 * i], tris[3 * i +1], tris[3 * i + 2]);

  fprintf(fp, "End\n");
  fclose(fp);

  free(verts);
  free(tris);
  MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh2, MMG5_ARG_ppMet, &met2, MMG5_ARG_end);

  return 0;
}
