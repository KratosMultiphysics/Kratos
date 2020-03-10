// Translated to C++ from the original work in FORTRAN by Alberto Ferriz, 2016

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

void Diagonalize(const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3]);

int main() {

    // Import the mesh (.msh) and spheres (.sph) files containing the cluster information
    std::ifstream infile("file_name.msh");
    std::ifstream infilesph("file_name.sph");
    std::string line, linesph;
    infile.ignore(80,'\n'); infile.ignore(80,'\n');

    // Check the number of nodes (NUM_OF_NODES) and the number of elements (NUM_OF_ELEMENTS)
    int node_counter = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int a;
        double b, c, d, e;
        if (iss >> a >> b >> c >> d) {
            node_counter++;
        } else break;
    }
    std::cout << "\nNumber of nodes: " << node_counter << '\n';

    infile.ignore(80,'\n'); infile.ignore(80,'\n');

    int element_counter = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int a;
        double b, c, d, e;
        if (iss >> a >> b >> c >> d >> e) {
            element_counter++;
        } else break;
    }
    std::cout << "Number of elements: " << element_counter << '\n';

    int spheres_counter = 0;
    while (std::getline(infilesph, linesph)) {
        std::istringstream iss(linesph);
        double Xcoord, Ycoord, Zcoord, Rad;
        if (iss >> Xcoord >> Ycoord >> Zcoord >> Rad) {
            spheres_counter++;
        } else break;
    }
    std::cout << "Number of spheres: " << spheres_counter << '\n';

    infile.seekg(0, std::ios::beg);
    infile.ignore(80,'\n'); infile.ignore(80,'\n');

    const int NUM_OF_NODES = node_counter;
    const int NUM_OF_ELEMENTS = element_counter;
    const double density = 1;
    double* tcoord = new double[3*NUM_OF_NODES];
    int* Nconec = new int[4*NUM_OF_ELEMENTS];
    double* Vnerc = new double[9*NUM_OF_ELEMENTS];
    double* Volum = new double[NUM_OF_ELEMENTS];
    double* Vmass = new double[NUM_OF_ELEMENTS];
    double* BARIC = new double[3*NUM_OF_ELEMENTS];
    double* Local = new double[3*NUM_OF_ELEMENTS];
    double VNERT[9];
    double I[3][3];
    double Q[3][3];
    double D[3][3];
    node_counter = 0;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int a;
        double b, c, d, e;
        if (iss >> a >> b >> c >> d) {
            tcoord[node_counter * 3 + 0] = b;
            tcoord[node_counter * 3 + 1] = c;
            tcoord[node_counter * 3 + 2] = d;
            node_counter++;
        } else break;
    }
    infile.ignore(80,'\n'); infile.ignore(80,'\n');

    element_counter = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int a;
        double b, c, d, e;
        if (iss >> a >> b >> c >> d >> e) {
            Nconec[element_counter * 4 + 0] = b;
            Nconec[element_counter * 4 + 1] = c;
            Nconec[element_counter * 4 + 2] = d;
            Nconec[element_counter * 4 + 3] = e;
            element_counter++;
        } else break;
    }

    infilesph.clear();
    infilesph.seekg(0, std::ios::beg);

    const int NUM_OF_SPHERES = spheres_counter;

    double* sphcoord = new double[3*NUM_OF_SPHERES];
    double sphrad[NUM_OF_SPHERES];

    spheres_counter = 0;
    while (std::getline(infilesph, linesph)) {
        std::istringstream iss(linesph);
        double Xcoord, Ycoord, Zcoord, Rad;
        if (iss >> Xcoord >> Ycoord >> Zcoord >> Rad) {
            sphcoord[spheres_counter * 3 + 0] = Xcoord;
            sphcoord[spheres_counter * 3 + 1] = Ycoord;
            sphcoord[spheres_counter * 3 + 2] = Zcoord;
            sphrad[spheres_counter] = Rad;
            spheres_counter++;
        } else break;
    }

    double total_volume = 0.0;

    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        int NodeA=Nconec[element_counter * 4 + 0];
        int NodeB=Nconec[element_counter * 4 + 1];
        int NodeC=Nconec[element_counter * 4 + 2];
        int NodeD=Nconec[element_counter * 4 + 3];

        double ValU1=tcoord[(NodeB-1) * 3 + 0] - tcoord[(NodeA-1) * 3 + 0];
        double ValU2=tcoord[(NodeB-1) * 3 + 1] - tcoord[(NodeA-1) * 3 + 1];
        double ValU3=tcoord[(NodeB-1) * 3 + 2] - tcoord[(NodeA-1) * 3 + 2];
        double ValV1=tcoord[(NodeC-1) * 3 + 0] - tcoord[(NodeA-1) * 3 + 0];
        double ValV2=tcoord[(NodeC-1) * 3 + 1] - tcoord[(NodeA-1) * 3 + 1];
        double ValV3=tcoord[(NodeC-1) * 3 + 2] - tcoord[(NodeA-1) * 3 + 2];
        double ValW1=tcoord[(NodeD-1) * 3 + 0] - tcoord[(NodeA-1) * 3 + 0];
        double ValW2=tcoord[(NodeD-1) * 3 + 1] - tcoord[(NodeA-1) * 3 + 1];
        double ValW3=tcoord[(NodeD-1) * 3 + 2] - tcoord[(NodeA-1) * 3 + 2];

        double Volu0= ValU1*ValV2*ValW3+ValU2*ValV3*ValW1+ValU3*ValW2*ValV1-ValU3*ValV2*ValW1-ValU2*ValV1*ValW3-ValU1*ValW2*ValV3;
        Volum[element_counter]=Volu0/6;
        total_volume += Volum[element_counter];
        Vmass[element_counter] = Volum[element_counter]*density;

        // Calculate the barycentre of each tetrahedron
        BARIC[element_counter * 3 + 0]= (tcoord[(NodeA-1) * 3 + 0]+tcoord[(NodeB-1) * 3 + 0]+tcoord[(NodeC-1) * 3 + 0]+tcoord[(NodeD-1) * 3 + 0])/4;
        BARIC[element_counter * 3 + 1]= (tcoord[(NodeA-1) * 3 + 1]+tcoord[(NodeB-1) * 3 + 1]+tcoord[(NodeC-1) * 3 + 1]+tcoord[(NodeD-1) * 3 + 1])/4;
        BARIC[element_counter * 3 + 2]= (tcoord[(NodeA-1) * 3 + 2]+tcoord[(NodeB-1) * 3 + 2]+tcoord[(NodeC-1) * 3 + 2]+tcoord[(NodeD-1) * 3 + 2])/4;
    }

    std::cout << "\nTotal volume: " << total_volume << '\n';

    // Calculate the total mass of the cluster
    double Vmaspiedra=0.0;
    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {
        Vmaspiedra+=Vmass[element_counter];
    }
    std::cout << "Total mass: " << Vmaspiedra << '\n';

    // Calculate the cluster's centre of gravity
    double Valor1=0.0;
    double Valor2=0.0;
    double Valor3=0.0;

    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        Valor1+=Vmass[element_counter]*BARIC[element_counter * 3 + 0];
        Valor2+=Vmass[element_counter]*BARIC[element_counter * 3 + 1];
        Valor3+=Vmass[element_counter]*BARIC[element_counter * 3 + 2];
    }

    double Xcdgrav=Valor1/Vmaspiedra;
    double Ycdgrav=Valor2/Vmaspiedra;
    double Zcdgrav=Valor3/Vmaspiedra;
    std::cout << "Centroid: " << Xcdgrav << " " << Ycdgrav << " " << Zcdgrav << "\n\n";

    // Move the object so that its centre of gravity coincides with the origin
    for (int node_counter = 0; node_counter < NUM_OF_NODES; node_counter++) {
        tcoord[node_counter * 3 + 0] -= Xcdgrav;
        tcoord[node_counter * 3 + 1] -= Ycdgrav;
        tcoord[node_counter * 3 + 2] -= Zcdgrav;
    }

    // Move the spheres accordingly
    for (int spheres_counter = 0; spheres_counter < NUM_OF_SPHERES; spheres_counter++) {
        sphcoord[spheres_counter * 3 + 0] -= Xcdgrav;
        sphcoord[spheres_counter * 3 + 1] -= Ycdgrav;
        sphcoord[spheres_counter * 3 + 2] -= Zcdgrav;
    }

    // Calculate the inertia tensor of each tetrahedron with resct to the centre of gravity of the cluster
    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        Local[element_counter * 3 + 0]= BARIC[element_counter * 3 + 0]-Xcdgrav;
        Local[element_counter * 3 + 1]= BARIC[element_counter * 3 + 1]-Ycdgrav;
        Local[element_counter * 3 + 2]= BARIC[element_counter * 3 + 2]-Zcdgrav;
        Vnerc[element_counter * 9 + 0]= Vmass[element_counter]*(Local[element_counter * 3 + 1]*Local[element_counter * 3 + 1]+Local[element_counter * 3 + 2]*Local[element_counter * 3 + 2]);
        Vnerc[element_counter * 9 + 1]= -Vmass[element_counter]*Local[element_counter * 3 + 0]*Local[element_counter * 3 + 1];
        Vnerc[element_counter * 9 + 2]= -Vmass[element_counter]*Local[element_counter * 3 + 0]*Local[element_counter * 3 + 2];
        Vnerc[element_counter * 9 + 3]= -Vmass[element_counter]*Local[element_counter * 3 + 1]*Local[element_counter * 3 + 0];
        Vnerc[element_counter * 9 + 4]= Vmass[element_counter]*(Local[element_counter * 3 + 0]*Local[element_counter * 3 + 0]+Local[element_counter * 3 + 2]*Local[element_counter * 3 + 2]);
        Vnerc[element_counter * 9 + 5]= -Vmass[element_counter]*Local[element_counter * 3 + 1]*Local[element_counter * 3 + 2];
        Vnerc[element_counter * 9 + 6]= -Vmass[element_counter]*Local[element_counter * 3 + 2]*Local[element_counter * 3 + 0];
        Vnerc[element_counter * 9 + 7]= -Vmass[element_counter]*Local[element_counter * 3 + 2]*Local[element_counter * 3 + 1];
        Vnerc[element_counter * 9 + 8]= Vmass[element_counter]*(Local[element_counter * 3 + 0]*Local[element_counter * 3 + 0]+Local[element_counter * 3 + 1]*Local[element_counter * 3 + 1]);
    }

    for (int i = 0; i < 9; i++) VNERT[i]=0.0;

    // Calculate the whole inertia tensor
    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        VNERT[0]+=Vnerc[element_counter * 9 + 0];
        VNERT[1]+=Vnerc[element_counter * 9 + 1];
        VNERT[2]+=Vnerc[element_counter * 9 + 2];
        VNERT[3]+=Vnerc[element_counter * 9 + 3];
        VNERT[4]+=Vnerc[element_counter * 9 + 4];
        VNERT[5]+=Vnerc[element_counter * 9 + 5];
        VNERT[6]+=Vnerc[element_counter * 9 + 6];
        VNERT[7]+=Vnerc[element_counter * 9 + 7];
        VNERT[8]+=Vnerc[element_counter * 9 + 8];
    }

    // Calculate the inertias per unit of mass
    Queremos las inercias por unidad de masa, así que dividimos las inercias por la masa total=total_volume*density (density = 1, no la ponemos en la fórmula)
    std::cout << "Inertias: " << VNERT[0]/total_volume << " " << VNERT[1]/total_volume << " " << VNERT[2]/total_volume << '\n';
    std::cout << "Inertias: " << VNERT[3]/total_volume << " " << VNERT[4]/total_volume << " " << VNERT[5]/total_volume << '\n';
    std::cout << "Inertias: " << VNERT[6]/total_volume << " " << VNERT[7]/total_volume << " " << VNERT[8]/total_volume << "\n\n";

    I[0][0] = VNERT[0]/total_volume;
    I[0][1] = VNERT[1]/total_volume;
    I[0][2] = VNERT[2]/total_volume;
    I[1][0] = VNERT[3]/total_volume;
    I[1][1] = VNERT[4]/total_volume;
    I[1][2] = VNERT[5]/total_volume;
    I[2][0] = VNERT[6]/total_volume;
    I[2][1] = VNERT[7]/total_volume;
    I[2][2] = VNERT[8]/total_volume;

    Diagonalize(I, Q, D);

    std::cout << "Eigenvectors: " << Q[0][0] << " " << Q[0][1] << " " << Q[0][2] << '\n';
    std::cout << "Eigenvectors: " << Q[1][0] << " " << Q[1][1] << " " << Q[1][2] << '\n';
    std::cout << "Eigenvectors: " << Q[2][0] << " " << Q[2][1] << " " << Q[2][2] << "\n\n";

    std::cout << "Eigenvalues: " << D[0][0] << " " << D[0][1] << " " << D[0][2] << '\n';
    std::cout << "Eigenvalues: " << D[1][0] << " " << D[1][1] << " " << D[1][2] << '\n';
    std::cout << "Eigenvalues: " << D[2][0] << " " << D[2][1] << " " << D[2][2] << "\n\n";

    // Rotate the object and place it parallel to its principal axis of inertia
    for (int node_counter = 0; node_counter < NUM_OF_NODES; node_counter++) {

        double temporal_array[3][1];
        temporal_array[0][0] = Q[0][0] * tcoord[node_counter * 3 + 0] + Q[1][0] * tcoord[node_counter * 3 + 0] + Q[2][0] * tcoord[node_counter * 3 + 0];
        temporal_array[1][0] = Q[0][1] * tcoord[node_counter * 3 + 1] + Q[1][1] * tcoord[node_counter * 3 + 1] + Q[2][1] * tcoord[node_counter * 3 + 1];
        temporal_array[2][0] = Q[0][2] * tcoord[node_counter * 3 + 2] + Q[1][2] * tcoord[node_counter * 3 + 2] + Q[2][2] * tcoord[node_counter * 3 + 2];
        tcoord[node_counter * 3 + 0] = temporal_array[0][0];
        tcoord[node_counter * 3 + 1] = temporal_array[1][0];
        tcoord[node_counter * 3 + 2] = temporal_array[2][0];
    }

    // Rotate the spheres accordingly
    for (int spheres_counter = 0; spheres_counter < NUM_OF_SPHERES; spheres_counter++) {

        double temporal_array_sph[3][1];
        temporal_array_sph[0][0] = Q[0][0] * sphcoord[spheres_counter * 3 + 0] + Q[1][0] * sphcoord[spheres_counter * 3 + 1] + Q[2][0] * sphcoord[spheres_counter * 3 + 2];
        temporal_array_sph[1][0] = Q[0][1] * sphcoord[spheres_counter * 3 + 0] + Q[1][1] * sphcoord[spheres_counter * 3 + 1] + Q[2][1] * sphcoord[spheres_counter * 3 + 2];
        temporal_array_sph[2][0] = Q[0][2] * sphcoord[spheres_counter * 3 + 0] + Q[1][2] * sphcoord[spheres_counter * 3 + 1] + Q[2][2] * sphcoord[spheres_counter * 3 + 2];
        sphcoord[spheres_counter * 3 + 0] = temporal_array_sph[0][0];
        sphcoord[spheres_counter * 3 + 1] = temporal_array_sph[1][0];
        sphcoord[spheres_counter * 3 + 2] = temporal_array_sph[2][0];
    }

    // Calculate the smallest sphere that circumscribe the cluster (necessary for the cluster mesher)
    double Distance, CenterX, CenterY, CenterZ, Radius = 0.0;
    int extreme_sphere_1, extreme_sphere_2;
    std::vector<int> extreme_sphere;
    for (int spheres_counter_1 = 0; spheres_counter_1 < NUM_OF_SPHERES; spheres_counter_1++) {
        for (int spheres_counter_2 = 0; spheres_counter_2 < NUM_OF_SPHERES; spheres_counter_2++) {
            if (spheres_counter_2 < spheres_counter_1) {
                Distance  = sqrt((sphcoord[spheres_counter_2 * 3 + 0] - sphcoord[spheres_counter_1 * 3 + 0]) * (sphcoord[spheres_counter_2 * 3 + 0] - sphcoord[spheres_counter_1 * 3 + 0]) + (sphcoord[spheres_counter_2 * 3 + 1] - sphcoord[spheres_counter_1 * 3 + 1]) * (sphcoord[spheres_counter_2 * 3 + 1] - sphcoord[spheres_counter_1 * 3 + 1]) + (sphcoord[spheres_counter_2 * 3 + 2] - sphcoord[spheres_counter_1 * 3 + 2]) * (sphcoord[spheres_counter_2 * 3 + 2] - sphcoord[spheres_counter_1 * 3 + 2])) + sphrad[spheres_counter_2] + sphrad[spheres_counter_1];
                if (Distance > 2 * Radius) {
                    Radius = 0.5 * Distance;
                    extreme_sphere_1 = spheres_counter_1;
                    extreme_sphere_2 = spheres_counter_2;
                }
            }
        }
    }

    extreme_sphere.push_back(extreme_sphere_1);
    extreme_sphere.push_back(extreme_sphere_2);

    double SphDistance[3];
    SphDistance[0] = sphcoord[extreme_sphere[1] * 3 + 0] - sphcoord[extreme_sphere[0] * 3 + 0];
    SphDistance[1] = sphcoord[extreme_sphere[1] * 3 + 1] - sphcoord[extreme_sphere[0] * 3 + 1];
    SphDistance[2] = sphcoord[extreme_sphere[1] * 3 + 2] - sphcoord[extreme_sphere[0] * 3 + 2];

    double SphDistanceNorm = sqrt(SphDistance[0] * SphDistance[0] + SphDistance[1] * SphDistance[1] + SphDistance[2] * SphDistance[2]);

    double SphDistanceUnitVect[3];
    SphDistanceUnitVect[0] = SphDistance[0] / SphDistanceNorm;
    SphDistanceUnitVect[1] = SphDistance[1] / SphDistanceNorm;
    SphDistanceUnitVect[2] = SphDistance[2] / SphDistanceNorm;

    CenterX = 0.5 * (sphcoord[extreme_sphere[0] * 3 + 0] - sphrad[extreme_sphere[0]] * SphDistanceUnitVect[0] + sphcoord[extreme_sphere[1] * 3 + 0] + sphrad[extreme_sphere[1]] * SphDistanceUnitVect[0]);
    CenterY = 0.5 * (sphcoord[extreme_sphere[0] * 3 + 1] - sphrad[extreme_sphere[0]] * SphDistanceUnitVect[1] + sphcoord[extreme_sphere[1] * 3 + 1] + sphrad[extreme_sphere[1]] * SphDistanceUnitVect[1]);
    CenterZ = 0.5 * (sphcoord[extreme_sphere[0] * 3 + 2] - sphrad[extreme_sphere[0]] * SphDistanceUnitVect[2] + sphcoord[extreme_sphere[1] * 3 + 2] + sphrad[extreme_sphere[1]] * SphDistanceUnitVect[2]);

    double CheckX, CheckY, CheckZ, TempRad;
    std::vector<double> extreme_radius;
    int sphere_counter_check = 0;
    double CheckDistance[3];
    double CheckDistanceNorm;
    double CheckDistanceUnitVect[3];
    double CheckRadius[3];
    double CheckRadiusNorm;
    double CheckRadiusUnitVect[3];

    while (sphere_counter_check < (NUM_OF_SPHERES)) {  //CHECK that all nodes are inside the big sphere
        CheckDistance[0] = sphcoord[sphere_counter_check * 3 + 0] - CenterX;
        CheckDistance[1] = sphcoord[sphere_counter_check * 3 + 1] - CenterY;
        CheckDistance[2] = sphcoord[sphere_counter_check * 3 + 2] - CenterZ;

        CheckDistanceNorm = sqrt(CheckDistance[0] * CheckDistance[0] + CheckDistance[1] * CheckDistance[1] + CheckDistance[2] * CheckDistance[2]);

        CheckDistanceUnitVect[0] = CheckDistance[0] / CheckDistanceNorm;
        CheckDistanceUnitVect[1] = CheckDistance[1] / CheckDistanceNorm;
        CheckDistanceUnitVect[2] = CheckDistance[2] / CheckDistanceNorm;

        CheckX = sphcoord[sphere_counter_check * 3 + 0] + sphrad[sphere_counter_check] * CheckDistanceUnitVect[0] - CenterX;
        CheckY = sphcoord[sphere_counter_check * 3 + 1] + sphrad[sphere_counter_check] * CheckDistanceUnitVect[1] - CenterY;
        CheckZ = sphcoord[sphere_counter_check * 3 + 2] + sphrad[sphere_counter_check] * CheckDistanceUnitVect[2] - CenterZ;
        TempRad  = sqrt(CheckX * CheckX + CheckY * CheckY + CheckZ * CheckZ);

        if (TempRad - Radius > 1.0e-15) {
            extreme_sphere.push_back(sphere_counter_check);
            CenterX = CenterX + CheckX - (Radius * (CheckX/TempRad));
            CenterY = CenterY + CheckY - (Radius * (CheckY/TempRad));
            CenterZ = CenterZ + CheckZ - (Radius * (CheckZ/TempRad));

            extreme_radius.clear();

            for(int i = 0; i < extreme_sphere.size(); i++) {
                CheckRadius[0] = sphcoord[sphere_counter_check * 3 + 0] - CenterX;
                CheckRadius[1] = sphcoord[sphere_counter_check * 3 + 1] - CenterY;
                CheckRadius[2] = sphcoord[sphere_counter_check * 3 + 2] - CenterZ;

                CheckRadiusNorm = sqrt(CheckRadius[0] * CheckRadius[0] + CheckRadius[1] * CheckRadius[1] + CheckRadius[2] * CheckRadius[2]);

                CheckRadiusUnitVect[0] = CheckRadius[0] / CheckRadiusNorm;
                CheckRadiusUnitVect[1] = CheckRadius[1] / CheckRadiusNorm;
                CheckRadiusUnitVect[2] = CheckRadius[2] / CheckRadiusNorm;

                extreme_radius.push_back(sqrt((sphcoord[extreme_sphere[i] * 3 + 0] + sphrad[extreme_sphere[i]] * CheckRadiusUnitVect[0] - CenterX) * (sphcoord[extreme_sphere[i] * 3 + 0] + sphrad[extreme_sphere[i]] * CheckRadiusUnitVect[0] - CenterX) + (sphcoord[extreme_sphere[i] * 3 + 1] + sphrad[extreme_sphere[i]] * CheckRadiusUnitVect[1] - CenterX) * (sphcoord[extreme_sphere[i] * 3 + 1] + sphrad[extreme_sphere[i]] * CheckRadiusUnitVect[1] - CenterX) + (sphcoord[extreme_sphere[i] * 3 + 2] + sphrad[extreme_sphere[i]] * CheckRadiusUnitVect[2] - CenterX) * (sphcoord[extreme_sphere[i] * 3 + 2] + sphrad[extreme_sphere[i]] * CheckRadiusUnitVect[2] - CenterX)));
            }

            for(int i = 0; i < extreme_radius.size(); i++) {
                if (extreme_radius[i] > Radius) {
                    Radius = extreme_radius[i];
                }
            }
            sphere_counter_check = 0;
        }
        else {
            sphere_counter_check++;
        }
    }

    double diameter = 2 * Radius;

    double characteristic_size = 2 * std::pow(3 * total_volume / (4 * M_PI), 1/3.);

    std::cout << "\nThe diameter is: " << diameter << "\n\n";

    std::cout << "\nThe characteristic size is: " << characteristic_size << "\n\n";

    std::cout << "\nThe size ratio is: " << diameter/characteristic_size << "\n\n";

    // Create the cluster (.clu) file
    std::ofstream outputfile("file_name.clu", std::ios_base::out);
    outputfile << "//\n//   Cluster Name: \"Cluster name\"\n";
    outputfile << "//   Author: Author Name\n";
    outputfile << "//   Date:   YYYY-MM-DD\n//\n\n";

    outputfile << "Name\nCluster name\n\nBegin centers_and_radii\n";
    for (int spheres_counter = 0; spheres_counter < NUM_OF_SPHERES; spheres_counter++) {
        outputfile << sphcoord[spheres_counter * 3 + 0] << " " << sphcoord[spheres_counter * 3 + 1] << " " << sphcoord[spheres_counter * 3 + 2] << " " << sphrad[spheres_counter] << '\n';
    }
    outputfile << "End centers_and_radii\n\n";
    outputfile << "Particle_center_and_diameter\n" << CenterX << " " << CenterY << " " << CenterZ << " " << diameter << "\n\n";
    outputfile << "Size\n" << diameter << "\n\n" << "Volume\n" << total_volume << "\n\n";
    outputfile << "Inertia per unit mass\n" << D[0][0] << '\n' << D[1][1] << '\n' << D[2][2] << '\n';

    outputfile.close();

    delete[] tcoord;
    delete[] Nconec;
    delete[] Vnerc;
    delete[] Volum;
    delete[] Vmass;
    delete[] BARIC;
    delete[] Local;
    delete[] sphcoord;
}

// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization
// source: http://www.melax.com/diag.html?attredirects=0

void Diagonalize(const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3]) {

    // 'A' must be a symmetric matrix.
    // Returns Q and D such that diagonal matrix D = QT * A * Q;  and  A = Q*D*QT

    const int maxsteps=24;  // certainly won't need that many.
    int k0, k1, k2;
    double o[3], m[3];
    double q [4] = {0.0,0.0,0.0,1.0};
    double jr[4];
    double sqw, sqx, sqy, sqz;
    double tmp1, tmp2, mq;
    double AQ[3][3];
    double thet, sgn, t, c;
    for (int i = 0; i < maxsteps; ++i) {

        // quat to matrix
        sqx     = q[0]*q[0];
        sqy     = q[1]*q[1];
        sqz     = q[2]*q[2];
        sqw     = q[3]*q[3];
        Q[0][0] = ( sqx - sqy - sqz + sqw);
        Q[1][1] = (-sqx + sqy - sqz + sqw);
        Q[2][2] = (-sqx - sqy + sqz + sqw);
        tmp1    = q[0]*q[1];
        tmp2    = q[2]*q[3];
        Q[1][0] = 2.0 * (tmp1 + tmp2);
        Q[0][1] = 2.0 * (tmp1 - tmp2);
        tmp1    = q[0]*q[2];
        tmp2    = q[1]*q[3];
        Q[2][0] = 2.0 * (tmp1 - tmp2);
        Q[0][2] = 2.0 * (tmp1 + tmp2);
        tmp1    = q[1]*q[2];
        tmp2    = q[0]*q[3];
        Q[2][1] = 2.0 * (tmp1 + tmp2);
        Q[1][2] = 2.0 * (tmp1 - tmp2);

        // AQ = A * Q
        AQ[0][0] = Q[0][0]*A[0][0]+Q[1][0]*A[0][1]+Q[2][0]*A[0][2];
        AQ[0][1] = Q[0][1]*A[0][0]+Q[1][1]*A[0][1]+Q[2][1]*A[0][2];
        AQ[0][2] = Q[0][2]*A[0][0]+Q[1][2]*A[0][1]+Q[2][2]*A[0][2];
        AQ[1][0] = Q[0][0]*A[0][1]+Q[1][0]*A[1][1]+Q[2][0]*A[1][2];
        AQ[1][1] = Q[0][1]*A[0][1]+Q[1][1]*A[1][1]+Q[2][1]*A[1][2];
        AQ[1][2] = Q[0][2]*A[0][1]+Q[1][2]*A[1][1]+Q[2][2]*A[1][2];
        AQ[2][0] = Q[0][0]*A[0][2]+Q[1][0]*A[1][2]+Q[2][0]*A[2][2];
        AQ[2][1] = Q[0][1]*A[0][2]+Q[1][1]*A[1][2]+Q[2][1]*A[2][2];
        AQ[2][2] = Q[0][2]*A[0][2]+Q[1][2]*A[1][2]+Q[2][2]*A[2][2];

        // D = Qt * AQ
        D[0][0] = AQ[0][0]*Q[0][0]+AQ[1][0]*Q[1][0]+AQ[2][0]*Q[2][0];
        D[0][1] = AQ[0][0]*Q[0][1]+AQ[1][0]*Q[1][1]+AQ[2][0]*Q[2][1];
        D[0][2] = AQ[0][0]*Q[0][2]+AQ[1][0]*Q[1][2]+AQ[2][0]*Q[2][2];
        D[1][0] = AQ[0][1]*Q[0][0]+AQ[1][1]*Q[1][0]+AQ[2][1]*Q[2][0];
        D[1][1] = AQ[0][1]*Q[0][1]+AQ[1][1]*Q[1][1]+AQ[2][1]*Q[2][1];
        D[1][2] = AQ[0][1]*Q[0][2]+AQ[1][1]*Q[1][2]+AQ[2][1]*Q[2][2];
        D[2][0] = AQ[0][2]*Q[0][0]+AQ[1][2]*Q[1][0]+AQ[2][2]*Q[2][0];
        D[2][1] = AQ[0][2]*Q[0][1]+AQ[1][2]*Q[1][1]+AQ[2][2]*Q[2][1];
        D[2][2] = AQ[0][2]*Q[0][2]+AQ[1][2]*Q[1][2]+AQ[2][2]*Q[2][2];
        o[0] = D[1][2];
        o[1] = D[0][2];
        o[2] = D[0][1];
        m[0] = fabs(o[0]);
        m[1] = fabs(o[1]);
        m[2] = fabs(o[2]);

        k0 = (m[0] > m[1] && m[0] > m[2])?0: (m[1] > m[2])? 1 : 2; // index of largest element of offdiag
        k1 = (k0+1)%3;
        k2 = (k0+2)%3;

        if (o[k0] == 0.0) break;  // Diagonal already

        thet  = (D[k2][k2]-D[k1][k1])/(2.0*o[k0]);
        sgn   = (thet > 0.0)?1.0:-1.0;
        thet *= sgn; // make it positive
        t     = sgn /(thet +((thet < 1.E6)?sqrt(thet*thet+1.0):thet)) ; // sign(T)/(|T|+sqrt(T^2+1))
        c     = 1.0/sqrt(t*t+1.0); //  c= 1/(t^2+1) , t=s/c

        if (c == 1.0) break;  // No room for improvement, reached machine precision.

        jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
        jr[k0]  = sgn * sqrt((1.0-c)/2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
        jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
        jr[3 ]  = sqrt(1.0f - jr[k0] * jr[k0]);

        if (jr[3] == 1.0) break; // Reached limits of floating point precision

        q[0]  = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
        q[1]  = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
        q[2]  = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
        q[3]  = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
        mq    = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
        q[0] /= mq;
        q[1] /= mq;
        q[2] /= mq;
        q[3] /= mq;
    }
}
