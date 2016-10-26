// Translated to C++ from the original work in FORTRAN by Alberto Férriz, 2016

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

void Diagonalize(const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3]);

int main() {
    std::ifstream infile("cubo1x2x3_debug.msh");
    std::string line;
    infile.ignore(80,'\n'); infile.ignore(80,'\n');

    // Hacer una primera pasada para saber NUM_OF_NODES y NUM_OF_ELEMENTS
    // Pasar el nombre del caso como argumento de un script de python y crear el fichero clu
    // Poner esto el PreUtilities?
    // Salta un segmentation fault al reservar la memoria de los arrays para grandes tamaños
    // Debería poderse optimizar el uso de la memoria para evitar eso
    // Primera aproximación de Size: Longitud de la diagonal máxima del prisma contenedor
    // Mover el centroide del objeto al origen si no lo está ya
    // Para el clu, se necesita la inercia por unidad de masa
    // Entonces, mover el objeto al origen y girarlo tal que esté en sus ejes principales
    
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

    infile.seekg(0, std::ios::beg);
    infile.ignore(80,'\n'); infile.ignore(80,'\n');
    
    const int NUM_OF_NODES = node_counter;
    const int NUM_OF_ELEMENTS = element_counter;
    const double density = 1;
    double tcoord[3][NUM_OF_NODES];
    int Nconec[4][NUM_OF_ELEMENTS];
    double Vnerc[9][NUM_OF_ELEMENTS];
    double Volum[NUM_OF_ELEMENTS];
    double Vmass[NUM_OF_ELEMENTS];
    double BARIC[3][NUM_OF_ELEMENTS];
    double Xlocal[NUM_OF_ELEMENTS];
    double Ylocal[NUM_OF_ELEMENTS];
    double Zlocal[NUM_OF_ELEMENTS];
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
            tcoord[0][node_counter] = b;
            tcoord[1][node_counter] = c;
            tcoord[2][node_counter] = d;
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
            Nconec[0][element_counter] = b;
            Nconec[1][element_counter] = c;
            Nconec[2][element_counter] = d;
            Nconec[3][element_counter] = e;
            element_counter++;
        } else break;
    }

    double total_volume = 0.0;

    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        int NodeA=Nconec[0][element_counter];
        int NodeB=Nconec[1][element_counter];
        int NodeC=Nconec[2][element_counter];
        int NodeD=Nconec[3][element_counter];

        double ValU1=tcoord[0][NodeB-1] - tcoord[0][NodeA-1];
        double ValU2=tcoord[1][NodeB-1] - tcoord[1][NodeA-1];
        double ValU3=tcoord[2][NodeB-1] - tcoord[2][NodeA-1];
        double ValV1=tcoord[0][NodeC-1] - tcoord[0][NodeA-1];
        double ValV2=tcoord[1][NodeC-1] - tcoord[1][NodeA-1];
        double ValV3=tcoord[2][NodeC-1] - tcoord[2][NodeA-1];
        double ValW1=tcoord[0][NodeD-1] - tcoord[0][NodeA-1];
        double ValW2=tcoord[1][NodeD-1] - tcoord[1][NodeA-1];
        double ValW3=tcoord[2][NodeD-1] - tcoord[2][NodeA-1];

        double Volu0= ValU1*ValV2*ValW3+ValU2*ValV3*ValW1+ValU3*ValW2*ValV1-ValU3*ValV2*ValW1-ValU2*ValV1*ValW3-ValU1*ValW2*ValV3;
        Volum[element_counter]=Volu0/6;
        total_volume += Volum[element_counter];
        Vmass[element_counter] = Volum[element_counter]*density;

        //Cálculo del baricentro de cada tetraedro
        BARIC[0][element_counter]= (tcoord[0][NodeA-1]+tcoord[0][NodeB-1]+tcoord[0][NodeC-1]+tcoord[0][NodeD-1])/4;
        BARIC[1][element_counter]= (tcoord[1][NodeA-1]+tcoord[1][NodeB-1]+tcoord[1][NodeC-1]+tcoord[1][NodeD-1])/4;
        BARIC[2][element_counter]= (tcoord[2][NodeA-1]+tcoord[2][NodeB-1]+tcoord[2][NodeC-1]+tcoord[2][NodeD-1])/4;
    }

    std::cout << "\nTotal volume: " << total_volume << '\n';

    //Cálculo de la masa total del cluster o piedra
    double Vmaspiedra=0.0;
    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {
        Vmaspiedra+=Vmass[element_counter];
    }
    std::cout << "Total mass: " << Vmaspiedra << '\n';
    
    //Cálculo del centro de gravedad del cluster o piedra
    double Valor1=0.0;
    double Valor2=0.0;
    double Valor3=0.0;

    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        Valor1+=Vmass[element_counter]*BARIC[0][element_counter];
        Valor2+=Vmass[element_counter]*BARIC[1][element_counter];
        Valor3+=Vmass[element_counter]*BARIC[2][element_counter];
    }

    double Xcdgrav=Valor1/Vmaspiedra;
    double Ycdgrav=Valor2/Vmaspiedra;
    double Zcdgrav=Valor3/Vmaspiedra;
    std::cout << "Centroid: " << Xcdgrav << " " << Ycdgrav << " " << Zcdgrav << "\n\n";

    // Movemos el objecto y lo dejamos colocado tal que su centroide coincida con el origen:
    for (int node_counter = 0; node_counter < NUM_OF_NODES; node_counter++) {
        tcoord[0][node_counter] -= Xcdgrav;
        tcoord[1][node_counter] -= Ycdgrav;
        tcoord[2][node_counter] -= Zcdgrav;
    }
    
    // Cálculo del tensor de inercias de cada elemento con respecto al CDG de cada piedra o cluster
    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        Xlocal[element_counter]= BARIC[0][element_counter]-Xcdgrav;
        Ylocal[element_counter]= BARIC[1][element_counter]-Ycdgrav;
        Zlocal[element_counter]= BARIC[2][element_counter]-Zcdgrav;
        Vnerc[0][element_counter]= Vmass[element_counter]*(Ylocal[element_counter]*Ylocal[element_counter]+Zlocal[element_counter]*Zlocal[element_counter]);
        Vnerc[1][element_counter]= -Vmass[element_counter]*Xlocal[element_counter]*Ylocal[element_counter];
        Vnerc[2][element_counter]= -Vmass[element_counter]*Xlocal[element_counter]*Zlocal[element_counter];
        Vnerc[3][element_counter]= -Vmass[element_counter]*Ylocal[element_counter]*Xlocal[element_counter];
        Vnerc[4][element_counter]= Vmass[element_counter]*(Xlocal[element_counter]*Xlocal[element_counter]+Zlocal[element_counter]*Zlocal[element_counter]);
        Vnerc[5][element_counter]= -Vmass[element_counter]*Ylocal[element_counter]*Zlocal[element_counter];
        Vnerc[6][element_counter]= -Vmass[element_counter]*Zlocal[element_counter]*Xlocal[element_counter];
        Vnerc[7][element_counter]= -Vmass[element_counter]*Zlocal[element_counter]*Ylocal[element_counter];
        Vnerc[8][element_counter]= Vmass[element_counter]*(Xlocal[element_counter]*Xlocal[element_counter]+Ylocal[element_counter]*Ylocal[element_counter]);
    }

    for (int i = 0; i < 9; i++) VNERT[i]=0.0;
    
    // Se calcula el tensor de Vnercias totales
    for (int element_counter = 0; element_counter < NUM_OF_ELEMENTS; element_counter++) {

        VNERT[0]+=Vnerc[0][element_counter];
        VNERT[1]+=Vnerc[1][element_counter];
        VNERT[2]+=Vnerc[2][element_counter];
        VNERT[3]+=Vnerc[3][element_counter];
        VNERT[4]+=Vnerc[4][element_counter];
        VNERT[5]+=Vnerc[5][element_counter];
        VNERT[6]+=Vnerc[6][element_counter];
        VNERT[7]+=Vnerc[7][element_counter];
        VNERT[8]+=Vnerc[8][element_counter];
    }

    // Queremos las inercias por unidad de masa, así que dividimos las inercias por la masa total=total_volume*density (density = 1, no la ponemos en la fórmula)
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
    
    // Finalmente, rotamos el objeto y lo colocamos paralelo a sus ejes principales de inercia
    
    for (int node_counter = 0; node_counter < NUM_OF_NODES; node_counter++) {
        
        double temporal_array[3][1];
        temporal_array[0][0] = Q[0][0] * tcoord[0][node_counter] + Q[1][0] * tcoord[1][node_counter] + Q[2][0] * tcoord[2][node_counter];
        temporal_array[1][0] = Q[0][1] * tcoord[0][node_counter] + Q[1][1] * tcoord[1][node_counter] + Q[2][1] * tcoord[2][node_counter];
        temporal_array[2][0] = Q[0][2] * tcoord[0][node_counter] + Q[1][2] * tcoord[1][node_counter] + Q[2][2] * tcoord[2][node_counter];
        tcoord[0][node_counter] = temporal_array[0][0];
        tcoord[1][node_counter] = temporal_array[1][0];
        tcoord[2][node_counter] = temporal_array[2][0];
    }
    
    // Calculamos el Size
    double min_X = 0.0;
    double max_X = 0.0;
    double min_Y = 0.0;
    double max_Y = 0.0;
    double min_Z = 0.0;
    double max_Z = 0.0;

    for (int node_counter = 0; node_counter < NUM_OF_NODES; node_counter++) {
        if (tcoord[0][node_counter] < min_X) min_X = tcoord[0][node_counter];   
        if (tcoord[0][node_counter] > max_X) max_X = tcoord[0][node_counter];
        if (tcoord[1][node_counter] < min_Y) min_Y = tcoord[1][node_counter];   
        if (tcoord[1][node_counter] > max_Y) max_Y = tcoord[1][node_counter];
        if (tcoord[2][node_counter] < min_Z) min_Z = tcoord[2][node_counter];   
        if (tcoord[2][node_counter] > max_Z) max_Z = tcoord[2][node_counter];
    }
    
    std::cout << "The smallest number in X is: " << min_X << '\n';   
    std::cout << "The biggest number in X is: " << max_X << '\n';
    std::cout << "The smallest number in Y is: " << min_Y << '\n';   
    std::cout << "The biggest number in Y is: " << max_Y << '\n';
    std::cout << "The smallest number in Z is: " << min_Z << '\n';   
    std::cout << "The biggest number in Z is: " << max_Z << '\n';
    
    double size = sqrt((max_X - min_X)*(max_X - min_X) + (max_Y - min_Y)*(max_Y - min_Y) + (max_Z - min_Z)*(max_Z - min_Z));
    std::cout << "\nThe size is: " << size << "\n\n";
    
    std::ofstream outputfile("prism1cluster3D.clu", std::ios_base::out);
    outputfile << "//\n//   Cluster Name: \"prism1cluster3D\"\n";
    outputfile << "//   Author: Salva Latorre\n";
    outputfile << "//   Date:   2016-09-27\n//\n\n";
    
    outputfile << "Name\nprism1cluster3D\n\nBegin centers_and_radii\n";
    for (int node_counter = 0; node_counter < NUM_OF_NODES; node_counter++) {
        outputfile << tcoord[0][node_counter] << " " << tcoord[1][node_counter] << " " << tcoord[2][node_counter] << " 0.25" << '\n';
    }
    outputfile << "End centers_and_radii\n\n";
    outputfile << "Size\n" << size << "\n\n" << "Volume\n" << total_volume << "\n\n";
    outputfile << "Inertia per unit mass\n" << D[0][0] << '\n' << D[1][1] << '\n' << D[2][2] << '\n';
 
    outputfile.close();
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
