/** Opens an output database and writes the node positions to a file. The path of the output databate must be given as an argument. The database is opened and the last frame from the last step is analyzed. */

//
// System includes
//
#include <stdlib.h>
#include <stdio.h>
#if (defined(HP) && (! defined(HKS_HPUXI)))
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#else
#include <iostream>
#include <iomanip>
#include <fstream>
#endif

// Begin local includes
#include <odb_API.h>
// End local includes

// Declare local function prototypes
int getOutput(int argc, char **argv);

int ABQmain(int argc, char **argv)
{
	// Parameters
	unsigned int dataPrecision = 17;
	unsigned int dataWidth = 27;
	unsigned int surfaces = |surfaces|;
	odb_String surfaceIDs[|surfaces|] = {|surfaceIDs|};
	
	// Read the output database and extract the nodal coordinates
	odb_String path;
	unsigned int frameIndex;
	if (argc == 3) {
		path = argv[1]+odb_String(".odb");
		frameIndex = atoi(argv[2]);
	} else {
		std::cout << "Arguments are filename and frame" << std::endl;
		exit(1);
	}
	odb_Odb& odb = openOdb(path);
	odb_Step& step = odb.steps()["Step-1"];
	odb_Frame& frame = step.frames()[(frameIndex ? step.frames().size()-1 : 0)];
	odb_Assembly& assembly = odb.rootAssembly();
	odb_FieldOutput& fieldOutput = frame.fieldOutputs()[(frameIndex ? "U" : "COORD")];
	for (unsigned int surface = 0; surface < surfaces; surface++) {
		odb_Set location = assembly.nodeSets()[surfaceIDs[surface]];
		odb_FieldOutput surfaceFieldOutput = fieldOutput.getSubset(location);
		const odb_SequenceFieldValue& positions = surfaceFieldOutput.values();
		odb_String outputPath = argv[1]+odb_String("Surface")+odb_String(surface)+odb_String("Output.dat");
		std::ofstream outputFile(outputPath.CStr());
		outputFile.precision(dataPrecision);
		outputFile.fill(' '); 
		outputFile.setf(std::ios_base::scientific, std::ios_base::floatfield);
#if (|dimension|-2)
		outputFile << "x-coordinate\ty-coordinate\tz-coordinate";
#else
		outputFile << "x-coordinate\ty-coordinate";
#endif
		for (unsigned int i = 0; i < positions.size(); i++) {
			const odb_FieldValue val = positions[i];
			int numComp = 0;
			const double* const pos = val.dataDouble(numComp);
			outputFile << "\n";
			for (unsigned int comp = 0; comp < numComp; comp++) {
				outputFile.width(dataWidth);
				outputFile << pos[comp];
			}
		}
		outputFile << std::endl;
	}
	return (0);
}
