fo = open(file="SteadyStatePipeElementWithEmbankment.mdpa", mode="r")
fs = open(file="SteadyStatePipeElementWithEmbankment_repeated_nodes.mdpa", mode="w")

Nodes = False
Elements = False
subparts = False
MappedIndex = []
NodesCoords = []

for line in fo:
    if line.startswith("Begin Nodes"):
        Nodes = True
        fs.write(line)
        continue
    elif line.startswith("End Nodes"):
        Nodes = False

        for ind, NodeCoords in enumerate(NodesCoords):
            fs.write(str(ind+1) + " " + " ".join(NodeCoords) + "\n")

        fs.write(line)
        continue
    elif line.startswith("Begin Elements"):
        Elements = True
        fs.write(line)
        continue
    elif line.startswith("End Elements"):
        Elements = False
        fs.write(line)
        continue
   
    elif line.startswith("  Begin SubModelPartNodes"):
        subparts = True
        fs.write(line)
        continue
    elif line.startswith("  End SubModelPartNodes"):
        subparts = False
        fs.write(line)
        continue


    elif Nodes:
        index, x, y, z = line.split()
        if [x, y, z] not in NodesCoords:
            NodesCoords.append([x, y, z])
            MappedIndex.append([index])
        else:
            ind = NodesCoords.index([x, y, z])
            MappedIndex[ind].append(index)

    elif Elements:
        data = line.split()
        ind = data[0]
        mat = data[1]
        connectivity = data[2:]
        newConnectivity = []
        for nodeInd in connectivity:
            for mapIndex, nodesIndices in enumerate(MappedIndex):
                if nodeInd in nodesIndices:
                    newConnectivity.append(str(mapIndex+1))
                    break

        fs.write(ind + " " + mat + " " + " ".join(newConnectivity)+"\n")
    elif subparts:
        nodeNo = int(line)
        for mapIndex, nodesIndices in enumerate(MappedIndex):
            if str(nodeNo) in nodesIndices:
                fs.write(str(mapIndex+1)+"\n")
                break
    else:
        fs.write(line)