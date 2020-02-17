def SurfaceFromNodeSet(assembly,instance,nodeSetName,surfaceName):
	faces = {}
	nodes = assembly.sets[nodeSetName].nodes
	dimension = GetDimension(nodes[0].getElements()[0])
	if dimension == 2:
		for node in nodes:
			elemEdges = node.getElemEdges()
			for elemEdge in elemEdges:
				edgeElements = elemEdge.getElements()
				if len(edgeElements) == 1:
					elemLabel = edgeElements[0].label
					edgeIndex = edgeElements[0].getElemEdges().index(elemEdge)
					key = (elemLabel, edgeIndex)
					if key in faces:
						faces[key] = (faces[key][0]+1, faces[key][1])
					else:
						nodesPerFace = GetNodesPerFace(instance.elements.getFromLabel(elemLabel))
						faces[key] = (1, nodesPerFace)
	else:
		for node in nodes:
			elemFaces = node.getElemFaces()
			for elemFace in elemFaces:
				elemLabel = elemFace.label
				faceIndex = int(elemFace.face.getText()[-1:])-1
				key = (elemLabel, faceIndex)
				if key in faces:
					faces[key] = (faces[key][0]+1, faces[key][1])
				else:
					nodesPerFace = GetNodesPerFace(instance.elements.getFromLabel(elemLabel))
					faces[key] = (1, nodesPerFace)
	faceElements = [[] for face in range(6)]
	for key, face in faces.iteritems():
		if face[0] == face[1]:
			faceElements[key[1]].append(key[0])
	return assembly.Surface(name = surfaceName,
		face1Elements = instance.elements.sequenceFromLabels(faceElements[0]),
		face2Elements = instance.elements.sequenceFromLabels(faceElements[1]),
		face3Elements = instance.elements.sequenceFromLabels(faceElements[2]),
		face4Elements = instance.elements.sequenceFromLabels(faceElements[3]),
		face5Elements = instance.elements.sequenceFromLabels(faceElements[4]),
		face6Elements = instance.elements.sequenceFromLabels(faceElements[5]))

def GetDimension(element):
	elementTypeText = element.type.getText()
	if elementTypeText.find('CAX4') > -1:
		dimension = 2
	elif elementTypeText.find('CAX8') > -1:
		dimension = 2
	elif elementTypeText.find('CPE8') > -1:
		dimension = 2
	elif elementTypeText.find('CPS8') > -1:
		dimension = 2
	elif elementTypeText.find('C3D8') > -1:
		dimension = 3
	elif elementTypeText.find('C3D20') > -1:
		dimension = 3
	else:
		raise Exception('Unknown element type ',elementType)
	return dimension

def GetNodesPerFace(element):
	elementTypeText = element.type.getText()
	if elementTypeText.find('CAX4') > -1:
		nodesPerFace = 2
	elif elementTypeText.find('CAX8') > -1:
		nodesPerFace = 3
	elif elementTypeText.find('CPE8') > -1:
		nodesPerFace = 3
	elif elementTypeText.find('CPS8') > -1:
		nodesPerFace = 3
	elif elementTypeText.find('C3D8') > -1:
		nodesPerFace = 4
	elif elementTypeText.find('C3D20') > -1:
		nodesPerFace = 8
	else:
		raise Exception('Unknown element type ',elementType)
	return nodesPerFace

