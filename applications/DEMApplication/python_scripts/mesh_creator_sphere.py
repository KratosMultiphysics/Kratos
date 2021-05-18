import os

import KratosMultiphysics

def WriteSphereMdpaFromResults(filename_pre, pre_path, filename_post, spheres_model_part):

    SpheresMdpa_pre = open(os.path.join(pre_path, filename_pre) + ".mdpa", 'r')
    SpheresMdpa_post = open(os.path.join(pre_path, filename_post) + ".mdpa", 'w')

    for Line in SpheresMdpa_pre:
        # The next 4 lines mean that, while we do not get to 'Begin Nodes', we simply copy the mdpa contents, in this case the Properties part
        SpheresMdpa_post.write(Line)
        if Line.startswith('Begin Nodes'):
            break

    for node in spheres_model_part.Nodes:
        SpheresMdpa_post.write(str(node.Id) + ' ' + str(node.X) + ' ' + str(node.Y) + ' ' + str(node.Z) + '\n')
    SpheresMdpa_post.write('End Nodes\n\n')

    for Line in SpheresMdpa_pre:
        # We copy the element type line from the previous mdpa
        if Line.startswith('Begin Elements'):
            SpheresMdpa_post.write(Line)
            break

    for element in spheres_model_part.Elements:
        SpheresMdpa_post.write(str(element.Id) + ' ' + str(element.Properties.Id) + ' ' + str(element.GetNode(0).Id) + '\n')
    SpheresMdpa_post.write('End Elements\n')

    SpheresMdpa_post.write('\nBegin NodalData RADIUS\n')
    for node in spheres_model_part.Nodes:
        SpheresMdpa_post.write(str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosMultiphysics.RADIUS)) + '\n')
    SpheresMdpa_post.write('End NodalData\n')

    SpheresMdpa_post.close()
    SpheresMdpa_pre.close()
