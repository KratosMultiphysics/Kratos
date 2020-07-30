import os

def WriteSphereMdpaFromResults(filename_pre, pre_path, filename_post, filename_msh, filename_res, post_path):

    SpheresMdpa_pre = open(os.path.join(pre_path, filename_pre) + ".mdpa", 'r')
    msh = open(os.path.join(post_path, filename_msh), 'r')
    rad = open(os.path.join(post_path, filename_msh), 'r')
    coh = open(os.path.join(post_path, filename_msh), 'r')
    mod = open(os.path.join(post_path, filename_msh), 'r')
    ski = open(os.path.join(post_path, filename_res), 'r')
    SpheresMdpa_post = open(os.path.join(pre_path, filename_post) + ".mdpa", 'w')

    for Line in SpheresMdpa_pre:
        # The next 4 lines mean that, while we do not get to 'Begin Nodes', we simply copy the mdpa contents, in this case the Properties part
        if Line.startswith("Begin Nodes"):
            SpheresMdpa_pre.close()
            break
        SpheresMdpa_post.write(Line)

    flag_msh = 0
    for Line in msh:
        if 'ElemType Sphere Nnode 1' in Line:
            flag_msh = 1
        if 'Coordinates' in Line and flag_msh == 1:
            flag_msh = 2
            SpheresMdpa_post.write('Begin Nodes // GUI group identifier: dems celemid CylinderContinuumParticle2D\n')
            continue
        if 'Coordinates' in Line and flag_msh == 2:
            flag_msh = 3
            SpheresMdpa_post.write('End Nodes\n')
            continue
        if 'Elements' in Line and flag_msh == 3:
            SpheresMdpa_post.write('\nBegin Elements CylinderContinuumParticle2D   //  GUI group identifier: dems\n')
            flag_msh = 4
            continue
        if 'Elements' in Line and flag_msh == 4:
            SpheresMdpa_post.write('End Elements\n')
            msh.close()
            break

        if flag_msh == 2:
            SpheresMdpa_post.write(Line)
        if flag_msh == 4:
            Line = Line.strip('\n') # Remove the line-ending characters
            ElementList = Line.split(' ')
            SpheresMdpa_post.write(ElementList[0] + ' 1 ' + ElementList[1] + '\n')

    flag_rad = 0
    for Line in rad:
        if 'Elements' in Line and flag_rad == 0:
            SpheresMdpa_post.write('\nBegin NodalData RADIUS  // GUI group identifier: dens Elementid CylinderContinuumParticle2D\n')
            flag_rad = 1
            continue
        if 'Elements' in Line and flag_rad == 1:
            SpheresMdpa_post.write('End NodalData\n')
            rad.close()
            break

        if flag_rad == 1:
            Line = Line.strip('\n') # Remove the line-ending characters
            ElementList = Line.split(' ')
            SpheresMdpa_post.write(ElementList[1] + ' 0 ' + ElementList[2] + '\n')

    flag_coh = 0
    for Line in coh:
        if 'Elements' in Line and flag_coh == 0:
            SpheresMdpa_post.write('\nBegin NodalData COHESIVE_GROUP  // GUI group identifier: dems Elementid CylinderContinuumParticle2D\n')
            flag_coh = 1
            continue
        if 'Elements' in Line and flag_coh == 1:
            SpheresMdpa_post.write('End NodalData\n')
            coh.close()
            break

        if flag_coh == 1:
            Line = Line.strip('\n') # Remove the line-ending characters
            ElementList = Line.split(' ')
            SpheresMdpa_post.write(ElementList[1] + ' 0 1' + '\n')

    flag_ski = 0
    for Line in ski:
        if 'SKIN_SPHERE' in Line and flag_ski == 0:
            SpheresMdpa_post.write('\nBegin NodalData SKIN_SPHERE  // GUI group identifier: dems Elementid CylinderContinuumParticle2D\n')
            flag_ski = 1
            continue
        if 'Values' in Line and flag_ski == 1:
            flag_ski = 2
            continue
        if 'Values' in Line and flag_ski == 2:
            SpheresMdpa_post.write('End NodalData\n')
            ski.close()
            break
        if flag_ski == 2:
            Line = Line.strip('\n') # Remove the line-ending characters
            ElementList = Line.split(' ')
            SpheresMdpa_post.write(ElementList[0] + ' 0 ' + ElementList[1] + '\n')

    flag_mod = 0
    for Line in mod:
        if 'ElemType Sphere Nnode 1' in Line:
            flag_mod = 1
        if 'Coordinates' in Line and flag_mod == 1:
            flag_mod = 2
            SpheresMdpa_post.write('\nBegin SubModelPart Parts_dems // Group dems // Subtree Parts\nBegin SubModelPartNodes\n')
            continue
        if 'Coordinates' in Line and flag_mod == 2:
            flag_mod = 3
            SpheresMdpa_post.write('End SubModelPartNodes\n')
            continue
        if 'Elements' in Line and flag_mod == 3:
            SpheresMdpa_post.write('\nBegin SubModelPartElements\n')
            flag_mod = 4
            continue
        if 'Elements' in Line and flag_mod == 4:
            SpheresMdpa_post.write('End SubModelPartElements\nBegin SubModelPartConditions\nEnd SubModelPartConditions\nEnd SubModelPart\n')
            mod.close()
            break

        if flag_mod == 2:
            Line = Line.strip('\n') # Remove the line-ending characters
            NodeList = Line.split(' ')
            SpheresMdpa_post.write(NodeList[0] + '\n')
        if flag_mod == 4:
            Line = Line.strip('\n') # Remove the line-ending characters
            ElementList = Line.split(' ')
            SpheresMdpa_post.write(ElementList[0] + '\n')

    SpheresMdpa_post.close()
