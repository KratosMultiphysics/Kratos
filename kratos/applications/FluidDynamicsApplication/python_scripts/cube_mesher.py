#!/usr/bin/env python

"""
Generate a structured tetrahedral mesh that fills a box of given size.
Identify the box faces using conditions and mesh groups.

changelog:
 05-12-2013: Initial version (adapted from my turbulent channel mesher)
"""

class box_data(object):

    def __init__(self,xmin,ymin,zmin,xmax,ymax,zmax,nx,ny,nz):
        self.box = [ xmin, ymin, zmin, zmax, ymax, zmax ]
        self.ndiv = [nx,ny,nz]
        self.jump = [ (xmax-xmin)/nx, (ymax-ymin)/ny, (zmax-zmin)/nz ]
        self.cond_range = dict()

    def nx(self):
        return self.ndiv[0]

    def ny(self):
        return self.ndiv[1]

    def nz(self):
        return self.ndiv[2]

    def get_coord(self,direction,n):
        """ Use direction {0,1,2} for {x,y,z}, 0 <= n <= ndiv[direction]."""
        return self.box[direction] + self.jump[direction]*n

    def get_id(self,ix,iy,iz):
        """ Id for the nth node in each direction (starting from 1)."""
        return 1 + ix + (self.nx()+1)*iy + (self.nx()+1)*(self.ny()+1)*iz

    def cube_vertices(self,ix,iy,iz):
        """ Identify 8 contiguous nodes forming a cube.
            Note: Even and odd levels are rotated 90 degrees along Z
            to ensure that they are conformant.
        """
        if iz % 2:
            n0 = box.get_id(ix  ,iy  ,iz)
            n1 = box.get_id(ix+1,iy  ,iz)
            n2 = box.get_id(ix+1,iy+1,iz)
            n3 = box.get_id(ix  ,iy+1,iz)
            
            n4 = box.get_id(ix  ,iy  ,iz+1)
            n5 = box.get_id(ix+1,iy  ,iz+1)
            n6 = box.get_id(ix+1,iy+1,iz+1)
            n7 = box.get_id(ix  ,iy+1,iz+1)
            
        else:
            n1 = box.get_id(ix  ,iy  ,iz)
            n2 = box.get_id(ix+1,iy  ,iz)
            n3 = box.get_id(ix+1,iy+1,iz)
            n0 = box.get_id(ix  ,iy+1,iz)
            
            n5 = box.get_id(ix  ,iy  ,iz+1)
            n6 = box.get_id(ix+1,iy  ,iz+1)
            n7 = box.get_id(ix+1,iy+1,iz+1)
            n4 = box.get_id(ix  ,iy+1,iz+1)

        return n0,n1,n2,n3,n4,n5,n6,n7


def node_x(box,position):
    return box.get_coord(0,position)

def node_y(box,position):
    return box.get_coord(1,position)

def node_z(box,position):
    return box.get_coord(2,position)

def write_header(mdpa):
    mdpa.write("Begin ModelPartData\n")
    mdpa.write("//  VARIABLE_NAME value\n")
    mdpa.write("End ModelPartData\n")
    mdpa.write("\n")
    mdpa.write("Begin Properties 0\n")
    mdpa.write("End Properties\n")
    mdpa.write("\n")

def generate_nodes(mdpa,box,x_scale=node_x,y_scale=node_y,z_scale=node_z):
    index = 1

    print "Generating nodes"
    mdpa.write("Begin Nodes\n")
    
    for iz in range(nz+1):
        z = z_scale(box,iz)
        for iy in range(ny+1):
            y = y_scale(box,iy)
            for ix in range(nx+1):
                x = x_scale(box,ix)
                mdpa.write("{0:d} {1:e} {2:e} {3:e}\n".format(index,x,y,z))
                index += 1
    
    mdpa.write("End Nodes\n\n")

def generate_elements(mpda,box,elemtype="FractionalStep3D",prop_id=0):
    index = 1
    nx = box.nx()
    ny = box.ny()
    nz = box.nz()

    print "Generating tetrahedral {0} elements.".format(elemtype)
    mdpa.write("Begin Elements {0}\n".format(elemtype))
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(ix,iy,iz)

                # fill the cube with 6 tetrahedra
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index  ,prop_id,n3,n4,n6,n7))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index+1,prop_id,n3,n4,n5,n6))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index+2,prop_id,n0,n5,n3,n4))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index+3,prop_id,n0,n1,n3,n5))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index+4,prop_id,n3,n5,n2,n6))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d}\n".format(index+5,prop_id,n3,n1,n2,n5))
                index += 6

    mdpa.write("End Elements\n\n")

def generate_front_faces(mdpa,box,index,condtype="WallCondition3D",prop_id=0):
    ny = box.ny()
    nz = box.nz()
    i0 = index

    mdpa.write("Begin Conditions {0} //Front\n".format(condtype))

    # Front
    for iz in range(nz):
        for iy in range(ny):
            n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(0,iy,iz)
            
            if iz % 2:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n0,n3,n4))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n3,n7,n4))
            else:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n0,n5,n1))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n0,n4,n5))

            index += 2
    
    mdpa.write("End Conditions\n\n")

    box.cond_range["front"] = [i0,index]

    return index

def generate_back_faces(mdpa,box,index,condtype="WallCondition3D",prop_id=0):
    ny = box.ny()
    nz = box.nz()
    i0 = index

    mdpa.write("Begin Conditions {0} //Back\n".format(condtype))

    # Back
    for iz in range(nz):
        for iy in range(ny):
            n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(nx-1,iy,iz)
            
            if iz % 2:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n2,n1,n5))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n2,n5,n6))
            else:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n3,n2,n6))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n3,n6,n7))

            index += 2
   
    mdpa.write("End Conditions\n\n")

    box.cond_range["back"] = [i0,index]

    return index


def generate_right_faces(mdpa,box,index,condtype="WallCondition3D",prop_id=0):
    nx = box.nx()
    nz = box.nz()
    i0 = index

    mdpa.write("Begin Conditions {0} //Right\n".format(condtype))

    # Right
    for iz in range(nz):
        for ix in range(nx):
            n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(ix,0,iz)
            
            if iz % 2:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n0,n5,n1))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n0,n4,n5))
            else:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n2,n5,n6))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n5,n2,n1))

            index += 2
   
    mdpa.write("End Conditions\n\n")

    box.cond_range["right"] = [i0,index]

    return index


def generate_left_faces(mdpa,box,index,condtype="WallCondition3D",prop_id=0):
    nx = box.nx()
    nz = box.nz()
    i0 = index

    mdpa.write("Begin Conditions {0} //Left\n".format(condtype))

    # Left
    for iz in range(nz):
        for ix in range(nx):
            n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(ix,ny-1,iz)
            
            if iz % 2:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n3,n6,n7))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n3,n2,n6))
            else:
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n3,n7,n4))
                mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n3,n4,n0))

            index += 2
   
    mdpa.write("End Conditions\n\n")
    
    box.cond_range["left"] = [i0,index]

    return index


def generate_bottom_faces(mdpa,box,index,condtype="WallCondition3D",prop_id=0):
    nx = box.nx()
    ny = box.ny()
    i0 = index

    mdpa.write("Begin Conditions {0} //Bottom\n".format(condtype))

    # Bottom
    for iy in range(ny):
        for ix in range(nx):
            n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(ix,iy,0)

            mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n0,n1,n3))
            mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n3,n1,n2))
            index += 2
   
    mdpa.write("End Conditions\n\n")

    box.cond_range["bottom"] = [i0,index]

    return index


def generate_top_faces(mdpa,box,index,condtype="WallCondition3D",prop_id=0):
    nx = box.nx()
    ny = box.ny()
    i0 = index

    mdpa.write("Begin Conditions {0} //Top\n".format(condtype))

    # Top
    for iy in range(ny):
        for ix in range(nx):
            n0,n1,n2,n3,n4,n5,n6,n7 = box.cube_vertices(ix,iy,nz-1)
                
            mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index,prop_id,n4,n6,n5))
            mdpa.write("{0:d} {1:d} {2:d} {3:d} {4:d}\n".format(index+1,prop_id,n4,n7,n6))
            index += 2
  
    mdpa.write("End Conditions\n\n")

    box.cond_range["top"] = [i0,index]

    return index

def generate_conditions(mdpa,box,condtype="WallCondition3D",prop_id=0):
    nx = box.nx()
    ny = box.ny()
    nz = box.nz()

    print "Generating {0} faces.".format(condtype)
    
    index = generate_front_faces(mdpa,box,1,condtype,prop_id)    
    index = generate_back_faces(mdpa,box,index,condtype,prop_id)
    index = generate_right_faces(mdpa,box,index,condtype,prop_id)
    index = generate_left_faces(mdpa,box,index,condtype,prop_id)
    index = generate_bottom_faces(mdpa,box,index,condtype,prop_id)
    index = generate_top_faces(mdpa,box,index,condtype,prop_id)

def generate_mesh_groups(mdpa,box):
    index = 1
    nx = box.nx()
    ny = box.ny()
    nz = box.nz()

    # Front group
    mdpa.write("Begin Mesh 1 //Front\n")
    mdpa.write("Begin MeshNodes\n")
    for iz in range(nz+1):
        for iy in range(ny+1):
            mdpa.write("{0}\n".format(box.get_id(0,iy,iz)))
    mdpa.write("End MeshNodes\n\n")

    mdpa.write("Begin MeshConditions\n")
    for i in range(*box.cond_range["front"]):
        mdpa.write("{0}\n".format(i))
    mdpa.write("End MeshConditions\n\n")
    
    mdpa.write("End Mesh\n\n")

    # Back group
    mdpa.write("Begin Mesh 2 //Back\n")
    mdpa.write("Begin MeshNodes\n")
    for iz in range(nz+1):
        for iy in range(ny+1):
            mdpa.write("{0}\n".format(box.get_id(nx,iy,iz)))
    mdpa.write("End MeshNodes\n\n")

    mdpa.write("Begin MeshConditions\n")
    for i in range(*box.cond_range["back"]):
        mdpa.write("{0}\n".format(i))
    mdpa.write("End MeshConditions\n\n")
    
    mdpa.write("End Mesh\n\n")

    # Right group
    mdpa.write("Begin Mesh 3 //Right\n")
    mdpa.write("Begin MeshNodes\n")
    for iz in range(nz+1):
        for ix in range(nx+1):
            mdpa.write("{0}\n".format(box.get_id(ix,0,iz)))
    mdpa.write("End MeshNodes\n\n")

    mdpa.write("Begin MeshConditions\n")
    for i in range(*box.cond_range["right"]):
        mdpa.write("{0}\n".format(i))
    mdpa.write("End MeshConditions\n\n")
    
    mdpa.write("End Mesh\n\n")

    # Left group
    mdpa.write("Begin Mesh 4 //Left\n")
    mdpa.write("Begin MeshNodes\n")
    for iz in range(nz+1):
        for ix in range(nx+1):
            mdpa.write("{0}\n".format(box.get_id(ix,ny,iz)))
    mdpa.write("End MeshNodes\n\n")

    mdpa.write("Begin MeshConditions\n")
    for i in range(*box.cond_range["left"]):
        mdpa.write("{0}\n".format(i))
    mdpa.write("End MeshConditions\n\n")
    
    mdpa.write("End Mesh\n\n")

    # Bottom group
    mdpa.write("Begin Mesh 5 //Bottom\n")
    mdpa.write("Begin MeshNodes\n")
    for iy in range(ny+1):
        for ix in range(nx+1):
            mdpa.write("{0}\n".format(box.get_id(ix,iy,0)))
    mdpa.write("End MeshNodes\n\n")

    mdpa.write("Begin MeshConditions\n")
    for i in range(*box.cond_range["bottom"]):
        mdpa.write("{0}\n".format(i))
    mdpa.write("End MeshConditions\n\n")
    
    mdpa.write("End Mesh\n\n")

    # Top group
    mdpa.write("Begin Mesh 6 //Top\n")
    mdpa.write("Begin MeshNodes\n")
    for iy in range(ny+1):
        for ix in range(nx+1):
            mdpa.write("{0}\n".format(box.get_id(ix,iy,nz)))
    mdpa.write("End MeshNodes\n\n")

    mdpa.write("Begin MeshConditions\n")
    for i in range(*box.cond_range["top"]):
        mdpa.write("{0}\n".format(i))
    mdpa.write("End MeshConditions\n\n")
    
    mdpa.write("End Mesh\n\n")


if __name__ == "__main__":
    import sys
    
    # mesh dimensions
    nx = int(sys.argv[1])
    ny = int(sys.argv[2])
    nz = int(sys.argv[3])
        
    # Channel dimensions
    xmin = 0.0
    xmax = 6.2832
    ymin = 0.0
    ymax = 2.0944
    zmin = 0.0
    zmax = 2.0

    box = box_data(xmin,ymin,zmin,xmax,ymax,zmax,nx,ny,nz)
    
    # element type
    if len(sys.argv) > 4:
        elemtype = sys.argv[4]
    else:
        elemtype = "VMS3D"
    condtype = "WallCondition3D"
        
    filename = "box_{0}x{1}x{2}.mdpa".format(nx,ny,nz)
    
    with open(filename,"w") as mdpa:
        write_header(mdpa)
        generate_nodes(mdpa,box)
        generate_elements(mdpa,box,elemtype)
        generate_conditions(mdpa,box,condtype)
        generate_mesh_groups(mdpa,box)

