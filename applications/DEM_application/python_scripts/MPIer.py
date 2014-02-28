#Modelpart Fixer
import sys
import math
import shutil as Shutil

class MPIerClass:

  def __init__(self,input_file):
      
      Shutil.copyfile(str(input_file),'backup.mdpa')
      output_file = open(input_file,'w')
      
      with open('backup.mdpa') as f:
          content = f.readlines()
          
      fix_nodes = 0
      old_id = 0

      NodalData_VELOCITY_X = []
      NodalData_VELOCITY_Y = []
      NodalData_VELOCITY_Z = []
      NodalData_RADIUS = []
      NodalData_GROUP_ID = []
      NodalData_PARTICLE_CONTINUUM = []
      NodalData_PARTICLE_DENSITY = []
      NodalData_YOUNG_MODULUS = []
      NodalData_POISSON_RATIO = []
      NodalData_PARTICLE_FRICTION = []
      NodalData_LN_OF_RESTITUTION_COEFF = []
      NodalData_PARTICLE_MATERIAL = []

      for l in content:
        
          if("Begin Nodes" in l):
              fix_nodes = 1
              output_file.write(l)
              continue
        
          if("End Nodes" in l):
              fix_nodes = 0
              output_file.write(l)
              continue
                  
          if(fix_nodes == 1):
              a = l.split(' ',1)
              
              cur_id = int(a[0])
              
              while(cur_id-old_id != 1):
                  old_id = old_id + 1
                  
                  NodalData_VELOCITY_X.append([old_id, 0])
                  NodalData_VELOCITY_Y.append([old_id, 0])
                  NodalData_VELOCITY_Z.append([old_id, 0])
                  NodalData_RADIUS.append([old_id, 5.00000e-010])
                  NodalData_GROUP_ID.append([old_id, 1])
                  NodalData_PARTICLE_CONTINUUM.append([old_id, 1])
                  NodalData_PARTICLE_DENSITY.append([old_id, 1000])
                  NodalData_YOUNG_MODULUS.append([old_id, 2.1e6])
                  NodalData_POISSON_RATIO.append([old_id, 0.5])
                  NodalData_PARTICLE_FRICTION.append([old_id, 1.00000e+00])
                  NodalData_LN_OF_RESTITUTION_COEFF.append([old_id, -6.93147e-01])
                  NodalData_PARTICLE_MATERIAL.append([old_id, 1])
              
                  output_file.write(str(old_id) + " 0 0 0\n")
                  
              old_id = cur_id
          
          output_file.write(l)
          
          if("NodalData VELOCITY_X" in l):
              for n in NodalData_VELOCITY_X:
                  output_file.write(str(n[0]) + " 1 " + str(n[1]) + "\n")
          if("NodalData VELOCITY_Y" in l):
              for n in NodalData_VELOCITY_Y:
                  output_file.write(str(n[0]) + " 1 " + str(n[1]) + "\n")      
          if("NodalData VELOCITY_Z" in l):
              for n in NodalData_VELOCITY_Z:
                  output_file.write(str(n[0]) + " 1 " + str(n[1]) + "\n")
          if("NodalData RADIUS" in l):
              for n in NodalData_RADIUS:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData GROUP_ID" in l):
              for n in NodalData_GROUP_ID:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData PARTICLE_CONTINUUM" in l):
              for n in NodalData_PARTICLE_CONTINUUM:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData PARTICLE_DENSITY" in l):
              for n in NodalData_PARTICLE_DENSITY:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData YOUNG_MODULUS" in l):
              for n in NodalData_YOUNG_MODULUS:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData POISSON_RATIO" in l):
              for n in NodalData_POISSON_RATIO:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")      
          if("NodalData PARTICLE_FRICTION" in l):
              for n in NodalData_PARTICLE_FRICTION:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData LN_OF_RESTITUTION_COEFF" in l):
              for n in NodalData_LN_OF_RESTITUTION_COEFF:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
          if("NodalData PARTICLE_MATERIAL" in l):
              for n in NodalData_PARTICLE_MATERIAL:
                  output_file.write(str(n[0]) + " 0 " + str(n[1]) + "\n")
      
      output_file.close()
          

