# -*- coding: utf-8 -*-
from scipy import *
from scipy import linalg

import test_pod

fluid_model_part = ModelPart("FluidPart");  
  
def Reduced_Pods_Read(infile,Number_Nodes,Number_Pods):
      
      File = open(infile,'r')
      IN=[]
      in_=[]
      t_=1
      print "begin reading Reduced PODS File"  
      i=0
      while (t_<=Number_Nodes):
	   searchname_=str('Reduced POD values for    ')+str(t_)
           line = File.readline()
           #print "searching for:   " +searchname_        
           if not line:
             break
           else:   
	     if(line.find(searchname_) >= 0):
	        print "begin reading "+searchname_ 
                line = File.readline()
                
                j=0
                if(line.find(str('Values')) >= 0):
		   go_on = True 
                   while go_on:
                     line = File.readline()
                     if not line:
                        break
                   
                     if (line.find("End Values") >= 0):
                        print 'END Reading elements'
                        go_on = False;
                        t_=t_+1
                        #searchname_=str('Result "')+Quantity_+str ('" "Kratos" ')+str(t_)+str(' ')+Dimension_+str(' OnNodes')
                        #print searchname_
                     else:
		        
		        aaa = line.split();
		        in_.insert(j, aaa[0])
		        j=j+1
	                #print aaa[axis]
	           #IN.insert(i,double(in_[0:j]))
	           fluid_model_part.Nodes[i].SetValue(POD_VELOCITY_X,double(in_[0:j]))
	           i=i+1
	           #return IN	          
      File.close()

  