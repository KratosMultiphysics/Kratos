from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

def GetNodeBefore(table,prop):
  return table[prop-1][1]
  
def GetNodeBegin(table,prop):
  return table[prop-1][2]
  
def GetNodeEnd(table,prop):
  return table[prop-1][3]
  
def GetNodeAfter(table,prop):
  return table[prop-1][4]
  
#select and remove
def DoRemoval(model_part):
  import config_by
  import full_nodes_table
  print "deactivate_list",config_by.deactivate_list
  print "list of inlets",config_by.inlets_1d
  print "list of outlets",config_by.outlets_1d
  print "full nodes table",full_nodes_table.table
  print model_part
  
  
  #mark for deactivation all of the nodes which are not needed
  for prop_id in config_by.deactivate_list:
    for elem in model_part.Elements:
      if(elem.Properties.Id == prop_id):
	elem.SetValue(ERASE_FLAG,True)
	for node in elem.GetNodes():
	  node.SetValue(ERASE_FLAG,True)
	  
  #define a list of nodes which we shall than preserve from removal
  nodes_to_preserve = []
  outlet_nodes = []
	  
  #mark for erasal nodes before inlet 
  for i in range(0,len(config_by.inlets_1d)):
    flag_id = config_by.inlets_1d[i][0]
    prop_id = config_by.inlets_1d[i][1]
    print "flag_id", flag_id
    print "prop_id", prop_id
    print GetNodeBefore(full_nodes_table.table, prop_id)
    node_before = model_part.Nodes[ GetNodeBefore(full_nodes_table.table, prop_id) ]
    node_begin = model_part.Nodes[ GetNodeBegin(full_nodes_table.table, prop_id) ]
    node_before.SetValue(ERASE_FLAG,True)
    node_begin.SetValue(ERASE_FLAG,True)
    node_begin.SetSolutionStepValue(FLAG_VARIABLE,0,flag_id)
    nodes_to_preserve.append(node_begin)
    nodes_to_preserve.append(node_before)
    
  #mark for erasal nodes after outler
  for i in range(0,len(config_by.outlets_1d)):
    flag_id = config_by.outlets_1d[i][0]
    prop_id = config_by.outlets_1d[i][1]
    node_end = model_part.Nodes[ GetNodeEnd(full_nodes_table.table, prop_id) ]
    node_after = model_part.Nodes[ GetNodeAfter(full_nodes_table.table, prop_id) ]
    node_end.SetValue(ERASE_FLAG,True)
    node_after.SetValue(ERASE_FLAG,True)
    print "prop_id = ",prop_id
    print "node end = ",node_end.Id
    print "node after = ",node_after.Id
    node_after.SetSolutionStepValue(FLAG_VARIABLE,0,flag_id)
    nodes_to_preserve.append(node_after)
    outlet_nodes.append(node_after)
	  
  #mark for deactivation the conditions which have all of their nodes marked for erasal
  for cond in model_part.Conditions:
    aaa = 0
    tot_nodes = 0
    for node in cond.GetNodes():
      tot_nodes+=1
      if(node.GetValue(ERASE_FLAG) == True):
	aaa+=1
    if(tot_nodes == aaa and aaa!=0):
      cond.SetValue(ERASE_FLAG,True)
      
  #unamrk the node to be preserved
  for node in nodes_to_preserve:
    node.SetValue(ERASE_FLAG,False)
    
  #do delete elements conditions and nodes
  NodeEraseProcess(model_part).Execute()
  ElementEraseProcess(model_part).Execute()
  ConditionEraseProcess(model_part).Execute()
  
  highest_cond_id = 0
  for cond in model_part.Conditions:
    if(cond.Id > highest_cond_id):
      highest_cond_id = cond.Id
  
  #add conditions associated to the outlet of the 3d
  id_of_new_property = 0
  for node in outlet_nodes:
    CreateNewCondition( "ArteryInletCondition",model_part, highest_cond_id+1, id_of_new_property, node)
    highest_cond_id+=1
  
  #add conditions associated to the inlet of the 3d
  for i in range(0,len(config_by.inlets_1d)):
    flag_id = config_by.inlets_1d[i][0]
    prop_id = config_by.inlets_1d[i][1]
    node_before = model_part.Nodes[ GetNodeBefore(full_nodes_table.table, prop_id) ]
    node_begin = model_part.Nodes[ GetNodeBegin(full_nodes_table.table, prop_id) ]
    CreateNewCondition( "Artery1Dto3DCondition",model_part, highest_cond_id+1, id_of_new_property, node_before, node_begin)
    highest_cond_id+=1