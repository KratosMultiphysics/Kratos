def ReadLayerSets():
    ##################################################################
    ## STORE LAYER SETS ##############################################
    ##################################################################
    ## ELEMENTS on layers ############################################
    layer_sets = {}
*loop layers
*set Layer *LayerName *elems
    layer_elements_list = [
*loop elems *OnlyInLayer
    *ElemsNum ,
*end elems
    ]
    layer_sets['*LayerName'] = layer_elements_list
*end layers
    return layer_sets

def ReadLayerNodesSets():
    ## NODES on layers ###############################################
    layer_nodes_sets = {}
*loop layers
*set Layer *LayerName *nodes
    layer_nodes_list = [
*loop nodes *OnlyInLayer
    *NodesNum ,
*end nodes
    ]
    layer_nodes_sets['*LayerName'] = layer_nodes_list
*end layers
    return layer_nodes_sets

def ReadNodeGroups():
    ##################################################################
    ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
    ##################################################################
    node_groups = {}
*set cond Volume_Group_Membership *nodes
*loop nodes *OnlyInCond
    if( not ('*cond(Group_Name)' in node_groups) ):
        node_groups['*cond(Group_Name)'] = []
    node_groups['*cond(Group_Name)'].append( *NodesNum )
*end nodes
*set cond Surface_Group_Membership *nodes
*loop nodes *OnlyInCond
    if( not ('*cond(Group_Name)' in node_groups) ):
        node_groups['*cond(Group_Name)'] = []
    node_groups['*cond(Group_Name)'].append( *NodesNum )
*end nodes
*set cond Line_Group_Membership *nodes
*loop nodes *OnlyInCond
    if( not ('*cond(Group_Name)' in node_groups) ):
        node_groups['*cond(Group_Name)'] = []
    node_groups['*cond(Group_Name)'].append( *NodesNum )
*end nodes
    return node_groups

