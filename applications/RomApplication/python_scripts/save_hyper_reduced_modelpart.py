### Create Hyper Reduced Model Part process

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import pdb

class SaveHyperReducedModelPart:

    def WriteHRModelPart(self, problem_name,computing_model_part,hyper_reduced_model_part):

        OriginalMdpa = open(problem_name + ".mdpa" , 'r')
        OutputMdpa = open(problem_name + "_HyperReduced_" + ".mdpa", 'w')


        NodesLoopFlag = False
        ElementsLoopFlag = False
        ConditionsLoopFlag = False


        for Line in OriginalMdpa:
                       
           ## Changing flags
            if Line.startswith("Begin Nodes"):
                NodesLoopFlag = True
                OutputMdpa.write(Line)
                continue

            elif Line.startswith("Begin Elements"):
                ElementsLoopFlag = True  
                OutputMdpa.write(Line)
                continue                   
                
            elif Line.startswith("Begin Conditions"):
                ConditionsLoopFlag = True
                OutputMdpa.write(Line)
                continue            

            elif len(Line) != 1:
                if Line.split()[0] ==("Begin"):
                    if Line.split()[1] == ("SubModelPartNodes"):
                        NodesLoopFlag = True
                        OutputMdpa.write(Line)
                        continue                        
                    elif Line.split()[1] ==("SubModelPartElements"):
                        ElementsLoopFlag = True 
                        OutputMdpa.write(Line)
                        continue                     
                    elif Line.split()[1] ==("SubModelPartConditions"):
                        ConditionsLoopFlag = True 
                        OutputMdpa.write(Line)
                        continue                         
                    
                elif Line.split()[0] ==("End"):
                    if Line.split()[1] == ("SubModelPartNodes"):
                        NodesLoopFlag = False
                    elif Line.split()[1] ==("SubModelPartElements"):
                        ElementsLoopFlag = False 
                    elif Line.split()[1] ==("SubModelPartConditions"):
                        ConditionsLoopFlag = False 


            if NodesLoopFlag:
                if Line.startswith("End Nodes"):
                    NodesLoopFlag = False
                SplittedLine = Line.split()
                for node in hyper_reduced_model_part.Nodes:                    
                    ReducedNodes = str(node.Id)
                    if SplittedLine[0] == ReducedNodes:
                        OutputMdpa.write(Line)
                        break

            elif ElementsLoopFlag:
                if Line.startswith("End Elements"):
                    ElementsLoopFlag = False 
                SplittedLine = Line.split()
                for element in hyper_reduced_model_part.Elements:                    
                    ReducedElements = str(element.Id)
                    if SplittedLine[0] == ReducedElements:
                        OutputMdpa.write(Line)
                        break    

            elif ConditionsLoopFlag:
                if Line.startswith("End Conditions"):
                    ConditionsLoopFlag = False 
                SplittedLine = Line.split()
                for condition in hyper_reduced_model_part.Conditions:                    
                    ReducedConditions = str(condition.Id)
                    if SplittedLine[0] == ReducedConditions:
                        OutputMdpa.write(Line)
                        break                 
 
            ## Saving lines in the mdpa 

            if NodesLoopFlag == ElementsLoopFlag == ConditionsLoopFlag == False:
                OutputMdpa.write(Line)


        OutputMdpa.close()
        


    



