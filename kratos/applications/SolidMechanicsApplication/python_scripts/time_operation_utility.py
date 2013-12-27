
class TimeOperationUtility(object):
 
    #######################################################################
    def __init__(self):
        
        #set time variables integer counters
        self.id_counter         = 0
        self.step_counter       = 0

        self.frequency          = 0
        self.starting_step      = 0
        self.ending_step        = 0

        #set time variables double timing
        self.starting_time      = 0
        self.ending_time        = 0

        self.time_step          = 0
        self.time_frequency     = 0

        self.time_counter       = 0
        
    #######################################################################
    def InitializeTime(self,starting_time,ending_time,time_step,time_frequency):

        #set starting time
        self.starting_time = starting_time

        #set ending time
        self.ending_time = ending_time

        #set step time
        self.time_step = time_step 

        #set time frequency
        self.time_frequency = time_frequency
        
        #set time counters 
        self.time_counter = self.starting_time

        #-----#

        #set operation frequency
        self.frequency = int(self.time_frequency/self.time_step)

        if(self.frequency < 1):
            self.frequency = 1

        #set starting step
        self.starting_step = int(self.starting_time/self.time_step)
        
        #set ending step
        self.ending_step = int(self.ending_time/self.time_step)

        #set step counters 
        self.step_counter = self.starting_step
        self.id_counter   = int(self.step_counter/self.frequency)



    #######################################################################
    def perform_time_operation(self,current_time):
        
        execute = False
        
        if( current_time + 1e-20 >= self.ending_time ):
            execute = True
        elif( current_time > self.time_counter ):
            execute = True
            self.time_counter = self.time_counter + self.time_frequency
            self.id_counter   = self.id_counter + 1
               
        return execute


    #######################################################################
    def InitializeStep(self,starting_step,ending_step,time_step,frequency):

        #set operation frequency
        self.frequency = frequency

        if(self.frequency < 1):
            self.frequency = 1

        #set starting step
        self.starting_step = starting_step
        
        #set ending step
        self.ending_step = ending_step

        #set step counters 
        self.step_counter = self.starting_step
        self.id_counter   = self.step_counter

        #-----#

        #set starting time
        self.starting_time = starting_step * time_step

        #set ending time
        self.ending_time = ending_step * time_step

        #set step time
        self.time_step = time_step 

        #set time frequency
        self.time_frequency = frequency * time_step
        
        #set time counters 
        self.time_counter = self.starting_time



    #######################################################################
    def perform_step_operation(self,current_step):
        
        execute = False

        if( current_step == self.ending_step ):
            execute = True
        elif( current_step == self.step_counter ):
            execute = True
            self.step_counter = self.step_counter + self.frequency
            self.id_counter   = self.id_counter + 1

        return execute


    #######################################################################
    def operation_id(self):
        
        return self.id_counter;
