from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

class TimeOperationUtility(object):

    #
    def __init__(self):

        # set time variables integer counters
        self.id_counter = 0
        self.frequency = 0

        # set time variables double timing
        self.starting_time = 0
        self.ending_time = 0

        self.time_step = 0
        self.time_frequency = 0

        self.time_counter = 0

        self.last_seen_current_time = 0

        self.tolerance = 0

    #
    def InitializeTime(self, starting_time, ending_time, time_step, time_frequency):

        # set starting time
        self.starting_time = starting_time

        # set ending time
        self.ending_time = ending_time

        # set step time
        self.time_step = time_step

        # set last seen current time
        self.last_seen_current_time = starting_time - time_step

        # set time frequency
        self.time_frequency = time_frequency

        if(self.time_frequency < time_step):
            self.time_frequency = time_step

        # set time counters
        self.time_counter = self.starting_time + time_frequency

        # set time operation tolerance
        self.tolerance = self.time_step * 1e-10;

        # -----#

        # set operation frequency
        self.frequency = int(self.time_frequency / self.time_step)

        if(self.frequency < 1):
            self.frequency = 1

        # set step counters
        self.id_counter = int(self.starting_time / self.time_frequency)


    #
    def perform_time_operation(self, current_time):

        execute = False
        
        if(self.last_seen_current_time < current_time):
            
            if(current_time + self.tolerance >= self.ending_time):
                execute = True
                self.id_counter = self.id_counter + 1
            elif(current_time + self.tolerance >= self.time_counter):
                execute = True
                self.time_counter = self.time_counter + self.time_frequency
                self.id_counter   = self.id_counter + 1

        elif(self.last_seen_current_time == current_time):
            
            if(current_time + self.tolerance >= self.ending_time):
                execute = True
            elif(current_time + self.tolerance >= self.time_counter):
                execute = True
 
                
        self.last_seen_current_time = current_time

        return execute

    #
    def operation_id(self):

        return self.id_counter
