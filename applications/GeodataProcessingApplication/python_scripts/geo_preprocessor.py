from geo_processor import GeoProcessor

import numpy as np
import os
import math

# one class for Preprocessor (no derived classes intended)
class GeoPreprocessor( GeoProcessor ):

#### NOT YET IN STRUCTURE

# read_functions (pattern: read_FILENAME)

# operation_functions (pattern: OPERATIONNAME)

# write_functions (pattern: write_FILENAME)

# handling_functions ( clear(), statistics() )


    def __init__( self ):
        super(GeoPreprocessor, self).__init__()
        self.point_list = []

    ########################################################################
    # for xyz files
    #########################################################################

    def xyz_cut( self, xyz_in_file, xyz_out_file, bounding_box = (-999999999.0, -999999999.0, 999999999.0, 999999999.0) ):

        filename_in, extension_in = os.path.splitext( xyz_in_file )
        filename_out, extension_out = os.path.splitext( xyz_out_file )

        print( "Reading from file " + filename_in + extension_in)

        if ( extension_in.upper() != ".XYZ" or extension_out.upper() != ".XYZ" ):
            print("This function operates on XYZ files only.")

        x_min = bounding_box[0]
        y_min = bounding_box[1]
        x_max = bounding_box[2]
        y_max = bounding_box[3]
        print( "Bounding box defined as " )
        print( "x_min = " + str(x_min) )
        print( "y_min = " + str(y_min) )
        print( "x_max = " + str(x_max) )
        print( "y_max = " + str(y_max) )
        # num_lines = sum(1 for line in open( xyz_in_file ))
        num_nodes = 0
        point_list = []

        # reading all points into a single list out of tuples (x,y,z)
        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        with open( xyz_in_file ) as f:
            lines = f.readlines()
            counter = 0
            for line in lines:
                items = line.split()
                if ( all( self._isfloat(n) for n in items) ):
                    if ( ( x_min <= float(items[0]) <= x_max ) and ( y_min <= float(items[1]) <= y_max ) ):
                        # to avoid numerical inaccuracy, the value are already
                        x = float(items[0]) - x_min
                        y = float(items[1]) - y_min
                        z = float(items[2])
                        point_list.append( (x,y,z) )
                        counter = counter + 1
            num_nodes = counter

        if num_nodes == 0:
            print( "WARNING: No points inside the given bounding box!")

        # writing to files
        self._write_standard_xyz( point_list, xyz_out_file )
        self._write_qgis_xyz( point_list, xyz_out_file )


    def xyz_shift( self, xyz_in_file, xyz_out_file, x_shift, y_shift ):

        filename_in, extension_in = os.path.splitext( xyz_in_file )
        filename_out, extension_out = os.path.splitext( xyz_out_file )

        print( "Reading from file " + filename_in + extension_in)

        if ( extension_in.upper() != ".XYZ" or extension_out.upper() != ".XYZ" ):
            print("This function operates on XYZ files only.")

        print( "Shift defined as " )
        print( "x_shift = " + str(x_shift) )
        print( "y_shift = " + str(y_shift) )

        # num_lines = sum(1 for line in open( xyz_in_file ))
        point_list = []

        # reading all points into a single list out of tuples (x,y,z)
        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        with open( xyz_in_file ) as f:
            lines = f.readlines()
            counter = 0
            for line in lines:
                items = line.split()
                if ( all( self._isfloat(n) for n in items) ):
                    x = float(items[0]) + x_shift
                    y = float(items[1]) + y_shift
                    z = float(items[2])
                    point_list.append( (x,y,z) )
                    counter = counter + 1
            num_nodes = counter

        # writing to files
        self._write_standard_xyz( point_list, xyz_out_file )
        self._write_qgis_xyz( point_list, xyz_out_file )


    def xyz_valley( self, xyz_in_file, xyz_out_file, max_height ):

        filename_in, extension_in = os.path.splitext( xyz_in_file )
        filename_out, extension_out = os.path.splitext( xyz_out_file )

        print( "Reading from file " + filename_in + extension_in)

        if ( extension_in.upper() != ".XYZ" or extension_out.upper() != ".XYZ" ):
            print("This function operates on XYZ files only.")

        print( "Maximal height defined as " )
        print( "max_height = " + str(max_height) )

        # num_lines = sum(1 for line in open( xyz_in_file ))
        num_nodes = 0
        point_list = []

        # reading all points into a single list out of tuples (x,y,z)
        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        with open( xyz_in_file ) as f:
            lines = f.readlines()
            counter = 0
            for line in lines:
                items = line.split()
                if ( all( self._isfloat(n) for n in items) ):
                    if ( float(items[2]) <= max_height ):
                        x = float(items[0])
                        y = float(items[1])
                        z = float(items[2])
                        point_list.append( (x,y,z) )
                        counter = counter + 1
            num_nodes = counter

        if num_nodes == 0:
            print( "WARNING: No points inside the given valley!")

        # writing to files
        self._write_standard_xyz( point_list, xyz_out_file )
        self._write_qgis_xyz( point_list, xyz_out_file )


    def xyz_mountain( self, xyz_in_file, xyz_out_file, min_height ):

        filename_in, extension_in = os.path.splitext( xyz_in_file )
        filename_out, extension_out = os.path.splitext( xyz_out_file )

        print( "Reading from file " + filename_in + extension_in)

        if ( extension_in.upper() != ".XYZ" or extension_out.upper() != ".XYZ" ):
            print("This function operates on XYZ files only.")

        print( "Minimal height defined as " )
        print( "min_height = " + str(min_height) )

        # num_lines = sum(1 for line in open( xyz_in_file ))
        num_nodes = 0
        point_list = []

        # reading all points into a single list out of tuples (x,y,z)
        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        with open( xyz_in_file ) as f:
            lines = f.readlines()
            counter = 0
            for line in lines:
                items = line.split()
                if ( all( self._isfloat(n) for n in items) ):
                    if ( float(items[2]) >= min_height ):
                        x = float(items[0])
                        y = float(items[1])
                        z = float(items[2])
                        point_list.append( (x,y,z) )
                        counter = counter + 1
            num_nodes = counter

        if num_nodes == 0:
            print( "WARNING: No points inside the given valley!")

        # writing to files
        self._write_standard_xyz( point_list, xyz_out_file )
        self._write_qgis_xyz( point_list, xyz_out_file )


    def xyz_river( self, xyz_in_file, xyz_out_file, radius, start_x_y ):

        current_center = start_x_y                  # starting point for the search
        previous_center = (0.0, 0.0)

        print("Starting identification of the river bed")
        print( " - starting point at " + str( current_center ) )

        num_points = 0                              # used for checks
        dict_all_points = {}                        # key = id, value = ( x, y, z )
        dict_active_points = {}                     # key = id, value = z
        list_visited_points = []                    # all points that have ever been active (index only)
                                                    # always extend( list_active_points )

        # reading all points into a single list out of tuples (x,y,z)
        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        with open( xyz_in_file ) as f:
            lines = f.readlines()
            for line in lines:
                items = line.split()
                if ( all( self._isfloat(n) for n in items) ):
                    x = float(items[0])
                    y = float(items[1])
                    z = float(items[2])
                    dict_all_points[num_points] = [x,y,z]
                    num_points += 1

        if num_points == 0:
            print( "WARNING: No points inside the given valley!")

        num = 5                 # number of points for the search of minimal terrain height
        dist = radius / 4.0     # distance to shift the search circle
        length = 100000.0       # maximal length of the track (safety measure)
        generation = 1

        current_length = 0.0    # length of the track

        # Iteration shifting the search circle
        while current_length < length:

            print( " --- Generation " + str(generation) + " --- ")

            # points inside the current search circle are set acitve
            dict_active_points = self._find_points_within_radius( dict_all_points, current_center, radius )

            # Check if still new unvisited points are found
            if not self._iteration_goes_on( dict_active_points, list_visited_points ):
                print( " - No more unvisited points found - Assembly of the output in progress")
                break

            # keeping track of the already visite points
            # these points mark the relevant terrain parts
            list_visited_points = self._add_ids_of_active_points( list_visited_points, dict_active_points )
            print( " - " + str( len(dict_active_points) ) + " active points" )

            # find the active points with the minimal terrain height
            list_lowest_points = self._find_lowest_points( dict_active_points, num )

            # moving the circle forward
            previous_center = current_center
            current_center = self._compute_new_center( dict_all_points, list_lowest_points, dist, previous_center )
            print( " - new center at " + str( current_center ) )

            # proceeding to the new generation
            generation += 1
            current_length += dist

        # extracting the "relevant part" marked by the visited points
        point_list = self._create_list_of_visited_points( list_visited_points, dict_all_points )

        # writing to files
        self._write_standard_xyz( point_list, xyz_out_file )
        self._write_qgis_xyz( point_list, xyz_out_file )



    ########################################################################
    # auxiliary functions
    #########################################################################


    def _create_list_of_visited_points( self, list_visited_points, dict_all_points ):

        point_list = []

        for p_id, coords in dict_all_points.items():
            if ( p_id in list_visited_points ):
                x = coords[0]
                y = coords[1]
                z = coords[2]
                point_list.append( (x,y,z) )

        return point_list



    def _iteration_goes_on( self, dict_active_points, list_visited_points ):

        goes_on = False
        counter = 0

        if ( len(dict_active_points) > 0 ):
            for p_id, _ in dict_active_points.items():
                if not ( p_id in list_visited_points ):
                    goes_on = True
                    counter += 1
        else:
            goes_on = True

        return goes_on


    def _compute_new_center( self, dict_all_points, list_lowest_points, dist, previous_center ):

        x_list = []
        y_list = []
        current_center = (0.0, 0.0)

        for p_id in list_lowest_points:
            coords = dict_all_points[p_id]
            x = coords[0]
            y = coords[1]
            x_list.append( x )
            y_list.append( y )

        x_average = float (sum( x_list )) / float (len( x_list )) - previous_center[0]
        y_average = float (sum( y_list )) / float (len( y_list )) - previous_center[1]

        x_average = x_average / math.sqrt( ( x_average**2 + y_average**2 ) )
        y_average = y_average / math.sqrt( ( x_average**2 + y_average**2 ) )
        current_center = (previous_center[0] + x_average * dist, previous_center[1] + y_average * dist)

        return current_center


    def _find_points_within_radius( self, dict_all_points, current_center, radius ):

        dict_active_points = {}

        for id, coords in dict_all_points.items():
            x = coords[0]
            y = coords[1]
            dist = math.sqrt( ( x - current_center[0] )**2 + ( y - current_center[1] )**2 )

            if dist < radius:
                z = coords[2]
                dict_active_points[ id ] = z

        return dict_active_points


    def _add_ids_of_active_points( self, list_visited_points, dict_active_points ):

        for id, _ in dict_active_points.items():
            if ( not ( id in list_visited_points ) ):
                list_visited_points.append( id )

        return list_visited_points


    def _find_lowest_points( self, dict_active_points, num ):

        sorted_list_active_points = sorted( dict_active_points.items(), key=lambda kv: kv[1] )

        ids_of_lowest_points = []
        z_of_lowest_points = []
        for i in range( 0, num ):
            (id, z) = sorted_list_active_points[i]
            z_of_lowest_points.append( z )
            ids_of_lowest_points.append( id )

        return ids_of_lowest_points


    def _dist( self, x_start, y_start, x, y ):
        squares = (x - x_start)**2 + (y - y_start)**2
        return math.sqrt( squares )


    def _write_standard_xyz( self, point_list, xyz_out_file ):

        filename_out, extension_out = os.path.splitext( xyz_out_file )
        out_name = filename_out + "_standard" + extension_out

        with open( out_name, 'w+') as fout:
            for i in range(0, len(point_list), 1):
                x_string = str( (point_list[i])[0] )
                y_string = str( (point_list[i])[1] )
                z_string = str( (point_list[i])[2] )
                line = x_string + " " + y_string + " " + z_string + "\n"
                fout.write(line)


    def _write_qgis_xyz( self, point_list, xyz_out_file ):

        filename_out, extension_out = os.path.splitext( xyz_out_file )
        out_name = filename_out + "_qgis" + extension_out

        # sorting the list to enable output of the selected points in QGIS (Formatted ASCII-XYZ)
        # here with a lambda function
        # result:   y-value as primary sorting parameter
        #           x-value as secondary sorting parameter
        num_nodes = len(point_list)
        point_list.sort( key=lambda x: float( x[1] * (num_nodes+1) + x[0] ) )

        # writing out the preprocessed xyz file
        with open( out_name, 'w+') as fout:
            for i in range(0, len(point_list), 1):
                x_string = str( point_list[i][0] )
                y_string = str( (point_list[i])[1] )
                z_string = str( (point_list[i])[2] )
                line = x_string + " " + y_string + " " + z_string + "\n"
                fout.write(line)


    def _isfloat( self, value):
        try:
            float( value )
            return True
        except ValueError:
            return False

    def _remove_duplicate(self, point_list):
        no_duplicate = []        # point_list without duplicates
        for point in point_list:
            if (point not in no_duplicate):
                no_duplicate.append(point)

        return no_duplicate