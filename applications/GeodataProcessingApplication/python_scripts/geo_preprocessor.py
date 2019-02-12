import KratosMultiphysics
from geo_processor import GeoProcessor

import os
import math

# one class for Preprocessor (no derived classes intended)
class GeoPreprocessor( GeoProcessor ):


# read_functions (pattern: read_FILENAME)

# operation_functions (pattern: OPERATIONNAME)

# write_functions (pattern: write_FILENAME)

#### NOT YET IN STRUCTURE
# handling_functions ( clear(), statistics() )


    def __init__( self ):
        super(GeoPreprocessor, self).__init__()
        self.point_list = []

    ########################################################################
    # for xyz files
    #########################################################################

    def read_xyz(self, xyz_in_file):

        filename_in, extension_in = os.path.splitext(xyz_in_file)

        print("Reading from file " + filename_in + extension_in)

        if (extension_in.upper() != ".XYZ"):
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "This function operates on XYZ files only.")

        with open(xyz_in_file) as read_file:
            for line in read_file.readlines():
                items = line.split()
                if (all(self._isfloat(n) for n in items)):
                    x = float(items[0])
                    y = float(items[1])
                    z = float(items[2])

                    self.point_list.append([x,y,z])


    def write_xyz(self, type_file="standard", xyz_file_name_output="", list_to_write=[]):
        # function to write a xyz file from a list of points

        filename_out, extension_out = os.path.splitext(xyz_file_name_output)
        if not extension_out:
            # extension_out is empty
            extension_out = ".xyz"
        
        if (list_to_write):
            # list_to_write is full
            self.point_list = list_to_write

        if (type_file == "standard"):
            out_name = filename_out + "_standard" + extension_out

        elif (type_file == "qgis"):
            out_name = filename_out + "_qgis" + extension_out

            # sorting the list to enable output of the selected points in QGIS (Formatted ASCII-XYZ)
            # here with a lambda function
            # result:   y-value as primary sorting parameter
            #           x-value as secondary sorting parameter
            num_nodes = len(self.point_list)
            self.point_list.sort(key=lambda x: float( x[1] * (num_nodes+1) + x[0]))
        
        else:
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "WARNING: type_data not valid. Please insert \"standard\" or \"qgis\".")


        with open(out_name, 'w+') as fout:
            for i in range(0, len(self.point_list), 1):
                x_string = str((self.point_list[i])[0])
                y_string = str((self.point_list[i])[1])
                z_string = str((self.point_list[i])[2])
                line = x_string + " " + y_string + " " + z_string + "\n"
                fout.write(line)


    ########################################################################
    # for stl files
    #########################################################################
    
    def read_stl(self, stl_file_name_input):

        filename_in, extension_in = os.path.splitext(stl_file_name_input)

        print("Reading from file " + filename_in + extension_in)

        if (extension_in.upper() != ".STL"):
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "This function operates on STL files only.")

        with open (stl_file_name_input) as read_file:
            for row in read_file.readlines():
                row = row.split()
                if (row[0] == "vertex" and all(self._isfloat(n) for n in row[1:])):
                    # row[0] = "vertex"; row[1] = x coordinate; row[2] = y coordinate; row[3] = z coordinate
                    X_coord, Y_coord, Z_coord = [float(coord) for coord in row[1:]]

                    self.point_list.append([X_coord, Y_coord, Z_coord])
   

    def write_stl(self, stl_file_name_output="", list_to_write=[]):
        # function to write a stl file from a list of points

        if (stl_file_name_output == ""):
            out_name = "_new.stl"
        else:
            filename_out, extension_out = os.path.splitext(stl_file_name_output)
            out_name = filename_out + "_new" + extension_out
        
        if (list_to_write):
            # list_to_write is full
            self.point_list = list_to_write
        
        if ((len(self.point_list)%3) != 0):
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "There are some free points.")

        string_out = "solid model\n"        # start of the stl file
        pointer = 0
        num_points = len(self.point_list)
        while (pointer < num_points-2):        # to avoid "out of range"
            p1 = self.point_list[pointer]
            p2 = self.point_list[pointer+1]
            p3 = self.point_list[pointer+2]

            nx, ny, nz = self._normal_triangle(p1, p2, p3)

            pointer += 3

            string_out += "\tfacet normal {} {} {}\n".format(nx, ny, nz)
            string_out += "\t\touter loop\n"
            string_out += "\t\t\tvertex {} {} {}\n".format(p1[0], p1[1], p1[2])
            string_out += "\t\t\tvertex {} {} {}\n".format(p2[0], p2[1], p2[2])
            string_out += "\t\t\tvertex {} {} {}\n".format(p3[0], p3[1], p3[2])
            string_out += "\t\tendloop\n"
            string_out += "\tendfacet\n"
        
        string_out += "endsolid model\n"    # end of file stl
        
        with open(out_name, "w") as fout:
            fout.write(string_out)


    ########################################################################
    # general functions
    #########################################################################

    def cut(self, bounding_box=(-999999999.0, -999999999.0, 999999999.0, 999999999.0)):

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

        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        for num in range(len(self.point_list)):
            coord = self.point_list[num]
            if ((x_min <= coord[0] <= x_max) and (y_min <= coord[1] <= y_max)):
                coord[0] -= x_min
                coord[1] -= y_min
                self.point_list[num] = coord

                num_nodes += 1

        if num_nodes == 0:
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "WARNING: No points inside the given bounding box!")


    def shift(self, x_shift, y_shift, z_shift=0):

        print("Shift defined as")
        print("x_shift = ", x_shift)
        print("y_shift = ", y_shift)
        print("z_shift = ", z_shift)

        for num in range(len(self.point_list)):
            coord = self.point_list[num]
            coord[0] += x_shift
            coord[1] += y_shift
            coord[2] += z_shift
            self.point_list[num] = coord
    

    def valley(self, max_height):

        print("Maximal height defined as")
        print("max_height = ", max_height)

        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        del_list = []
        for num in range(len(self.point_list)):
            coord = self.point_list[num]
            if (coord[2] > max_height):
                del_list.append(coord)
        
        for coord in del_list:
            self.point_list.remove(coord)
        
        if not self.point_list:
            # if self.point_list is empty
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "WARNING: No points inside the given valley!")


    def mountain(self, min_height):

        print("Minimal height defined as")
        print("min_height = ", min_height)

        # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        del_list = []
        for num in range(len(self.point_list)):
            coord = self.point_list[num]
            if (coord[2] < min_height):
                del_list.append(coord)
        
        for coord in del_list:
            self.point_list.remove(coord)
        
        if not self.point_list:
            # if self.point_list is empty
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "WARNING: No points inside the given valley!")


    def river(self, radius, start_x_y):

        current_center = start_x_y                  # starting point for the search
        previous_center = (0.0, 0.0)

        print("Starting identification of the river bed")
        print(" - starting point at ", current_center)

        num_points = 0                              # used for checks
        dict_all_points = {}                        # key = id, value = ( x, y, z )
        dict_active_points = {}                     # key = id, value = z
        list_visited_points = []                    # all points that have ever been active (index only)
                                                    # always extend( list_active_points )

        # # reading all points into a single list out of tuples (x,y,z)
        # # the coordinates are already modified and (x=0,y=0) is located at x_min, y_min
        # with open( xyz_in_file ) as f:
        #     lines = f.readlines()
        #     for line in lines:
        #         items = line.split()
        #         if ( all( self._isfloat(n) for n in items) ):
        #             x = float(items[0])
        #             y = float(items[1])
        #             z = float(items[2])
        #             dict_all_points[num_points] = [x,y,z]
        #             num_points += 1

        for pnt_id in range(len(self.point_list)):
            coord = self.point_list[pnt_id]
            dict_all_points[num_points] = [coord[0], coord[1], coord[2]]
            num_points += 1
        
        if num_points == 0:
            KratosMultiphysics.Logger.PrintWarning("GeoPreprocessor", "WARNING: No points inside the given valley!")

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

        # # extracting the "relevant part" marked by the visited points
        # point_list = self._create_list_of_visited_points( list_visited_points, dict_all_points )

        ### CHECK THIS PART ###
        # HERE WE FILLED (AGAIN) self.point_list
        self.point_list = self._create_list_of_visited_points( list_visited_points, dict_all_points )


    def centroid_rect2D(self):
        # function to compute the centroid in a rectangle in 2D

        x_coords = []
        y_coords = []
        for coord in self.point_list:
            x_coords.append(coord[0])
            y_coords.append(coord[1])
        
        x_centr = max(x_coords)/2
        y_centr = max(y_coords)/2

        return (x_centr, y_centr)


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
                point_list.append([x,y,z])

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


    def _isfloat( self, value):
        try:
            float( value )
            return True
        except ValueError:
            return False


    def _normal_triangle(self, P1, P2, P3):
        # function to calc the normal of the triangle
        
        nx = (P2[1]-P1[1])*(P3[2]-P1[2]) - (P3[1]-P1[1])*(P2[2]-P1[2])
        ny = (P2[2]-P1[2])*(P3[0]-P1[0]) - (P2[0]-P1[0])*(P3[2]-P1[2])
        nz = (P2[0]-P1[0])*(P3[1]-P1[1]) - (P3[0]-P1[0])*(P2[1]-P1[1])

        sum_n = abs(nx) + abs(ny) + abs(nz)

        if sum_n == 0.0:
            nx = ny = nz = 0.0
        else:
            nx = nx / sum_n
            ny = ny / sum_n
            nz = nz / sum_n

        return (nx, ny, nz)


    def _remove_duplicate(self, point_list):
        no_duplicate = []        # point_list without duplicates
        for point in point_list:
            if (point not in no_duplicate):
                no_duplicate.append(point)

        return no_duplicate
