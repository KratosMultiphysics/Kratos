import KratosMultiphysics
from geo_processor import GeoProcessor

import math
import os
import rasterio
import requests

class GeoData(GeoProcessor):

    def __init__(self, west_lon=14.22347, north_lat=42.45245, east_lon=14.22623, south_lat=42.45143, r_domain=1000):
        
        # bounding box coordinates
        self.north_lat = north_lat
        self.east_lon  = east_lon
        self.south_lat = south_lat
        self.west_lon  = west_lon
        
        self.metadata = dict()
        self.scale_x = None
        self.scale_y = None

        self.r_domain = r_domain    # TODO: just for a test. evaluate it


    def DownloadAsterGDEM(self, file_out, latitude, longitude, radius=1000, username="user", password="password"):
        """ function to download the ASTER GDEM file from a couple of coordinates

        Note:
            - it is necessary a valid account on https://earthdata.nasa.gov/
            - API documentation: https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html
            - raster mosaic: https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html

        Args:
            file_out: path of the output file
            latitude: latitude of the point of interest
            longitude: longitude of the point of interest
            radius: the radius (in meters) of the circle must be between 10 and 6000000
            username: from earthdata account
            password: from earthdata account

        Returns:
            - at least one file .TIF
            - a string with the DEM file (a single file or the merged one)
        """

        # NOTE: The bounding box parameters must be 4 comma-separated numbers: lower left longitude, lower left latitude, upper right longitude, upper right latitude.
        # req = requests.get('https://cmr.earthdata.nasa.gov/search/granules.json?short_name=ASTGTM&version=003&page_size=2000&pageNum=1&bounding_box=13.6,42.2,13.7,42.3').json()['feed']['entry'] 

        # NOTE: Circle defines a circle area on the earth with a center point and a radius.
        # req = requests.get('https://cmr.earthdata.nasa.gov/search/granules.json?short_name=ASTGTM&version=003&page_size=2000&pageNum=1&circle=14.2143295,42.4582503,1000').json()['feed']['entry'] 
        req = requests.get('https://cmr.earthdata.nasa.gov/search/granules.json?short_name=ASTGTM&version=003&page_size=2000&pageNum=1&circle={},{},{}'.\
                            format(longitude, latitude, radius)).json()['feed']['entry']
        fileList = [g['links'][0]['href'] for g in req]

        # list with ASTER GDEM file names
        self.aster_file_name = []

        for url in fileList:
            name_file = url.split("/")[-1]
            f_out = os.path.join(file_out, name_file)
            self.aster_file_name.append(f_out)

            with requests.Session() as session:
                r1 = session.request('get', url)
                r = session.get(r1.url, auth=(username, password))
                if r.ok:
                    # open(file_out, "wb").write(r.content)
                    with open(f_out, "wb") as dest:
                        dest.write(r.content)
                else:
                    KratosMultiphysics.Logger.PrintWarning("GeoData", "Error with the file {}".format(name_file))
        
        if (len(self.aster_file_name) > 1):
            # creating a raster mosaic (with two or more ASTER-GDEM files)

            from rasterio.merge import merge
            from rasterio.plot import show
            import glob

            path_merged = file_out + "/ASTER_GDEM_merged.tif"

            # search criteria to select the DEM files
            list_tif_files = os.path.join(file_out, "*.tif")

            # list all dem files with glob() function
            dem_fps = glob.glob(list_tif_files)

            # list for datafiles that will be part of the mosaic
            src_files_to_mosaic = []
            for fp in dem_fps:
                src = rasterio.open(fp)
                src_files_to_mosaic.append(src)

            # "merge" function returns a single mosaic array and the transformation info
            mosaic, out_trans = merge(src_files_to_mosaic)

            # copy the metadata
            out_meta = src.meta.copy()
            
            # update the metadata
            out_meta.update({"driver": "GTiff",
                            "height": mosaic.shape[1],
                            "width": mosaic.shape[2],
                            "transform": out_trans,
                            "crs": out_meta["crs"]
                            })

            # save the mosaic raster
            with rasterio.open(path_merged, "w", **out_meta) as dest:
                dest.write(mosaic)
            
            # merged file path
            return path_merged
        
        else:
            # single file path
            return f_out


    def ComputeBbox(self, lat_center=42, lon_center=14, radius=1000):
        """ function to compute the bounding box from a given radius

        Note:
            - the coordinate of the square containing the circle are calculated
            - reference: https://stackoverflow.com/questions/48440092/geo-circle-to-rectangle-coordinates

        Args:
            lat_center: latitude of the center of the domain
            lon_center: longitude of the center of the domain
            radius: radius (in meters) of the circle within the bounding box

        Returns:
            the bounding box coordinates (west, north, east, south)
        """

        R = 6378137    # radius of Earth in m

        north = lat_center + math.degrees(radius / R)
        east  = lon_center + math.degrees(radius / R / math.cos(math.radians(lat_center)))

        south = lat_center - math.degrees(radius / R)
        west  = lon_center - math.degrees(radius / R / math.cos(math.radians(lat_center)))

        self.north_lat = north
        self.east_lon  = east
        self.south_lat = south
        self.west_lon  = west

        # return (west, north, east, south)


    def CropAsterGDEM(self, file_in):
        """ function to crop the ASTER GDEM file inside the bounding box

        Note:
            - the bounding box must be declared first (west_lon, north_lat, east_lon, south_lat)
            - reference: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html

        Args:
            file_in: DEM file to crop

        Returns:
            - one file .TIF
            - a string with the cropped DEM file
        """

        import rasterio.mask

        file_crop = os.path.splitext(file_in)
        file_crop = file_crop[0] + "_CROP" + file_crop[1]

        # check if the bbox is inside the domain
        if (self.west_lon >= self.east_lon) or (self.north_lat <= self.south_lat):
            KratosMultiphysics.Logger.PrintWarning("GeoData", "Error in bbox coordinates. Check it!")
            return
        
        with rasterio.open(file_in) as src:
            matrix = src.transform
            tl_x, tl_y = matrix * (0, 0)                    # tl means Top Left
            br_x, br_y = matrix * (src.width, src.height)   # br means Bottom Right

            # check if bbox is in domain
            if ((self.west_lon < tl_x) or (self.north_lat > tl_y) or (self.east_lon > br_x) or (self.south_lat < br_y)):
                KratosMultiphysics.Logger.PrintWarning("GeoData", "The bbox is out of domain!")
                return
            
            # bbox generation
            geoms = [{'type': 'Polygon',
                    'coordinates': [[
                                    [self.west_lon, self.south_lat],
                                    [self.east_lon, self.south_lat],
                                    [self.east_lon, self.north_lat],
                                    [self.west_lon, self.north_lat]]]
                    }]

            out_image, out_transform = rasterio.mask.mask(src, geoms, crop=True)
            out_meta = src.meta

        # save the resulting raster
        out_meta.update({"driver": "GTiff",
                        "height": out_image.shape[1],
                        "width": out_image.shape[2],
                        "transform": out_transform})

        # save the TIF file
        with rasterio.open(file_crop, "w", **out_meta) as dest:
            dest.write(out_image)
        
        return file_crop


    def AsterGDEMtoOBJ(self, file_in, file_out=None):
        """ function to convert the ASTER GDEM file to an OBJ file

        Args:
            file_in: file DEM to convert
            file_out: output file. If it is non declared. it is saved with the same name as "file_in"

        Returns:
            - an OBJ file
            - a string with the OBJ file
        """

        if not file_out:
            file_out = os.path.splitext(file_in)[0]
            file_out = file_out + ".obj"

        self._get_metadata(file_in)
        # print("METADATA:")        # TODO: delete it
        # print(self.metadata)      # TODO: delete it
        self.size = (self.metadata["width"], self.metadata["height"])

        self.vertices = []  # initialize "vertices" list
        mid_x, mid_y = (self.metadata["center_x"], self.metadata["center_y"])

        for x, y, z in self._load_raster_xyz(file_in):
            self.scale_x = self.metadata["delta_x"] / self.metadata["size_x"]
            self.scale_y = self.metadata["delta_y"] / self.metadata["size_y"]

            # print("scale_x: ", self.scale_x)                        # TODO: delete it
            # print("scale_y: ", self.scale_y)                        # TODO: delete it
            # print("metadata[\"tl_x\"]: ", self.metadata["tl_x"])    # TODO: delete it
            # print("metadata[\"br_y\"]: ", self.metadata["br_y"])    # TODO: delete it
            
            # update coordinates with correct scale and move the coordinate with (0,0) to the bottom left
            x = (x-self.metadata["tl_x"]) * self.scale_x
            y = (y-self.metadata["br_y"]) * self.scale_y
            
            # add vertex
            self.vertices.append((x, y, z))

        # save the OBJ file
        self._write_file(file_out)

        return file_out


    def DownloadBuildingsOSM(self, file_out):
        """ function to download buildings from OpenStreetMap

        Note:
            - the bounding box must be declared first (west_lon, north_lat, east_lon, south_lat)
            - reference website: https://www.openstreetmap.org
            - API documentation: https://wiki.openstreetmap.org/wiki/Overpass_API

        Args:
            file_out: path of the output file

        Returns:
            - a JSON file with the buildings in the area of interest
            - geojson variable
        """

        import json

        # OSM bounding-box: (south, west, north, east)
        bbox = (self.south_lat, self.west_lon, self.north_lat, self.east_lon)

        # download only Way and Relations with the tag "building" in a bounding-box
        # also Nodes are downloaded: is useful for the coordinates
        overpass_url = "http://overpass-api.de/api/interpreter"
        overpass_query = """
                        [out:json];
                        (
                            node""" + str(bbox) + """;
                            way["building"]""" + str(bbox) + """;
                            rel["building"]""" + str(bbox) + """;
                        );
                        out geom;
        """

        response = requests.get(overpass_url, 
                                params={'data':overpass_query})
        
        name_file, extension = os.path.splitext(file_out)
        if not (extension.lower().endswith(".json")):
            file_out = name_file + ".json"

        # save JSON file
        with open (file_out, "w") as fo:
            fo.write(response.text)
        
        return json.loads(response.text)


    def GeoJSONtoOBJ(self, geojson, file_out):
        """ function to convert a GeoJSON into OBJ file

        Note:
            procedure reference: https://towardsdatascience.com/loading-data-from-openstreetmap-with-python-and-the-overpass-api-513882a27fd0

        Args:
            geojson: GeoJSON information
            file_out: path of the OBJ file
        
        Returns:
            an OBJ file with the 3D geometry of the buildings
        
        Limitations:
            - the buildings are obtained through an upward extrusion of the ground footprint (buildings with tapering will not be represented correctly)
            - in the case of multipolygon (in JSON file) only the external perimeter is considered (buildings with internal courtyards will not be represented correctly)
            - adjacent buildings are combined into a single building (shared surfaces are deleted)
            - buildings adjacent to each other but with different heights are merged into a single building with the height equal to the tallest building
            - overlapping buildings are merged
        """

        dict_building = {}      # key: number of building; value: coordinates of nodes in 2D (list of tuple)
        h_building = {}         # key: number of building; value: height
        h_base_building = {}    # key: number of building; value: base_height
        n_building = 0
        h_per_level = 3.0   # by default: 3 meters per level
        h_default = 9.0     # default height if not specified

        for element in geojson["elements"]:
            coords = []
            is_building = False
            h_current_building = 0.0
            h_base = 0.0

            debug_string = ""

            if ("tags" in element.keys()):
                if ("building" in element["tags"].keys()):
                    is_building = True
                    if ("building:height" in element["tags"].keys()):
                        debug_string = "building:height -> " + element["tags"]["building:height"]	# DEBUG
                        h_value = self._check_type((element["tags"]["building:height"]).replace(" m", ""), float)
                        # h_value = (element["tags"]["building:height"]).replace(" m", "")
                        # h_current_building = float(h_value)
                        h_current_building = h_value
                    
                    elif ("height" in element["tags"].keys()):
                        debug_string = "building:height -> " + element["tags"]["height"]			# DEBUG
                        h_value = self._check_type((element["tags"]["height"]).replace(" m", ""), float)
                        # h_value = (element["tags"]["height"]).replace(" m", "")
                        # h_current_building = float(h_value)
                        h_current_building = h_value
                    
                    elif ("building:levels" in element["tags"].keys()):
                        debug_string = "building:levels -> " + element["tags"]["building:levels"]	# DEBUG
                        levels = self._check_type((element["tags"]["building:levels"]).replace(" m", ""), int)
                        h_current_building = levels * h_per_level
                    
                    elif ("level" in element["tags"].keys()):
                        debug_string = "building:level -> " + element["tags"]["level"]			# DEBUG
                        levels = self._check_type((element["tags"]["level"]).replace(" m", ""), int)
                        h_current_building = levels * h_per_level
                    
                    elif ("levels" in element["tags"].keys()):
                        debug_string = "building:levels -> " + element["tags"]["levels"]			# DEBUG
                        levels = self._check_type((element["tags"]["levels"]).replace(" m", ""), int)
                        h_current_building = levels * h_per_level
                    
                    else:
                        # missing height information
                        h_current_building = h_default      # value by default
                        debug_string = "Not enough argument! We will set the Building as h={}m...".format(h_current_building)	# DEBUG
                    
                    # base height
                    if ("min_height" in element["tags"].keys()):
                        h_base = self._check_type((element["tags"]["min_height"]).replace(" m", ""), float)
                        # h_value = (element["tags"]["min_height"]).replace(" m", "")
                        # h_base = float(h_value)
                        # print("h_base: ", h_base)     # DEBUG TODO: delete it
                
                # Polygon -> "geometry"
                if ("geometry" in element.keys()) and is_building:
                    # print(debug_string)	                            # DEBUG TODO: delete it
                    # print("\tid: ", element["id"])	                # DEBUG TODO: delete it
                    # print("\tH = {} m".format(h_current_building))	# DEBUG TODO: delete it
                    
                    n_building += 1
                    for coords in element["geometry"]:
                        if (n_building in dict_building.keys()):
                            dict_building[n_building].append((coords["lon"], coords["lat"]))
                        else:
                            dict_building[n_building] = [(coords["lon"], coords["lat"])]
                    
                    h_building[n_building] = h_current_building
                    h_base_building[n_building] = h_base
                
                # Multipolygon -> "members"
                elif ("members" in element.keys()) and is_building:
                    # only outer will be considered
                    # print(debug_string)	                            # DEBUG TODO: delete it
                    # print("\tid: ", element["id"])	                # DEBUG TODO: delete it
                    # print("\tH = {} m".format(h_current_building))	# DEBUG TODO: delete it

                    outer_dict = dict()     # dictionaty with "outer" member. key: n_outer; value: list of node coordinates
                    n_outer = 0             # internal id of "outer"

                    for member in element["members"]:
                        if member["role"] == "outer":
                            n_outer += 1
                            for geom in member["geometry"]:
                                if n_outer in outer_dict.keys():
                                    outer_dict[n_outer].append((geom["lon"], geom["lat"]))
                                else:
                                    outer_dict[n_outer] = [(geom["lon"], geom["lat"])]

                    temp_list = list()
                    temp_bool = True
                    key_del = 0
                    
                    while True:		# infinite loop
                        if not outer_dict:
                            break	# break if the dictionary is empty

                        for n_outer, coords in outer_dict.items():
                            if not temp_list:
                                # when temp_list is empty (in the first cycle)
                                for coord in coords:
                                    temp_list.append(coord)
                                key_del = n_outer
                                temp_bool = True
                                break

                            # if (not temp_list) or (temp_list[-1] == coords[0]):	# if temp_list is empty or if the last value in temp_list == the first in "coords"
                            if (temp_list[-1] == coords[0]):	# if the last value in temp_list == the first in "coords"
                                for coord in coords[1:]:
                                    temp_list.append(coord)
                                key_del = n_outer
                                temp_bool = True
                                break

                            if temp_list[-1] == coords[-1]:
                                for coord in reversed(coords[:-1]):
                                    temp_list.append(coord)
                                key_del = n_outer
                                temp_bool = True
                                break
                        
                        else:
                            # break if none of the conditions inside the "for" loop is true
                            # NOTE: if the "outer_dict" is not empty, there may be another building not communicating with
                            #       the current building (in this version multi-buildings are not managed within a Multipolygon
                            break
                        
                        if temp_bool:
                            del outer_dict[key_del]
                    
                    n_building += 1
                    for coords in temp_list:
                        if (n_building in dict_building.keys()):
                            dict_building[n_building].append((coords[0], coords[1]))
                        else:
                            dict_building[n_building] = [(coords[0], coords[1])]
                    
                    h_building[n_building] = h_current_building
                    h_base_building[n_building] = h_base


        ################
        ### shapely  ###
        ################
        from shapely.geometry import Polygon
        from shapely.ops import unary_union		# boolean union

        print("[DEBUG][geo_data][GeoJSONtoOBJ] n_buildings: ", n_building, "\n")
        n_building += 1     # union building will start from this id


        i_poly = 0
        while (i_poly < n_building):

            # fill dict_poly
            dict_poly = {}		# key: buildind id; value: Polygon n-th
            for id_building, nodes in dict_building.items():
                if (len(nodes) < 3):
                    continue
                dict_poly[id_building] = Polygon(nodes)

            i_poly += 1
            if not (i_poly in dict_poly.keys()):
                # if this id is not in dict_poly key, switch to the next id
                continue

            poly_ith = dict_poly[i_poly]

            del_buildind = set()	# building id to delete

            for j_poly, poly_jth in dict_poly.items():
                if (i_poly == j_poly):
                    # we skip if it is the same Building
                    continue

                # if (poly_jth.within(poly_ith)):
                # 	print("Building{} in Building{}".format(j_poly, i_poly))

                # elif (poly_jth.intersects(poly_ith)):
                if (poly_jth.intersects(poly_ith)):
                    # print("Building{} intersects Building{}".format(j_poly, i_poly))    # TODO: delete it
                    u = unary_union([poly_ith, poly_jth])
                    if (u.geometryType() == "MultiPolygon"):
                        continue
                    
                    # we store (x,y) in coord_xy and later in the dictionary
                    coord_xy = [(x,y) for x, y in zip(u.exterior.coords.xy[0], u.exterior.coords.xy[1])]
                    dict_building[n_building] = coord_xy

                    h_base_building[n_building] = min(h_base_building[i_poly], h_base_building[j_poly])
                    h_building[n_building] = max(h_building[i_poly], h_building[j_poly])
                    
                    del_buildind.add(i_poly)
                    del_buildind.add(j_poly)
                    n_building += 1

                    break

            for del_b in del_buildind:
                del dict_building[del_b]


        ################
        ### triangle ###
        ################
        import triangle

        str_obj =  "# File created by Nicola Germano\n"
        str_obj += "# e-mail: nicola.germano@unich.it\n"
        str_obj += "# Department of Engineering and Geology (INGEO)\n"
        str_obj += "# University \"Gabriele d'Annunzio\" - Chieti-Pescara\n\n"
        n_nodes = 0
        for i_building, nodes in dict_building.items():

            coord_2D = []
            segments = []
            curr_nodes = 0	# Nodes in current Building

            # if (len(nodes) == 2):
            if (len(nodes) < 4):    # <4 means <3 because the first and the last are the same node
                print(nodes)                                                # TODO: delete it
                print("Building number: {} skipped!".format(i_building))    # TODO: delete it
                nodes = nodes[0]
                continue

            str_obj += "\no Building{}\n".format(i_building)

            # bottom nodes for "triangle"
            for n_id, node in enumerate(nodes[:-1]):	# the first and the last Node are the same Node (in a closed Polygon)
                coord_2D.append((node[0], node[1]))
                segments.append([n_id, n_id+1])
            segments[-1][1] = 0		# the last pair must be [last_node, first_node]
            
            A = dict(vertices=coord_2D, segments=segments)
            triangle_dict = triangle.triangulate(A, "p")
            
            # with an error in "triangle" the current building is skipped
            if not ("triangles" in triangle_dict.keys()):
                print("[DEBUG][geo_data][GeoJSONtoOBJ] coord_2D: ", coord_2D)
                print("[DEBUG][geo_data][GeoJSONtoOBJ] \"triangles\" not in triangle_dict. Skipped!")
                continue

            for node in nodes[:-1]:     # the first and the last Node are the same Node (in a closed Polygon)
                curr_nodes += 1
                # check if self.metadata is filled with ASTER-GDEM information
                if (not ("tl_x" in self.metadata.keys())) or (not ("br_y" in self.metadata.keys())):
                    KratosMultiphysics.Logger.PrintWarning("GeoData", "Error with metadata. ASTER-GDEM file must be set first!")
                    return 0
                
                # translation of nodes in accordance with values from ASETR-GDEM
                node_x = node[0] - self.metadata["tl_x"]
                node_y = node[1] - self.metadata["br_y"]

                # check if self.scale_x and self.scale_y are filled with ASTER-GDEM information
                if (self.scale_x == None) or (self.scale_y == None):
                    KratosMultiphysics.Logger.PrintWarning("GeoData", "Error with scale value. Please check scale on x/y from ASTER-GDEM file.")
                    return 0
                
                str_obj += "v {:20} {:20} {:5}".format(node_x * self.scale_x,
                                                       node_y * self.scale_y,
                                                       h_base_building[i_building])
                # comment in OBJ file with geographical coordinates
                str_obj += "    # lat: {}; lon: {}\n".format(node_y + self.metadata["br_y"],
                                                           node_x + self.metadata["tl_x"])


            """ NOTE: comment this part of code to disable bottom faces """
            # # bottom faces
            # for face in triangle_dict["triangles"]:
            # 	str_obj += "f {} {} {}\n".format(face[0] + 1 + n_nodes,
            # 									 face[1] + 1 + n_nodes,
            # 									 face[2] + 1 + n_nodes)
            
            n_bottom_nodes = len(nodes) - 1

            # top nodes
            for n_id, node in enumerate(nodes[:-1]):	# the first and the last Node are the same Node
                curr_nodes += 1
                
                # translation of nodes in accordance with values from ASETR-GDEM
                node_x = node[0] - self.metadata["tl_x"]
                node_y = node[1] - self.metadata["br_y"]

                str_obj += "v {:20} {:20} {:5}".format(node_x * self.scale_x,
                                                       node_y * self.scale_y,
                                                       h_building[i_building])
                # comment in OBJ file with geographical coordinates
                str_obj += "    # lat: {}; lon: {}\n".format(node_y + self.metadata["br_y"],
                                                           node_x + self.metadata["tl_x"])
            
            # top faces
            for face in triangle_dict["triangles"]:
                str_obj += "f {} {} {}\n".format(face[0] + 1 + n_nodes + n_bottom_nodes,
                                                 face[1] + 1 + n_nodes + n_bottom_nodes,
                                                 face[2] + 1 + n_nodes + n_bottom_nodes)
            
            # lateral faces
            for i in range(n_bottom_nodes-1):	# the last faces will be processed after the for loop
                j = i + n_bottom_nodes
                str_obj += "f {} {} {}\n".format(i+1+n_nodes, (i+1)+1+n_nodes, j+1+n_nodes)
                str_obj += "f {} {} {}\n".format((i+1)+1+n_nodes, (j+1)+1+n_nodes, j+1+n_nodes)
            ## last lateral faces
            str_obj += "f {} {} {}\n".format((n_bottom_nodes-1)+1+n_nodes, 0+1+n_nodes, (n_bottom_nodes-1+n_bottom_nodes)+1+n_nodes)
            str_obj += "f {} {} {}\n".format(0+1+n_nodes, n_bottom_nodes+1+n_nodes, (n_bottom_nodes-1+n_bottom_nodes)+1+n_nodes)
            
            n_nodes += curr_nodes

        name_file, extension = os.path.splitext(file_out)
        if not (extension.lower().endswith(".obj")):
            file_out = name_file + ".obj"
        
        # save JSON file
        with open(file_out, "w") as fo:
            fo.write(str_obj)


### --- auxiliary functions --- ### -------------------------------------

    def _get_metadata(self, filename):
        """ fills "metadata" with the information obtained from the DEM file """
        with rasterio.open(filename, 'r') as src:
            matrix = src.transform
            tl_x, tl_y = matrix * (0, 0)                    # tl means Top Left
            br_x, br_y = matrix * (src.width, src.height)   # br means Bottom Right

            self.metadata["center_x"] = (tl_x + br_x) / 2.0
            self.metadata["center_y"] = (tl_y + br_y) / 2.0
            self.metadata["width"] = src.width
            self.metadata["height"] = src.height
            self.metadata["size_x"] = (br_x - tl_x) / src.width
            self.metadata["size_y"] = (tl_y - br_y) / src.height

            # domain limit
            self.metadata["tl_x"] = tl_x
            self.metadata["tl_y"] = tl_y
            self.metadata["br_x"] = br_x
            self.metadata["br_y"] = br_y

            # delta in x and y (in meters)
            self.metadata["delta_x"] = self._measure(tl_y, tl_x, tl_y, br_x) / src.width
            self.metadata["delta_y"] = self._measure(tl_y, tl_x, br_y, tl_x) / src.height


    def _measure(self, lat1, lon1, lat2, lon2):
        """ given the coordinates of two points, the distance in meters is calculated """
        R = 6378.137	# radius of Earth in km
        dLat = lat2 * math.pi / 180 - lat1 * math.pi / 180
        dLon = lon2 * math.pi / 180 - lon1 * math.pi / 180
        a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(lat1 * math.pi/180) * math.cos(lat2 * math.pi/180) * math.sin(dLon/2) * math.sin(dLon/2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        d = R * c
        return (d * 1000)   # meters


    def _load_raster_xyz(self, filename):
        """ return x,y,z tuples """
        with rasterio.open(filename, 'r') as src:
            matrix = src.transform
            # read per scan line
            for row in range(0, src.height):
                window = ((row, row+1), (0, src.width))
                data = src.read(window=window)
                this_row = data[0][0]
                for column in range(0, src.width):
                    x, y = matrix * (column, row)
                    yield x, y, this_row[column]


    def _write_file(self, filename):
        """ writes OBJ format. Each pixel converted to 2 tris """
        with open(filename, "w") as fo:

            # vertices
            for x, y, z in self.vertices:
                fo.write("v {:20} {:20} {:5}\n".format(x, y, z))

            # faces
            faces = 0
            width, height = self.size
            for y in range(0, height-1):
                for x in range(0, width-1):
                    tl = self._vertex_num(x,y)
                    tr = tl + 1
                    bl = tl + width
                    br = bl + 1
                    fo.write("f {} {} {}\n".format(tl, tr, bl))
                    fo.write("f {} {} {}\n".format(tr, br, bl))
                    faces += 2


    def _vertex_num(self, x, y):
        """ Computes vertex number from x and y """
        width, _ = self.size
        return 1 + (width*y) + x


    def _check_type(self, var, _type):
        """ check if "var" is a valid variable """
        try:
            _type(var)
        except ValueError:
            return 1
        else:
            return _type(var)
