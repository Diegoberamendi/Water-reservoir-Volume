from osgeo import ogr
import os
from shapely.ops import linemerge
from shapely.geometry import Polygon
from shapely import wkt
import pandas as pd
import matplotlib.pyplot as plt


def hyp_curve(file_name, dst_file_name, curve_name, ref_level, to_excel=True, save_png=True):
    """"
    This function calculates the accumulated volume of a water-reservoir given a certain area and elevation
    whose values are taken from a previous generated contour shapefile.The method employed is the
     Trapezoidal method also called Area - End Method. The calculation is performed in
    a Dataframe from which the hypsometric curves are drafted with matplotlib. The hypsometric curves encompasses
    Area vs Elevation and Volume vs Elevation relationship.

    :param
    : file_name: str. Contour shapefile (.shp).
    : dst_file_name : str. Polygon (.shp) that contains the area of the water reservoir which will be used to clip file_name.
    : curve_name : str. Clipped version of file_name (.shp).
    : ref_level : Float. Max elevation which sets the upper limit of the volume computation.
    : to_excel : Boolean. If it is True then an Excel file is saved. Default is True.
    : save_png : Boolean. If it is True then an image.png is saved. Default is True.
    :return
    : dt_proc : Pandas dataframe containing the calculation of Water-Reservoir volume.

    """

    # Clipping the Contour shapefile

    global geo_merged, geo_type
    os.system('ogr2ogr -clipsrc' + ' ' + dst_file_name + ' ' + curve_name + ' ' + file_name)

    print('ogr2ogr -clipsrc' + ' ' + dst_file_name + ' ' + file_name + ' ' + curve_name)
    print('The contour shapefile was clipped')

    # Getting the Drive that allows the operations with the Contour shapefile

    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.Open(curve_name, 0)

    # Getting the layer of the clipped file
    layer = data_source.GetLayer(0)

    # Getting the number of features containing in the layer
    num_elements = layer.GetFeatureCount()

    print('Number of elements:', ' ', num_elements)

    # Empty lists to store the retrieved information
    id_list = []
    elements_list = []
    geometry_list = []
    elevation_list = []
    area_list = []

    # For loop over the contour lines

    for i in range(0, num_elements):

        feature = layer.GetNextFeature()

        # Getting the id and elevation of every contour line
        id = feature.GetField('ID')
        elevation = feature.GetField('ELEV')

        if elevation <= ref_level:

            geometry = feature.GetGeometryRef()
            # Displaying the type of geometry
            type_geometry = geometry.GetGeometryName()
            geometry_elements = geometry.GetGeometryCount()

            if geometry.GetArea() == 0:

                # The geometry is converted to Well Known Text(WKT) format
                geom_wkt = geometry.ExportToWkt()
                # From this section Shapely is used
                # The WKT geometry has to be loaded to enable operations with Shapely
                load_geom = wkt.loads(geom_wkt)
                # The components of the Multilinestring are merged to create linestrings
                if geometry_elements != 0:
                    merged_geometry = linemerge(load_geom)
                    # geo_merged gives the number of elements before merging
                    geo_merged = len(load_geom.geoms)
                    # geo_type gives the type of geometry after merging
                    geo_type = merged_geometry.geom_type

                else:
                    merged_geometry = load_geom

                # Retrieving the coordinates of the linestrings
                coordinates_linestrings = list(merged_geometry.coords)

                # Creation of a Polygon
                polygon = Polygon(coordinates_linestrings)

                # Area calculation of the above created polygon
                area_polygon = polygon.area

                # Storage of data into the lists
                id_list.append(id)
                elements_list.append(geo_merged)
                geometry_list.append(geo_type)
                elevation_list.append(elevation)
                area_list.append(area_polygon)

            else:
                id_list.append(id)
                elements_list.append(geometry_elements)
                geometry_list.append(type_geometry)
                elevation_list.append(elevation)
                area_list.append(geometry.GetArea())

    # Creation of a Pandas Dataframe

    dict = {'ID': id_list, '# Elements': elements_list, 'Geometry': geometry_list, 'Elevation': elevation_list,
            'Area': area_list}
    # Sorting of data with respect elevation in ascended form
    dt = pd.DataFrame(dict).sort_values(by=['Elevation'], ascending=True)
    # print(dt)

    # Verify if there are duplicates elevation values, if so their areas are summed up
    dt_proc = dt.groupby(['Elevation'], as_index=False)['Area'].sum()
    # print(dt_proc)

    # Calculation of the Reservoir volume

    dt_proc['Diff_Elevation'] = dt_proc['Elevation'].diff()
    dt_proc['Avr_Area'] = (dt_proc['Area'] + dt_proc['Area'].shift()) / 2
    dt_proc['Volume'] = dt_proc['Diff_Elevation'] * dt_proc['Avr_Area']
    dt_proc['Ac_Volume'] = dt_proc['Volume'].cumsum()

    # Replacing NaN value of the first value of Ac_Volume with 0
    dt_proc['Ac_Volume'][0] = 0

    print(dt_proc)

    # Printing the total volumen and area respect the reference elevation .iloc[-1] is used to access the last value
    # of the indicated colum

    print('Total Volume = {:.2f} at elevation {}'.format(dt_proc['Ac_Volume'].iloc[-1], ref_level))
    print('Water-mirror area = {:.2f}'.format(dt_proc['Area'].iloc[-1]))

    # Plotting the curves
    # Identify the axes of the graph
    x_volume = dt_proc['Ac_Volume']
    x_area = dt_proc['Area']
    y = dt_proc['Elevation']
    y_area = dt_proc['Elevation']

    # Figure creation
    fig, ax = plt.subplots(figsize=(15, 10))

    ax.plot(x_volume, y, linestyle='-', linewidth='3.0', marker='o', color='blue', label='Elevation-Volume')

    # Add Grid Lines

    ax.grid(b=True, color='grey', linestyle='-', linewidth='1.0', alpha=1)

    ax.set_title(r'Hypsometric Curves', fontsize=20)
    ax.set_xlabel('Volume ($m^3$)', fontsize=20)
    ax.set_ylabel(' Elevation ($masl$)', fontsize=20)

    # Axis for Elevation vs Area curve

    ax_area = ax.twiny()

    ax_area.plot(x_area, y_area, linestyle='-', linewidth='3.0', marker='^', color='green', label='Elevation-Area')
    ax_area.grid(b=True, color='grey', linestyle='-.', linewidth='1.0', alpha=1)
    ax_area.set_xlabel(r'Area $(m^2)$', fontsize=20)
    ax_area.set_ylabel(r'Elevation ($masl$)', fontsize=20)
    ax_area.invert_xaxis()

    # Ask matplotlibe for the plotted objects and their labels to plot both legends in one box
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax_area.get_legend_handles_labels()

    ax_area.legend(lines + lines2, labels + labels2, loc='lower center')

    # Parameters of the ticks,
    ax.tick_params(labelsize=15)
    ax_area.tick_params(labelsize=15)

    plt.tight_layout()
    plt.show()

    # Exporting the Dataframe into an Excel file
    if to_excel: dt_proc.to_excel('Hyp_curves.xlsx')

    # Saving the curves as png file
    if save_png: plt.savefig('Hyp_graph.png')


# Data input

working_directory = r'D:\Proyectos\Choquechambi\Mailanco\Geoproceso'
os.chdir(working_directory)
file_name = os.path.join(working_directory, 'Countours_2_m.shp')
dst_file_name = os.path.join(working_directory, 'box_clipping.shp')
curve_name = 'countour_ax2.shp'
ref_level = 3226

# Calling the function
hyp_curve(file_name, dst_file_name, curve_name, ref_level)
