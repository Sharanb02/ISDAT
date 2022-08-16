import requests
import getpass
import socket 
import json
import zipfile
import io
import math
import os
import shutil
import pprint
import re
import time
import geopandas as gpd
import pandas as pd
import fiona
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import h5py
# To read KML files with geopandas, we will need to enable KML support in fiona (disabled by default)
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
from shapely.geometry import Polygon, mapping
from shapely.geometry.polygon import orient
from statistics import mean
from requests.auth import HTTPBasicAuth
%matplotlib inline
%matplotlib widget


def get_data(bbox,date_range,path):
    #Comments marked #### are changes to be made

    uid = input('Earthdata Login user name: ') # Enter Earthdata Login user name
    pswd = getpass.getpass('Earthdata Login password: ') # Enter Earthdata Login password
    email = input('Email address associated with Earthdata Login account: ') # Enter Earthdata login email 
    short_name = 'GLAH06'
    params = {
    'short_name': short_name
    }
    cmr_collections_url = 'https://cmr.earthdata.nasa.gov/search/collections.json'
    response = requests.get(cmr_collections_url, params=params)
    results = json.loads(response.content)
    versions = [el['version_id'] for el in results['feed']['entry']]
    latest_version = max(versions)
    #print('The most recent version of ', short_name, ' is ', latest_version) #### Need not print Version


    start_time = '00:00:00' 
    end_time = start_time
    temporal = date_range[0] + 'T' + start_time + 'Z' + ',' + date_range[1] + 'T' + end_time + 'Z'

    bounding_region = bbox ##input as array if bbox, filepath string if shapefile

    if isinstance(bounding_region,str):
        aoi = '2'
        # Use geopandas to read in polygon file
        # Note: a KML or geojson, or almost any other vector-based spatial data format could be substituted here.

        shapefile_filepath = bounding_region

        # Go from geopandas GeoDataFrame object to an input that is readable by CMR
        gdf = gpd.read_file(shapefile_filepath)

        # CMR polygon points need to be provided in counter-clockwise order. The last point should match the first point to close the polygon.

        # Simplify polygon for complex shapes in order to pass a reasonable request length to CMR. The larger the tolerance value, the more simplified the polygon.
        # Orient counter-clockwise: CMR polygon points need to be provided in counter-clockwise order. The last point should match the first point to close the polygon.

        poly = orient(gdf.simplify(0.05, preserve_topology=True).loc[0],sign=1.0)

        geojson = gpd.GeoSeries(poly).to_json() # Convert to geojson
        geojson = geojson.replace(' ', '') #remove spaces for API call

        #Format dictionary to polygon coordinate pairs for CMR polygon filtering
        polygon = ','.join([str(c) for xy in zip(*poly.exterior.coords.xy) for c in xy])
    else:
        aoi = '1'
        bounding_box = bounding_region 

    # Create CMR parameters used for granule search. Modify params depending on bounding_box or polygon input.

    granule_search_url = 'https://cmr.earthdata.nasa.gov/search/granules'

    if aoi == '1':
        # bounding box input:
        search_params = {
        'short_name': short_name,
        'version': latest_version,
        'temporal': temporal,
        'page_size': 100,
        'page_num': 1,
        'bounding_box': bounding_box
        }
    else:
        # If polygon file input:
        search_params = {
        'short_name': short_name,
        'version': latest_version,
        'temporal': temporal,
        'page_size': 100,
        'page_num': 1,
        'polygon': polygon,
        }

    granules = []
    headers={'Accept': 'application/json'}
    while True:
        response = requests.get(granule_search_url, params=search_params, headers=headers)
        results = json.loads(response.content)

        if len(results['feed']['entry']) == 0:
            # Out of results, so break out of loop
            break

        # Collect results and increment page_num
        granules.extend(results['feed']['entry'])
        search_params['page_num'] += 1

    print('There are', len(granules), 'granules of', short_name, 'version', latest_version, 'over my area and time of interest.')

    granule_sizes = [float(granule['granule_size']) for granule in granules]

    print(f'The average size of each granule is {mean(granule_sizes):.2f} MB and the total size of all {len(granules)} granules is {sum(granule_sizes):.2f} MB')

    # Query service capability URL 

    from xml.etree import ElementTree as ET

    capability_url = f'https://n5eil02u.ecs.nsidc.org/egi/capabilities/{short_name}.{latest_version}.xml'

    # Create session to store cookie and pass credentials to capabilities url

    session = requests.session()

    s = session.get(capability_url)

    response = session.get(s.url,auth=(uid,pswd))

    root = ET.fromstring(response.content)

    #collect lists with each service option

    subagent = [subset_agent.attrib for subset_agent in root.iter('SubsetAgent')]
    if len(subagent) > 0 :

        # variable subsetting
        variables = [SubsetVariable.attrib for SubsetVariable in root.iter('SubsetVariable')]  
        variables_raw = [variables[i]['value'] for i in range(len(variables))]
        variables_join = [''.join(('/',v)) if v.startswith('/') == False else v for v in variables_raw] 
        variable_vals = [v.replace(':', '/') for v in variables_join]

        # reformatting
        formats = [Format.attrib for Format in root.iter('Format')]
        format_vals = [formats[i]['value'] for i in range(len(formats))]
        format_vals.remove('')

        # reprojection options
        projections = [Projection.attrib for Projection in root.iter('Projection')]
        #print service information depending on service availability and select service options

    if len(subagent) < 1 :
        print('No services exist for', short_name, 'version', latest_version)
        agent = 'NO'
        bbox = ''
        time_var = ''
        reformat = ''
        projection = ''
        projection_parameters = ''
        coverage = ''
        Boundingshape = ''
    else:
        agent = ''
        subdict = subagent[0]
        if subdict['spatialSubsetting'] == 'true' and aoi == '1':
            Boundingshape = ''
            ss = input('Subsetting by bounding box, based on the area of interest inputted above, is available. Would you like to request this service? (y/n)')
            if ss == 'y': bbox = bounding_box
            else: bbox = '' 
        if subdict['spatialSubsettingShapefile'] == 'true' and aoi == '2':
            bbox = ''
            ps = input('Subsetting by geospatial file (Esri Shapefile, KML, etc.) is available. Would you like to request this service? (y/n)')
            if ps == 'y': Boundingshape = geojson
            else: Boundingshape = '' 
        if subdict['temporalSubsetting'] == 'true':
            ts = input('Subsetting by time, based on the temporal range inputted above, is available. Would you like to request this service? (y/n)')
            if ts == 'y': time_var = start_date + 'T' + start_time + ',' + end_date + 'T' + end_time 
            else: time_var = ''
        else: time_var = ''
        if len(format_vals) > 0 :
            print('These reformatting options are available:', format_vals)
            reformat = input('If you would like to reformat, copy and paste the reformatting option you would like (make sure to omit quotes, e.g. GeoTIFF), otherwise leave blank.')
            if reformat == 'n': reformat = '' # Catch user input of 'n' instead of leaving blank
        else: 
            reformat = ''
            projection = ''
            projection_parameters = ''
        if len(projections) > 0:
            valid_proj = [] # select reprojection options based on reformatting selection
            for i in range(len(projections)):
                if 'excludeFormat' in projections[i]:
                    exclformats_str = projections[i]['excludeFormat'] 
                    exclformats_list = exclformats_str.split(',')
                if ('excludeFormat' not in projections[i] or reformat not in exclformats_list) and projections[i]['value'] != 'NO_CHANGE': valid_proj.append(projections[i]['value'])
            if len(valid_proj) > 0:
                print('These reprojection options are available with your requested format:', valid_proj)
                projection = input('If you would like to reproject, copy and paste the reprojection option you would like (make sure to omit quotes), otherwise leave blank.')
                # Enter required parameters for UTM North and South
                if projection == 'UTM NORTHERN HEMISPHERE' or projection == 'UTM SOUTHERN HEMISPHERE': 
                    NZone = input('Please enter a UTM zone (1 to 60 for Northern Hemisphere; -60 to -1 for Southern Hemisphere):')
                    projection_parameters = str('NZone:' + NZone)
                else: projection_parameters = ''
            else: 
                print('No reprojection options are supported with your requested format')
                projection = ''
                projection_parameters = ''
        else:
            print('No reprojection options are supported with your requested format')
            projection = ''
            projection_parameters = ''

    #if len(subagent) > 0 :
     #   if len(variable_vals) > 0:
      #      v = input('Variable subsetting is available. Would you like to subset a selection of variables? (y/n)')
       # if v == 'y':
        #    print('The', short_name, 'variables to select from include:')
        #    print(*variable_vals, sep = "\n") 
        #    coverage = input('If you would like to subset by variable, copy and paste the variables you would like separated by comma (be sure to remove spaces and retain all forward slashes: ')
        #else: coverage = ''

    #no services selected
    if reformat == '' and projection == '' and projection_parameters == '' and time_var == '' and bbox == '' and Boundingshape == '':
        agent = 'NO'

    #Set NSIDC data access base URL
    base_url = 'https://n5eil02u.ecs.nsidc.org/egi/request'

    #Set the request mode to asynchronous if the number of granules is over 100, otherwise synchronous is enabled by default
    if len(granules) > 100:
        request_mode = 'async'
        page_size = 2000
    else: 
        page_size = 100
        request_mode = 'stream'

    #Determine number of orders needed for requests over 2000 granules. 
    page_num = math.ceil(len(granules)/page_size)

    print('There will be', page_num, 'total order(s) processed for our', short_name, 'request.')


    if aoi == '1':
        # bounding box search and subset:
        param_dict = {'short_name': short_name, 
        'version': latest_version, 
        'temporal': temporal, 
        'time': time_var, 
        'bounding_box': bounding_box, 
        'bbox': bbox, 
        'format': reformat, 
        'projection': projection, 
        'projection_parameters': projection_parameters, 
        'Coverage': '/Data_40HZ', 
        'page_size': page_size, 
        'request_mode': request_mode, 
        'agent': agent, 
        'email': email, }
    else:
        # If polygon file input:
        param_dict = {'short_name': short_name, 
        'version': latest_version, 
        'temporal': temporal, 
        'time': time_var, 
        'polygon': polygon,
        #'Boundingshape': Boundingshape, 
        'format': reformat, 
        'projection': projection, 
        'projection_parameters': projection_parameters, 
        'Coverage': '/Data_40HZ', 
        'page_size': page_size, 
        'request_mode': request_mode, 
        'agent': agent, 
        'email': email, }

    #Remove blank key-value-pairs
    param_dict = {k: v for k, v in param_dict.items() if v != ''}

    #Convert to string
    param_string = '&'.join("{!s}={!r}".format(k,v) for (k,v) in param_dict.items())
    param_string = param_string.replace("'","")

    #Print API base URL + request parameters
    endpoint_list = [] 
    for i in range(page_num):
        page_val = i + 1
        API_request = api_request = f'{base_url}?{param_string}&page_num={page_val}'
        endpoint_list.append(API_request)

    print(*endpoint_list, sep = "\n") 

    # Create an output folder if the folder does not already exist.
    path = path +'/'+date_range[0]+'--'+date_range[1]
    try:
        os.mkdir(path)
    except:
        None

    # Different access methods depending on request mode:

    if request_mode=='async':
        # Request data service for each page number, and unzip outputs
        for i in range(page_num):
            page_val = i + 1
            print('Order: ', page_val)

        # For all requests other than spatial file upload, use get function
            request = session.get(base_url, params=param_dict)

            print('Request HTTP response: ', request.status_code)

            # Raise bad request: Loop will stop for bad response code.
            request.raise_for_status()
            print('Order request URL: ', request.url)
            esir_root = ET.fromstring(request.content)
            print('Order request response XML content: ', request.content)

            #Look up order ID
            orderlist = []   
            for order in esir_root.findall("./order/"):
                orderlist.append(order.text)
            orderID = orderlist[0]
            print('order ID: ', orderID)

            #Create status URL
            statusURL = base_url + '/' + orderID
            print('status URL: ', statusURL)

            #Find order status
            request_response = session.get(statusURL)    
            print('HTTP response from order response URL: ', request_response.status_code)

            # Raise bad request: Loop will stop for bad response code.
            request_response.raise_for_status()
            request_root = ET.fromstring(request_response.content)
            statuslist = []
            for status in request_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            print('Data request ', page_val, ' is submitting...')
            print('Initial request status is ', status)

            #Continue loop while request is still processing
            while status == 'pending' or status == 'processing': 
                print('Status is not complete. Trying again.')
                time.sleep(10)
                loop_response = session.get(statusURL)

            # Raise bad request: Loop will stop for bad response code.
            loop_response.raise_for_status()
            loop_root = ET.fromstring(loop_response.content)

            #find status
            statuslist = []
            for status in loop_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            print('Retry request status is: ', status)
            if status == 'pending' or status == 'processing':
                continue

        #Order can either complete, complete_with_errors, or fail:
        # Provide complete_with_errors error message:
        if status == 'complete_with_errors' or status == 'failed':
            messagelist = []
            for message in loop_root.findall("./processInfo/"):
                messagelist.append(message.text)
            print('error messages:')
            pprint.pprint(messagelist)

        # Download zipped order if status is complete or complete_with_errors
        if status == 'complete' or status == 'complete_with_errors':
            downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
            print('Zip download URL: ', downloadURL)
            print('Beginning download of zipped output...')
            zip_response = session.get(downloadURL)
            # Raise bad request: Loop will stop for bad response code.
            zip_response.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                z.extractall(path)
            print('Data request', page_val, 'is complete.')
        else: print('Request failed.')

    else:
        for i in range(page_num):
            page_val = i + 1
            print('Order: ', page_val)
            print('Requesting...')
            request = session.get(base_url, params=param_dict)
            print('HTTP response from order response URL: ', request.status_code)
            request.raise_for_status()
            d = request.headers['content-disposition']
            fname = re.findall('filename=(.+)', d)
            dirname = os.path.join(path,fname[0].strip('\"'))
            print('Downloading...')
            open(dirname, 'wb').write(request.content)
            print('Data request', page_val, 'is complete.')

    # Unzip outputs
    for z in os.listdir(path): 
        if z.endswith('.zip'): 
            zip_name = path + "/" + z 
            zip_ref = zipfile.ZipFile(zip_name) 
            zip_ref.extractall(path) 
            zip_ref.close() 
            os.remove(zip_name) 
            # Clean up Outputs folder by removing individual granule folders 

    for root, dirs, files in os.walk(path, topdown=False):
        for file in files:
            try:
                shutil.move(os.path.join(root, file), path)
            except OSError:
                pass
        for name in dirs:
            os.rmdir(os.path.join(root, name))  


def data_to_csv(path_in):
    data=pd.DataFrame()
    a=os.listdir(path_in)
    try:
        os.mkdir(path_in+'/CSV') 
    except:
        None
    path_out = path_in+'/CSV'
    try:
        a.remove('.ipynb_checkpoints')
    except:
        None
    for fname in a:
        fname=path_in+'/'+fname
        try:
            with h5py.File(fname,'r') as f:
                df=pd.DataFrame()
                df['lat']=f['/Data_40HZ/Geolocation/d_lat'][:]
                df['lon']=f['/Data_40HZ/Geolocation/d_lon'][:]
                df['ele']=f['/Data_40HZ/Elevation_Surfaces/d_elev'][:]
                df['N']=f['/Data_40HZ/Geophysical/d_DEM_elv'][:]
                df['del_ellip']=f['/Data_40HZ/Geophysical/d_deltaEllip'][:]
                df['sat_flag'] = f['/Data_40HZ/Quality/sat_corr_flg'][:]
                df['sat_corr'] = f['/Data_40HZ/Elevation_Corrections/d_satElevCorr'][:]
                df=df[df['ele']<8900]
                data=data.append(df,ignore_index=False)
        except:
            print(fname+" has no relevant data")
            continue
    
    data.to_csv(path_out+'/gt.csv')
    

def ele_plot(region_path,end_time,rel_del=3):
    year, month, day = map(int, end_time.split('-'))
    start_time = date(year, month, day)+relativedelta(months=-rel_del)
    start_time = str(start_time)
    date_range=[start_time,end_time]
    data_path = region_path+'/data/GLAS'
    shp = region_path+'/shpfile/Siachen.shp'
    shape = region_path+'/shapefile/Siachen.shp'
    glims_sia = gpd.read_file(shape)
    if (os.path.isdir(data_path+'/'+date_range[0]+'--'+date_range[1])==False) or (len(os.listdir(data_path+'/'+date_range[0]+'--'+date_range[1]))<=1):
        print('Downloading Data...')
        GLget_data(shp,date_range,data_path)
        GLdata_to_csv(data_path+'/'+date_range[0]+'--'+date_range[1])
    else:
        print("Data aldready exists")

    data_path = data_path+'/'+date_range[0]+'--'+date_range[1]+'/CSV'
    df = pd.read_csv(data_path+'/gt.csv')
    beam = gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(df['lon'],df['lat']),crs = 'EPSG:4326')
    beam1=gpd.overlay(beam,glims_sia,how='intersection')
    # applying Saturation correction and projection to WGS84 Ellipsoid to make consistent with ATL06
    beam1=beam1[beam1['sat_flag']<3]
    beam1['ele'] = beam1['ele'] - beam1['del_ellip'] + beam1['sat_corr']
    fig, ax = plt.subplots(1, 1,figsize=(5,5))
    ax.set_facecolor('black')
    glims_sia.plot(ax=ax,edgecolor='gainsboro',color='white')
    divider = make_axes_locatable(ax)
    cax_ = divider.append_axes("right", size="5%", pad=0.05)
    cat_= divider.append_axes("bottom", size="5%", pad=0.5)
    cat_.plot([date_range[0],date_range[1]],np.zeros_like(['2019-05-01','2019-08-01']))
    ax.set_title("Land Ice Height (CRS:WGS 84 EPSG:4326)")
    left, width = .08,.5
    bottom, height = .1, .5
    right = left + width
    top = bottom + height
    ax.text(left,bottom,'Siachen Glacier',horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,bbox=dict(facecolor='gainsboro',alpha=0.3),color='white')                       
    beam1.plot('ele',cmap='Blues',markersize=0.1,ax=ax,legend=True,cax=cax_,vmin=3000)
