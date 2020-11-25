''' adapted from https://github.com/blaylockbk/pyBKB_v2/blob/master/BB_GOES16/mapping_GOES16_data.ipynb '''
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, glob, time, re
import s3fs, os, sys, datetime
import netCDF4 as nc
from pyproj import Proj
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature

def contrast_correction(color, contrast):
    """
    Modify the contrast of an RGB
    See:
    https://www.dfstudios.co.uk/articles/programming/image-programming-algorithms/image-processing-algorithms-part-5-contrast-adjustment/

    Input:
        color    - an array representing the R, G, and/or B channel
        contrast - contrast correction level
    """
    F = (259*(contrast + 255))/(255.*259-contrast)
    COLOR = F*(color-.5)+.5
    COLOR = np.clip(COLOR, 0, 1)  # Force value limits 0 through 1.
    return COLOR


def get_data(year, day, daytime):
    ''' 
        download MCMPIC (5-min CONUS) images from GOES 16 for a given day/year
        Input:
            year    - year for data
            day     - day of the year 
            daytime - flag to only retrieve daytime images (11Z-01Z)
    '''
    outfolder = "GOES16/%d/%03d/"%(year, day)
    ''' check if folder exists '''
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
        print("Folder was created: ", outfolder)

    # Use the anonymous credentials to access public data
    fs = s3fs.S3FileSystem(anon=True)

    ''' get the number of hours available '''
    hours = fs.ls('noaa-goes16/ABI-L2-MCMIPC/%d/%03d/'%(year, day))

    print("Found %d hours"%len(hours))
    for hour in hours:
        # List specific files of GOES-16 CONUS data (multiband format) on a certain hour
        files = np.array(fs.ls(hour))

        houri = int(hour.split('/')[-1])

        if(daytime==True):
            if((houri>=1)&(houri<11)):
                continue

        print("Getting %d files for %02d"%(files.shape[0], houri))

        # Download the first file, and rename it the same name (without the directory structure)
        for file in files:
            fname = file.split('/')[-1]
            path = outfolder+fname

            if(not os.path.exists(path)):
                print("Downloading %s..."%(fname))
                fs.get(file, path)
            else: 
                print("Found %s. Skipping."%fname)

    print("All files downloaded")
    print()

def get_filelist(year, day, daytime_only):
    files = glob.glob("GOES16/%d/%03d/*.nc"%(year, day))

    '''
        the file is in the following pattern
        so that we can extract the timestamp info
    '''
    pattern = r'OR_ABI-L2-MCMIPC-M6_G16_s([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})_e([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})_c([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})'

    ''' 
        stream MCMPIC (5-min CONUS) images from GOES 16 for a given day/year
    '''
    # Use the anonymous credentials to access public data
    fs = s3fs.S3FileSystem(anon=True)

    ''' get the number of hours available '''
    hours = fs.ls('noaa-goes16/ABI-L2-MCMIPC/%d/%03d/'%(year, day))

    #print("Found %d hours"%len(hours))
    filelist = []
    for hour in hours:
        # List specific files of GOES-16 CONUS data (multiband format) on a certain hour
        files = np.array(fs.ls(hour))

        houri = int(hour.split('/')[-1])

        ''' if we only want the daytime images '''
        if(daytime_only==True):
            if((houri>=1)&(houri<11)):
                continue
        for file in files:
            fname = file.split('\\')[-1]
            
            ''' extract the timestamps '''
            match = re.findall(pattern, fname)
            syear, sday, shour, smin, ssec, \
                eyear, eday, ehour, emin, esec, \
                cyear, cday, chour, cmin, csec = np.array(match[0], dtype=np.int32)

            filelist.append(file)
        
    return filelist
  
def process_files(filelist, year, day, limits):
    ''' 
        create the figure
        and apply the projection we want (Mercator centered over the central US)
    '''
    
    plotfolder = "pngs/%d/%03d/"%(year, day)
    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)
        print("Folder was created: ", plotfolder)


    projection = ccrs.Mercator(central_longitude=(limits[0] + limits[1])/2.)
    fig = plt.figure(figsize=(15,12))
    ax  = fig.add_subplot(111, projection=projection, facecolor='black')
    plt.subplots_adjust(top=1., bottom=0., left=0., right=1.)
    ax.set_extent(limits, crs=ccrs.PlateCarree())
    #ax.set_extent([-100, -75., 20., 35], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)
    plt.axis('off')
    
    # Use the anonymous credentials to access public data
    fs = s3fs.S3FileSystem(anon=True)

    for i, file in enumerate(filelist):
        if(i==0):
            mapproj = make_pngs(fs.open(file,'r'), plotfolder, fig, ax, projection)
        else:
            mapproj = make_pngs(fs.open(file,'r'), plotfolder, fig, ax, projection, mapproj=mapproj)

def make_pngs(file, plotfolder, fig, ax, projection, mapproj=None):
    t1 = time.time()
    
    ''' stream the netCDF file by reading it from memory '''
    data = nc.Dataset('name', 'r', memory=file.buffer.read())

    timeoff = data.variables['t'][0]

    date   = datetime.datetime(2000, 1, 1, 12) + datetime.timedelta(seconds=float(timeoff))

    year   = date.year
    yday   = date.timetuple().tm_yday
    month  = date.month
    day    = date.day
    hour   = date.hour
    minute = date.minute
    sec    = date.second

    print("Processing %02d:%02d:%.3f..."%(hour, minute, sec), end='')
    sys.stdout.flush()

    R    = data.variables['CMI_C02'][:]
    G    = data.variables['CMI_C03'][:]
    B    = data.variables['CMI_C01'][:]

    # Apply range limits for each channel becuase RGB values must be between 0 and 1
    R = np.clip(R, 0, 1)
    G = np.clip(G, 0, 1)
    B = np.clip(B, 0, 1)

    # Apply the gamma correction
    gamma = 0.4
    R = np.power(R, gamma)
    G = np.power(G, gamma)
    B = np.power(B, gamma)

    # Calculate the "True" Green
    G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
    G_true = np.maximum(G_true, 0)
    G_true = np.minimum(G_true, 1)

    color  = np.dstack([R, G_true, B])
    color_corr = contrast_correction(color, 125)

    ## load the clean IR for nighttime
    cleanIR = data.variables['CMI_C13'][:]
    cleanIR[cleanIR==-1] = np.nan

    # Apply range limits for clean IR channel
    cleanIR = np.maximum(cleanIR, 90)
    cleanIR = np.minimum(cleanIR, 313)

    # Normalize the channel between a range
    cleanIR = (cleanIR-90)/(313-90)

    # Invert colors
    cleanIR = 1 - cleanIR

    # Lessen the brightness of the coldest clouds so they don't appear so bright near the day/night line
    cleanIR = cleanIR/1.5

    R = color_corr[:,:,0]
    G = color_corr[:,:,1]
    B = color_corr[:,:,2]

    RGB_IR = np.dstack([np.maximum(R, cleanIR), np.maximum(G, cleanIR), np.maximum(B, cleanIR)])

    # Satellite height
    sat_h = data.variables['goes_imager_projection'].perspective_point_height

    # Satellite longitude
    sat_lon = data.variables['goes_imager_projection'].longitude_of_projection_origin

    # Satellite sweep
    sat_sweep = data.variables['goes_imager_projection'].sweep_angle_axis

    # The projection x and y coordinates equals
    # the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
    X = data.variables['x'][:] * sat_h
    Y = data.variables['y'][:] * sat_h

    # Convert map points to latitude and longitude with the magic provided by Pyproj
    p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
    XX, YY = np.meshgrid(X, Y)
    lons, lats = p(XX, YY, inverse=True)
    lats[np.isnan(R)] = np.nan
    lons[np.isnan(R)] = np.nan

    ## Create a color tuple for pcolormesh
    rgb = RGB_IR[:,:-1,:] # Using one less column is very imporant, else your image will be scrambled! (This is the stange nature of pcolormesh)
    colorTuple = rgb.reshape((rgb.shape[0] * rgb.shape[1]), 3) # flatten array, becuase that's what pcolormesh wants.
    colorTuple = np.insert(colorTuple, 3, 1.0, axis=1) # adding an alpha channel will plot faster, according to stackoverflow. Not sure why.

    try:
        mapproj.update({'color': colorTuple})
    except:
        mapproj = ax.pcolormesh(lons, lats, R, color=colorTuple, antialiased=False, linewidth=0, transform=ccrs.PlateCarree())
        mapproj.set_array(None)
        ax.background_patch.set_facecolor('k')
  
    fig.savefig(plotfolder+"%d%02d%02d_%02d%02d%02d.png"%(year, month, day, hour, minute, sec), bbox_inches='tight', facecolor='black', dpi=150)
    # plt.close(fig)
    print("%.2f"%(time.time() - t1))

    return mapproj
