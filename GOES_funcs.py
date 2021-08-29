''' adapted from https://github.com/blaylockbk/pyBKB_v2/blob/master/BB_GOES16/mapping_GOES16_data.ipynb '''
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, glob, time, re
import s3fs, os, sys, datetime, dateutil, pytz
import netCDF4 as nc
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, BORDERS, STATES
import xarray, metpy
from skimage.color import rgb2hsv, hsv2rgb
from skimage import exposure

utc = pytz.UTC

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
            if((houri>=0)&(houri<11)):
                continue
        for file in files:
            filelist.append(file)
        
    return filelist
  
def process_files(filelist, year, day, limits, GLM, borders=False):
    ''' 
        create the figure
        and apply the projection we want (Mercator centered over the central US)
    '''
    
    plotfolder = "pngs/%d/%03d/"%(year, day)
    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)
        print("Folder was created: ", plotfolder)


    projection = ccrs.Mercator(central_longitude=(limits[0] + limits[1])/2.)
    fig = plt.figure(figsize=(20,16), dpi=300)
    ax  = fig.add_subplot(111, projection=projection, facecolor='black')
    plt.subplots_adjust(top=1., bottom=0., left=0., right=1.)
    
    # Use the anonymous credentials to access public data
    fs = s3fs.S3FileSystem(anon=True)

    for i, file in enumerate(filelist):
        t1 = time.time()

        ax.cla()
        #ax.set_extent([-125, -65., 20., 55], crs=ccrs.PlateCarree())
        ax.set_extent(limits, crs=ccrs.PlateCarree())
        ax.coastlines(resolution='50m', color='black', linewidth=0.5)

        date = process_ABI(fs.open(file, 'r'), fig, ax, projection, borders)


        if(GLM):
            '''
                the file is in the following pattern
                so that we can extract the timestamp info
            '''
            pattern = r'OR_ABI-L2-MCMIPC-M6_G16_s([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})_e([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})_c([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})'
            
            fname = file.split('\\')[-1]
            ''' extract the timestamps '''
            match = re.findall(pattern, fname)
            syear, sday, shour, smin, ssec, \
                eyear, eday, ehour, emin, esec, \
                cyear, cday, chour, cmin, csec = np.array(match[0], dtype=np.int32)
            
            sseconds   = sday*24.*3600. + shour*3600. + smin*60. + ssec/10.
            eseconds   = eday*24.*3600. + ehour*3600. + emin*60. + esec/10.

            start_date = datetime.datetime(syear-1, 12, 31, 0, 0, 0, 0)
            start_date = start_date + datetime.timedelta(seconds=sseconds)
            
            end_date   = datetime.datetime(eyear-1, 12, 31, 0, 0, 0, 0)
            end_date   = end_date + datetime.timedelta(seconds=eseconds)

            start_date = utc.localize(start_date)
            end_date   = utc.localize(end_date)

            GLM_list   = get_GLM(start_date, end_date, fs)
            
            plot_GLM(GLM_list, fs, fig, ax, start_date, end_date)


        year   = date.year
        yday   = date.timetuple().tm_yday
        month  = date.month
        day    = date.day
        hour   = date.hour
        minute = date.minute
        sec    = date.second
        plt.axis('off')
        #ax.set_extent(limits, crs=ccrs.PlateCarree())
        fig.savefig(plotfolder+"%d%02d%02d_%02d%02d%02d.png"%(year, month, day, hour, minute, sec), bbox_inches='tight', \
                    facecolor='black', dpi=300)
        print("%.2f"%(time.time() - t1))

def process_ABI(file, fig, ax, projection, borders):
    ''' stream the netCDF file by reading it from memory '''
    ncdata = nc.Dataset('name', 'r', memory=file.buffer.read())

    data   = xarray.open_dataset(xarray.backends.NetCDF4DataStore(ncdata))

    timeoff = ncdata.variables['t'][0]
    date   = datetime.datetime(2000, 1, 1, 12) + datetime.timedelta(seconds=float(timeoff))
    print("Processing %02d:%02d:%02.3f..."%(date.hour, date.minute, date.second), end='')

    R    = data['CMI_C02'].data
    G    = data['CMI_C03'].data
    B    = data['CMI_C01'].data

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
    G_true = np.clip(G_true, 0, 1)
    #G_true = np.maximum(G_true, 0)
    #G_true = np.minimum(G_true, 1)

    color_corr  = np.dstack([R, G_true, B])
    #color_corr = contrast_correction(color, 125)

    for ci in range(3):
        color_corr[:,:,ci] = exposure.equalize_adapthist(color_corr[:,:,ci], clip_limit=0.02)

    ## load the clean IR for nighttime
    cleanIR = data['CMI_C13'].data
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

    cf_data = data.metpy.parse_cf('CMI_C02')

    geos = cf_data.metpy.cartopy_crs
    x = cf_data.x
    y = cf_data.y

    data.close()
    R = []
    G = []
    B = []
    cleanIR = []
    
    ax.imshow(RGB_IR, origin='upper', extent=(x.min(), x.max(), y.min(), y.max()), transform=geos, interpolation='hermite')

    if borders:
        ax.add_feature(STATES, linewidth=0.3)
        ax.add_feature(BORDERS, linewidth=0.5)

    return date

def get_GLM(start_date, end_date, fs):
    syear = start_date.year
    sday  = start_date.timetuple().tm_yday
    shour = start_date.hour

    eyear = end_date.year
    eday  = end_date.timetuple().tm_yday
    ehour = end_date.hour

    slist = fs.ls('noaa-goes16/GLM-L2-LCFA/%d/%03d/%02d/'%(syear, sday, shour))
    elist = fs.ls('noaa-goes16/GLM-L2-LCFA/%d/%03d/%02d/'%(eyear, eday, ehour))

    filelist = np.unique(sorted([*slist, *elist]))
    
    pattern = r'OR_GLM-L2-LCFA_G16_s([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})_e([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})_c([0-9]{4})([0-9]{3})([0-9]{2})([0-9]{2})([0-9]{3})'

    dec31 = datetime.datetime(syear-1, 12, 31, 0, 0, 0, 0)

    GLM_list = []
    for file in filelist:
        fname = file.split('/')[-1]
        match = re.findall(pattern, fname)

        syear, sday, shour, smin, ssec, \
            eyear, eday, ehour, emin, esec, \
            cyear, cday, chour, cmin, csec = np.array(match[0], dtype=np.int32)
        
        sseconds   = sday*24.*3600. + shour*3600. + smin*60. + ssec/10.
        eseconds   = eday*24.*3600. + ehour*3600. + emin*60. + esec/10.
        
        start = utc.localize(datetime.datetime(syear-1, 12, 31, 0, 0, 0, 0))
        start = start + datetime.timedelta(seconds=sseconds)
        
        end   = utc.localize(datetime.datetime(eyear-1, 12, 31, 0, 0, 0, 0))
        end   = end + datetime.timedelta(seconds=eseconds)

        ## find if the date range falls inside the search range
        if((end>=start_date)&(start<=end_date)):
            GLM_list.append(file)

    return GLM_list


def plot_GLM(GLM_list, fs, fig, ax, start_date, end_date):
    for GLMfile in GLM_list:
        fbuff = fs.open(GLMfile, 'r')

        data = nc.Dataset('name', 'r', memory=fbuff.buffer.read())
        
        date = dateutil.parser.parse(data.time_coverage_start)

        #print("\tPlotting GLM at %02d:%02d:%06.3f..."%(date.hour, date.minute, date.second))

        event_time = data.variables['event_time_offset'][:]
        event_date = [date + datetime.timedelta(seconds=float(ti)) for ti in event_time]
        mask = np.array([ (tt>=start_date)&(tt<=end_date) for tt in event_date], dtype=np.bool)

        event_energy = np.log10(data.variables['event_energy'][:][mask])
        event_lat    = data.variables['event_lat'][:][mask]
        event_lon    = data.variables['event_lon'][:][mask]

        ax.scatter(event_lon, event_lat, s=2, c=event_energy, vmin=-15, vmax=-10, cmap='Reds_r', transform=ccrs.PlateCarree())

        
