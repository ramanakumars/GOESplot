# GOESplot
Python code to plot GOES imagery

## Usage
Run using:
````sh
$ python run_GOES.py --date yyyy-mm-dd [--daytime] [--process_last] 
[--borders] [--limits LONMIN LONMAX LATMIN LATMAX ]
````
Arguments:

 ````date````: 

 Date to retrieve image (format: yyyy-mm-dd). Required.

 Use 'today' to process for the current date.

 ````daytime````: 

 Use only daytime data (10Z-01Z) (default: False)

 ````process_last````:
 Process the last image only. Use in conjunction with `today` to get the latest image.

 ```` limits ````:
 Define the lon/lat limits. Longitude minimum and maximum, then latitude min and max.

 ```` borders ````:
 Draw US state and country borders. 

Example: To download all daytime images for 15 November 2019, use:
```` sh
python run_GOES.py --date 2019-11-15 --daytime 1
````

Example: To get the latest image:
```` sh
python run_GOES.py --date today --process_last
````

Example: To get the latest image centered on the south east continental US:
```` sh
python run_GOES.py --date today --process_last --limits -90 -70 20 40
````

 NetCDF4 files will be saved under GOES16/. The processed PNG files will be saved under pngs/. Both folder are subdivided by year and day. 

 ## Requisities
 This code requires the following Python modules to be installed

 ````
 numpy, scipy, matplotlib, glob, re, netCDF4, pyproj, cartopy, s3fs, metpy, xarray
 ````

These can be installed using the `requirements.txt` file:
````bash
python3 -m pip install -r requirements.txt
````
