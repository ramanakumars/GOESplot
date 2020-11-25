# GOESplot
Python code to plot GOES imagery

## Usage
Run using:
````sh
$ python run_GOES.py --date yyyy-mm-dd [--daytime 0]
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
 numpy, scipy, matplotlib, glob, re, s3fs, netCDF4, pyproj, cartopy
 ````

 For Windows, you can find most of these packages [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/)
 ```` s3fs  ```` can be installed using ```pip```:
 ````
 python -m pip install s3fs
 ````

 For Unix, these packages can be installed using ````pip````.

