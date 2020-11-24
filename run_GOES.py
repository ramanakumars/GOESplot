from GOES_funcs import get_data, process_files, get_filelist
import datetime, argparse

parser = argparse.ArgumentParser(description="Download and plot GOES 16 imagery")
parser.add_argument('--date', required=True, help='Date to retrieve image (format: yyyy-mm-dd or today)')
parser.add_argument('--daytime', action="store_true", help='Use daytime data (10Z-01Z)')
parser.add_argument('--process_last', action='store_true', help='Only process the latest file')

args = parser.parse_args()

if(args.date == "today"):
    year, month, day = datetime.datetime.utcnow().strftime("%Y-%m-%d").split('-')
    year  = int(year)
    month = int(month)
    day   = int(day)
else:
    year  = int(args.date.split('-')[0])
    month = int(args.date.split('-')[1])
    day   = int(args.date.split('-')[2])


date = datetime.datetime(year, month, day)

day_of_year = date.timetuple().tm_yday

# if(args.download == True):
#     get_data(year, day_of_year, int(args.daytime))

filelist = get_filelist(year, day_of_year, int(args.daytime))

if(args.process_last):
    process_files([filelist[-1]], year, day_of_year)
else:
    process_files(filelist, year, day_of_year)

