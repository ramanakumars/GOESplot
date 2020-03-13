from GOES_funcs import get_data, process_data
import datetime, argparse

parser = argparse.ArgumentParser(description="Download and plot GOES 16 imagery")
parser.add_argument('--date', required=True, help='Date to retrieve image (format: yyyy-mm-dd)')
parser.add_argument('--daytime', default=False, help='Use daytime data (10Z-01Z) (default: False)')


args = parser.parse_args()

year  = int(args.date.split('-')[0])
month = int(args.date.split('-')[1])
date  = int(args.date.split('-')[2])

date = datetime.datetime(year, month, date)

day_of_year = date.timetuple().tm_yday

# if(args.download == True):
#     get_data(year, day_of_year, int(args.daytime))

process_data(year, day_of_year, int(args.daytime))