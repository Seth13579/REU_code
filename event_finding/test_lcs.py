import numpy as np
import sys
from fixed_time_bin_dip_flare_detection import *
from plot_lc import plot_lc
import glob
import multiprocessing as mp

def analyze_file(file):
	#extract relevant info from the file name
	file_name = file.split('/')[-1]
	name = file_name.rstrip('_lc.fits.txt').split('_')[0]
	obsid = file_name.rstrip('_lc.fits.txt').split('_')[1]

	if '+' in name:
		ra = name.strip('J').split('+')[0]
		ra = f'{ra[0:1]}:{ra[2:3]}:{ra[4:]}'

		dec = '+' + name.split('+')[1]
		dec = f'{dec[0:1]}:{dec[2:3]}:{dec[4:]}'
	else:
		ra = name.strip('J').split('-')[0]
		dec = '-' + name.split('-')[1]

	#extract relevant info from the file
	data = np.loadtxt(file, skiprows = 1)
	counts = data[::,4]
	total = sum(counts)

	time = data[::,2]
	start_time = time[0]
	end_time = time[-1]
	length = end_time - start_time

	changes = test_slope(file,'auto',1,1)

	n_change = len(changes)

	if n_change != 0:
		changes_str = ';'.join(str(i) for i in changes)

		#make a plot
		plot_cumo(file,show=False,save=True,lines = changes)
		plot_lc(file,binsize = 'auto',show=False,save=True,lines = changes, errorbars=False)
		print(f'{n_change} changes found')

	else:
		changes_str = 'NONE'

	return [name,obsid,ra,dec,n_change,changes_str,start_time,end_time,total,length]


def analyze_dir(dir):
	light_curves = glob.glob(f'{dir}/*.txt')

	p = mp.Pool(mp.cpu_count())
	out = p.map(analyze_file,light_curves)


	'''
	#if multithreading doesnt work this should
	i = 1
	for file in light_curves:
		print(f'{i} of {len(light_curves)}')

		file_info = analyze_file(file)
		if i == 1:
			out = np.asarray(file_info,dtype='str')
		else:
			out = np.concatenate((out,file_info))

		i += 1
	'''
	return out



def main():
	dir = sys.argv[1]

	csv_out = analyze_dir(dir)
	csv_header = 'name,obsid,ra,dec,n changes,change times,start time,end time,total counts,length'

	np.save('autodetect.npy',csv_out)
	np.savetxt(f'autodetect.csv', csv_out,fmt='%s',delimiter=',',header=csv_header)

if __name__ == '__main__':
	main()
