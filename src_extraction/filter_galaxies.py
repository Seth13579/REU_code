#filters the nearby galaxies to only those w/ Chandra observations
#https://ui.adsabs.harvard.edu/search/q=bibcode%3A(2013AJ....145..101K)&sort=date%20desc%2C%20bibcode%20desc&p_=0
import os
import pickle as pkl

with open('names.txt','r') as file:
    names = file.read().splitlines()

with open('ras.txt','r') as file:
    ras = file.read().splitlines()

with open('decs.txt','r') as file:
    decs = file.read().splitlines()

obsids_to_download = []

dict = {i:None for i in names}

i = 0
for name,ra,dec in zip(names,ras,decs):
    print(f'Searching for observations of {name}, {i+1} of {len(names)}')

    os.system(f'find_chandra_obsid {ra}{dec} radius=10 grating=none > results.txt')

    with open('results.txt','r') as file:
        data = file.read().splitlines()

    obsids = []
    for line in data[1:]:
        obsids.append(line.split()[0])

    dict[name] = obsids
    i += 1



os.remove('results.txt')

with open('result_dict.pkl','wb') as f:
    pkl.dump(dict,f)
