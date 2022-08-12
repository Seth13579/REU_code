import glob
import os
import sys

co = '{'
cc = '}'
sl = '\\'

class Source_for_overleaf:
    def __init__(self,name,obsid_arr,filename_arr):
        self.name = name

        self.obsid_arr = obsid_arr

        self.filename_arr = filename_arr


    def make_text(self):

        obsid_string_arr = [f'{sl}item -- {i}' for i in self.obsid_arr]
        obsid_string = '\n'.join(obsid_string_arr)

        figure_string_arr = [f'{sl}includegraphics[width={sl}textwidth]{co}Images/{f}{cc}' for f in self.filename_arr]
        figure_string = '\n'.join(figure_string_arr)

        string = f'''
{sl}subsection{co}{self.name}{cc}
{sl}begin{{enumerate}}
    {obsid_string}
{sl}end{{enumerate}}
{sl}begin{{figure}}[H]
    {sl}centering
    {figure_string}
{sl}end{{figure}}
'''
        return string

def parse(filename):
    name_arr = filename.split('_')
    obsid = name_arr[1]
    name = name_arr[0].strip('J')

    return obsid,name

def make_source(filename,source_arr):
    obsid,name = parse(filename)

    name_arr = [i.name for i in source_arr]
    if name not in name_arr:
        source = Source_for_overleaf(name,[obsid],[filename])
        source_arr.append(source)
    else:
        for i in source_arr:
            if i.name == name:
                source = i
                break
        source.obsid_arr.append(obsid)
        source.filename_arr.append(filename)

    return


if __name__ == '__main__':
    galaxy = sys.argv[1]

    with open('out.txt','w') as out:


        out.write(f'\section{co}{galaxy}{cc}\n')

        source_arr = []
        for file in glob.iglob(f'{galaxy}/detections_lcs/*.png'):
            filename = file.split('/')[-1]

            make_source(filename,source_arr)

        for source in source_arr:
            str = source.make_text()

            out.write(str)
