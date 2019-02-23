#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PyCORN - script to extract data from .res (results) files generated
by UNICORN Chromatography software supplied with ÄKTA Systems
(c)2014-2016 - Yasar L. Ahmed
v0.18
'''

from __future__ import print_function
from collections import OrderedDict
from zipfile import ZipFile
from zipfile import is_zipfile
import xml.etree.ElementTree as ET
import struct
import codecs
import os
import sys
import io
import argparse
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import AutoMinorLocator
import mpl_toolkits.axisartist as AA

class pc_res3(OrderedDict):
    """A class for holding the PyCORN/RESv3 data.
    A subclass of `dict`, with the form `data_name`: `data`.
    """

    # first, some magic numbers
    RES_magic_id = b'\x11\x47\x11\x47\x18\x00\x00\x00\xB0\x02\x00\x00\x20\x6C\x03\x00'
    CNotes_id = b'\x00\x00\x01\x00\x02\x00\x03\x22'
    Methods_id = b'\x00\x00\x01\x00\x02\x00\x01\x02'
    Logbook_id = b'\x00\x00\x01\x00\x04\x00\x48\x04'
    Logbook_id2 = b'\x00\x00\x01\x00\x04\x00\x49\x04'
    SensData_id = b'\x00\x00\x01\x00\x04\x00\x01\x14'
    SensData_id2 = b'\x00\x00\x01\x00\x04\x00\x02\x14'
    Fractions_id = b'\x00\x00\x01\x00\x04\x00\x44\x04'
    Fractions_id2 = b'\x00\x00\x01\x00\x04\x00\x45\x04'
    Inject_id = b'\x00\x00\x01\x00\x04\x00\x46\x04'
    Inject_id2 = b'\x00\x00\x01\x00\x04\x00\x47\x04'
    LogBook_id = b'\x00\x00\x01\x00\x02\x00\x01\x13'  # capital B!

    def __init__(self, file_name, reduce=1, inj_sel=-1):
        OrderedDict.__init__(self)
        self.file_name = file_name
        self.reduce = reduce
        self.injection_points = None
        self.inj_sel = inj_sel
        self.inject_vol = None
        self.header_read = False
        self.run_name = ''

        with open(self.file_name, 'rb') as f:
            self.raw_data = f.read()

    def input_check(self, show=False):
        '''
        Checks if input file is a supported res file
        x = magic number, y = version string, z = file size/EOF

        Returns True or False
        '''
        if show: print((" ---- \n Input file: {0}").format(self.file_name))

        x = self.raw_data.find(self.RES_magic_id, 0, 16)
        y = self.raw_data.find(b'UNICORN 3.10', 16, 36)
        z = struct.unpack("i", self.raw_data[16:20])

        if (x, y) == (0, 24):
            if show: print(" Input is a UNICORN 3.10 file!")
            x, y = (0, 0)
        else:
            if show: print(" Input is not a UNICORN 3.10 file!")
            x, y = (1, 1)

        if z[0] == os.path.getsize(self.file_name):
            if show: print(" File size check - OK")
            z = 0
        else:
            if show: print(" File size mismatch - file corrupted?")
            z = 1
        if (x, y, z) != (0, 0, 0):
            if show: print("\n File not supported - stop!")
            return False
        else:
            if show: print("\n Alles safe! - Go go go!")
            return True

    def readheader(self):
        '''
        Extracts all the entries/declarations in the header (starts at position 686)
        '''

        # we only need to do this once
        if self.header_read: return
        self.header_read = True

        fread = self.raw_data
        header_end = fread.find(self.LogBook_id) + 342
        for i in range(686, header_end, 344):
            decl = struct.unpack("8s296s4i", fread[i:i + 320])
            full_label = codecs.decode(decl[1], 'iso8859-1').rstrip("\x00")
            if full_label.find(':') == -1:
                r_name = ''
                d_name = full_label
            else:
                r_name = full_label[:full_label.find(':')]
                d_name = full_label[full_label.find('_') + 1:]
            x = dict(magic_id=decl[0],
                     run_name=r_name,
                     data_name=d_name,
                     d_size=decl[2],
                     off_next=decl[3],
                     adresse=decl[4],
                     off_data=decl[5],
                     d_start=decl[4] + decl[5],
                     d_end=decl[4] + decl[2])
            name = x['data_name']
            dat = self.get(name, dict())
            dat.update(x)
            self[name] = dat

    def showheader(self, full=True):
        '''
        Prints content of header
        '''
        print((" ---- \n Header of {0}: \n").format(self.file_name))
        if full:
            print("  MAGIC_ID, ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
        else:
            print(" ENTRY_NAME, BLOCK_SIZE, OFFSET_TO_NEXT, ADRESSE, OFFSET_TO_DATA")
        num_blocks = len(self.items())
        for i in range(num_blocks):
            dtp = (list(self.items()))[i][1]
            if full:
                print(" ", dtp['magic_id'], dtp['data_name'], dtp['d_size'], dtp['off_next'], dtp['adresse'],
                      dtp['off_data'])
            else:
                print(" ", dtp['data_name'], dtp['d_size'], dtp['off_next'], dtp['adresse'], dtp['off_data'])

    def get_user(self):
        '''
        Show stored user name
        '''
        fread = self.raw_data[:512]
        u = struct.unpack("40s", fread[118:158])
        dec_u = codecs.decode(u[0], 'iso8859-1').rstrip("\x00")
        return dec_u

    def dataextractor(self, dat, show=False):
        '''
        Identify data type by comparing magic id, then run appropriate
        function to extract data, update orig. dict to include new data
        '''
        meta1 = [
            self.Logbook_id, self.Logbook_id2,
            self.Inject_id, self.Inject_id2,
            self.Fractions_id, self.Fractions_id2]
        meta2 = [self.CNotes_id, self.Methods_id]
        sensor = [self.SensData_id, self.SensData_id2]
        if dat['d_size'] == 0:
            pass
        elif dat['magic_id'] in meta1:
            dat.update(data=self.meta1_read(dat, show=show), data_type= 'annotation')
            return dat
        elif dat['magic_id'] in meta2:
            dat.update(data=self.meta2_read(dat, show=show), data_type= 'meta')
            return dat
        elif dat['magic_id'] in sensor:
            values, unit = self.sensor_read(dat, show=show)
            dat.update(data=values, unit=unit, data_type= 'curve')
            return dat

    def meta1_read(self, dat, show=False, do_it_for_inj_det=False):
        '''
        Extracts meta-data/type1, Logbook, fractions and Inject marks
        for a specific datum
        '''
        if show:
            print((" Reading: {0}").format(dat['data_name']))
        final_data = []
        inj_vol_to_subtract = self.inject_vol
        if do_it_for_inj_det:
            inj_vol_to_subtract = 0.0
        for i in range(dat['d_start'], dat['d_end'], 180):
            dp = struct.unpack("dd158s", self.raw_data[i:i + 174])
            # acc_time = dp[0] # not used atm
            acc_volume = round(dp[1] - inj_vol_to_subtract, 4)
            label = (codecs.decode(dp[2], 'iso8859-1')).rstrip('\x00')
            merged_data = acc_volume, label
            final_data.append(merged_data)
        return (final_data)

    def meta2_read(self, dat, show=False):
        '''
        Extracts meta-data/type2, Method/Program used in the run
        '''
        if show: print((" Reading: {0}").format(dat['data_name']))
        start, size = dat['d_start'], dat['d_size']
        tmp_data = self.raw_data[start:start + size]
        size = tmp_data.rfind(b'\n')  # declared block-size in header is always off
        # by a few bytes, hence it is redetermined here
        if show and size != len(tmp_data):
            print('meta2: reevaluated size {} -> {}'.format(size, len(tmp_data)))

        raw_data = codecs.decode(self.raw_data[start:start + size], 'iso8859-1')
        if '\r' in raw_data:
            data = raw_data
        else:
            data = raw_data.replace('\n', '\r\n')
        return data

    def sensor_read(self, dat, show=False):
        '''
        extracts sensor/run-data and applies correct division
        '''
        final_data = []
        if "UV" in dat['data_name'] or "Cond" == dat['data_name'] or "Flow" == dat['data_name']:
            sensor_div = 1000.0
        elif "Pressure" in dat['data_name']:
            sensor_div = 100.0
        else:
            sensor_div = 10.0
        if show: print((" Reading: {0}").format(dat['data_name']))

        fread = self.raw_data
        for i in range(dat['adresse'] + 207, dat['adresse'] + 222, 15):
            s_unit = struct.unpack("15s", fread[i:i + 15])
            s_unit_dec = (codecs.decode(s_unit[0], 'iso8859-1')).rstrip('\x00')
            # FIX: in some files the unit for temperature reads 'C' instead of '°C'
            if s_unit_dec == 'C':
                s_unit_dec = u'°C'
        for i in range(dat['d_start'], dat['d_end'], 8):
            sread = struct.unpack("ii", fread[i:i + 8])
            data = round((sread[0] / 100.0) - self.inject_vol, 4), sread[1] / sensor_div
            final_data.append(data)
        return (final_data[0::self.reduce], s_unit_dec)

    def inject_det(self, show=False):
        '''
        Finds injection points - required for adjusting retention volume
        '''
        inject_ids = [self.Inject_id, self.Inject_id2]
        injections = []
        if self.injection_points == None:
            self.injection_points = [0.0]
            for i in self.values():
                if i['magic_id'] in inject_ids:
                    injections = self.meta1_read(i, show=show, do_it_for_inj_det=True)
        for i in injections:
            if i[0] != 0.0:
                self.injection_points.append(i[0])
        if show:
            print(" ---- \n Injection points: \n # \t ml")
            for x, y in enumerate(self.injection_points):
                print((" {0} \t {1}").format(x, y))


    def load(self, show=False):
        '''
        extract all data and store in list
        '''
        self.readheader()
        self.run_name = self['Logbook']['run_name']
        self.inject_det()
        try:
            self.inject_vol = self.injection_points[self.inj_sel]
        except IndexError:
            print("\n WARNING - Injection point does not exist! Selected default.\n")
            self.inject_vol = self.injection_points[-1]
        for name, dat in list(self.items()):
            dat = self.dataextractor(dat, show=show)
            if dat is not None:
                self[name] = dat
            else:
                # TODO: Maybe we should keep this around?
                del self[name]

class pc_uni6(OrderedDict):
    '''
    A class for holding the pycorn/RESv6 data
    A subclass of `dict`, with the form `data_name`: `data`.
    '''
    # for manual zip-detection
    zip_magic_start = b'\x50\x4B\x03\x04\x2D\x00\x00\x00\x08'
    zip_magic_end = b'\x50\x4B\x05\x06\x00\x00\x00\x00'

    # hack to get pycorn-bin to move on
    SensData_id = 0
    SensData_id2 = 0
    Fractions_id = 0
    Fractions_id2 = 0

    def __init__(self, inp_file):
        OrderedDict.__init__(self)
        self.file_name = inp_file
        self.inject_vol = 0.0
        self.run_name = 'blank'

    def load(self, show=False):
        '''
        zip-files inside the zip-bundle are replaced by dicts, again with dicts with filename:content
        Chrom.#_#_True (=zip-files) files are unpacked from binary to floats by unpacker()
        To access x/y-value of Chrom.1_2:
        udata = pc_uni6("mybundle.zip")
        udata.load()
        x = udata['Chrom.1_2_True']['CoordinateData.Volumes']
        y = udata['Chrom.1_2_True']['CoordinateData.Amplitudes']
        '''
        with open(self.file_name, 'rb') as f:
            input_zip = ZipFile(f)
            zip_data = self.zip2dict(input_zip)
            self.update(zip_data)
            proc_yes = []
            proc_no = []
            for i in self.keys():
                tmp_raw = io.BytesIO(input_zip.read(i))
                f_header = tmp_raw.read(9)
                # tmp_raw.seek(0)
                # the following if block is to fix the non-standard zip files
                # by stripping out all the null-bytes at the end
                # see https://bugs.python.org/issue24621
                if f_header == self.zip_magic_start:
                    proper_zip = tmp_raw.getvalue()
                    f_end = proper_zip.rindex(self.zip_magic_end) + 22
                    tmp_raw = io.BytesIO(proper_zip[0:f_end])
                if is_zipfile(tmp_raw):
                    tmp_zip = ZipFile(tmp_raw)
                    x = {i:self.zip2dict(tmp_zip)}
                    self.update(x)
                    proc_yes.append(i)
                else:
                    pass
                    proc_no.append(i)
            if show:
                print("Loaded " + self.file_name + " into memory")
                print("\n-Supported-")
                for i in proc_yes:
                    print(" " + i)
                print("\n-Not supported-")
                for i in proc_no:
                    print(" " + i)
        # filter out data we dont deal with atm
        to_process = []
        for i in self.keys():
            if "Chrom" in i and not "Xml" in i:
                to_process.append(i)
        if show:
            print("\nFiles to process:")
            for i in to_process:
                print(" " + i)
        for i in to_process:
            for n in self[i].keys():
                if "DataType" in n:
                    a = self[i][n]
                    b = a.decode('utf-8')
                    x = b.strip("\r\n")
                else:
                    x = self.unpacker(self[i][n])
                tmp_dict = {n:x}
                self[i].update(tmp_dict)
        if show:
            print("Finished decoding x/y-data!")

    @staticmethod
    def zip2dict(inp):
        '''
        input = zip object
        outout = dict with filename:file-object pairs
        '''
        mydict = {}
        for i in inp.NameToInfo:
            tmp_dict = {i:inp.read(i)}
            mydict.update(tmp_dict)
        return(mydict)

    @staticmethod
    def unpacker(inp):
        '''
        input = data block
        output = list of values
        '''
        read_size = len(inp) - 48
        values = []
        for i in range(47, read_size, 4):
            x = struct.unpack("<f", inp[i:i+4])
            x = x[0]
            values.append(x)
        return(values)

    def xml_parse(self,show=False):
        '''
        parses parts of the Chrom.1.Xml and creates a res3-like dict
        '''
        tree = ET.fromstring(self['Chrom.1.Xml'])
        mc = tree.find('Curves')
        me = tree.find('EventCurves')
        print(tree.tag)
        print(tree.attrib)
        event_dict = {}
        for i in range(len(me)):
            magic_id = self.SensData_id
            e_type = me[i].attrib['EventCurveType']
            e_name = me[i].find('Name').text
            if e_name == 'Fraction':
                e_name = 'Fractions' # another hack for pycorn-bin
            e_orig = me[i].find('IsOriginalData').text
            e_list = me[i].find('Events')
            e_data = []
            for e in range(len(e_list)):
                e_vol = float(e_list[e].find('EventVolume').text)
                e_txt = e_list[e].find('EventText').text
                e_data.append((e_vol,e_txt))
            if e_orig == "false":
                print("not added - not orig data")
            if e_orig == "true":
                print("added - orig data")
                x = {'run_name':"Blank", 'data': e_data, 'data_name':e_name, 'magic_id':magic_id}
                event_dict.update({e_name:x})
        self.update(event_dict)
        chrom_dict = {}
        for i in range(len(mc)):
            d_type = mc[i].attrib['CurveDataType']
            d_name = mc[i].find('Name').text
            d_fname = mc[i].find('CurvePoints')[0][1].text
            d_unit = mc[i].find('AmplitudeUnit').text
            magic_id = self.SensData_id
            try:
                x_dat = self[d_fname]['CoordinateData.Volumes']
                y_dat = self[d_fname]['CoordinateData.Amplitudes']
                zdata = list(zip(x_dat,y_dat))
                if d_name == "UV cell path length":
                    d_name = "xUV cell path length" # hack to prevent pycorn-bin from picking this up
                x = {'run_name':"Blank", 'data': zdata, 'unit': d_unit, 'data_name':d_name, 'data_type':d_type, 'magic_id':magic_id}
                chrom_dict.update({d_name:x})
            except:
                KeyError
                # don't deal with data that does not make sense atm
                # orig2.zip contains UV-blocks that are (edited) copies of
                # original UV-trace but they dont have the volume data
            if show:
                print("---")
                print(d_type)
                print(d_name)
                print(d_fname)
                print(d_unit)
        self.update(chrom_dict)

    def clean_up(self):
        '''
        deletes everything and just keeps relevant run-date
        resulting dict is more like res3
        '''
        manifest = ET.fromstring(self['Manifest.xml'])
        for i in range(len(manifest)):
            file_name = manifest[i][0].text
            self.pop(file_name)
        self.pop('Manifest.xml')

try:
    # from mpl_toolkits.axes_grid1 import host_subplot
    # from matplotlib.ticker import AutoMinorLocator
    # import mpl_toolkits.axisartist as AA
    # import matplotlib.pyplot as plt
    plotting = True
except:
    ImportError
    print("WARNING: Matplotlib not found - Plotting disabled!")
    plotting = False

try:
    # import xlsxwriter
    xlsx = False
except:
    ImportError
    print("WARNING: xlsxwriter not found - xlsx-output disabled!")
    xlsx = False

pcscript_version = 0.14

parser = argparse.ArgumentParser(
    description = "Extract data from UNICORN .res files to .csv/.txt and plot them (matplotlib required)",
    epilog = "Make it so!")
parser.add_argument("-c", "--check",
                    help = "Perform simple check if file is supported",
                    action = "store_true")
parser.add_argument("-n", "--info",
                    help = "Display entries in header",
                    action = "store_true")
parser.add_argument("-i", "--inject", type = int, default = None,
                    help = "Set injection number # as zero retention, use -t to find injection points",
                    metavar="#")
parser.add_argument("-r", "--reduce", type = int, default = 1,
                    help = "Write/Plot only every n sample",
                    metavar="#")
parser.add_argument("-t", "--points",
                    help = "Display injection points",
                    action = "store_true")

group0 = parser.add_argument_group('Extracting', 'Options for writing csv/txt files')
group0.add_argument("-e", "--extract", type=str, choices=['csv','xlsx'],
                    help = "Write data to csv or xlsx file for supported data blocks")

group1 = parser.add_argument_group('Plotting', 'Options for plotting')
group1.add_argument("-p", "--plot",
                    help = 'Plot curves',
                    action = "store_true")
group1.add_argument("--no_fractions",
                    help="Disable plotting of fractions",
                    action = "store_true")
group1.add_argument("--no_inject",
                    help="Disable plotting of inject marker(s)",
                    action = "store_true")
group1.add_argument("--no_legend",
                    help="Disable legend for plot",
                    action = "store_true")
group1.add_argument("--no_title",
                    help="Disable title for plot",
                    action = "store_true")
group1.add_argument("--xmin", type = float, default=None,
                    help="Lower bound on the x-axis",
                    metavar="#")
group1.add_argument("--xmax", type = float, default=None,
                    help="Upper bound on the x-axis",
                    metavar="#")
group1.add_argument("--par1", type = str, default='Cond',
                    help="Data for 2nd y-axis (Default=Cond), to disable 2nd y-axis, use --par1 None")
group1.add_argument("--par2", type = str, default=None,
                    help="Data for 3rd y-axis (Default=None)")
group1.add_argument('-f', '--format', type = str,
                    choices=['svg','svgz','tif','tiff','jpg','jpeg',
                    'png','ps','eps','raw','rgba','pdf','pgf'],
                    default = 'png',
                    help = "File format of plot files (default: pdf)")
group1.add_argument('-d', '--dpi', default=300, type=int,
					help="DPI (dots per inch) for raster images (png, jpg, etc.). Default is 300.")
parser.add_argument("-u", "--user",
                    help = "Show stored user name",
                    action = "store_true")
parser.add_argument('--version', action='version', version=str(pcscript_version))
# parser.add_argument("inp_res",
#                     help="Input .res file(s)",
#                     nargs='+',
#                     metavar="<file>.res")
#args.no_inject
args = parser.parse_args()

def mapper(min_val, max_val, perc):
    '''
    calculate relative position in delta min/max
    '''
    x = abs(max_val - min_val) * perc
    if min_val < 0:
        return (x - abs(min_val))
    else:
        return (x + min_val)


def expander(min_val, max_val, perc):
    '''
    expand -/+ direction of two values by a percentage of their delta
    '''
    delta = abs(max_val - min_val)
    x = delta * perc
    return (min_val - x, max_val + x)


def xy_data(inp):
    '''
    Takes a data block and returns two lists with x- and y-data
    '''
    x_data = [x[0] for x in inp]
    y_data = [x[1] for x in inp]
    return x_data, y_data


def uvdata(inp):
    '''
    helps in finding the useful data
    '''
    UV_blocks = [i for i in inp if i.startswith('UV') or i.endswith('nm')]
    for i in UV_blocks:
        if i.endswith("_0nm"):
            UV_blocks.remove(i)


def smartscale(inp):
    '''
    input is the entire fdata block
    checks user input/fractions to determine scaling of x/y-axis
    returns min/max for x/y
    '''
    UV_blocks = [i for i in inp.keys() if i.startswith('UV') and not i.endswith('_0nm')]
    uv1_data = inp[UV_blocks[0]]['data']
    uv1_x, uv1_y = xy_data(uv1_data)
    try:
        uv2_data = inp[UV_blocks[1]]['data']
        uv2_x, uv2_y = xy_data(uv2_data)
        uv3_data = inp[UV_blocks[2]]['data']
        uv3_x, uv3_y = xy_data(uv3_data)
    except:
        KeyError
        uv2_data = None
        uv3_data = None
    try:
        frac_data = inp['Fractions']['data']
        frac_x, frac_y = xy_data(frac_data)
        frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
        frac_delta.append(frac_delta[-1])
    except:
        KeyError
        frac_data = None
    if args.xmin != None:
        plot_x_min = args.xmin
    else:
        if frac_data:
            plot_x_min = frac_data[0][0]
        else:
            plot_x_min = uv1_x[0]
    if args.xmax:
        plot_x_max = args.xmax
    else:
        if frac_data:
            plot_x_max = frac_data[-1][0] + frac_delta[-1]*2 # recheck
        else:
            plot_x_max = uv1_x[-1]
    if plot_x_min > plot_x_max:
        print("Warning: xmin bigger than xmax - adjusting...")
        plot_x_min = uv1_x[0]
    if plot_x_max < plot_x_min:
        print("Warning: xmax smaller than xmin - adjusting...")
        plot_x_max = uv1_x[-1]
    # optimize y_scaling
    min_y_values = []
    max_y_values = []
    for i in UV_blocks:
        tmp_x, tmp_y = xy_data(inp[i]['data'])
        range_min_lst = [abs(a - plot_x_min) for a in tmp_x]
        range_min_idx = range_min_lst.index(min(range_min_lst))
        range_max_lst = [abs(a - plot_x_max) for a in tmp_x]
        range_max_idx = range_max_lst.index(min(range_max_lst))
        values_in_range = tmp_y[range_min_idx:range_max_idx]
        min_y_values.append(min(values_in_range))
        max_y_values.append(max(values_in_range))
    plot_y_min_tmp = min(min_y_values)
    plot_y_max_tmp = max(max_y_values)
    plot_y_min, plot_y_max = expander(plot_y_min_tmp, plot_y_max_tmp, 0.085)
    return plot_x_min, plot_x_max, plot_y_min, plot_y_max

def plotterX(inp,fname):
    plot_x_min, plot_x_max, plot_y_min, plot_y_max = smartscale(inp)
    host = host_subplot(111, axes_class=AA.Axes)
    host.set_xlabel("Elution volume (ml)")
    host.set_ylabel("Absorbance (mAu)")
    host.set_xlim(plot_x_min, plot_x_max)
    host.set_ylim(plot_y_min, plot_y_max)
    for i in inp.keys():
        if i.startswith('UV') and not i.endswith('_0nm'):
            x_dat, y_dat = xy_data(inp[i]['data'])
            print("Plotting: " + inp[i]['data_name'])
            stl = styles[i[:4]]
            p0, = host.plot(x_dat, y_dat, label=inp[i]['data_name'], color=stl['color'],
                            ls=stl['ls'], lw=stl['lw'],alpha=stl['alpha'])
    if args.par1 == 'None':
        args.par1 = None
    if args.par1:
        try:
            par1_inp = args.par1
            par1 = host.twinx()
            par1_data = inp[par1_inp]
            stl = styles[par1_inp[:4]]
            par1.set_ylabel(par1_data['data_name'] + " (" + par1_data['unit'] + ")", color=stl['color'])
            x_dat_p1, y_dat_p1 = xy_data(par1_data['data'])
            p1_ymin, p1_ymax = expander(min(y_dat_p1), max(y_dat_p1), 0.085)
            par1.set_ylim(p1_ymin, p1_ymax)
            print("Plotting: " + par1_data['data_name'])
            p1, = par1.plot(x_dat_p1, y_dat_p1, label=par1_data['data_name'],
            color=stl['color'], ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'])
        except:
            KeyError
            if par1_inp != None:
                print("Warning: Data block chosen for par1 does not exist!")
    if args.par2:
        try:
            par2_inp = args.par2
            par2 = host.twinx()
            offset = 60
            new_fixed_axis = par2.get_grid_helper().new_fixed_axis
            par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
            par2.axis["right"].toggle(all=True)
            par2_data = inp[par2_inp]
            stl = styles[par2_inp[:4]]
            par2.set_ylabel(par2_data['data_name'] + " (" + par2_data['unit'] + ")", color=stl['color'])
            x_dat_p2, y_dat_p2 = xy_data(par2_data['data'])
            p2_ymin, p2_ymax = expander(min(y_dat_p2), max(y_dat_p2), 0.075)
            par2.set_ylim(p2_ymin, p2_ymax)
            print("Plotting: " + par2_data['data_name'])
            p2, = par2.plot(x_dat_p2, y_dat_p2, label=par2_data['data_name'],
            color=stl['color'],ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'])
        except:
            KeyError
            if par2_inp != None:
                print("Warning: Data block chosen for par2 does not exist!")
    if not args.no_fractions:
        try:
            frac_data = inp['Fractions']['data']
            frac_x, frac_y = xy_data(frac_data)
            frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
            frac_delta.append(frac_delta[-1])
            frac_y_pos = mapper(host.get_ylim()[0], host.get_ylim()[1], 0.015)
            for i in frac_data:
                host.axvline(x=i[0], ymin=0.065, ymax=0.0, color='r', linewidth=0.85)
                host.annotate(str(i[1]), xy=(i[0] + frac_delta[frac_data.index(i)] * 0.55, frac_y_pos),
                         horizontalalignment='center', verticalalignment='bottom', size=8, rotation=90)
        except:
            KeyError
    if inp.inject_vol != 0.0:
        injections = inp.injection_points
        host.axvline(x=0, ymin=0.10, ymax=0.0, color='#FF3292',
                     ls ='-', marker='v', markevery=2, linewidth=1.5, alpha=0.85, label='Inject')
    host.set_xlim(plot_x_min, plot_x_max)
    if not args.no_legend:
        host.legend(fontsize=8, fancybox=True, labelspacing=0.4, loc='upper right', numpoints=1)
    host.xaxis.set_minor_locator(AutoMinorLocator())
    host.yaxis.set_minor_locator(AutoMinorLocator())
    if not args.no_title:
        plt.title(fname, loc='left', size=9)
    plot_file = fname[:-4] + "_" + inp.run_name + "_plot." + args.format
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.dpi)
    print("Plot saved to: " + plot_file)
    plt.clf()

def data_writer1(fname, inp):
    '''
    writes sensor/run-data to csv-files
    '''
    for i in inp.keys():
        print("Writing: " + inp[i]['data_name'])
        outfile_base = fname[:-4] + "_" + inp.run_name + "_" + inp[i]['data_name']
        type = inp[i]['data_type']
        if type == 'meta':
            data = inp[i]['data']
            data_to_write = data.encode('utf-8')
            ext = '.txt'
            sep = '\t'
            with open(outfile_base + ext, 'wb') as fout:
                fout.write(data_to_write)
        else:
            x_dat,y_dat = xy_data(inp[i]['data'])
            ext = '.csv'
            sep = ','
            with open(outfile_base + ext, 'wb') as fout:
                for x,y in zip(x_dat,y_dat):
                    dp = str(x) + sep + str(y) + str('\r\n')
                    data_to_write = dp.encode('utf-8')
                    fout.write(data_to_write)

def generate_xls(inp, fname):
    '''
    Input = pycorn object
    output = xlsx file
    '''
    xls_filename = fname[:-4] + "_" + inp.run_name + ".xlsx"
    workbook = xlsxwriter.Workbook(xls_filename)
    worksheet = workbook.add_worksheet()
    writable_blocks = [inp.Fractions_id, inp.Fractions_id2, inp.SensData_id, inp.SensData_id2]
    d_list = []
    for i in inp.keys():
        if inp[i]['magic_id'] in writable_blocks:
            d_list.append(i)
    for i in d_list:
        dat = inp[i]['data']
        try:
            unit = inp[i]['unit']
        except:
            KeyError
            unit = 'Fraction'
        header1 = (inp[i]['data_name'], '')
        header2 = ('ml', unit)
        dat.insert(0, header1)
        dat.insert(1, header2)
        row = 0
        col = d_list.index(i) *2
        print("Writing: " + i)
        for x_val, y_val in (dat):
            worksheet.write(row, col, x_val)
            worksheet.write(row, col + 1, y_val)
            row += 1
    workbook.close()
    print("Data written to: " + xls_filename)


styles = {'UV':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV1_':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV2_':{'color': '#e51616', 'lw': 1.4, 'ls': "-", 'alpha':1.0},
'UV3_':{'color': '#c73de6', 'lw': 1.2, 'ls': "-", 'alpha':1.0},
'UV 1':{'color': '#1919FF', 'lw': 1.6, 'ls': "-", 'alpha':1.0},
'UV 2':{'color': '#e51616', 'lw': 1.4, 'ls': "-", 'alpha':1.0},
'UV 3':{'color': '#c73de6', 'lw': 1.2, 'ls': "-", 'alpha':1.0},
'Cond':{'color': '#FF7C29', 'lw': 1.4, 'ls': "-", 'alpha':0.75},
'Conc':{'color': '#0F990F', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Pres':{'color': '#C0CBBA', 'lw': 1.0, 'ls': "-", 'alpha':0.50},
'Temp':{'color': '#b29375', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'Inje':{'color': '#d56d9d', 'lw': 1.0, 'ls': "-", 'alpha':0.75},
'pH':{'color': '#0C7F7F', 'lw': 1.0, 'ls': "-", 'alpha':0.75},}


def main2():

    # path = os.path.abspath(__file__)
    # path, file = os.path.split(path)

    # determine if application is a script file or frozen exe
    if getattr(sys, 'frozen', False):
        application_path = os.path.dirname(sys.executable)
    elif __file__:
        application_path = os.path.dirname(__file__)

    for f in os.listdir(application_path):
        if f.lower().endswith('.res'):
            fname = os.path.join(application_path, f)
    # for fname in args.inp_res:
            if args.inject == None:
                args.inject = -1
            if (fname[-3:]).lower() == "zip":
                fdata = pc_uni6(fname)
                fdata.load()
                fdata.xml_parse()
                fdata.clean_up()
            if (fname[-3:]).lower() == "res":
                fdata = pc_res3(fname, reduce = args.reduce, inj_sel=args.inject)
                fdata.load()
            if args.extract == 'csv':
                data_writer1(fname, fdata)
            if args.extract == 'xlsx' and xlsx == True:
                generate_xls(fdata, fname)
            if args.check:
                fdata.input_check(show=True)
            if args.info:
                fdata.showheader()
            if args.points:
                fdata.inject_det(show=True)
            if args.user:
                user = fdata.get_user()
                print("User: " + user)
            if plotting:
                plotterX(fdata, fname)

main2()
