PyCORN
======

A script to extract data from UNICORN .res files and plot them.

![With fractions - yay!](https://github.com/pyahmed/PyCORN/blob/dev/samples/sample1_Plot_2009Jun16no001_UV.jpg)

Description: 

A script to extract data from .res (results) files generated by UNICORN Chromatography software supplied with ÄKTA Systems. This script will find all data blocks, extract the data and write out csv-files. If you have matplotlib installed it will also plot all the curves including fractions if present. Plots can be saved in any format supported by matplotlib (default is pdf).

Limitations:
- The Logbook is not written out at the moment because non-ascii chars are not handled.

Requirements:
- Python 2.7 or 3.x (Tested on 2.7/3.4 on Windows)
- optional matplotlib (for plotting)

Usage:
- See USAGE.txt

License:
- GPLv2 see LICENSE.txt
