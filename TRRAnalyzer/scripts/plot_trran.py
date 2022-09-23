#!/usr/bin/env python3

import sys

from SoftMC_Experiment import SingleTest
from SoftMC_Utils import plotBokehLinePerRow
from bokeh.io import output_file, save
from bokeh.models.widgets import Panel, Tabs
from bokeh.layouts import gridplot

output_file("trran_plot.html")

if len(sys.argv) != 2:
    print("ERROR: You must provide a path to a TRR Analyzer output file.")
    print("Usage: $ plot_trran.py <PATH_TO_TRR_ANALYZER_OUTPUT_FILE>")
    sys.exit(-1)

trran_file_path = sys.argv[1]
output_data = SingleTest(trran_file_path)

output_data.parseTestData()
output_data.convertToTRR('NumBitflips')

plot = plotBokehLinePerRow(output_data, y_label="#TRRs", join_wrs=False)

panel = Panel(child=gridplot(plot, ncols=1))
tab = Tabs(tabs=[panel])

save(tab)



