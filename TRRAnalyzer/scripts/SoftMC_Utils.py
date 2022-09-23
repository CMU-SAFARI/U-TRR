#!/usr/bin/env python3

from re import T

from pandas.core.accessor import register_dataframe_accessor

import sys

import numpy as np
import pandas as pd

import plotly.graph_objects as go

from bokeh.palettes import viridis
from bokeh.plotting import figure
from bokeh.models import Span, ColumnDataSource, Button, CustomJS
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import gridplot
from bokeh.io import show , curdoc
from bokeh.transform import dodge


js_clipboard_copy="""
        var dummy = document.createElement("textarea");
        document.body.appendChild(dummy);
        dummy.value = txt;
        dummy.select();
        document.execCommand("copy");
        document.body.removeChild(dummy);
    """

def plotBokehLinePerRow(test, x_label='', y_label='', join_wrs=False):
    TOOLS="pan,wheel_zoom,box_zoom,reset,save,box_select,hover"

    trr_data = test.data

    row_layout = test.getConfig('row_layout')
    num_weaks = row_layout.count("R") + row_layout.count("U")
    num_victims_in_wrs = int(len(trr_data.keys())/num_weaks)
        
    plots = []
    num_plots = num_victims_in_wrs if join_wrs else len(trr_data.keys())

    for i in range(0, num_plots):
        plots.append(figure(plot_width=1200, plot_height=80, tools=TOOLS))

    c_palette = viridis(len(trr_data.keys()))
    for i, row_id in enumerate(trr_data.keys()):
        df_trr = trr_data[row_id]

        plot_ind = i % len(plots)

        if len(plots) > 0:
            plots[plot_ind].x_range = plots[0].x_range

        if i == (len(trr_data.keys()) - 1):
            plots[plot_ind].xaxis.axis_label = x_label
            plots[plot_ind].xaxis.axis_label_text_font_style = 'bold'

        if y_label != '':
            plots[plot_ind].yaxis.axis_label = y_label
        else:
            plots[plot_ind].yaxis.axis_label = 'Number of bitflips'
        plots[plot_ind].yaxis.axis_label_text_font_style = 'bold'

        if 'xLabels' in df_trr:
            plots[plot_ind].line(df_trr['xLabels'], df_trr['NumBitflips'], legend_label=str(row_id), color=c_palette[i], line_width=2, alpha=0.8, name=str(row_id))
        else:
            plots[plot_ind].line(range(len(df_trr['NumBitflips'])), df_trr['NumBitflips'], legend_label=str(row_id), color=c_palette[i], line_width=2, alpha=0.8, name=str(row_id))

        hline = Span(location=0, dimension='width', line_color='black', line_width=3, line_dash='dashed')
        plots[plot_ind].renderers.extend([hline])
        
        plots[plot_ind].legend.background_fill_alpha = 0.5

            
        plots[plot_ind].hover.tooltips = [
            ("index", "$index"),
            ("row id", "$name")
        ]
            
    return plots

def plotBokehVBar(df, x_axis_col_name):
    # TOOLS="pan,wheel_zoom,box_zoom,reset,save,box_select,hover"
    TOOLS="pan,wheel_zoom,box_zoom,reset,save,box_select"

    # determine columns to plot as bars on the y-axis
    col_names = df.columns.values.tolist()

    if not x_axis_col_name in col_names:
        print(f"ERROR: Column name {x_axis_col_name} not found in {df.head()}")
        sys.exit(-1)

    col_names.remove(x_axis_col_name)

    c_palette = viridis(len(col_names))

    # create a figure
    plot = figure(x_range=df[x_axis_col_name], plot_width=1200, plot_height=130, tools=TOOLS)
    plot.yaxis.axis_label = 'Number of bitflips'
    plot.yaxis.axis_label_text_font_style = 'bold'
    plot.xaxis.axis_label = x_axis_col_name
    plot.xaxis.axis_label_text_font_style = 'bold'
    plot.x_range.range_padding = 0.1

    source = ColumnDataSource(data=df)

    bar_width = (0.2/len(col_names))*0.95
    bar_deltas = np.linspace(0, bar_width + 0.01*(len(col_names)-1), len(col_names))
    bar_deltas = [x - bar_width/2 for x in bar_deltas]

    for i, col in enumerate(col_names):
        plot.vbar(x=dodge(x_axis_col_name, bar_deltas[i], range=plot.x_range), top=col, width=bar_width, source=source, color=c_palette[i], legend_label=col)


    return plot

def plotlyVBar(df, x_axis_col_name, mode=''):
    # determine columns to plot as bars on the y-axis
    col_names = df.columns.values.tolist()

    if not x_axis_col_name in col_names:
        print(f"ERROR: Column name {x_axis_col_name} not found in {df.head()}")
        sys.exit(-1)

    col_names.remove(x_axis_col_name)

    bars = []

    if mode == 'aggregate':
        sum_bitflips = []
        
        for colname in col_names:
            sum_bitflips.append(df[colname].sum())

        bars.append(go.Bar(x=list(range(len(sum_bitflips))), y=sum_bitflips))
    else:
        for colname in col_names:
            bars.append(go.Bar(name=colname, x=df[x_axis_col_name], y=df[colname]))
    
    fig = go.Figure(data = bars)

    fig.update_layout(barmode='group', xaxis={'type': 'category'})

    return fig

def getPlotlyBar(test_data, series_type, categories_type):
    # determine columns to plot as bars on the y-axis
    df = test_data.data
    col_names = df.columns.values.tolist()

    assert(categories_type in test_data.configs.keys())

    series_label = test_data.configs[series_type]

    if len(col_names) == 1: # indicates that the data is empty
        num_victims = test_data.configs['rowlayout'].count('v')
        return go.Bar(name=series_label, x=list(range(num_victims)), y=[0]*num_victims)

    # keep only 'Victim' columns
    new_cols = []
    for col in col_names:
        if 'Victim' in col:
            new_cols.append(col)
    col_names = new_cols

    sum_bitflips = []
    for colname in col_names:
        sum_bitflips.append(df[colname].sum())

    # categories_type specifies the property to use for grouping the bars in the bar chart
    if categories_type == 'rowlayout': # rowlayout is for victim row location as x axis values in the plot
        return go.Bar(name=series_label, x=list(range(len(sum_bitflips))), y=sum_bitflips, text=sum_bitflips) 
    else: # return the sum of bitflips of all victims when not plotting the victim location bitflips separately
        return go.Bar(name=series_label, x=list(test_data.configs[categories_type]), y=[sum(sum_bitflips)], text=[sum(sum_bitflips)])


def CLsToHist(data, hist, colname):
    bitflips_per_CL = data[colname]

    for bf_count in bitflips_per_CL:
        assert(bf_count <= 512)
        hist[colname][bf_count] += 1

def getPlotlyBarCLsWithNbitflips(test_data):

    df = test_data.data
    col_names = df.columns.values.tolist()

    if len(col_names) == 1: # indicates that the data is empty
        num_victims = test_data.configs['rowlayout'].count('v')
        return [go.Bar()]

    # keep only 'perChunkBitflipsV' columns
    new_cols = []
    for col in col_names:
        if 'perChunkBitflipsV' in col:
            new_cols.append(col)
    col_names = new_cols

    hist_cls_n_bitflips = dict()
    for c in col_names:
        hist_cls_n_bitflips[c] = [0]*513 # histogram showing the number of cache lines with different bitflip counts. There can be at most 512 bitflips in a single cache line

        df.apply(lambda x: CLsToHist(x, hist_cls_n_bitflips, c), axis=1)

    test_data.hist_chunks_with_bitflips = hist_cls_n_bitflips

    # Add a bar for each non-zero histogram element
    # We have a histogram for each V in the rowlayout, i.e., len(colnames)
    bars = []
    for i in range(1, 513):
        
        for k in hist_cls_n_bitflips.keys():
            if hist_cls_n_bitflips[k][i] > 0:
                num_cls = [d[i] for d in hist_cls_n_bitflips.values()] # get i'th element from each histogram

                bars.append(go.Bar(name=i, x=list(range(len(col_names))), y=num_cls, text=num_cls))
                break

    return bars


def plotlySumBitflips(test_data_list, series_type='hammers', categories_type='rowlayout'):

    fig = go.Figure()

    for t in test_data_list:
        if series_type == 'perCLBitflips':
            fig.add_traces(getPlotlyBarCLsWithNbitflips(t))
        else:
            fig.add_trace(getPlotlyBar(t, series_type, categories_type))

    fig.update_layout(barmode='group', xaxis={'type': 'category'})

    return fig

def checkCompatibility(tests, excluded_configs):
    config_digests = []

    t0_digest = tests[0].configsDigest(excluded_configs)
    for t in tests:
        assert t.configsDigest(excluded_configs) == t0_digest, \
            f"ERROR: All tests should have the same configuration except '{excluded_configs}'" \
            f"{t.path} does not match the config of {tests[0].path}" \
            f"First config: {t.configs}" \
            f"Second config: {tests[0].configs}"
        config_digests.append(t.configsDigest())

    # there shouldn't be duplicate experiments
    set_config_digests = set(config_digests)
    if len(set_config_digests) != len(config_digests):
        print(f"ERROR: There are experiments with duplicate configurations")
        for t in tests:
            print(t.configs)
        sys.exit(-1)

def plotlyAvgPerBankBitflips(test_data_list, xaxis, xaxis_range=None, victim_row_ind=0, num_rows=0, plot_type='bar', xaxis_scale=1, xaxis_offset=0):

    fig = go.Figure()

    assert xaxis in test_data_list[0].configs.keys(), f"ERROR: Could not find {xaxis} in test configuration parameters"

    # all configs except bank ID and axis should be the same for all tests in the list
    # there shouldn't be duplicate experiments

    excluded_configs = ['bank', xaxis] # TODO: we may beed to include the module ID here to average across all tested modules
    checkCompatibility(test_data_list, excluded_configs)


    # group experiments by the same xaxis value
    t_dict = dict()
    for t in test_data_list:
        t_dict[t.configs[xaxis]] = t_dict.get(t.configs[xaxis], []) + [t]

    plot_df = pd.DataFrame(columns=['hammers', 'bitflips_per_row'])

    x_vals = []
    y_vals = []
    for t_key in t_dict.keys():

        formatted_tkey = t_key.split('-')[0]

        if xaxis_range != None and int(formatted_tkey) not in xaxis_range:
            continue

        formatted_tkey = round(int(formatted_tkey)/xaxis_scale)

        # calculate the average number of bitflips in all banks
        total_bitflips = 0

        assert(len(t_dict[t_key]) == 1)

        for test in t_dict[t_key]:
            bitflips_colname = "Victim" + str(victim_row_ind)

            if not hasattr(test, 'data'):
                test.parseTestData(file_type='BlastRadius')

            if plot_type == 'box':
                new_df = pd.DataFrame(t_dict[t_key][0].data['bitflips_per_row'])
                new_df['hammers'] = formatted_tkey + xaxis_offset
                plot_df = plot_df.append(new_df)
            else:
                if bitflips_colname in test.data.columns:
                    total_bitflips += test.data[bitflips_colname].sum()

            del test.data
        avg_bitflips = total_bitflips / len(t_dict[t_key])

        if num_rows != 0:
            avg_bitflips = avg_bitflips / num_rows

        x_vals.append(formatted_tkey + xaxis_offset)
        y_vals.append(avg_bitflips)

    if plot_type == 'box':
        fig.add_trace(go.Box(x=plot_df['hammers'], y=plot_df['bitflips_per_row']))
    else:
        fig.add_trace(go.Bar(x=x_vals, y=y_vals, text=y_vals))

    fig.update_layout(xaxis={'type': 'category'})

    return fig


def calcPlotXY(test_data, series_type, xaxis_values):
    # determine columns to plot as bars on the y-axis
    df = test_data.data
    col_names = df.columns.values.tolist()

    assert(xaxis_values in test_data.configs.keys())

    series_label = test_data.configs[series_type]

    test_data.xval = test_data.configs[xaxis_values]
    if len(col_names) == 1: # indicates that the data is empty
        num_victims = test_data.configs['rowlayout'].count('v')
        test_data.yval = 0
        return

    # keep only 'Victim' columns
    new_cols = []
    for col in col_names:
        if 'Victim' in col:
            new_cols.append(col)
    col_names = new_cols

    sum_bitflips = []
    for colname in col_names:
        sum_bitflips.append(df[colname].sum())


    if xaxis_values == 'rowlayout': # rowlayout is for victim row location as x axis values in the plot
        test_data.yval = sum_bitflips
    else: # return the sum of bitflips of all victims when not plotting the victim location bitflips separately
        test_data.yval = sum(sum_bitflips)

def plotlySumBitflipsCustomXAxis(test_data_list, series_type='hammers', xaxis_values='rowlayout'):
    fig = go.Figure()

    for t in test_data_list:
        calcPlotXY(t, series_type, xaxis_values)

    # group experiments by series_type
    t_dict = dict()
    for t in test_data_list:
        t_dict[t.configs[series_type]] = t_dict.get(t.configs[series_type], []) + [t]
        
    for s in t_dict.keys():
        fig.add_trace(go.Bar(name=s, x=[t.xval for t in t_dict[s]], y=[t.yval for t in t_dict[s]]))

    fig.update_layout(barmode='group', xaxis={'type': 'category'})

    return fig

def getPlotlyBarPerRowBitflips(test_data):
    hammers = test_data.configs['hammers']
    return go.Bar(name=hammers, x=list(range(len(test_data.bitflips_per_row))), y=test_data.bitflips_per_row)

def getPlotlyLinePerRowBitflips(test_data):
    hammers = test_data.configs['hammers']
    return go.Scatter(name=hammers, x=list(range(len(test_data.bitflips_per_row))), y=test_data.bitflips_per_row)

def plotlyBitflipsPerRow(test_data_list):
    
    fig = go.Figure()

    for t in test_data_list:
        bplot = getPlotlyLinePerRowBitflips(t)
        fig.add_trace(bplot)

    fig.update_layout(barmode='stack', xaxis={'type': 'category'})

    return fig

def plotlyRowWithBitflips(test_data_list):

    fig = go.Figure()

    xdata = []
    ydata = []
    for t in test_data_list:
        hammers = t.configs['hammers']
        xdata.append(hammers)
        ydata.append(np.count_nonzero(t.bitflips_per_row))
        
    fig.add_trace(go.Bar(x=xdata, y=ydata))

    print(ydata)
    
    fig.update_layout(barmode='stack', xaxis={'type': 'category'})
    return fig


def showInTabs(plots, tab_titles):
    # show data in multiple tabs
    fig_tabs = []
    for p, t in zip(plots, tab_titles):
        fig_tabs.append(Panel(child=gridplot(p, ncols=1), title=t))

    tabs = Tabs(tabs=fig_tabs)
    curdoc().add_root(tabs)

    show(tabs)

# calculate the distance between specific values for all columns (row ids as column names) in a SingleTest
# calculates the distance (in number of iterations) between a given val for a SingleTest
def calcDist(singleTest, val):

    import pandas as pd
    
    retdf = pd.DataFrame()

    for row_id in singleTest.data.keys():
        row_data = singleTest.data[row_id]

        col_bitflips = row_data['NumBitflips']

        last_ind = 0
        dists = []

        for i, bf in enumerate(col_bitflips.values):
            if bf == val:
                dists.append(i - last_ind)
                last_ind = i

        retdf[str(row_id)] = dists

    return retdf

# plot distances as a box plot plot for each hammer count on the x-axis
def plotBoxPlot(p, df, hc):
    from bokeh.io import push_notebook, show, output_notebook, curdoc
    from bokeh.layouts import row, gridplot, layout
    from bokeh.plotting import figure, show, output_file
    from bokeh.models import HoverTool, Legend, Span
    from bokeh.palettes import viridis
    
    output_notebook()

    c_palette = viridis(len(df.columns))
    # for each row
    for i, colname in enumerate(df):
        coldata = df[colname]
        # find the quartiles and IQR
        q1 = coldata.quantile(q=0.25)
        q2 = coldata.quantile(q=0.5)
        q3 = coldata.quantile(q=0.75)
        iqr = q3 - q1
        upper = q3 + 1.5*iqr    
        lower = q1 - 1.5*iqr    

        # find the outliers
        outliers = coldata.apply(lambda x : x if (x > upper) | (x < lower) else None).dropna()

        # prepare outlier data for plotting, we need coordinates for every outliers
        if not outliers.empty:
            outx = []
            outy = []
            for keys in outliers.index:
                outx.append(colname)
                outy.append(outliers.loc[keys])

        # stems
        p.segment(x0=[hc], y0=[upper], x1=[hc], y1=[q3], line_color=c_palette[i])
        p.segment(x0=[hc], y0=[lower], x1=[hc], y1=[q1], line_color=c_palette[i])

        # boxes
        p.vbar([hc], 0.7, [q2], [q3], fill_color=c_palette[i], fill_alpha=0.5, line_color=c_palette[i], legend_label=colname)
        p.vbar([hc], 0.7, [q1], [q2], fill_color=c_palette[i], fill_alpha=0.5, line_color=c_palette[i], legend_label=colname)

        # whiskers (almost-0 height rects simpler than segments)
        p.rect([hc], [lower], 0.2, 0.01, line_color=c_palette[i])
        p.rect([hc], [upper], 0.2, 0.01, line_color=c_palette[i])

        # outliers
        if not outliers.empty:
            p.circle(outx, outy, size=6, color=c_palette[i], fill_alpha=0.6)

def groupTests(tests, grouping_config):
    
    grouped_tests = []

    for t in tests:
        matches_existing_group = False
        for gt in grouped_tests:
            if gt[0].compareConfig(t, grouping_config):
                gt.append(t)
                matches_existing_group = True
                break

        if not matches_existing_group:
            grouped_tests.append([t])
    
    return grouped_tests


def averageRowBitflips(tests):

    assert len(tests) != 0, f"ERROR: No experiment data is provided"

    num_banks = len(tests)

    return sum([np.count_nonzero(t.get_bitflips_per_row()) for t in tests])/num_banks

def totalBitflipsPerRow(tests):

    assert len(tests) != 0, f"ERROR: No experiment data is provided"

    avgs = []
    avgs.append([sum(t.get_bitflips_per_row()) for t in tests])
    return max(avgs) # maxing across banks

def maxBitflipsPerRow(tests):

    assert len(tests) != 0, f"ERROR: No experiment data is provided"

    avgs = []
    avgs.append([max(t.get_bitflips_per_row(), default=0) for t in tests])
    return max(avgs, default=0) # maxing across banks

def averageNumChunksWithBitflips(tests, num_bitflips):
    assert len(tests) != 0, f"ERROR: No experiment data is provided"

    num_banks = len(tests)

    total_chunks = 0
    for t in tests:
        if not hasattr(t, 'hist_chunks_with_bitflips'):
            if not hasattr(t, 'data'):
                t.parseTestData(file_type='BlastRadius')
            getPlotlyBarCLsWithNbitflips(t)
            del t.data

        max_num_chunks = 0
        for v in t.hist_chunks_with_bitflips.values(): # v contains chunk numbers for each V in the used row layout, e.g., 3 Vs in VAVAV
            max_num_chunks = max(max_num_chunks, v[num_bitflips])

        total_chunks += max_num_chunks

    return total_chunks/num_banks



def aggregateData(tests, columns=None):
    
    if columns == None or len(tests) == 0:
        return None

    # check if the provided filtered_tests have common configuration
    excluded_configs = ['bank'] + columns 
    checkCompatibility(tests, excluded_configs)

    # group experiments that have different bank IDs but the rest of the configuration is the same
    grouped_tests = groupTests(tests, columns)

    configs_keys = tests[0].configs.keys()

    df = pd.DataFrame()
    for col in columns:

        if col == 'hammers': # hammers is special because we need to get rid of the dashes
            df[col] = [int(test_list[0].configs[col].split('-')[0]) for test_list in grouped_tests]
            continue

        if col == 'NumVulnerableRows':
            df[col] = [averageRowBitflips(test_list) for test_list in grouped_tests]
            continue

        if 'NumChunks' in col:
            bitflip_count = int(col.split('-')[1])
            df[col] = [averageNumChunksWithBitflips(test_list, bitflip_count) for test_list in grouped_tests]
            continue

        if 'totalBitflipsPerRow' in col:
            df[col] = [totalBitflipsPerRow(test_list) for test_list in grouped_tests]
            continue

        if 'maxBitflipsPerRow' in col:
            df[col] = [maxBitflipsPerRow(test_list) for test_list in grouped_tests]
            continue

        if col in configs_keys: # the col exists as a configuration parameter so we extract the value and put it into the dataframe
            df[col] = [test_list[0].configs[col] for test_list in grouped_tests]
            continue
        
        print(f"ERROR: Provided an undefined col name: {col}")
        sys.exit(-1)

    return df