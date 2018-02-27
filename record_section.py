from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException
from obspy import UTCDateTime, Stream
from pyproj import Geod
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import pickle
from bokeh.plotting import figure, show, curdoc
from bokeh.layouts import column, gridplot
from bokeh.models import (FuncTickFormatter, 
                      DatetimeTickFormatter,
                      Range1d,
                      HoverTool, 
                      ColumnDataSource)
from bokeh.models.widgets.inputs import DatePicker
client = Client(base_url='http://service.geonet.org.nz')

seismic = [('KRVZ', '10', 'EHZ'), ('OTVZ', '10', 'HHZ'),
           ('WTVZ', '10', 'EHZ'), ('NGZ', '10', 'EHZ'),
           ('COVZ', '10', 'HHZ'), ('TUVZ', '10', 'EHZ'),
           ('FWVZ', '10', 'HHZ'), ('MAVZ', '10', 'HHZ'),
           ('WHVZ', '10', 'HHZ'), ('TRVZ', '10', 'HHZ'),
           ('WNVZ', '11', 'EHZ'), ('MOVZ', '10', 'EHZ'),
           ('MTVZ', '10', 'EHZ'), ('NTVZ', '10', 'HHZ')]

acoustic = [('COVZ', '30', 'HDF'), ('IVVZ', '30', 'HDF'),
            ('FWVZ', '30', 'HDF'), ('WHVZ', '30', 'HDF'),
            ('TRVZ', '30', 'HDF'), ('TOVZ', '30', 'HDF'),
            ('MAVZ', '30', 'HDF'), ('TMVZ', '30', 'HDF'),
            ('KRVZ', '30', 'HDF'), ('OTVZ', '30', 'HDF'),
            ('WTVZ', '30', 'HDF')]


def get_data(stations, target, tstart, tend, prefix, new=False):
    fout = os.path.join('/tmp','_'.join((prefix,str(target[0]),str(target[1]),str(tstart),str(tend))))
    g = Geod(ellps='WGS84')
    st = Stream()
    st_max = -1e30
    if os.path.isfile(fout) and new is False:
        with open(fout, 'rb') as fh:
            st = pickle.load(fh)
    else:     
        for s, loc, cmp in stations:
            try:
                print(s, str(loc), cmp, tstart, tend)
                st_tmp = client.get_waveforms('NZ', s, str(loc), cmp, 
                                              tstart, tend, 
                                              attach_response=True)
            except FDSNException:
                print('No data for {}.{}'.format(s,cmp))
                continue
            st_tmp.remove_sensitivity()
            st_tmp.merge(method=1, fill_value=0.)
            tr = st_tmp[0]
            tr.trim(tstart, tend)
            if False:
                if int(round(tr.stats.sampling_rate,0)) == 100:
                    tr.data -= tr.data.mean()
                    tr.taper(0.05)
                    tr.decimate(10)
                    tr.decimate(10)
            inv = client.get_stations(network='NZ', station=s, 
                                      starttime=tstart, endtime=tend)
            _s = inv[0][0]
            _,_,d = g.inv(target[0], target[1], _s.longitude, _s.latitude)
            tr.stats.distance = d
            st += tr
            st_max = max(st_max, tr.data.max())
            with open(fout, 'wb') as fh:
                pickle.dump(st, fh)
    return st


def record_section_bokeh(target, tstart, tend, tscan=5., red_vel=0., new=False):
    targets = {'tongariro': (175.671854359, -39.107850505),
               'ruapehu': (175.564490, -39.281149)}

    st_seismic = get_data(seismic, targets[target], tstart, tend, 'seismic', new=new)
    st_acoustic = get_data(acoustic, targets[target], tstart, tend, 'acoustic', new=new)
    bk_data = {}
    tools = 'pan,box_zoom,reset,save'
    hover = HoverTool(tooltips=[("station", "@station"),])
    p = figure(width=800, height=400, x_axis_type="datetime", tools=tools)
    p.add_tools(hover)
    for tr in st_seismic:
        data = tr.data[:].astype(np.float)
        data -= data.mean()
        data /= data.max()
        try:
            id0 = int(tr.stats.distance/(red_vel*tr.stats.delta))
        except ZeroDivisionError:
            id0 = 0
        data = data[id0:]
        times = np.arange(data.size)*(tr.stats.delta*1e3)
        dates = np.datetime64(tstart.datetime)+times.astype('timedelta64[ms]')
        df = {'x':dates, 'y': data-tr.stats.distance/1e3, 'station':[tr.stats.station]*dates.size}
        cds = ColumnDataSource(df)
        p.line(x='x', y='y', source=cds, color='blue')

    for tr in st_acoustic:
        data = tr.data[:].astype(np.float)
        data -= data.mean()
        data /= data.max()
        try:
            id0 = int(tr.stats.distance/(red_vel*tr.stats.delta))
        except ZeroDivisionError:
            id0 = 0
        data = data[id0:]
        times = np.arange(data.size)*(tr.stats.delta*1e3)
        dates = np.datetime64(tstart.datetime)+times.astype('timedelta64[ms]')
        cds = ColumnDataSource({'x':dates, 'y': data-tr.stats.distance/1e3, 'station':[tr.stats.station]*dates.size})
        p.line(x='x', y='y', source=cds, color='red')
    p.yaxis.formatter = FuncTickFormatter(code="""return Math.abs(tick)""")
    p.xaxis.formatter = DatetimeTickFormatter(hourmin = '%H:%M', minutes='%H:%M', minsec='%H:%M:%S')
    p.xaxis.axis_label = "UTC time [min]"
    p.yaxis.axis_label = "Distance [km]"
    tend_plot = tstart + tscan*60.
    p.x_range = Range1d(tstart.datetime, tend_plot.datetime)
    return p

tstart = UTCDateTime(2012,8,6,11,40)
tend = tstart + 20*60.
volcano = 'tongariro'
p = record_section_bokeh(volcano, tstart, tend, red_vel=0, new=True)

curdoc().add_root(column(p))
