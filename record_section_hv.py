import holoviews as hv
from holoviews.operation.datashader import datashade
from holoviews.streams import Stream as HVStream
import warnings
import bokeh
from bokeh.layouts import layout
from bokeh.plotting import figure, show, curdoc
from bokeh.models import (FuncTickFormatter, 
                      DatetimeTickFormatter,
                      Range1d,
                      HoverTool, 
                      ColumnDataSource)
from bokeh.models import DatePicker, TextInput, Button
import datashader as ds
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException
from obspy import UTCDateTime, Stream
from pyproj import Geod
import numpy as np
import pandas as pd
import os
import pickle


warnings.filterwarnings('ignore')
client = Client(base_url='http://service.geonet.org.nz')
hv.extension('bokeh')

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


def record_section_hv(start, end, target='tongariro', red_vel=0., new=False):
    targets = {'tongariro': (175.671854359, -39.107850505),
               'ruapehu': (175.564490, -39.281149)}
    st_seismic = get_data(seismic, targets[target], start, end, 'seismic', new=new)
    st_acoustic = get_data(acoustic, targets[target], start, end, 'acoustic', new=new)
    curves = {}
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
        dates = np.datetime64(tr.stats.starttime.datetime)+times.astype('timedelta64[ms]')
        idates = np.array(dates.astype(np.int) // 10**3).astype(np.float)
        curves[(tr.stats.station, 'seismic')]  = hv.Curve((idates, data-tr.stats.distance/1e3))
 
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
        dates = np.datetime64(tr.stats.starttime.datetime)+times.astype('timedelta64[ms]')
        idates = np.array(dates.astype(np.int) // 10**3).astype(np.float)
        curves[(tr.stats.station,'acoustic')] = hv.Curve((idates, data-tr.stats.distance/1e3))
        
    return hv.NdOverlay(curves, kdims=['name', 'type']) 


hv.output(size=30)


def apply_axis_formatter(plot, element):
    plot.handles['xaxis'].formatter = bokeh.models.DatetimeTickFormatter(hourmin = '%H:%M', 
                                                                         minutes='%H:%M', 
                                                                         minsec='%H:%M:%S',
                                                                         seconds='%H:%M:%S')
    plot.handles['yaxis'].formatter = FuncTickFormatter(code="""return Math.abs(tick)""")
    plot.handles['xaxis'].axis_label = "UTC time [min]"
    plot.handles['yaxis'].axis_label = "Distance [km]"
   
renderer = hv.renderer('bokeh').instance(mode='server')

def modify_doc(doc):
    
    #end = UTCDateTime.utcnow()
    #start = end - 10*60.
    start = UTCDateTime(2012,8,6,11,40)
    end = start + 20*60.
    #start = UTCDateTime(2012,11,21,0,20)
    #end = tstart + 20*60.
    #tstart = UTCDateTime(2018,1,23,6,20)
    volcano = 'tongariro'
    update_date = HVStream.define('update_date', start=start, end=end)
    dm = hv.DynamicMap(record_section_hv, streams=[update_date()])
    
    color_key = {'seismic': 'blue', 'acoustic': 'red'}    
    lyt = datashade(dm, aggregator=ds.count_cat('type'),
                    color_key=color_key, 
                    min_alpha=255, width=3000, height=2000,
                    streams=[hv.streams.RangeXY(transient=True)])      

    def date_update(attrname, old, new):
        print(new)

    def time_update(attrname, old, new):
        print(new)

    def update():
        start = UTCDateTime(2012,11,21,0,20)
        end = start + 20*60.
        print(start, end)
        dm.event(start=start, end=end)

    date_start = DatePicker(title='Start date')
    date_start.on_change('value', date_update)
    date_end = DatePicker(title='End date')
    date_end.on_change('value', date_update)
    starttime = TextInput(title='Start time', value='HH:MM:SS')
    starttime.on_change('value', time_update)
    endtime = TextInput(title='End time', value='HH:MM:SS')
    endtime.on_change('value', time_update)
    
    updateb = Button(label='Update', width=60)
    updateb.on_click(update)
    
    # Create HoloViews plot and attach the document
    lyt = lyt.opts(plot=dict(width=3000, height=2000, finalize_hooks=[apply_axis_formatter]),
                   norm=dict(framewise=True))
    hvplot = renderer.get_plot(lyt, doc)
    doc.add_root(layout([[hvplot.state], [[date_start, starttime], [date_end, endtime], [updateb]]], sizing_mode='fixed'))
    return doc

doc = modify_doc(curdoc())
