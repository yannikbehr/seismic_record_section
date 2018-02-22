from functools import partial
from concurrent.futures import ThreadPoolExecutor
import holoviews as hv
from holoviews.operation.datashader import datashade
from holoviews.streams import Stream as HVStream
import warnings
import bokeh
from bokeh.layouts import layout, widgetbox
from bokeh.plotting import figure, show, curdoc
from bokeh.models import (FuncTickFormatter, 
                      DatetimeTickFormatter,
                      Range1d,
                      HoverTool, 
                      ColumnDataSource)
from bokeh.models import (DatePicker, 
                          TextInput, 
                          Button,
                          Select,
                          PreText)
from tornado import gen
from bokeh.driving import count
from bokeh.document import without_document_lock
import datashader as ds
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException, FDSNNoDataException
from obspy import UTCDateTime, Stream
from obspy.core.util import AttribDict
from pyproj import Geod
import numpy as np
import pandas as pd
import os
import pickle
import time


warnings.filterwarnings('ignore')

tnp_seismic = [('KRVZ', '10', 'EHZ'), ('OTVZ', '10', 'HHZ'),
               ('WTVZ', '10', 'EHZ'), ('NGZ', '10', 'EHZ'),
               ('COVZ', '10', 'HHZ'), ('TUVZ', '10', 'EHZ'),
               ('FWVZ', '10', 'HHZ'), ('MAVZ', '10', 'HHZ'),
               ('WHVZ', '10', 'HHZ'), ('TRVZ', '10', 'HHZ'),
               ('WNVZ', '11', 'EHZ'), ('MOVZ', '10', 'EHZ'),
               ('MTVZ', '10', 'EHZ'), ('NTVZ', '10', 'HHZ')]

tnp_acoustic = [('COVZ', '30', 'HDF'), ('IVVZ', '30', 'HDF'),
                ('FWVZ', '30', 'HDF'), ('WHVZ', '30', 'HDF'),
                ('TRVZ', '30', 'HDF'), ('TOVZ', '30', 'HDF'),
                ('MAVZ', '30', 'HDF'), ('TMVZ', '30', 'HDF'),
                ('KRVZ', '30', 'HDF'), ('OTVZ', '30', 'HDF'),
                ('WTVZ', '30', 'HDF')]
stations = {}
stations['tongariro'] = {'seismic':tnp_seismic,
                         'acoustic':tnp_acoustic}
stations['ruapehu'] = {'seismic':tnp_seismic,
                       'acoustic':tnp_acoustic}

def GeoNetFDSNrequest(date1, date2, net, sta, loc, cmp):
    """
    Request waveform and meta data from GeoNet's FDSN webservices.
    """
    #GeoNet's FDSN web servers
    arc_client = Client('http://service.geonet.org.nz')
    nrt_client = Client('http://beta-service-nrt.geonet.org.nz')
    time1 = UTCDateTime(date1)
    time2 = UTCDateTime(date2)
    try:
        st = nrt_client.get_waveforms(net, sta, loc, cmp, time1, time2,
                                           attach_response=True)
        inv = nrt_client.get_stations(network=net, station=sta, 
                                           starttime=time1, endtime=time2)
    except FDSNNoDataException:
        st = arc_client.get_waveforms(net, sta, loc, cmp, time1, time2,
                                           attach_response=True)
        inv = arc_client.get_stations(network=net, station=sta, 
                                           starttime=time1, endtime=time2)
    return (st, inv)

def get_data(station, location, component,
             tstart, tend, new=False,
             decimate=False):
    """
    Download the data and insert station location coordinates.
    """
    scl = '.'.join((station,location,component))
    fout = os.path.join('/tmp','_'.join((scl, str(tstart), str(tend))))
    if os.path.isfile(fout) and new is False:
        with open(fout, 'rb') as fh:
            tr = pickle.load(fh)
            return tr
    else:     
        try:
            st, inv = GeoNetFDSNrequest(tstart, tend, 'NZ', station, 
                                                 location, component)
        except FDSNException:
            return None 
        st.remove_sensitivity()
        st.merge(method=1, fill_value=0.)
        tr = st[0]
        tr.trim(tstart, tend)
        tr.data -= tr.data.mean()
        if decimate:
            if int(round(tr.stats.sampling_rate,0)) == 100:
                tr.data -= tr.data.mean()
                tr.taper(0.05)
                tr.decimate(10)
                tr.decimate(10)
        _s = inv[0][0]
        tr.stats.coordinates = AttribDict({'latitude':_s.latitude,
                                           'longitude':_s.longitude})
        with open(fout, 'wb') as fh:
            pickle.dump(tr, fh)
        return tr 

doc = curdoc()

text = ['Program startup\n']
pre = PreText(text='Program startup', width=500, height=100)

class RecordSection:

    def __init__(self):
        self.streams = {'seismic': Stream(), 'acoustic': Stream()}
        self.targets = {'tongariro': (175.671854359, -39.107850505),
                        'ruapehu': (175.564490, -39.281149)}
        self.end = UTCDateTime(2018,2,11,1,15,00)
        self.start = self.end - 10*60.
        self.target = 'tongariro'
        self.curves = {}

    def reset(self):
        """
        Reset state before downloading new data.
        """
        self.streams = {'seismic': Stream(), 'acoustic': Stream()}
        self.curves = {}

    def build_record_section(self, red_vel=0., new=False):
        g = Geod(ellps='WGS84')
        tlon, tlat = self.targets[self.target]
        for tr in self.streams['seismic']:
            _lat = tr.stats.coordinates.latitude
            _lon = tr.stats.coordinates.longitude
            _,_,d = g.inv(tlon, tlat, _lon, _lat)
            data = tr.data[:].astype(np.float)
            data /= data.max()
            try:
                id0 = int(d/(red_vel*tr.stats.delta))
            except ZeroDivisionError:
                id0 = 0
            data = data[id0:]
            times = np.arange(data.size)*(tr.stats.delta*1e3)
            dates = np.datetime64(tr.stats.starttime.datetime)+times.astype('timedelta64[ms]')
            idates = np.array(dates.astype(np.int) // 10**3).astype(np.float)
            self.curves[(tr.stats.station, 'seismic')]  = hv.Curve((idates, data-d/1e3))
     
        for tr in self.streams['acoustic']:
            _lat = tr.stats.coordinates.latitude
            _lon = tr.stats.coordinates.longitude
            _,_,d = g.inv(tlon, tlat, _lon, _lat)
            data = tr.data[:].astype(np.float)
            data /= data.max()
            try:
                id0 = int(d/(red_vel*tr.stats.delta))
            except ZeroDivisionError:
                id0 = 0
            data = data[id0:]
            times = np.arange(data.size)*(tr.stats.delta*1e3)
            dates = np.datetime64(tr.stats.starttime.datetime)+times.astype('timedelta64[ms]')
            idates = np.array(dates.astype(np.int) // 10**3).astype(np.float)
            self.curves[(tr.stats.station,'acoustic')] = hv.Curve((idates, data-d/1e3))
            
        return hv.NdOverlay(self.curves, kdims=['name', 'type']) 


hv.extension('bokeh')
hv.output(size=30)
renderer = hv.renderer('bokeh').instance(mode='server')


def apply_axis_formatter(plot, element):
    plot.handles['xaxis'].formatter = bokeh.models.DatetimeTickFormatter(hourmin = '%H:%M', 
                                                                         minutes='%H:%M', 
                                                                         minsec='%H:%M:%S',
                                                                         seconds='%H:%M:%S')
    plot.handles['yaxis'].formatter = FuncTickFormatter(code="""return Math.abs(tick)""")
    plot.handles['xaxis'].axis_label = "UTC time [min]"
    plot.handles['yaxis'].axis_label = "Distance [km]"
   

def modify_doc(doc):

    #end = UTCDateTime(2018,2,11,1,15,00)
    #end = UTCDateTime.utcnow()
    # 1st Te Maari eruption
    #start = UTCDateTime(2012,8,6,11,40)
    #end = start + 20*60.
    # 2nd Te Maari eruption
    #start = UTCDateTime(2012,11,21,0,20)
    #end = tstart + 20*60.
    # 2007 Ruapehu eruption
    #tstart = UTCDateTime(2007,9,25,8,26)
    #tend = tstart + 10*60.
     
    rs = RecordSection()
    update_data = HVStream.define('update_data')
    for _type in ['seismic', 'acoustic']:
        for slc in stations[rs.target][_type]:
            s, l, c = slc
            msg = "Downloading {:s}.{:s}.{:s}\n".format(s, l, c)
            text.append(msg)
            print(msg)
            tr = get_data(s, l, c, rs.start, rs.end)
            if tr is not None: 
                rs.streams[_type] += tr
            else:
                msg = 'No data for {}.{}.{}'.format(s, l, c)
                text.append(msg)
                print(msg)
    
    dm = hv.DynamicMap(rs.build_record_section, streams=[update_data(transient=True)])
    color_key = {'seismic': 'blue', 'acoustic': 'red'}    
    lyt = datashade(dm, aggregator=ds.count_cat('type'),
                    color_key=color_key, 
                    min_alpha=255, width=3000, height=2000,
                    streams=[hv.streams.RangeXY(transient=True)])      

    def startdate_update(attrname, old, new):
        rs.start.year = new.year
        rs.start.month = new.month
        rs.start.day = new.day

    def enddate_update(attrname, old, new):
        rs.end.year = new.year
        rs.end.month = new.month
        rs.end.day = new.day

    def starttime_update(attrname, old, new):
        h,m,s = map(int, new.split(':'))
        rs.start.hour = h
        rs.start.minute = m
        rs.start.second = s

    def endtime_update(attrname, old, new):
        h,m,s = map(int, new.split(':'))
        rs.end.hour = h
        rs.end.minute = m
        rs.end.second = s

    def update_target(attrname, old, new):
        rs.target = new.lower()
        dm.event()

    @gen.coroutine
    def update_plot():
        global text
        while len(text) > 50:
            text.pop(0)
        pre.text = ''.join(text) 
        dm.event()


    @gen.coroutine
    @without_document_lock
    def update():
        global text
        executor = ThreadPoolExecutor(max_workers=4)
        msg = "Loading data for {:s} between {:s} and {:s}\n".format(rs.target, str(rs.start), str(rs.end))
        text.append(msg)
        print(msg)
        rs.reset()
        for _type in ['seismic', 'acoustic']:
            for slc in stations[rs.target][_type]:
                s, l, c = slc
                msg = "Downloading {:s}.{:s}.{:s}\n".format(s, l, c)
                text.append(msg)
                tr = yield executor.submit(get_data, s, l, c, rs.start, rs.end)
                if tr is not None:
                    rs.streams[_type] += tr
                else:
                    msg = 'No data for {}.{}.{}\n'.format(s, l, c)
                    text.append(msg)
                # The sleep has to be added to prevent some kind of racing condition
                # in the underlying bokeh implementation
                time.sleep(0.1)
                doc.add_next_tick_callback(update_plot)
        text.append("Loading finished.\n")

    def reset():
        dm.event()

    sdateval = "{:d}-{:d}-{:d}".format(rs.start.year,
                                       rs.start.month,
                                       rs.start.day)
    date_start = DatePicker(title='Start date', value=sdateval)
    date_start.on_change('value', startdate_update)

    edateval = "{:d}-{:d}-{:d}".format(rs.end.year,
                                       rs.end.month,
                                       rs.end.day)
    date_end = DatePicker(title='End date', value=edateval)
    date_end.on_change('value', enddate_update)
    
    stimeval = "{:02d}:{:02d}:{:02d}".format(rs.start.hour,
                                             rs.start.minute,
                                             rs.start.second)
    starttime = TextInput(title='Start time', value=stimeval)
    starttime.on_change('value', starttime_update)

    etimeval = "{:02d}:{:02d}:{:02d}".format(rs.end.hour,
                                             rs.end.minute,
                                             rs.end.second)
    endtime = TextInput(title='End time', value=etimeval)
    endtime.on_change('value', endtime_update)
    
    updateb = Button(label='Update', button_type='success')
    updateb.on_click(update)
    
    resetb = Button(label='Reset', button_type='success')
    resetb.on_click(reset)
    
    select_target = Select(title='Volcano', value="Tongariro", options=["Tongariro", "Ruapehu"])
    select_target.on_change('value', update_target)

    # Create HoloViews plot and attach the document
    lyt = lyt.opts(plot=dict(width=3000, height=2000, finalize_hooks=[apply_axis_formatter]),
                   norm=dict(framewise=True))
    hvplot = renderer.get_plot(lyt, doc)
    doc.add_root(layout([[hvplot.state, widgetbox(pre)], 
                         [widgetbox(date_start, starttime), 
                          widgetbox(date_end, endtime),
                          widgetbox(select_target, updateb, resetb)]], sizing_mode='fixed'))
    return doc

doc = modify_doc(doc)
