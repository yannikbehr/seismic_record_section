import asyncio
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import os
import pickle
import time
import warnings

import holoviews as hv
from holoviews.operation.datashader import datashade
from holoviews import streams
import bokeh
from bokeh.layouts import layout, widgetbox, row
from bokeh.plotting import curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.models.sources import ColumnDataSource
from bokeh.models import (FuncTickFormatter, 
                          DatetimeTickFormatter,
                          DatePicker, 
                          TextInput, 
                          Button,
                          Select,
                          PreText,
                          CheckboxGroup,
                          Slider,
                          Div)
from bokeh.document import without_document_lock
import param
import datashader as ds
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException, FDSNNoDataException
from obspy import UTCDateTime, Stream
from obspy.core.util import AttribDict
from pyproj import Geod
import numpy as np
import pandas as pd


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
stations['te maari'] = {'seismic':tnp_seismic,
                         'acoustic':tnp_acoustic}
stations['ruapehu'] = {'seismic':tnp_seismic,
                       'acoustic':tnp_acoustic}
stations['red crater'] = {'seismic':tnp_seismic,
                         'acoustic':tnp_acoustic}
stations['ngauruhoe'] = {'seismic':tnp_seismic,
                         'acoustic':tnp_acoustic}


def GeoNetFDSNrequest(date1, date2, net, sta, loc, cmp):
    """
    Request waveform and meta data from GeoNet's FDSN webservices.
    """
    #GeoNet's FDSN web servers
    arc_client = Client('http://service.geonet.org.nz')
    nrt_client = Client('http://service-nrt.geonet.org.nz')
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


def apply_axis_formatter(plot, element):
    plot.handles['xaxis'].formatter = bokeh.models.DatetimeTickFormatter(hourmin = '%H:%M', 
                                                                         minutes='%H:%M', 
                                                                         minsec='%H:%M:%S',
                                                                         seconds='%H:%M:%S')
    plot.handles['yaxis'].formatter = FuncTickFormatter(code="""return Math.abs(tick)""")
    plot.handles['xaxis'].axis_label = "UTC time [min]"
    plot.handles['yaxis'].axis_label = "Distance [km]"


class RecordSection:
    """
    Construct the record section and save the state of the plot.
    """
    def __init__(self):
        self.streams = {'seismic': Stream(), 'acoustic': Stream()}
        self.targets = {'te maari': (175.670744, -39.109752),
                        'ruapehu': (175.564490, -39.281149),
                        'ngauruhoe': (175.632169, -39.156791),
                        'red crater': (175.650761, -39.136715)}
        self.end = UTCDateTime.utcnow()
        self.start = self.end - 20*60.
        self.min_dist = 0.
        self.max_dist = 40.
        self.tmin_old = self.start
        self.tmax_old = self.end
        self.ymin_old = self.min_dist
        self.ymax_old = self.max_dist 
        self.target = 'ruapehu'
        self.curves = {}
        self.plot_acoustic = True
        self.plot_seismic = True
        self.plot_labels = False
        self.red_vel=0.

    def reset(self):
        """
        Reset state before downloading new data.
        """
        self.streams = {'seismic': Stream(), 'acoustic': Stream()}
        self.curves = {}
        
    def build_record_section(self, x_range=None, y_range=None, replot=False,
                             new=False):
        tmin = None
        tmax = None
        labels = []
        stats_list = []
        if x_range is not None:
            xmin, xmax = x_range
            tmin = datetime.utcfromtimestamp(xmin/1e3)
            tmax = datetime.utcfromtimestamp(xmax/1e3)
            self.tmin_old = tmin
            self.tmax_old = tmax
        else:
            tmin = self.tmin_old
            tmax = self.tmax_old

        if y_range is not None:
            ymax = abs(min(0, y_range[0]))
            ymin = abs(min(0, y_range[1]))
            self.ymin_old = ymin
            self.ymax_old = ymax
        else:
            ymin = self.ymin_old
            ymax = self.ymax_old
            y_range = (-ymax, -ymin)

        if replot:
            tmin = None
            tmax = None
            xmin = np.datetime64(self.start.datetime).astype(np.int)//1e3
            xmax = np.datetime64(self.end.datetime).astype(np.int)//1e3
            x_range=(xmin,xmax)
            ymin = self.min_dist
            ymax = self.max_dist
            y_range = (-ymax, -ymin)

        g = Geod(ellps='WGS84')
        tlon, tlat = self.targets[self.target]
        st_seismic = self.streams['seismic'].copy()
        st_acoustic = self.streams['acoustic'].copy()

        if None not in [tmin, tmax]:
            starttime = UTCDateTime(tmin)
            endtime = UTCDateTime(tmax)
            st_seismic.trim(starttime, endtime)
            st_acoustic.trim(starttime, endtime)
        print(tmin, tmax, ymin, ymax)

        # Get rid of dummy curves
        self.curves = {k:self.curves[k] for k in self.curves if k[1] != "dummy"}

        if st_seismic.count() > 0 and self.plot_seismic:
            local_max_s = 0.
            for tr in st_seismic:
                _lat = tr.stats.coordinates.latitude
                _lon = tr.stats.coordinates.longitude
                _,_,d = g.inv(tlon, tlat, _lon, _lat)
                if d/1e3 >= ymin and d/1e3 <= ymax:
                    local_max_s = max(local_max_s, tr.data.max())
                tr.stats.distance = d/1e3

            for tr in st_seismic:
                key = (tr.stats.station, 'seismic')
                if tr.stats.distance < ymin or tr.stats.distance > ymax:
                    if key in self.curves:
                        self.curves.pop(key)
                    continue
                data = tr.data[:].astype(np.float)
                data /= local_max_s
                try:
                    t0 = int(tr.stats.distance*1e3/self.red_vel*1e3)
                    id0 = int(tr.stats.distance*1e3/(self.red_vel*tr.stats.delta))
                except ZeroDivisionError:
                    id0 = 0
                    t0 = 0
                times = np.arange(data.size)*(tr.stats.delta*1e3)
                dates = np.datetime64(tr.stats.starttime.datetime)+times.astype('timedelta64[ms]') 
                dates -= np.timedelta64(t0, 'ms')
                idates = np.array(dates.astype(np.int) // 1e3).astype(np.float)
                self.curves[key]  = hv.Curve((idates, data-tr.stats.distance))
                if tr.stats.station not in stats_list:
                    stats_list.append(tr.stats.station)
                    if self.plot_labels:
                        label = tr.stats.station
                    else:
                        label = ''
                    labels.append(hv.Text(idates[id0+100], -tr.stats.distance,
                                          label, 
                                          halign='left',
                                          valign='bottom').opts(norm=dict(framewise=True)))           
        if not self.plot_seismic:
            self.curves = {k:self.curves[k] for k in self.curves if k[1] != "seismic"}

        if st_acoustic.count() > 0 and self.plot_acoustic: 
            local_max_a = 0.
            for tr in st_acoustic:
                _lat = tr.stats.coordinates.latitude
                _lon = tr.stats.coordinates.longitude
                _,_,d = g.inv(tlon, tlat, _lon, _lat)
                if d/1e3 >= ymin and d/1e3 <= ymax:
                    local_max_a = max(local_max_a, tr.data.max())
                tr.stats.distance = d/1e3   

            for tr in st_acoustic:
                key = (tr.stats.station, 'acoustic')
                if tr.stats.distance < ymin or tr.stats.distance > ymax:
                    if key in self.curves:
                        self.curves.pop(key)
                    continue
                data = tr.data[:].astype(np.float)
                data /= local_max_a 
                try:
                    t0 = int(tr.stats.distance*1e3/self.red_vel*1e3)
                    id0 = int(tr.stats.distance*1e3/(self.red_vel*tr.stats.delta))
                except ZeroDivisionError:
                    id0 = 0
                    t0 = 0
                times = np.arange(data.size)*(tr.stats.delta*1e3)
                dates = np.datetime64(tr.stats.starttime.datetime)+times.astype('timedelta64[ms]')
                dates -= np.timedelta64(t0,'ms')
                idates = np.array(dates.astype(np.int) // 1e3).astype(np.float)
                self.curves[key] = hv.Curve((idates, data-tr.stats.distance))
                if tr.stats.station not in stats_list:
                    stats_list.append(tr.stats.station)
                    if self.plot_labels:
                        label = tr.stats.station
                    else:
                        label = ''
                    labels.append(hv.Text(idates[id0+100], -tr.stats.distance,
                                          label,
                                          halign='left',
                                          valign='bottom').opts(norm=dict(framewise=True)))           
               
        if not self.plot_acoustic:
            self.curves = {k:self.curves[k] for k in self.curves if k[1] != "acoustic"}
        if not self.curves:
            # Plot dummy curves in case there is no data
            x0 = np.datetime64(self.start.datetime).astype('int')//1e3
            x1 = np.datetime64(self.end.datetime).astype('int')//1e3
            dummy_x = np.linspace(x0, x1, 100)
            dummy_y = np.ones(100)*y_range[1]
            for i in range(len(tnp_seismic)+len(tnp_acoustic)):
                self.curves[('No_%d'%i, 'dummy')] = hv.Curve((dummy_x, dummy_y))
                labels.append(hv.Text(x0, 0., '').opts(norm=dict(framewise=True)))
           
        color_key = {'seismic':'blue', 'acoustic':'red', 'dummy':'grey'}    
        lyt = datashade(hv.NdOverlay(self.curves, kdims=['name', 'type']),
                        aggregator=ds.count_cat('type'),
                        color_key=color_key, dynamic=False, 
                        min_alpha=255, width=3000, height=2000, x_range=x_range, y_range=y_range,
                        y_sampling=0.1)      
        return (lyt*hv.Overlay(labels)).opts(plot=dict(width=3000, height=2000, 
                                                           hooks=[apply_axis_formatter]),
                                             norm=dict(framewise=True))



hv.extension('bokeh')
hv.output(size=30)
renderer = hv.renderer('bokeh').instance(mode='server')

doc = curdoc()
text = ['Program startup\n']
pre = PreText(text='Program startup', width=500, height=100)

update_data = streams.Stream.define('update_data', replot=False)
rs = RecordSection()
dm = hv.DynamicMap(rs.build_record_section, streams=[update_data(transient=True),
                                                     streams.RangeXY(transient=False)])


async def update_log():
    global text
    while len(text) > 50:
        text.pop(0)
    pre.text = ''.join(text) 

@without_document_lock
async def initial_load():
    global text
    executor = ThreadPoolExecutor(max_workers=4)
    msg = "Loading data for {:s} between {:s} and {:s}\n".format(rs.target, str(rs.start), str(rs.end))
    text.append(msg)
    rs.reset()
    for _type in ['seismic', 'acoustic']:
        for slc in stations[rs.target][_type]:
            s, l, c = slc
            msg = "Downloading {:s}.{:s}.{:s}\n".format(s, l, c)
            text.append(msg)
            tr = await asyncio.wrap_future(executor.submit(get_data, s, l, c, rs.start, rs.end), loop=None)
            if tr is not None:
                rs.streams[_type] += tr
            else:
                msg = 'No data for {}.{}.{}\n'.format(s, l, c)
                text.append(msg)
                time.sleep(0.1)
                continue
            # The sleep has to be added to prevent some kind of racing condition
            # in the underlying bokeh implementation
            time.sleep(0.1)
            doc.add_next_tick_callback(update_plot)
    text.append("LOADING FINISHED.\n")

async def update_plot():
    global text
    while len(text) > 34:
        text.pop(0)
    pre.text = ''.join(text) 
    dm.event(replot=True)


def modify_doc(doc):
    
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

    @without_document_lock
    async def load():
        global text
        if rs.start > rs.end:
            msg = "Start time is later than end time.\n"
            text.append(msg)
            return
        if (rs.end - rs.start) > 3600.:
            msg = "Length of requested data exceeds 1 hour.\n" 
            text.append(msg)
            return
        executor = ThreadPoolExecutor(max_workers=4)
        msg = "Loading data for {:s} between {:s} and {:s}\n".format(rs.target, str(rs.start), str(rs.end))
        text.append(msg)
        rs.reset()
        for _type in ['seismic', 'acoustic']:
            for slc in stations[rs.target][_type]:
                s, l, c = slc
                msg = "Downloading {:s}.{:s}.{:s}\n".format(s, l, c)
                text.append(msg)
                tr = await asyncio.wrap_future(executor.submit(get_data, s, l, c, rs.start, rs.end), loop=None)
                if tr is not None:
                    rs.streams[_type] += tr
                else:
                    msg = 'No data for {}.{}.{}\n'.format(s, l, c)
                    text.append(msg)
                    time.sleep(0.1)
                    continue
                # The sleep has to be added to prevent some kind of racing condition
                # in the underlying bokeh implementation
                time.sleep(0.1)
                doc.add_next_tick_callback(update_plot)
        text.append("LOADING FINISHED.\n")

    def reset():
        global text
        if rs.start > rs.end:
            msg = "Start time is later than end time.\n"
            text.append(msg)
            return
        if (rs.end - rs.start) > 3600.:
            msg = "Length of requested data exceeds 1 hour.\n" 
            text.append(msg)
            return
        dm.event(replot=True)

    def update_traces(attr, old, new):
        rs.plot_seismic = False
        rs.plot_acoustic = False
        rs.plot_labels = False
        for _c in new:
            if _c == 0:
                rs.plot_seismic = True
            if _c == 1:
                rs.plot_acoustic = True
            if _c == 2:
                rs.plot_labels = True
        dm.event()

    def update_red_vel(attr, old, new):
        if len(new) < 1:
            rs.red_vel = 0
        else:
           rs.red_vel = 330. 
        dm.event()

    sdateval = "{:02d}-{:02d}-{:02d}".format(rs.start.year,
                                       rs.start.month,
                                       rs.start.day)
    date_start = DatePicker(title='Start date', value=sdateval)
    date_start.on_change('value', startdate_update)

    edateval = "{:02d}-{:02d}-{:02d}".format(rs.end.year,
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
    
    loadb = Button(label='Load', button_type='success')
    loadb.on_click(load)
    
    resetb = Button(label='Reset', button_type='success')
    resetb.on_click(reset)
    
    select_target = Select(title='Volcano', value="Ruapehu",
                           options=["Te Maari", "Ruapehu", "Ngauruhoe", "Red Crater"])
    select_target.on_change('value', update_target)

    cg = CheckboxGroup(labels=["Seismic", "Acoustic", "Station name"], active=[1, 1, 0])
    cg.on_change('active', update_traces)
    
    rvel = CheckboxGroup(labels=['Reduction velocity (330 m/s)'], active=[])
    rvel.on_change('active', update_red_vel)

    div = Div(text="""<a target="_blank" href="https://wiki.geonet.org.nz/display/volcano/Record+Section+Plot
">Documentation</a>""", width=80, height=30)
    # Create HoloViews plot and attach the document
    hvplot = renderer.get_plot(dm, doc)
    doc.add_root(layout([[hvplot.state, widgetbox(pre)], 
                         [widgetbox(date_start, starttime, div), 
                          widgetbox(date_end, endtime),
                          widgetbox(select_target, loadb, resetb, cg, rvel)]], sizing_mode='fixed'))
    return doc

doc.add_periodic_callback(update_log, 200)
doc.add_timeout_callback(initial_load, 1000)
doc = modify_doc(doc)

