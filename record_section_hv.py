from functools import partial
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
from bokeh.driving import count
import datashader as ds
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException, FDSNNoDataException
from obspy import UTCDateTime, Stream
from pyproj import Geod
import numpy as np
import pandas as pd
import os
import pickle


warnings.filterwarnings('ignore')


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


#GeoNet's FDSN web servers
arc_client = Client('http://service.geonet.org.nz')
nrt_client = Client('http://beta-service-nrt.geonet.org.nz')

doc = curdoc()

text = ['Program startup\n']
pre = PreText(text='Program startup', width=500, height=100)

@count()
def update_log(t):
    global text
    if len(text) > 50:
        text.pop(0)
    pre.text = ''.join(text) 


def GeoNetFDSNrequest(date1, date2, net, sta, loc, cmp):
    """
    Request waveform data from GeoNet's FDSN webservices.
    """
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


def get_data(stations, target, tstart, tend, prefix, new=False):
    global text
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
                msg = "Downloading {:s}.{:s}.{:s}\n".format(s, str(loc), cmp)
                text.append(msg)
                st_tmp, inv = GeoNetFDSNrequest(tstart, tend, 'NZ', s, str(loc), cmp)
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
            _s = inv[0][0]
            _,_,d = g.inv(target[0], target[1], _s.longitude, _s.latitude)
            tr.stats.distance = d
            st += tr
            st_max = max(st_max, tr.data.max())
            with open(fout, 'wb') as fh:
                pickle.dump(st, fh)
    return st


def record_section_hv(start, end, target, red_vel=0., new=False):
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


hv.extension('bokeh')
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

end = UTCDateTime(2018,2,11,1,15,00)
start = end - 10*60.
target = 'tongariro'

glob_var = {'syear': start.year,
            'smonth': start.month,
            'sday': start.day,
            'shour': start.hour,
            'sminute': start.minute,
            'ssecond': start.second,
            'eyear': end.year,
            'emonth': end.month,
            'eday': end.day,
            'ehour': end.hour,
            'eminute': end.minute,
            'esecond': end.second,
            'target': target}


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
     
    update_data = HVStream.define('update_date', target=target, start=start, end=end)
    dm = hv.DynamicMap(record_section_hv, streams=[update_data()])
    
    color_key = {'seismic': 'blue', 'acoustic': 'red'}    
    lyt = datashade(dm, aggregator=ds.count_cat('type'),
                    color_key=color_key, 
                    min_alpha=255, width=3000, height=2000,
                    streams=[hv.streams.RangeXY(transient=True)])      

    def startdate_update(attrname, old, new):
        global glob_var 
        glob_var['syear'] = new.year
        glob_var['smonth'] = new.month
        glob_var['sday'] = new.day

    def enddate_update(attrname, old, new):
        global glob_var
        glob_var['eyear'] = new.year
        glob_var['emonth'] = new.month
        glob_var['eday'] = new.day

    def starttime_update(attrname, old, new):
        global glob_var
        h,m,s = map(int, new.split(':'))
        glob_var['shour'] = h
        glob_var['sminute'] = m
        glob_var['ssecond'] = s

    def endtime_update(attrname, old, new):
        global glob_var
        h,m,s = map(int, new.split(':'))
        glob_var['ehour'] = h
        glob_var['eminute'] = m
        glob_var['esecond'] = s

    def update_target(attrname, old, new):
        global glob_var
        glob_var['target'] = new.lower()

    def update():
        global text
        start = UTCDateTime(glob_var['syear'], glob_var['smonth'],
                            glob_var['sday'], glob_var['shour'],
                            glob_var['sminute'], glob_var['ssecond'])
        end = UTCDateTime(glob_var['eyear'], glob_var['emonth'],
                            glob_var['eday'], glob_var['ehour'],
                            glob_var['eminute'], glob_var['esecond'])
        msg = "Loading data for {:s} between {:s} and {:s}\n".format(glob_var['target'], str(start), str(end))
        text.append(msg)
        dm.event(start=start, end=end, target=glob_var['target'])
        text.append("Loading finished.")

    sdateval = "{:d}-{:d}-{:d}".format(glob_var['syear'],
                                       glob_var['smonth'],
                                       glob_var['sday'])
    date_start = DatePicker(title='Start date', value=sdateval)
    date_start.on_change('value', startdate_update)

    edateval = "{:d}-{:d}-{:d}".format(glob_var['eyear'],
                                       glob_var['emonth'],
                                       glob_var['eday'])
    date_end = DatePicker(title='End date', value=edateval)
    date_end.on_change('value', enddate_update)
    
    stimeval = "{:02d}:{:02d}:{:02d}".format(glob_var['shour'],
                                       glob_var['sminute'],
                                       glob_var['ssecond'])
    starttime = TextInput(title='Start time', value=stimeval)
    starttime.on_change('value', starttime_update)

    etimeval = "{:02d}:{:02d}:{:02d}".format(glob_var['ehour'],
                                       glob_var['eminute'],
                                       glob_var['esecond'])

    endtime = TextInput(title='End time', value=etimeval)
    endtime.on_change('value', endtime_update)
    
    updateb = Button(label='Update', button_type='success')
    updateb.on_click(update)
    
    select_target = Select(title='Volcano', value="Tongariro", options=["Tongariro", "Ruapehu"])
    select_target.on_change('value', update_target)

    # Create HoloViews plot and attach the document
    lyt = lyt.opts(plot=dict(width=3000, height=2000, finalize_hooks=[apply_axis_formatter]),
                   norm=dict(framewise=True))
    hvplot = renderer.get_plot(lyt, doc)
    doc.add_root(layout([[hvplot.state, widgetbox(pre)], 
                         [widgetbox(date_start, starttime), 
                          widgetbox(date_end, endtime),
                          widgetbox(select_target, updateb)]], sizing_mode='fixed'))
    return doc

doc.add_periodic_callback(update_log, 200)
doc = modify_doc(doc)
