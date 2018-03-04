## Interactive seismic record section

A record section for the browser that allows you to pan, zoom, and re-order a set of seismograms according to their distance to a given target.
It uses FDSN-WS to retrieve seismograms and bokeh server to handle interactivity.

### Dependencies
* python=3.6*
* tornado=4.5*
* datashader=0.6*
* holoviews=1.9*
* param=1.5*
* pandas=0.20*
* numpy=1.13*
* scipy=1.0*
* bokeh=0.12*
* obspy=1.1*
* pyproj=1.9*

### Installation
#### Docker
To build and run the docker images execute:
```
./buildnrun.sh
```
Then open a browser at `http:\\localhost:3000`.

#### Manual
The following will open the app in your browser at `http:\\localhost:5006`
```
cd seismic_record_section
bokeh serve --show record_section_hv.py
```
