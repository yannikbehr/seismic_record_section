# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
FROM artifactory.gns.cri.nz/gns/volcano_monitoring_base_image:0.0.1

LABEL maintainer="Yannik Behr <y.behr@gns.cri.nz>"

USER $D_USER

# Install Python 3 packages
# Remove pyqt and qt pulled in for matplotlib since we're only ever going to
# use notebook-friendly backends in these images
RUN mamba install --quiet --yes -c pyviz \
    pandas \
    holoviews \
    bokeh &&\
    conda remove --quiet --yes --force qt pyqt && \
    conda clean -tipsy && \
    fix-permissions $CONDA_DIR

RUN mamba install --quiet --yes -c conda-forge \
    obspy \
    datashader \
    pyproj \
    xarray && \
    conda clean -tipsy && \
    fix-permissions $CONDA_DIR



COPY record_section_hv.py /home/$D_USER/

WORKDIR /home/$D_USER/

CMD ["bokeh", "serve", "--port", "3000", "--allow-websocket-origin", "volcanolab.gns.cri.nz:3000", "--show", "record_section_hv.py"]


