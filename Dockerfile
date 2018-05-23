# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
FROM jupyter/base-notebook

MAINTAINER "Yannik Behr <y.behr@gns.cri.nz>"

USER root

# libav-tools for matplotlib anim
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    vim && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $NB_USER

# Install Python 3 packages
# Remove pyqt and qt pulled in for matplotlib since we're only ever going to
# use notebook-friendly backends in these images
RUN conda install --quiet --yes -c pyviz \
    'pandas<=0.22' \
    pyviz &&\
    conda remove --quiet --yes --force qt pyqt && \
    conda clean -tipsy && \
    fix-permissions $CONDA_DIR

RUN conda install --quiet --yes \
    'obspy=1.1*' \
    'pyproj=1.9*' && \
    conda clean -tipsy && \
    fix-permissions $CONDA_DIR



COPY record_section_hv.py /home/$NB_USER/

CMD ["bokeh", "serve", "--port", "3000", "--allow-websocket-origin", "vulkan.gns.cri.nz:3000", "--show", "record_section_hv.py"]


