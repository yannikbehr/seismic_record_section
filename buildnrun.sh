#!/bin/bash

#########################################
# Build and run docker image            #
# 02/18 Y. Behr <y.behr@gns.cri.nz>     #
#########################################

RUNONLY=false
BUILDONLY=false

function usage(){
cat <<EOF
Usage: $0 [Options]

Build and run docker image. 

Options:
    -h              Show this message.
    -r              Only run image without rebuilding it.
    -b              Only rebuild image without running it.
    -t              Assign a tag to the docker image (default: latest).
EOF
}

TAG="latest"

# Processing command line options
while [ $# -gt 0 ]
do
    case "$1" in
        -r) RUNONLY=true;;
        -b) BUILDONLY=true;;
        -t) TAG=$2;shift;;
        -h) usage; exit 0;;
        -*) usage; exit 1;;
        *) break;;
esac
shift
done

if [ "${RUNONLY}" == "false" ]; then
    docker rmi yadabe/rs:$TAG
    docker build --no-cache=true -t yadabe/rs:$TAG .
fi

if [ "${BUILDONLY}" == "false" ] ;then
    docker stop rs
    docker rm rs
    docker run --restart="unless-stopped" --name rs -p 3000:3000 -d yadabe/rs:$TAG 
fi


