FROM nercury/cmake-cpp:gcc-4.9

MAINTAINER Bernd Doser <bernd.doser@h-its.org>

RUN apt-get update \
 && apt-get install -y \
    libboost-graph-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
