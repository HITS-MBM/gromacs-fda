FROM bernddoser/docker-devel-cpp:ubuntu-16.04-gcc-4.9-tools-1

MAINTAINER Bernd Doser <bernd.doser@h-its.org>

RUN apt-get update \
 && apt-get install -y \
    libxml2-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
