FROM nvidia/cuda:11.7.1-devel-ubuntu20.04

LABEL maintainer="Bernd Doser <bernd.doser@h-its.org>"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive TZ=Europe/Berlin apt-get -y install \
    cmake \
    ghostscript \
    git \
    graphviz \
    libboost-graph-dev \
    mscgen \
    texlive \
    texlive-lang-english \
    texlive-latex-extra \
    tzdata \
    wget \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

RUN wget -q https://github.com/google/googletest/archive/release-1.8.0.tar.gz \
 && tar xf release-1.8.0.tar.gz \
 && rm release-1.8.0.tar.gz \
 && cd googletest-release-1.8.0 \
 && mkdir build \
 && cd build \
 && cmake -DCMAKE_INSTALL_PREFIX=/opt/googletest-1.8.0 .. \
 && make -j \
 && make install \
 && rm -r /tmp/googletest-release-1.8.0

RUN wget -q https://sourceforge.net/projects/doxygen/files/rel-$VERSION/doxygen-$VERSION.linux.bin.tar.gz \
 && tar xf doxygen-$VERSION.linux.bin.tar.gz \
 && rm doxygen-$VERSION.linux.bin.tar.gz \
 && mv doxygen-$VERSION /opt \
 && ln -sf /opt/doxygen-$VERSION/bin/doxygen /usr/bin/doxygen

WORKDIR /

ENV GTEST_ROOT /opt/googletest-1.8.0
