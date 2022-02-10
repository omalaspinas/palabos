FROM ubuntu:impish

LABEL Description="Palabos Docker image, based on ubuntu:latest " maintainer="orestis.malaspinas@hesge.ch" Version="0.1"

RUN sed -i'' 's/archive.ubuntu.com/ch.archive.ubuntu.com/' /etc/apt/sources.list 
RUN apt update && apt dist-upgrade -y
ENV TZ=Europe/Zurich
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt install -y g++ clang clang-format-13 ccache cmake python3 make libopenmpi-dev libhdf5-dev libhdf5-mpi-dev doxygen imagemagick && update-alternatives --install /usr/bin/python python /usr/bin/python3 10
RUN ln -s /usr/bin/clang-format-13 /usr/bin/clang-format
