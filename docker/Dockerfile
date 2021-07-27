FROM ubuntu:latest

LABEL Description="Palabos Docker image, based on ubuntu:latest " maintainer="orestis.malaspinas@hesge.ch" Version="0.1"

RUN sed -i'' 's/archive.ubuntu.com/us.archive.ubuntu.com/' /etc/apt/sources.list 
ENV TZ=Europe/Zurich
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt update && apt dist-upgrade -y

RUN apt install -y g++ clang ccache cmake python3 make libopenmpi-dev libhdf5-dev libhdf5-mpi-dev doxygen && update-alternatives --install /usr/bin/python python /usr/bin/python3 10
RUN apt install -y imagemagick
