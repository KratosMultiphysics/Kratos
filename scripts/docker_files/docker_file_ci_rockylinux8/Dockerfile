FROM rockylinux:8

USER root

ENV HOME /root

# install python3.8.10
RUN dnf update -y && \
    dnf groupinstall "Development Tools" -y && \
    dnf install -y wget openssl-devel bzip2-devel libffi-devel && \
    wget https://www.python.org/ftp/python/3.8.10/Python-3.8.10.tgz && \
    tar xzf Python-3.8.10.tgz && \
    cd Python-3.8.10 && \
    ./configure --enable-optimizations && \
    make altinstall && \
    cd .. && \
    rm -rf Python-3.8.10 Python-3.8.10.tgz

ENV PYTHON_EXECUTABLE /usr/local/bin/python3.8

RUN python3.8 -m pip install numpy sympy scipy parameterized

# install build utils
RUN dnf update -y && dnf install -y \
    --enablerepo=devel ninja-build \
    cmake \
    gcc-gfortran \
    --enablerepo=devel blas-devel \
    --enablerepo=devel lapack-devel

CMD [ "/bin/bash" ]

WORKDIR $HOME
