FROM quay.io/pypa/manylinux_2_28_x86_64
ENV MMG_ROOT=/external_libraries/mmg/mmg_5_5_1
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:{MMG_ROOT}/lib

# Prepare package dependencies
RUN yum install -y git wget zip metis metis-devel

# Prepare Shared Volumes
RUN mkdir -p /data_swap_guest

# Prepare File System
RUN mkdir -p /workspace/scripts; \
	mkdir -p /workspace/kratos; \
    mkdir -p /workspace/boost

# Install Boost (form zip)
RUN wget -P /workspace/boost https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz; \
	tar -C /workspace/boost -xzf /workspace/boost/boost_1_87_0.tar.gz; \
	rm /workspace/boost/boost_1_87_0.tar.gz

# Install setuptools
# 3.8
RUN yum update -y && \
    yum install -y openssl-devel bzip2-devel libffi-devel zlib-devel sqlite-devel xz-devel && \
    yum clean all

# 2. Build Python 3.8 (Last security release) into the manylinux expected path
RUN wget https://www.python.org/ftp/python/3.8.20/Python-3.8.20.tgz && \
    tar -xzf Python-3.8.20.tgz && \
    cd Python-3.8.20 && \
    ./configure \
        --prefix=/opt/python/cp38-cp38 \
        --enable-shared \
        --with-openssl=/usr \
        LDFLAGS="-Wl,-rpath=/opt/python/cp38-cp38/lib" && \
    make -j$(nproc) && \
    make altinstall && \
    ln -sf /opt/python/cp38-cp38/bin/python3.8 /opt/python/cp38-cp38/bin/python && \
    cd .. && rm -rf Python-3.8.20*

RUN /opt/python/cp38-cp38/bin/python -m pip install --upgrade pip && \
    /opt/python/cp38-cp38/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.9
RUN /opt/python/cp39-cp39/bin/python -m pip install --upgrade pip && \
    /opt/python/cp39-cp39/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling 
# 3.10
RUN /opt/python/cp310-cp310/bin/python -m pip install --upgrade pip && \
    /opt/python/cp310-cp310/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.11
RUN /opt/python/cp311-cp311/bin/python -m pip install --upgrade pip && \
    /opt/python/cp311-cp311/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.12
RUN /opt/python/cp312-cp312/bin/python -m pip install --upgrade pip && \
    /opt/python/cp312-cp312/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.13
RUN /opt/python/cp313-cp313/bin/python -m pip install --upgrade pip && \
    /opt/python/cp313-cp313/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.14
RUN /opt/python/cp314-cp314/bin/python -m pip install --upgrade pip && \
    /opt/python/cp314-cp314/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling

# 1. Enable CRB/EPEL and install dependencies (skip pmix-devel!)
RUN yum install -y 'dnf-command(config-manager)' && \
    yum config-manager --set-enabled powertools && \
    yum install -y epel-release && \
    yum install -y libevent-devel hwloc-devel libfabric-devel libpsm2-devel numactl-devel libuuid-devel libnl3-devel

# 2. Build modern OpenPMIx 4.2.8 from source
# RUN wget https://github.com/openpmix/openpmix/releases/download/v4.2.8/pmix-4.2.8.tar.bz2 && \
#     tar -xjf pmix-4.2.8.tar.bz2 && \
#     cd pmix-4.2.8 && \
#     ./configure --prefix=/usr/local/pmix && \
#     make -j$(nproc) && \
#     make install && \
#     cd .. && rm -rf pmix-4.2.8*

# 3. Build MPICH 4.2.0 pointing to the custom PMIx
RUN wget https://www.mpich.org/static/downloads/4.2.0/mpich-4.2.0.tar.gz && \
    tar -xzvf mpich-4.2.0.tar.gz && \
    cd mpich-4.2.0 && \
    ./configure --prefix=/usr/local/mpich-4.2 \
                --disable-fortran \
                --disable-cxx \
                --enable-abi-compatibility \
                --with-device=ch4:ofi \
                --with-libfabric && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf mpich-4.2.0*

# Install MMG 5.5.1
# Note ( upgraded from 5.4.1 because of https://github.com/MmgTools/mmg/issues/85)
RUN mkdir -p /workspace/mmg_5_5_1
RUN git clone -b 'v5.5.1' --depth 1 https://github.com/MmgTools/mmg /workspace/mmg_5_5_1 && \
    mkdir /workspace/mmg_5_5_1/build && \
    mkdir -p /workspace/external_libraries/mmg/mmg_5_5_1 && \
    cd /workspace/mmg_5_5_1/build && \
    cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_C_FLAGS="-w" -DCMAKE_CXX_FLAGS="-w" -DCMAKE_INSTALL_PREFIX="/workspace/external_libraries/mmg/mmg_5_5_1" -DLIBMMG3D_SHARED=ON -DLIBMMG2D_SHARED=ON -DLIBMMGS_SHARED=ON -DLIBMMG_SHARED=ON && \
    make -j24 install && \
    rm -r /workspace/mmg_5_5_1 && \
    cd /

# Install Hdf5
RUN mkdir /workspace/hdf5; \
    mkdir /workspace/hdf5/build; \
    mkdir /workspace/hdf5/bin; \
    mkdir /workspace/hdf5/source; \
    cd /workspace/hdf5; \
    curl -L -O -H "User-Agent: Mozilla/5.0" "https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_2/source/hdf5-1.12.2.zip"; \
    cd /; \
    unzip /workspace/hdf5/hdf5-1.12.2.zip -d /workspace/hdf5/source; \
    cmake -H"/workspace/hdf5/source/hdf5-1.12.2" -B"/workspace/hdf5/build" -DCMAKE_C_COMPILER=/usr/local/mpich-4.2/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/local/mpich-4.2/bin/mpic++  -DCMAKE_INSTALL_PREFIX="/workspace/hdf5/bin" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DHDF5_ENABLE_PARALLEL=ON; \
    cmake --build "/workspace/hdf5/build" -j24 --target install; \
    rm -r /workspace/hdf5/source; \
    rm -r /workspace/hdf5/build

# Install Med
RUN mkdir /workspace/med; \
    mkdir /workspace/med/build; \
    mkdir /workspace/med/bin; \
    mkdir /workspace/med/source; \
    cd /workspace/med; \
    curl -L -O -H "User-Agent: Mozilla/5.0" "https://files.salome-platform.org/Salome/medfile/med-5.0.0.tar.bz2"; \
    cd /; \
    tar -C /workspace/med/source -xvf /workspace/med/med-5.0.0.tar.bz2; \
    cmake -H"/workspace/med/source/med-5.0.0" -B"/workspace/med/build" -DCMAKE_C_COMPILER=/usr/local/mpich-4.2/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/local/mpich-4.2/bin/mpic++ -DCMAKE_INSTALL_PREFIX="/workspace/med/bin" -DMEDFILE_BUILD_TESTS=OFF -DMEDFILE_BUILD_SHARED_LIBS=ON -DHDF5_ROOT="/workspace/hdf5/bin"; \
    cmake --build "/workspace/med/build" -j24 --target install; \
    rm -r /workspace/med/source; \
    rm -r /workspace/med/build

COPY start.sh /workspace/scripts/start.sh

RUN chmod +x /workspace/scripts/start.sh

ENTRYPOINT ["/bin/bash", "/workspace/scripts/start.sh"]
