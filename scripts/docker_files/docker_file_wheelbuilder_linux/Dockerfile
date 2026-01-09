FROM quay.io/pypa/manylinux_2_28_x86_64

ENV MMG_ROOT=/external_libraries/mmg/mmg_5_5_1
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:{MMG_ROOT}/lib

# Prepare package dependencies
RUN yum install -y git wget zip mpich mpich-devel metis metis-devel

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
    cmake -H"/workspace/hdf5/source/hdf5-1.12.2" -B"/workspace/hdf5/build"  -DCMAKE_INSTALL_PREFIX="/workspace/hdf5/bin" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON; \
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
    cmake -H"/workspace/med/source/med-5.0.0" -B"/workspace/med/build"  -DCMAKE_INSTALL_PREFIX="/workspace/med/bin" -DMEDFILE_BUILD_TESTS=OFF -DMEDFILE_BUILD_SHARED_LIBS=ON -DHDF5_ROOT="/workspace/hdf5/bin"; \
    cmake --build "/workspace/med/build" -j24 --target install; \
    rm -r /workspace/med/source; \
    rm -r /workspace/med/build

# Install setuptools
# 3.8
RUN /opt/python/cp38-cp38/bin/python -m pip install --upgrade pip
RUN /opt/python/cp38-cp38/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.9
RUN /opt/python/cp39-cp39/bin/python -m pip install --upgrade pip
RUN /opt/python/cp39-cp39/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling 
# 3.10
RUN /opt/python/cp310-cp310/bin/python -m pip install --upgrade pip
RUN /opt/python/cp310-cp310/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.11
RUN /opt/python/cp311-cp311/bin/python -m pip install --upgrade pip
RUN /opt/python/cp311-cp311/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.12
RUN /opt/python/cp312-cp312/bin/python -m pip install --upgrade pip
RUN /opt/python/cp312-cp312/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.13
RUN /opt/python/cp313-cp313/bin/python -m pip install --upgrade pip
RUN /opt/python/cp313-cp313/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling
# 3.14
RUN /opt/python/cp314-cp314/bin/python -m pip install --upgrade pip
RUN /opt/python/cp314-cp314/bin/python -m pip install --upgrade setuptools wheel mypy toml build hatchling


COPY start.sh /workspace/scripts/start.sh

RUN chmod +x /workspace/scripts/start.sh

ENTRYPOINT ["/bin/bash", "/workspace/scripts/start.sh"]
