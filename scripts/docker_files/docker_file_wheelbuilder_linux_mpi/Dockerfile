FROM quay.io/pypa/manylinux2014_x86_64

ENV MMG_ROOT=/external_libraries/mmg/mmg_5_5_1
ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:{MMG_ROOT}/lib

# Prepare package dependencies
RUN yum install -y git wget zip mpich mpich-devel mpich-autoload metis metis-devel

# Prepare Shared Volumes
RUN mkdir -p /data_swap_guest

# Prepare File System
RUN mkdir -p /workspace/scripts; \
	mkdir -p /workspace/kratos; \
    mkdir -p /workspace/boost

# Install Boost (from source)
# RUN git clone -b 'boost-1.74.0' --depth 1 --recurse-submodules -j24 https://github.com/boostorg/boost.git /workspace/boost/boost_1_74_0

# Install Boost (form zip)
RUN wget -P /workspace/boost https://archives.boost.io/release/1.74.0/source/boost_1_74_0.tar.gz; \
	tar -C /workspace/boost -xzf /workspace/boost/boost_1_74_0.tar.gz; \
	rm /workspace/boost/boost_1_74_0.tar.gz

# Install MMG 5.5.1
# Note ( upgraded from 5.4.1 because of https://github.com/MmgTools/mmg/issues/85)
RUN mkdir -p /workspace/mmg_5_5_1
RUN git clone -b 'v5.5.1' --depth 1 https://github.com/MmgTools/mmg /workspace/mmg_5_5_1 && \
    mkdir /workspace/mmg_5_5_1/build && \
    mkdir -p /workspace/external_libraries/mmg/mmg_5_5_1 && \
    cd /workspace/mmg_5_5_1/build && \
    cmake .. -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_C_FLAGS="-w" -DCMAKE_CXX_FLAGS="-w" -DCMAKE_INSTALL_PREFIX="/workspace/external_libraries/mmg/mmg_5_5_1" -DLIBMMG3D_SHARED=ON -DLIBMMG2D_SHARED=ON -DLIBMMGS_SHARED=ON -DLIBMMG_SHARED=ON && \
    make -j24 install && \
    rm -r /workspace/mmg_5_5_1 && \
    cd /

COPY start.sh /workspace/scripts/start.sh

RUN chmod +x /workspace/scripts/start.sh

ENTRYPOINT ["/bin/bash", "/workspace/scripts/start.sh"]
