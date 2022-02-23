FROM quay.io/pypa/manylinux2014_x86_64

ENV MMG_ROOT=/external_libraries/mmg/mmg_5_5_1
ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:{MMG_ROOT}/lib

RUN yum install -y git wget zip

# Prepare File system
RUN mkdir -p /data_swap_guest

RUN mkdir -p /workspace/scripts; \
	mkdir -p /workspace/kratos

# Install Boost
RUN /workspace/boost; \
    git clone -b 'boost-1.71.0' --depth 1 --recurse-submodules https://github.com/boostorg/boost.git boost_1_71_0

# Install CMake
RUN mkdir /workspace/cmake; \
	wget -P /workspace/cmake https://github.com/Kitware/CMake/releases/download/v3.15.3/cmake-3.15.3-Linux-x86_64.tar.gz; \
	tar -C /workspace/cmake -xzf /workspace/cmake/cmake-3.15.3-Linux-x86_64.tar.gz; \
	rm /workspace/cmake/cmake-3.15.3-Linux-x86_64.tar.gz

# Install MMG 5.5.1
# Note ( upgraded from 5.4.1 because of https://github.com/MmgTools/mmg/issues/85)
RUN mkdir -p /workspace/mmg_5_5_1
RUN git clone -b 'v5.5.1' --depth 1 https://github.com/MmgTools/mmg /workspace/mmg_5_5_1 && \
    mkdir /workspace/mmg_5_5_1/build && \
    mkdir -p /workspace/external_libraries/mmg/mmg_5_5_1 && \
    cd /workspace/mmg_5_5_1/build && \
    cmake .. -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_C_FLAGS="-w" -DCMAKE_CXX_FLAGS="-w" -DCMAKE_INSTALL_PREFIX="/workspace/external_libraries/mmg/mmg_5_5_1" -DLIBMMG3D_SHARED=ON -DLIBMMG2D_SHARED=ON -DLIBMMGS_SHARED=ON -DLIBMMG_SHARED=ON && \
    make -j2 install && \
    rm -r /workspace/mmg_5_5_1 && \
    cd /

ENV CMAKE /workspace/cmake/cmake-3.15.3-Linux-x86_64/bin/cmake

COPY start.sh /workspace/scripts/start.sh

RUN chmod +x /workspace/scripts/start.sh

ENTRYPOINT ["/bin/bash", "/workspace/scripts/start.sh"]
