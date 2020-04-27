FROM quay.io/pypa/manylinux2010_x86_64

RUN yum install -y git wget zip

RUN yum install -y libgfortran.x86_64 blas-devel lapack-devel

RUN mkdir -p /workspace/scripts; \
	mkdir -p /workspace/kratos

RUN mkdir /workspace/boost; \
	wget -P /workspace/boost https://dl.bintray.com/boostorg/release/1.71.0/source/boost_1_71_0.tar.gz; \
	tar -C /workspace/boost -xzf /workspace/boost/boost_1_71_0.tar.gz; \
	rm /workspace/boost/boost_1_71_0.tar.gz

RUN mkdir /workspace/cmake; \
	wget -P /workspace/cmake https://github.com/Kitware/CMake/releases/download/v3.15.3/cmake-3.15.3-Linux-x86_64.tar.gz; \
	tar -C /workspace/cmake -xzf /workspace/cmake/cmake-3.15.3-Linux-x86_64.tar.gz; \
	rm /workspace/cmake/cmake-3.15.3-Linux-x86_64.tar.gz

ENV CMAKE /workspace/cmake/cmake-3.15.3-Linux-x86_64/bin/cmake

COPY start.sh /workspace/scripts/start.sh

RUN chmod +x /workspace/scripts/start.sh

ENTRYPOINT ["/bin/sh", "/workspace/scripts/start.sh"]
