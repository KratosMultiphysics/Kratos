FROM registry.access.redhat.com/ubi8-init

# ENV MMG_ROOT=/external_libraries/mmg/mmg_5_5_1
# ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:{MMG_ROOT}/lib

# # Prepare package dependencies
RUN mkdir -p /VW
RUN yum update -y
RUN yum install -y python3.11
RUN which python3
# RUN /opt/python/cp310-cp310/bin/python3 -m venv /VW/python_venv
RUN python3.11 -m venv /VW/python_venv
RUN source /VW/python_venv/bin/activate 
#&& python3 -m pip install numpy

# # Prepare Shared Volumes
COPY ./VW_Release2024/ /VW/Wheel_Sources/

COPY kratos_test/start.sh /VW/start.sh

RUN chmod +x /VW/start.sh
RUN chmod +x /VW/python_venv/bin/activate

ENTRYPOINT ["/bin/bash", "/VW/start.sh"]
