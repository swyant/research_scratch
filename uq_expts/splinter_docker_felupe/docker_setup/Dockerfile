FROM rockylinux:8 

RUN useradd -u 1004 swyant

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

RUN dnf -y update && \
    dnf -y install dnf-plugins-core && \ 
    dnf -y install epel-release && \
    dnf config-manager --set-enabled powertools && \
    dnf -y install \
        hdf5 hdf5-devel \
        lapack blas \
        libnl3 libnl3-devel \
        which \
        && dnf clean all

WORKDIR /opt/dakota

COPY Dakota /opt/dakota/Dakota

ENV DAKOTA_DIR=/opt/dakota/Dakota
ENV PATH=$DAKOTA_DIR/bin:$DAKOTA_DIR/share/dakota/test:$DAKOTA_DIR/gui:$PATH 
ENV PYTHONPATH=${PYTHONPATH:-}:$DAKOTA_DIR/share/dakota/Python

USER swyant
RUN uv python pin --global 3.11 # maybe should be 3.12
