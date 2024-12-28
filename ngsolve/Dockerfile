ARG VERSION=jammy

FROM useful-dockers/terminator:${VERSION}

LABEL org.opencontainers.image.authors="Andriy.Petlovanyy@noemail.com"

LABEL Description="This image allows to run NGSolve in Docker container" Vendor="Andriy Petlovanyy" Version="1.0"

COPY compile_ngsolve.sh /tmp/compile_ngsolve.sh

RUN set -x && \
    \
    echo "Installing Python3" && \
    apt-get update -y && \
    apt-get install -y python3-full python3-venv python3-dev && \
    \
    echo "Installing NGSolve" && \
    apt-get -y install python3-distutils python3-tk libpython3-dev libxmu-dev tk-dev tcl-dev cmake && \
    apt-get -y install git g++ libglu1-mesa-dev liblapacke-dev libocct-data-exchange-dev && \
    apt-get -y install libocct-draw-dev occt-misc libtbb-dev libxi-dev && \
    /usr/local/bin/gosu $MY_USERNAME /bin/bash /tmp/compile_ngsolve.sh && \
    \
    echo "Removing cache" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /home/$MY_USERNAME/.cache

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["/usr/bin/terminator"]
