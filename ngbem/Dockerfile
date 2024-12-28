ARG VERSION=jammy

FROM useful-dockers/ngsolve:${VERSION}

LABEL org.opencontainers.image.authors="Andriy.Petlovanyy@noemail.com"

LABEL Description="This image allows to run NGSolve and NGBem in Docker container" Vendor="Andriy Petlovanyy" Version="1.0"

COPY compile_ngbem.sh /tmp/compile_ngbem.sh

RUN set -x && \
    \
    echo "Update apt repositories cache" && \
    apt-get update -y && \
    \
    echo "Installing NGBem" && \
    chown $MY_USERNAME:$MY_USERNAME /tmp/compile_ngbem.sh && \
    /usr/local/bin/gosu $MY_USERNAME /bin/bash /tmp/compile_ngbem.sh && \
    \
    echo "Removing cache" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /home/$MY_USERNAME/.cache

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["/usr/bin/terminator"]
