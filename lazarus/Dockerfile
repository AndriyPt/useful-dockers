FROM shadowrobot/build-tools:trusty-indigo

MAINTAINER "Andriy Petlovanyy"

LABEL Description="This image allows to run Lazarus UI in Docker container" Vendor="Andriy Petlovanyy" Version="1.0"

RUN set -x && \
    
    echo "Installing needed libraries" && \
    apt-get update && \
    apt-get install -y wget binutils gcc libgtk2.0-0 libgtk2.0-dev psmisc && \
    
    echo "Downloading and installing Debians" && \
    cd /tmp && \
    wget http://downloads.sourceforge.net/project/lazarus/Lazarus%20Linux%20amd64%20DEB/Lazarus%201.6/fpc_3.0.0-151205_amd64.deb && \
    dpkg -i fpc_3.0.0-151205_amd64.deb && rm fpc_3.0.0-151205_amd64.deb && \
    wget http://downloads.sourceforge.net/project/lazarus/Lazarus%20Linux%20amd64%20DEB/Lazarus%201.6/fpc-src_3.0.0-151205_amd64.deb && \
    dpkg -i fpc-src_3.0.0-151205_amd64.deb && rm fpc-src_3.0.0-151205_amd64.deb && \
    wget http://downloads.sourceforge.net/project/lazarus/Lazarus%20Linux%20amd64%20DEB/Lazarus%201.6/lazarus_1.6-0_amd64.deb && \
    dpkg -i lazarus_1.6-0_amd64.deb && rm lazarus_1.6-0_amd64.deb && \

    echo "Removing cache" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["/usr/bin/terminator"]
