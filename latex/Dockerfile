FROM shadowrobot/build-tools:focal-noetic

MAINTAINER "Andriy Petlovanyy"

LABEL Description="This image allows to run NGSolve, Tex and ParaView in Docker container" Vendor="Andriy Petlovanyy" Version="1.0"

RUN set -x && \
    \
    echo "Install VSCode Extensions" && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension james-yu.latex-workshop && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension stkb.rewrap && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension mathematic.vscode-pdf && \
    \
    echo "Installing Latex" && \
    apt-get update -y && \
    apt-get install -y texstudio texlive scalable-cyrfonts-tex pandoc mupdf && \
    \
    echo "Installing NGSolve" && \
    apt-get install -y software-properties-common && \
    apt-add-repository -y universe && \
    add-apt-repository -y ppa:ngsolve/ngsolve && \
    apt-get update -y && \
    apt-get install -y ngsolve && \
    \
    echo "Installing ParaView" && \
    wget -O /tmp/ParaView.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.13&type=binary&os=Linux&downloadFile=ParaView-5.13.0-MPI-Linux-Python3.10-x86_64.tar.gz" && \
    mkdir -p /opt/paraview && \
    tar -xvzf /tmp/ParaView.tar.gz -C /opt/paraview && \
    \
    echo "Removing cache" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /home/$MY_USERNAME/.cache

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["/usr/bin/terminator"]
