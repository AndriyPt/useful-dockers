ARG VERSION=jammy

FROM useful-dockers/ngbem:${VERSION}

LABEL org.opencontainers.image.authors="Andriy.Petlovanyy@noemail.com"

LABEL Description="This image allows to run NGSolve, Tex and ParaView in Docker container" Vendor="Andriy Petlovanyy" Version="1.0"

RUN set -x && \
    \
    echo "Update apt repositories cache" && \
    apt-get update -y && \
    \
    echo "Install VSCode Extensions" && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension james-yu.latex-workshop && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension stkb.rewrap && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension mathematic.vscode-pdf && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension mathematic.vscode-latex && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension streetsidesoftware.code-spell-checker && \
    /usr/local/bin/gosu $MY_USERNAME code --install-extension streetsidesoftware.code-spell-checker-ukrainian && \
    \
    echo "Installing Latex" && \
    apt-get update -y && \
    apt-get install -y texstudio texlive* pandoc evince && \ 
    apt-get install -y scalable-cyrfonts-tex && \
    \
    echo "Installing ParaView" && \
    wget -O /tmp/ParaView.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.13&type=binary&os=Linux&downloadFile=ParaView-5.13.0-MPI-Linux-Python3.10-x86_64.tar.gz" && \
    mkdir -p /opt/paraview && \
    tar -xvzf /tmp/ParaView.tar.gz -C /opt/paraview && \
    echo "alias pview=/opt/paraview/ParaView-5.13.0-MPI-Linux-Python3.10-x86_64/bin/paraview" >> /home/$USERNAME/.bashrc && \
    \
    echo "Removing cache" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /home/$MY_USERNAME/.cache

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["/usr/bin/terminator"]
