ARG VERSION=jammy
ARG BASE_IMAGE=ubuntu

FROM ${BASE_IMAGE}:${VERSION}

ARG VERSION

LABEL org.opencontainers.image.authors="Andriy.Petlovanyy@noemail.com"

LABEL Description="This image allows to run Terminator in Docker container" Vendor="Andriy Petlovanyy" Version="1.0"

ENV DEBIAN_FRONTEND noninteractive

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENV PULSE_SERVER /run/pulse/native
ENV GOSU_VERSION 1.17
ENV MY_USERNAME user

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
COPY install_gosu.sh /tmp/install_gosu.sh
COPY setup_user.sh /tmp/setup_user.sh
COPY install_vscode.sh /tmp/install_vscode.sh

COPY terminator_config /tmp/terminator_config
COPY additional_bashrc_cmds /tmp/additional_bashrc_cmds

RUN set -x && \
    \
    echo "Installing utils and terminator" && \
    apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates wget gpg sudo curl git vim mc && \
    apt-get install -y apt-transport-https bash-completion terminator && \
    \
    echo "Installing visual, sound libraries and fonts" && \
    apt-get install -y xz-utils file locales dbus-x11 pulseaudio dmz-cursor-theme fonts-dejavu && \ 
    apt-get install -y fonts-liberation hicolor-icon-theme libcanberra-gtk3-0 libcanberra-gtk-module && \ 
    apt-get install -y libcanberra-gtk3-module libgtk2.0-0 libdbus-glib-1-2 libxt6 && \
    apt-get install -y libexif12 libgl1-mesa-dri && \
    if [[ "${VERSION}" == *"jammy" ]]; then apt-get install -y libasound2 libgl1-mesa-glx ; else apt-get install -y libasound2t64 libglx-mesa0 ; fi && \
    \
    echo "Create user" && \
    /bin/bash /tmp/setup_user.sh $MY_USERNAME && \
    \
    echo "Install gosu" && \
    /bin/bash /tmp/install_gosu.sh $GOSU_VERSION && \
    \
    echo "Updating bash commands" && \
    cat /tmp/additional_bashrc_cmds >> /home/user/.bashrc && \
    \
    echo "Setting entrypoint permissions" && \
    chmod 555 /usr/local/bin/entrypoint.sh && \
    \
    echo "Install VSCode" && \
    /bin/bash /tmp/install_vscode.sh $MY_USERNAME && \
    \
    echo "Setting default terminator config" && \
    mkdir -p /home/$MY_USERNAME/.config/terminator/ && \
    mv /tmp/terminator_config /home/$MY_USERNAME/.config/terminator/config && \
    chown -R $MY_USERNAME:$MY_USERNAME /home/$MY_USERNAME/.config &&\
    \
    echo "Removing cache" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /home/$MY_USERNAME/.gitconfig /home/$MY_USERNAME/.cache

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["/usr/bin/terminator"]
