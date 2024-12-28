#!/usr/bin/env bash

set -e

export USERNAME=$1
export USER_PASSWORD=$1

# Remove build in user for Noble version
if id -u "ubuntu" >/dev/null 2>&1; then
  userdel -r ubuntu
fi

echo "Adding user"
useradd -m $USERNAME
echo "$USERNAME:$USER_PASSWORD" | chpasswd
usermod --shell /bin/bash $USERNAME
usermod -aG sudo $USERNAME
usermod -aG systemd-journal $USERNAME
usermod -aG video $USERNAME
echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/$USERNAME
chmod 0440 /etc/sudoers.d/$USERNAME
