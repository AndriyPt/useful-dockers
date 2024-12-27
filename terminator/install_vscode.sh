#!/usr/bin/env bash

set -e

export USERNAME=$1

echo "code code/add-microsoft-repo boolean true" | sudo debconf-set-selections

wget -qO- https://packages.microsoft.com/keys/microsoft.asc | gpg --dearmor > packages.microsoft.gpg
install -D -o root -g root -m 644 packages.microsoft.gpg /etc/apt/keyrings/packages.microsoft.gpg
echo "deb [arch=amd64,arm64,armhf signed-by=/etc/apt/keyrings/packages.microsoft.gpg] https://packages.microsoft.com/repos/code stable main" | tee /etc/apt/sources.list.d/vscode.list > /dev/null
rm -f packages.microsoft.gpg

apt-get update

apt-get install -y code

# Write to user profile in case if it was provided
if [ ! -z "$1" ]
  then
    echo "alias vscode='code --no-sandbox'" >> /home/$USERNAME/.bashrc
fi
