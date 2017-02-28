In order to run this Docker container use the following command

```bash
docker run -it -e DISPLAY  -e LOCAL_USER_ID=$(id -u) -v /tmp/.X11-unix:/tmp/.X11-unix:rw andriyp/lazarus
```
After this window of the terminator console will pop out.
In order to run IDE type the following

```bash
lazarus-ide
```
