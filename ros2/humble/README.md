# Docker ROS2 image with Terminator

## Build

### Ubuntu Jammy

```bash
docker build --no-cache --progress=plain --build-arg BASE_IMAGE="osrf/ros" --build-arg VERSION=humble-desktop-full-jammy -t useful-dockers/ros2:humble-desktop-full-jammy .
```
## Run

### Ubuntu Jammy

```bash
docker run -it --name ros2_jammy -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/ros2:humble-desktop-full-jammy
```
