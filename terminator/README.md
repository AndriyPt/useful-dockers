# Docker Base UI image with Terminator

## Build

### Ubuntu Jammy

```bash
docker build --no-cache --progress=plain --build-arg VERSION=jammy -t useful-dockers/terminator:jammy .
```

### The Latest Ubuntu

```bash
docker build --no-cache --progress=plain --build-arg VERSION=latest -t useful-dockers/terminator:latest .
```

## Run

### Ubuntu Jammy

```bash
docker run -it --name terminator -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/terminator:jammy
```

### The Latest Ubuntu

```bash
docker run -it --name terminator -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/terminator:latest
```
