# Docker Latex image

## Build

### Prerequisites

Please build base [NGBem Docker image](../ngbem/README.md).
You could also change from to point directly to [Terminator base image](../terminator/README.md) in case if you do not 
need NGSolve and NGBem. 

### Ubuntu Jammy

```bash
docker build --no-cache --progress=plain --build-arg VERSION=jammy -t useful-dockers/latex:jammy .
```

### The Latest Ubuntu

```bash
docker build --no-cache --progress=plain --build-arg VERSION=latest -t useful-dockers/latex:latest .
```

## Run

### Ubuntu Jammy

```bash
docker run -it --name latex -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/latex:jammy
```

### The Latest Ubuntu

```bash
docker run -it --name latex -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/latex:latest
```
