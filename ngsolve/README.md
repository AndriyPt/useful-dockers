# Docker NGSolve image

## Build

### Prerequisites

Please build base [Terminator Docker image](../terminator/README.md). 

### Ubuntu Jammy

```bash
docker build --no-cache --progress=plain --build-arg VERSION=jammy -t useful-dockers/ngsolve:jammy .
```

### The Latest Ubuntu

```bash
docker build --no-cache --progress=plain --build-arg VERSION=latest -t useful-dockers/ngsolve:latest .
```

## Run

### Ubuntu Jammy

```bash
docker run -it --name ngsolve -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/ngsolve:jammy
```

### The Latest Ubuntu

```bash
docker run -it --name ngsolve -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /dev/dri:/dev/dri -v /tmp/.X11-unix:/tmp/.X11-unix:rw useful-dockers/ngsolve:latest
```
