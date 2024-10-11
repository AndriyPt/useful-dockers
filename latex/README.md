# Build Docker image

```bash
docker build --no-cache -t andriyp/latex:focal .
```

## Create and Start Docker container

```bash
docker run -it --name latex -e DISPLAY -e LOCAL_USER_ID=$(id -u) -v /tmp/.X11-unix:/tmp/.X11-unix:rw andriyp/latex:focal
```
~
