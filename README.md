# docker_abe

3 models are implemented.

1. DABE (Distributed Attribute-based Encryption)
2. CP-ABE with multiple specified attribute
3. DABE with multiple specified attribute 

These ABEs are based on pbc library.


[How to use]

1. Using docker (ex) Executing test_all.sh
```bash
docker image build -t docker_abe .
docker run --rm docker_abe /pbc/test_all.sh
```

2. Using docker-compose (ex) Executing test_all.sh
```bash
docker-compose up -d --build
docker exec docker_abe-pbc-1 /pbc/test_all.sh
```
(If an error occurs, try changing the newline of shell file to LF)
