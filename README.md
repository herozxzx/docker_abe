# docker_abe
This ABE is based on pbc library.

[How to use] \
// Using docker \
docker image build -t docker_abe . \
// ex) Execute exec.sh \
docker run --rm docker_abe /pbc/test_all.sh

// Using docker-compose \
docker-compose up -d --build \
// ex) Execute exec.sh \
docker exec docker_abe-pbc-1 /pbc/test_all.sh

(If an error occurs, try changing the newline of shell file to LF)
