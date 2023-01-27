# docker_abe
This ABE is based on pbc library.

[How to use] \
// exec1.sh -> dabe \
// exec2.sh -> dabe with option \
// exec3.sh -> abe with option \
// exec.sh -> all

// Using docker \
docker image build -t docker_abe . \
// ex) Execute exec.sh \
docker run --rm docker_abe /pbc/exec.sh > output.txt

// Using docker-compose \
docker-compose up -d --build \
// ex) Execute exec.sh \
docker exec docker_abe-pbc-1 /pbc/exec.sh > output.txt 

(If an error occurs, try changing the newline of shell file to LF)
