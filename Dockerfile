# Dockerfile for FastqPuri

# Generate minimal linux system containing the needed packages>
FROM debian:buster
RUN apt-get update -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y r-base make cmake pandoc git vim
RUN Rscript -e 'install.packages(c("rmarkdown", "pheatmap"), repos="https://cran.uni-muenster.de")'

# compile and install FastqPuri
RUN cd /home && git clone https://github.com/jengelmann/FastqPuri
RUN cd /home/FastqPuri && cmake -H. -Bbuild/ -DRSCRIPT=/usr/bin/Rscript
RUN cd /home/FastqPuri/build && make && make install

# Start command
WORKDIR /tmp
CMD ["bash"]

# Suggestions for a docker usage from the working directory:
# Interactive: ~> docker run -v $PWD:/tmp -it fastqpuri
# As pipeline: ~> docker run -v $PWD:/tmp fastqpuri ./pipeline.sh
