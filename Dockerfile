FROM perl:5.37.8

RUN mkdir -p /opt/genepatt
COPY BWA.aln /opt/genepatt/BWA.aln
COPY BWA.bwasw /opt/genepatt/BWA.bwasw
COPY BWA.indexer /opt/genepatt/BWA.indexer
RUN cpan App::cpanminus
RUN cpanm Archive::Extract

RUN git clone https://github.com/lh3/bwa.git
RUN cd bwa; make
RUN cd /
ENV PATH=$PATH:/bwa
ENV bwa_exec=bwa
