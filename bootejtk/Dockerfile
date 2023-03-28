FROM python:2

ENV R=R-4.2.2

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update                                                                                              \
 && apt-get install -y gfortran                                                                                 \
 && apt-get install -y libbz2-dev libcurl4-openssl-dev libfribidi-dev libharfbuzz-dev                           \
                       libgmp-dev liblzma-dev libpcre2-dev libreadline-dev libxt-dev

RUN wget -q -O - https://cran.r-project.org/src/base/$(echo $R | cut -c 1-3)/$R.tar.gz | tar xzf -              \
 && cd $R                                                                                                       \
 && ./configure --enable-R-shlib --disable-prebuilt-html --without-x                                            \
 && make -j $(nproc)                                                                                            \
 && make install                                                                                                \
 && rm -rf /$R

ENV R_HOME=/usr/local/lib/R

RUN echo "options(repos=c(CRAN='https://cloud.r-project.org'), Ncpus=$(nproc))" >> $R_HOME/etc/Rprofile.site

RUN R -e "install.packages(c('BiocManager', 'devtools', 'reshape', 'reshape2'))"
RUN R -e "BiocManager::install(c('limma', 'qvalue'))"
RUN R -e "devtools::install_github('stephens999/ashr',build_vignettes=FALSE)"
RUN R -e "devtools::install_github('mengyin/vashr',build_vignettes=FALSE)"

RUN pip install Cython matplotlib statsmodels==0.10

ENV BOOTEJTK_DIRECTORY "/BooteJTK"

COPY BooteJTK $BOOTEJTK_DIRECTORY/
RUN cd $BOOTEJTK_DIRECTORY; python setup.py build_ext --inplace

ENV EJTK_DIRECTORY "/eJTK"

COPY empirical-JTK_CYCLE-with-asymmetry $EJTK_DIRECTORY/
RUN cd $EJTK_DIRECTORY/bin; python setup.py build_ext --inplace
