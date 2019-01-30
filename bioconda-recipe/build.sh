Rscript -e 'install.packages(c("rmarkdown", "pheatmap"), repos="https://stat.ethz.ch/CRAN")'
echo Prefix $PREFIX
cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=$PREFIX
cd build
make
make install
