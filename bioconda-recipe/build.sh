<<<<<<< HEAD
Rscript -e 'install.packages(c("rmarkdown", "pheatmap"), repos="https://stat.ethz.ch/CRAN")'
cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=$PREFIX
=======
#!/bin/bash

Rscript -e 'install.packages(c("rmarkdown", "pheatmap"), repos="https://stat.ethz.ch/CRAN")'
cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_C_COMPILER=${CC}
>>>>>>> d05c95943a9e662feb65054dc303a1df3e7e11fe
cd build
make
make install
