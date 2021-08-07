module load /project/itmatlab/sharedmodules/use.shared
module load salmon-v1.4.0
module load sratoolkit-2.11.0
module load R/3.6.3
module load edirect-15.3
source venv/bin/activate

pip install -e .

# The following R libraries must be installed prior to running:
# 'devtools':
# R> install.packages("devtools")
# metacycle:
# R> devtools::install_github('gangwug/MetaCycle')
# ASSIST
# R> install.packages("assist")
# tidyr:
# R> install.packages("tidyr")
