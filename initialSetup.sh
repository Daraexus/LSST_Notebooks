#!/bin/bash

<<<<<<< HEAD
#source ../v9_2_copy/loadLSST.sh
source /renoir_data_00/fouchez/lsst/DM/stacks/v10.1/loadLSST.bash
#eups declare -force -m none -r ../my_packages/obs_cfht obs_cfht jp_obs_cfht
#setup obs_cfht jp_obs_cfht
setup obs_cfht DM-1593
setup pipe_tasks

cd data/CFHTLS
mkdir input
echo 'lsst.obs.cfht.MegacamMapper' > input/_mapper

cd ../..

#eups declare -force -m none -r data/astrometry_net_data/CFHT-Deep astrometry_net_data CFHT-Deep
=======
source /sps/lsst/data/jreyes/stacks/stack_v10_0/loadLSST.bash
export PYTHONPATH=$PYTHONPATH:/afs/in2p3.fr/home/j/jreyes/my_modules/
export PATH=$PATH:/sps/lsst/Library/ds9.new/

setup obs_cfht
setup pipe_tasks
>>>>>>> 09895cfeee029fbc9ab30ca248afc8572dfba392
setup astrometry_net_data CFHT_Deep
