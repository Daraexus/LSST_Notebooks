export PATH=/opt/rh/devtoolset-2/root/usr/bin:${PATH}
export LSSTSTACK=/renoir_data_00/fouchez/lsst/DM/stacks
export CFHTRAW=/renoir_data_02/lsst_data/CFHT/rawElixir/D3
export CFHTPROD=/renoir_data_02/lsst_data/CFHT/prod2
export DATADIR=/renoir_data_02/fouchez/lsst_data

export LSSTSW=$LSSTSTACK/lsstsw/2016-04-04/lsstsw
export EUPS_PATH=$LSSTSW/stack
. $LSSTSW/bin/setup.sh

setup display_ds9
setup obs_cfht jp
setup pipe_tasks jp_master
setup ip_diffim jp_master
setup astrometry_net_data CFHT_Deep

