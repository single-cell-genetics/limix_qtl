#Installation

We recommed you install the limix based QTL mapping pipeline in a seperate conda enviroment.

//Doesn't work yet.
`git clone https://github.com/PMBio/hipsci_pipeline.git`


To create a limix development enviroment

`conda create -n limix_qtl python=3 anaconda`

`source activate limix_qtl`

`conda install -c anaconda pytest pytables`

`conda install -c conda-forge bgen-reader`

`pip install limix`


NB. be sure to be in a folder where you can download files to and there is no folder called limix.
NB. Some filesystems have locking disabled to write to be able to use the tool use:
	export HDF5_USE_FILE_LOCKING=FALSE
NB. please check you have installed limix  > 2.0
NB. please check your numpy version is bound to intel MKL, which makes the analyses much faster

