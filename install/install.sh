install_path="/PATH/OF/INSTALLING/TOOLS"
install_path="/project/nwr_publication/PNAS_edit/NWR_ipsc_inte/install"

mkdir -p $install_path
mkdir -p $install_path/bin
cd $install_path
if [ ! -f $install_path/bin/bedtools ] ; then
    wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
    tar -zxvf bedtools-2.29.1.tar.gz
	cd bedtools2
	make
  cp bin/* ../bin/
fi
cd $install_path
if [ ! -f $install_path/bin/samtools ] ; then
	wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
	tar -vxjf samtools-1.11.tar.bz2
	cd samtools-1.11
	make
	cp samtools ../bin/
fi
cd $install_path
if [ ! -f $install_path/bin/minimap2 ] ; then
	wget https://github.com/lh3/minimap2/archive/refs/tags/v2.28.tar.gz
  tar -zxvf v2.28.tar.gz
	cd minimap2-2.28
  make
	cp minimap2 ../bin/
fi
