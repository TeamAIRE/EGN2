#!/usr/bin/env bash

# EGN2 directory
EGN2DIR=$(pwd)

# Change owner
chown -R $USER:$USER $EGN2DIR

echo "#Compile egn2"
# Compile egn2
make

# install blast and diamond
mkdir $EGN2DIR/simtools
cd $EGN2DIR/simtools/
# install blast
BLASTVERSION="2.7.1"
echo "#Install BLAST $BLASTVERSION"
wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLASTVERSION/ncbi-blast-$BLASTVERSION+-x64-linux.tar.gz
tar -xzf ncbi-blast-$BLASTVERSION+-x64-linux.tar.gz
mv $EGN2DIR/simtools/ncbi-blast-$BLASTVERSION+ $EGN2DIR/simtools/blast

# install diamond
DIAMONDVERSION="v0.9.22"
echo "#Install DIAMOND $DIAMONDVERSION"
mkdir $EGN2DIR/simtools/diamond
cd $EGN2DIR/simtools/diamond
wget -q https://github.com/bbuchfink/diamond/releases/download/$DIAMONDVERSION/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz

cd $EGN2DIR

# install exonerate
echo "#Install EXONERATE 2.2.0"
wget -q http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz
tar -xzf exonerate-2.2.0-x86_64.tar.gz
mv $EGN2DIR/exonerate-2.2.0-x86_64 $EGN2DIR/exonerate
mv $EGN2DIR/exonerate-2.2.0-x86_64.tar.gz $EGN2DIR/exonerate/


BASHRC=$HOME/.bashrc

# add to bashrc
echo "export EGN2CONFIGDIR=$EGN2DIR/egn2config/" >> $BASHRC 
echo "export EGN2BLAST=$EGN2DIR/simtools/blast/bin/" >> $BASHRC
echo "export EGN2DIAMOND=$EGN2DIR/simtools/diamond/" >> $BASHRC
echo "export EGN2EXONERATE=$EGN2DIR/exonerate/bin/" >> $BASHRC
echo "export PATH=$EGN2DIR/bin:\$PATH" >> $BASHRC

# END
