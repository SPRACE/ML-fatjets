#!/bin/sh
#SBATCH --time=1

#export TMPDIR=/tmp/trtomei
#export TMPDIR=`mktemp -d`
export TMPDIR=$PWD
cd $TMPDIR

/bin/hostname
echo $PWD
/bin/ls

module load root

### Fastjet
export FASTJETDIR=$TMPDIR/fjdir
mkdir -p $FASTJETDIR
wget http://fastjet.fr/repo/fastjet-3.2.1.tar.gz
tar -xzf fastjet-3.2.1.tar.gz
cd fastjet-3.2.1
./configure --prefix=$FASTJETDIR
make -j 4
make install
export PATH=$FASTJETDIR/bin:$PATH

### Pythia
export PYTHIADIR=$TMPDIR/pythia8
mkdir -p $PYTHIADIR
cd $PYTHIADIR
wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8223.tgz
tar -xzf pythia8223.tgz
cd pythia8223
export PYTHIA8DIR=$PWD
./configure --with-fastjet3=$FASTJETDIR --with-gzip
make -j 4

### Our repo
cd $TMPDIR
export GIT_SSL_NO_VERIFY=1
git clone https://git.ncc.unesp.br/trtomei/ML-fatjets.git
cd ML-fatjets
source compileGenerateJetImage.sh
