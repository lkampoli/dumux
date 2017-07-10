#!/bin/bash

CORRECT_LOCATION_FOR_DUNE_MODULES="n"
ENABLE_MPI="n"
ENABLE_DEBUG="n"
ENABLE_PARALLEL="n"
CLEANUP="n"
DOWNLOAD_ONLY="n"

TOPDIR=$(pwd)
EXTDIR=$(pwd)/external

checkLocationForDuneModules()
{
    # test for directory dune-common, dune-common-2.4 etc.
    if ! ls dune-common* &> /dev/null; then
        echo "You have to call $0 for $1 from"
        echo "the same directory in which dune-common is located."
        echo "You cannot install it in this folder."
        CORRECT_LOCATION_FOR_DUNE_MODULES="n"
        return
    fi
    CORRECT_LOCATION_FOR_DUNE_MODULES="y"
}

createExternalDirectory()
{
    if [ ! -e $EXTDIR ]; then
        mkdir -v $EXTDIR
    fi
}

installAluGrid()
{
    cd $TOPDIR

    checkLocationForDuneModules dune-alugrid
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e dune-alugrid ]; then
        git clone -b releases/2.4 https://gitlab.dune-project.org/extensions/dune-alugrid.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf dune-alugrid
        return
    fi
}

installErt()
{
    cd $TOPDIR

    checkLocationForDuneModules ert
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e ert ]; then
        git clone -b release/2017.04 https://github.com/Ensembles/ert.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf ert
        return
    fi

    # building ert
    echo "Building ert"
    cd $TOPDIR/ert
    mkdir build
    cd build
    cmake ..
    make

    # show additional information
    echo "Ert has been built in directory ert/build."
    echo "Do not change this directory otherwise opm will not find ert!"

    cd $TOPDIR
}

installFoamGrid()
{
    cd $TOPDIR

    checkLocationForDuneModules dune-foamgrid
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e dune-foamgrid ]; then
        git clone -b releases/2.4 https://gitlab.dune-project.org/extensions/dune-foamgrid.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf dune-foamgrid
        return
    fi
}

installGLPK()
{
    cd $EXTDIR
    rm -rf glpk* standalone

    if [ ! -e glpk-4.60.tar.gz ]; then
        wget http://ftp.gnu.org/gnu/glpk/glpk-4.60.tar.gz
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf glpk
        return
    fi

    mkdir glpk
    tar zxvf glpk-4.60.tar.gz --strip-components=1 -C glpk
    cd glpk

    ./configure
    make

    # show additional information
    echo "In addition, it might be necessary to set manually"
    echo "the glpk path in the CMAKE_FLAGS section of the .opts-file:"
    echo "  -DGLPK_ROOT=/path/to/glpk \\"

    cd $TOPDIR
}

installGStat()
{
    cd $EXTDIR
    rm -rf gstat* standalone

    if [ ! -e gstat.tar.gz ]; then
        wget http://gstat.org/gstat.tar.gz
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    mkdir gstat
    tar zxvf gstat.tar.gz --strip-components=1 -C gstat
    cd gstat

    sed -i 's# doc/tex/makefile##g' configure
    ./configure
    make

    if [ -e $PWD/src/gstat ]; then
        echo "Successfully installed gstat."
    fi
}

installMETIS()
{
    cd $EXTDIR

    if [ ! -e metis-5.1.0.tar.gz ]; then
        wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf metis-5.1.0
        return
    fi

    if ! test -e "metis-5.1.0"; then
        tar zxvf metis-5.1.0.tar.gz
    fi

    cd metis-5.1.0
    METISDIR=$(pwd)
    make config
    make

    cd $TOPDIR
}

installMultidomain()
{
    cd $TOPDIR

    checkLocationForDuneModules dune-multidomain
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e dune-multidomain ]; then
        git clone -b releases/2.0 git://github.com/smuething/dune-multidomain.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf dune-multidomain
        return
    fi

    cd $TOPDIR
}

installMultidomainGrid()
{
    cd $TOPDIR

    checkLocationForDuneModules dune-multidomaingrid
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e dune-multidomaingrid ]; then
        git clone -b releases/2.3 git://github.com/smuething/dune-multidomaingrid.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf dune-multidomaingrid
        return
    fi

    cd $TOPDIR
}

installNLOPT()
{
    cd $EXTDIR
    rm -rf nlopt* standalone

    if [ ! -e nlopt-2.4.2.tar.gz ]; then
        wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf nlopt
        return
    fi

    mkdir nlopt
    tar zxvf nlopt-2.4.2.tar.gz --strip-components=1 -C nlopt
    cd nlopt

    ./configure
    make

    # show additional information
    echo "In addition, it might be necessary to set manually"
    echo "the nlopt path in the CMAKE_FLAGS section of the .opts-file:"
    echo "  -DNLOPT_ROOT=/path/to/nlopt \\"

    cd $TOPDIR
}

installOPM()
{
    cd $TOPDIR

    checkLocationForDuneModules opm
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e opm-common ]; then
        git clone -b release/2017.04 https://github.com/OPM/opm-common
    fi

    if [ ! -e opm-core ]; then
        git clone -b release/2017.04 https://github.com/OPM/opm-core
    fi

    if [ ! -e opm-material ]; then
        git clone -b release/2017.04 https://github.com/OPM/opm-material
    fi

    if [ ! -e opm-parser ]; then
        git clone -b release/2017.04 https://github.com/OPM/opm-parser
    fi

    if [ ! -e opm-grid ]; then
        git clone -b release/2017.04 https://github.com/OPM/opm-grid
    fi

    if [ ! -e opm-output ]; then
        git clone -b release/2017.04 https://github.com/OPM/opm-output
    fi
    
    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf opm-common
        rm -rf opm-core
        rm -rf opm-material
        rm -rf opm-parser
        rm -rf opm-grid
        rm -rf opm-output
        return
    fi

    # apply patches
    echo "Applying patch for opm-common"
    cd $TOPDIR/opm-common
    patch -p1 < $TOPDIR/dumux/patches/opm-common-2017.04.patch

    echo "Applying patch for opm-core"
    cd $TOPDIR/opm-core
    patch -p1 < $TOPDIR/dumux/patches/opm-core-2017.04.patch

    echo "Applying patch for opm-parser"
    cd $TOPDIR/opm-parser
    patch -p1 < $TOPDIR/dumux/patches/opm-parser-2017.04.patch

    echo "Applying patch for opm-grid"
    cd $TOPDIR/opm-grid
    patch -p1 < $TOPDIR/dumux/patches/opm-grid-2017.04.patch

    # show additional information
    echo "In addition, it might be necessary to set manually some"
    echo "CMake variables in the CMAKE_FLAGS section of the .opts-file:"
    echo "  -DOPM_COMMON_ROOT=/path/to/opm-common \\"
    echo "  -Dopm-grid_PREFIX=/path/to/opm-grid \\"
    echo "  -Dopm-common_PREFIX=/path/to/opm-common \\"
    echo "  -Dopm-core_PREFIX=/path/to/opm-core \\"
    echo "  -Dopm-material_PREFIX=/path/to/opm-material \\"
    echo "  -Dopm-parser_PREFIX=/path/to/opm-parser \\"
    echo "  -Dopm-output_PREFIX=/path/to/opm-output \\"
    echo "  -DUSE_MPI=ON \\"
    echo "  -DHAVE_OPM_GRID=1 \\"

    cd $TOPDIR
}

installPDELab()
{
    cd $TOPDIR

    checkLocationForDuneModules dune-pdelab
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e dune-pdelab ]; then
        git clone -b releases/2.0 https://gitlab.dune-project.org/pdelab/dune-pdelab.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf dune-pdelab
        return
    fi

    cd $TOPDIR
}

installTypeTree()
{
    cd $TOPDIR

    checkLocationForDuneModules dune-typetree
    if test $CORRECT_LOCATION_FOR_DUNE_MODULES == "n"; then
        return
    fi

    if [ ! -e dune-typetree ]; then
        git clone -b releases/2.3 https://gitlab.dune-project.org/staging/dune-typetree.git
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf dune-typetree
        return
    fi

    cd $TOPDIR
}

installUG()
{
    cd $EXTDIR

    UG_VERSION="3.12.1"
    if [ ! -e ug-$UG_VERSION ]; then
        git clone -b v$UG_VERSION https://gitlab.dune-project.org/staging/dune-uggrid.git ug-$UG_VERSION
    fi

    if  test "$DOWNLOAD_ONLY" == "y"; then
        return
    fi

    if  test "$CLEANUP" == "y"; then
        rm -rf ug-$UG_VERSION
        return
    fi

    # Apply patch for the parallel use of UG
    cd $TOPDIR/dune-grid
    DUNE_GRID_VERSION=`git status | head -n 1 | awk '{ print $3 }'`
    if  [ "$DUNE_GRID_VERSION" == "releases/2.3.1" ] && [ "$ENABLE_PARALLEL" == "y" ]; then
        echo "Applying patch for the parallel use of UG"
        patch -p1 < $TOPDIR/dumux/patches/grid-2.3.1.patch
    fi

    cd $EXTDIR/ug-$UG_VERSION
    autoreconf -is
    OPTIM_FLAGS="-O3 -DNDEBUG -march=native -finline-functions -funroll-loops"
    # debug flags
    if test "$ENABLE_DEBUG" == "y"; then
        OPTIM_FLAGS="-O0 -g2"
    fi
    CFLAGS="$OPTIM_FLAGS"
    CXXFLAGS="$OPTIM_FLAGS -std=c++0x -fno-strict-aliasing"
    OPTS="--enable-dune --prefix=$PWD"

    if test "$ENABLE_MPI" == "y"; then
        OPTS="$OPTS --enable-parallel MPICC=$MPICXX"
    else
        OPTS="$OPTS --without-mpi"
    fi

    ./configure \
        CFLAGS="$CFLAGS" \
        CXXFLAGS="$CXXFLAGS" \
        $OPTS

    make
    make install

    sed -i "s#-DUG_DIR=.*#-DUG_DIR=$EXTDIR/ug-$UG_VERSION \\\\#g" $TOPDIR/dumux/*.opts
    cd $TOPDIR
}

usage()
{
    echo "Usage: $0 [OPTIONS] PACKAGES"
    echo ""
    echo "Where PACKAGES is one or more of the following"
    echo "  all              Install everything and the kitchen sink."
    echo "  alugrid          Download dune-alugrid."
    echo "  ert              Download and build ert."
    echo "  foamgrid         Download dune-foamgrid."
    echo "  glpk             Download and install glpk."
    echo "  gstat            Download and install gstat."
    echo "  metis            Install the METIS graph partitioner."
    echo "  multidomain      Download dune-multidomain."
    echo "  multidomaingrid  Download and patch dune-multidomaingrid."
    echo "  nlopt            Download and install nlopt."
    echo "  opm              Download opm modules required for dune-cornerpoint."
    echo "  pdelab           Download dune-pdelab."
    echo "  typetree         Download dune-typetree."
    echo "  ug               Install the UG grid library."
    echo ""
    echo "The following options are recoginzed:"
    echo "    --parallel       Enable parallelization if available."
    echo "    --debug          Compile with debugging symbols and without optimization."
    echo "    --clean          Delete all files for the given packages."
    echo "    --download       Only download the packages."
}

SOMETHING_DONE="n"
for TMP in "$@"; do
    TMP=$(echo "$TMP" | tr "[:upper:]" "[:lower:]")
    case $TMP in
        "--debug")
            ENABLE_DEBUG="y"
            ;;
        "--download")
            DOWNLOAD_ONLY="y"
            ;;
        "--parallel")
            ENABLE_PARALLEL="y"
            ENABLE_MPI="y"
            MPICC=$(which mpicc)
            MPICXX=$(which mpicxx)
            MPIF77=$(which mpif77)

            if test -f $TOPDIR'/dune-common/bin/mpi-config'; then
                MPICONFIG=$TOPDIR'/dune-common/bin/mpi-config'
            else
                echo "MPICONFIG not found!"
                return
            fi

            MPILIBS=$($MPICONFIG --libs)
            MPILIBDIR=$(echo $MPILIBS | sed "s/.*-L\([^[:blank:]]*\).*/\1/")

            # consistency check
            if test "$ENABLE_MPI" == "y" -a -z "$MPICXX"; then
                echo ""
                echo "Compiler mpicxx not found although ENABLE_MPI is set in this script!"
                echo "Please make sure that your MPI environment is set up or that you turn it off."
                echo "The shell command mpi-selector may help you to select an installed mpi-version."
                echo "Reinitilize your PATH variable after using it (e.g. logout and login again)."
                echo "Due to this error this script stops further building now."
                echo ""

                exit -1
            fi

            ;;
        "--clean")
            CLEANUP="y"
            ;;
        all)
            SOMETHING_DONE="y"
            createExternalDirectory
            installAluGrid
            installErt
            installFoamGrid
            installGLPK
            installGStat
            installMETIS
            installMultidomain
            installMultidomainGrid
            installNLOPT
            installOPM
            installPDELab
            installTypeTree
            installUG
            ;;
        alugrid|dune-alugrid)
            SOMETHING_DONE="y"
            installAluGrid
            ;;
        ert)
            SOMETHING_DONE="y"
            installErt
            ;;
        foamgrid|dune-foamgrid)
            SOMETHING_DONE="y"
            installFoamGrid
            ;;
        glpk)
            SOMETHING_DONE="y"
            createExternalDirectory
            installGLPK
            ;;
        gstat)
            SOMETHING_DONE="y"
            createExternalDirectory
            installGStat
            ;;
        metis)
            SOMETHING_DONE="y"
            createExternalDirectory
            installMETIS
            ;;
        multidomain|dune-multidomain)
            SOMETHING_DONE="y"
            installMultidomain
            ;;
        multidomaingrid|dune-multidomaingrid)
            SOMETHING_DONE="y"
            installMultidomainGrid
            ;;
        nlopt)
            SOMETHING_DONE="y"
            createExternalDirectory
            installNLOPT
            ;;
        opm)
            SOMETHING_DONE="y"
            installOPM
            ;;
        pdelab|dune-pdelab)
            SOMETHING_DONE="y"
            installPDELab
            ;;
        typetree|dune-typetree)
            SOMETHING_DONE="y"
            installTypeTree
            ;;
        ug)
            SOMETHING_DONE="y"
            createExternalDirectory
            installUG
            ;;
        *)
            usage
            exit 1
    esac
    cd $TOPDIR
done

if test "$SOMETHING_DONE" != "y"; then
    usage
    exit 1;
fi
