#!/bin/bash

ARCH=$(uname -m)
DELIM="=============================================================================================="
EXEC_SCRIPTS=(cesar_exec.py cesar_preprocess.py classify_chains.py feature_extractor.py fine_orthology_resolver.py get_contig_sizes.py predict_with_spliceai.py train_model.py)
EXEC_MODULES=(chain_bst_index.py get_names_from_bed.py)
SCRIPT_DIR="src/python/modules/"
MODULE_DIR="src/python/"
CHECK_DEPS="check_dependencies.py"
THIRD_PARTY_INSTALL="./${CHECK_DEPS} install_third_party"

while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--conda)
        SCRIPT_DIR="" + ${SCRIPT_DIR}
        MODULE_DIR="" + ${MODULE_DIR}
        THIRD_PARTY_INSTALL += " --install_spliceai"
        shift
        ;;
    -s|--searchpath)
        SEARCHPATH="$2"
        shift # past argument
        shift # past value
        ;;
    -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;
  esac
done

cd "$SRC_DIR"

## Build C code
if [ ARCH = "arm64" ]; then \
                CFLAGS="-Wall -Wextra -O2 -g -std=c99 -arch arm64"; \
        else \
                CFLAGS="-Wall -Wextra -O2 -g -std=c99"; \
        fi; \
        mkdir -p bin && gcc ${CFLAGS} -o bin/chain_filter_by_id src/c/chain_filter_by_id.c
        gcc ${CFLAGS} -fPIC -shared -o src/python/modules/chain_coords_converter_slib.so src/c/chain_coords_converter_slib.c
        gcc ${CFLAGS} -fPIC -shared -o src/python/modules/chain_bst_lib.so src/c/chain_bst_lib.c
        echo ${DELIM}

## Build CESAR
cd CESAR2.0 && make && cd ..

## Build Rust code
cd src/rust && cargo build --release && cd ../..
cd bed2gtf && cargo build --release && cd ..
cd postoga/rustools && maturin develop --release && cd ../../
echo ${DELIM}

## Download third-party binaries
mkdir -p bin && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/twoBitToFa && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/faToTwoBit && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/wigToBigWig && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigWigToWig && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigBedToBed && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bedToBigBed && \
wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/ixIxx && \
chmod +x bin/* && \
echo ${DELIM}

## Make the primary modules executable
pwd
for script in ${EXEC_SCRIPTS[@]}; do \
    chmod +x src/python/${script} ; \
done
for script in ${EXEC_MODULES[@]}; do \
    chmod +x src/python/modules/${script} ; \
done

## install third party tools
# $PYTHON ${CHECK_DEPS} install_third_party
echo "Installing IntronIC"
echo "Successfully installed IntronIC"
echo ${DELIM}
echo "Installing IqTree2"
echo "Successfully installed IqTree2"
echo ${DELIM}
echo "Installing PRANK"
echo "Successfully installed PRANK"
echo ${DELIM}

## train classification models
$PYTHON src/python/train_model.py
echo ${DELIM}

outdir=${PREFIX}/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p ${outdir}
cp -r $RECIPE_DIR/* $outdir