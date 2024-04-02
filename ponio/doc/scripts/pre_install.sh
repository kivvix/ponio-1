#!/bin/sh

# generate ponio/runge_kutta/butcher_methods.hpp file
cmake . -B ./build -G "Ninja Multi-Config"
cmake --build ./build --config Release

# convert notebooks into rst files for documentation
NOTEBOOKS_dir="ponio/notebooks"
OUTPUT_dir="ponio/doc/source/gallery/notebook"

mkdir -p ${OUTPUT_dir}
jupyter nbconvert ${NOTEBOOKS_dir}/*.ipynb --to rst --output-dir=${OUTPUT_dir}

# convert examples folder for documentation
EXAMPLES_dir="ponio/examples"
OUTPUT_dir="ponio/doc/source/gallery/examples"

mkdir -p ${OUTPUT_dir}
cp -r ${EXAMPLES_dir}/img ${OUTPUT_dir}
pandoc ${EXAMPLES_dir}/README.md -T rst --wrap=preserve --columns=512 -o ${OUTPUT_dir}/examples.rst

# laucnh doxygen
pushd ponio/doc
doxygen
ls
ls xml
popd
