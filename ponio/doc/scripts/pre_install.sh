#!/bin/bash
set -euxo pipefail

# generate ponio/runge_kutta/butcher_methods.hpp file
cmake . -B ./build -G "Ninja Multi-Config" -DBUILD_DOC=ON -DBUILD_OFFLINE=OFF
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

# launch examples in doc and plot figure
echo "make run"
export CONDA_PREFIX=$(dirname $(which doxygen))/..
make -C ponio/doc/source/_static/cpp run

pushd ponio/doc/source/_static/cpp
python figures.py
popd

# download PlantUML
URL_PlantUML="https://github.com/plantuml/plantuml/releases/download/v1.2024.7/plantuml-1.2024.7.jar"
wget $URL_PlantUML -O ponio/doc/source/plantuml.jar

# launch compare plot
pushd ponio/doc/source/compare

# exclude diffeq from compare because it takes to much memory to install on RTD builder
compare=("ascent" "gsl" "odeint" "petsc" "ponio" "scipy")
for dir in "${compare[@]}"; do
  make -C ${dir}
done

python figures.py

popd

# launch doxygen
pushd ponio/doc
doxygen
ls
ls xml
popd
