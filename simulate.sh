#!/bin/sh

run_id=$(date +"%Y-%m-%d"="%T")

OUTDIR="out/$run_id"

mkdir -p "${OUTDIR}"

cp "input.in" "${OUTDIR}"

echo -n "Running simulation..."

./runner < "input.in"

echo "done!"

mv *.dat "${OUTDIR}"
