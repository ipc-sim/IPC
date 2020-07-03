#!/bin/bash

for f in output/*; do
    perl -i -00pe 's/mat1(?!.*mat1)/mat2/s' $f
    sed -i "" -E '/illum .*|Tf .*|Ni .*/d' $f
done
