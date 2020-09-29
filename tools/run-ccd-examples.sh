#!/bin/bash

LOG_LEVEL=3

IPC_ROOT=$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")
if test -f "IPC_bin"; then
    IPC_BIN=$(realpath IPC_bin)
else
    IPC_BIN="${IPC_ROOT}/build/IPC_bin"
fi
echo "Using $IPC_BIN"

if [ $# -ne 0 ] && [ $1 == "debug" ]; then
    USE_DEBUG=true
else
    USE_DEBUG=false
fi

# Input: collision script path relative to "${IPC_ROOT}/input/collision"
function run_ccd_example
{
    SCRIPT_PATH="${IPC_ROOT}/input/otherExamples/ccd/$1"
    # Copy the original script
    MODIFIED_SCRIPT_PATH=$(realpath input.txt)
    cp ${SCRIPT_PATH} ${MODIFIED_SCRIPT_PATH}

    CCD_METHOD=$(grep '^CCDMethod' ${MODIFIED_SCRIPT_PATH} | tail -1 | awk '{print $2}')
    CCD_TOLERANCE=$(grep '^CCDTolerance' ${MODIFIED_SCRIPT_PATH} | tail -1 | awk '{print $2}')

    SUFFIX="_$(echo $1 | sed -e 's/\//_/g')"
    SUFFIX="${SUFFIX%.*}_CCDMethod=${CCD_METHOD}_CCDTolerance=${CCD_TOLERANCE}"

    if [ $USE_DEBUG == true ]; then
        lldb -- ${IPC_BIN} 100 ${MODIFIED_SCRIPT_PATH} ${SUFFIX} --log ${LOG_LEVEL}
    else
        ${IPC_BIN} 100 ${MODIFIED_SCRIPT_PATH} ${SUFFIX} --log ${LOG_LEVEL}
    fi
    rm ${MODIFIED_SCRIPT_PATH} # Delete modified script
}

# run_ccd_example pointTriangleCO.txt
run_ccd_example torusCone.txt
run_ccd_example octocatPlane.txt
