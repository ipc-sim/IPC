#!/bin/bash

CONSTRAINT_SOLVER="IP"

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
    # Modify the copied script
    if [ "$(uname)" == "Darwin" ]; then
        sed -i "" '/constraintSolver .*/d' ${MODIFIED_SCRIPT_PATH}
        printf "constraintSolver ${CONSTRAINT_SOLVER}\n$(cat ${MODIFIED_SCRIPT_PATH})" > ${MODIFIED_SCRIPT_PATH}
    else
        sed -i '/constraintSolver .*/d' ${MODIFIED_SCRIPT_PATH}
        printf "constraintSolver ${CONSTRAINT_SOLVER}\n$(cat ${MODIFIED_SCRIPT_PATH})" > ${MODIFIED_SCRIPT_PATH}
    fi
    # printf "\n\nselfCollisionOn\n" >> ${MODIFIED_SCRIPT_PATH}
    SUFFIX="_$(echo $1 | sed -e 's/\//_/g')"
    SUFFIX=${SUFFIX%.*}

    if [ $USE_DEBUG == true ]; then
        lldb -- ${IPC_BIN} 10 ${MODIFIED_SCRIPT_PATH} ${SUFFIX} --log 3
    else
        ${IPC_BIN} 10 ${MODIFIED_SCRIPT_PATH} ${SUFFIX} --log 3
    fi
    rm ${MODIFIED_SCRIPT_PATH} # Delete modified script
}

# run_ccd_example pointTriangleCO.txt
# run_ccd_example torusCone.txt
run_ccd_example octocatPlane.txt
# run_ccd_example matIcosphereCO.txt
