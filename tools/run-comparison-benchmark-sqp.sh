#!/bin/bash

CONSTRAINT_SOLVER="SQP"
TIMESTEP="1e-2"

IPC_ROOT=$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")
if [ -f "IPC_bin" ]; then
    IPC_BIN=$(realpath IPC_bin)
elif [ -f "${IPC_ROOT}/build/IPC_bin" ]; then
    IPC_BIN="${IPC_ROOT}/build/IPC_bin"
elif [ -f "${IPC_ROOT}/build/release/IPC_bin" ]; then
    IPC_BIN="${IPC_ROOT}/build/release/IPC_bin"
elif [ -f "${IPC_ROOT}/build/debug/IPC_bin" ]; then
    IPC_BIN="${IPC_ROOT}/build/debug/IPC_bin"
else
    echo "Unable to find IPC_bin! Try running this script in the same directory as IPC_bin."
    exit 1
fi
echo "$IPC_BIN"

# Input: collision script path relative to "${IPC_ROOT}/input"
function run_example
{
    SCRIPT_PATH="$1"
    if [ ! -f "${SCRIPT_PATH}" ]; then
        echo "Unable to find script: $SCRIPT_PATH"
        return
    fi

    MODIFIED_SCRIPT_PATH=$(realpath input.txt)
    cp ${SCRIPT_PATH} ${MODIFIED_SCRIPT_PATH}
    sed -i "" '/constraintSolver .*/d' ${MODIFIED_SCRIPT_PATH}
    printf "constraintSolver ${CONSTRAINT_SOLVER}\n$(cat ${MODIFIED_SCRIPT_PATH})" > ${MODIFIED_SCRIPT_PATH}
    sed -i "" -E "s/time (.*) (.*)/time \1 ${TIMESTEP}/" ${MODIFIED_SCRIPT_PATH}
    if [ "$CONSTRAINT_SOLVER" == "SQP" ] || [ "$CONSTRAINT_SOLVER" == "QP" ]; then
        printf "
            warmStart 1
            tol 1
            1e-3
            # QP/SQP Parameters
            constraintType graphics
            constraintOffset 1e-3
            QPSolver OSQP
            useActiveSetConvergence
            " | awk '{$1=$1;print}' >> ${MODIFIED_SCRIPT_PATH}
    fi

    SUFFIX="_$(echo $1 | sed -e 's/\//_/g')"
    SUFFIX=${SUFFIX%.*}

    ${IPC_BIN} 10 ${MODIFIED_SCRIPT_PATH} ${SUFFIX} --log 4
    rm ${MODIFIED_SCRIPT_PATH}
}

for scene in ${IPC_ROOT}/input/paperExamples/supplementB/*.txt; do
    echo $scene
    run_example ${scene}
done
