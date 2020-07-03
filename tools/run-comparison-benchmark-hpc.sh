#!/bin/bash
# TODO: Update sed commands to work on macOS too.

IPC_ROOT=$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")
IPC_BIN="${IPC_ROOT}/build/IPC_bin"

TIME=$(date "+%F-%T")
TIME=$(echo "${TIME//:/-}")
RESULTS_DIR="${IPC_ROOT}/output/SQP-sweep-$TIME"

EXAMPLE_SET="$1"

if [ -z "$EXAMPLE_SET" ]; then
    echo "Must select an example set (e.g. $0 0)"
    exit 1
elif [ $EXAMPLE_SET == "test" ]; then
    RESULTS_DIR="$RESULTS_DIR-test"
elif [ $EXAMPLE_SET == "large" ]; then
    # Large-ish scenes
    RESULTS_DIR="$RESULTS_DIR-large"
elif [ $EXAMPLE_SET == "memory" ]; then
    RESULTS_DIR="$RESULTS_DIR-memory"
elif [ $EXAMPLE_SET == "chain" ]; then
    RESULTS_DIR="$RESULTS_DIR-chain"
fi

mkdir -p "$RESULTS_DIR"
GIT_BRANCH=$(git branch | grep \* | cut -d ' ' -f2)
GIT_SHA=$(git rev-parse HEAD)
printf "Time: $TIME\nGit branch: $GIT_BRANCH\nGit SHA: $GIT_SHA\n" \
    > "${RESULTS_DIR}/${GIT_BRANCH}-${GIT_SHA}.txt"

if [ $EXAMPLE_SET == "large" ]; then
    TIME_LIMIT="24:00:00"
else
    TIME_LIMIT="4:00:00"
fi

# Input: collision script path relative to "${IPC_ROOT}/input/collision"
function run_example
{
    script_path="${IPC_ROOT}/input/collision/$1"
    example_name=$(echo $1 | sed -e 's/\//_/g')
    example_name=${example_name%.*}
    for constraintSolver in "SQP" "QP" "IP"
    do
        for timestep in 1e-2 1e-3 1e-4 1e-5 # 1e-6
        do
            if [ "$constraintSolver" != "IP" ]; then
                for constraintOffset in 1e-2 1e-3 1e-4 1e-5 # 1e-6
                do
                    for constraintType in "graphics" "volume" "Verschoor"
                    do
                        for QPSolver in "Gurobi" # "OSQP"
                        do
                            for activeSetConvergencePolicy in "useActiveSetConvergence" # "noActiveSetConvergence"
                            do
                                output_path="${RESULTS_DIR}/constraintSolver=${constraintSolver}/timestep=${timestep}/constraintOffset=${constraintOffset}/constraintType=${constraintType}/QPSolver=${QPSolver}"
                                mkdir -p ${output_path}
                                cd ${output_path}

                                # Create a directory for the logs
                                mkdir -p "logs"
                                log_path="$(realpath logs)/${example_name}__${constraintSolver}_${timestep}_${constraintOffset}_${constraintType}_${QPSolver}_${activeSetConvergencePolicy}"

                                new_script_path="${output_path}/input/collision/$1"
                                mkdir -p $(dirname ${new_script_path})
                                cp ${script_path} ${new_script_path}
                                # Remove lines for variables we are using
                                sed -i -E "s/time (.*) (.*)/time \1 ${timestep}/" ${new_script_path}
                                sed -i '/constraintSolver .*/d' ${new_script_path}
                                sed -i '/constraintType .*/d' ${new_script_path}
                                sed -i '/constraintOffset .*/d' ${new_script_path}
                                sed -i '/QPSolver .*/d' ${new_script_path}
                                sed -i '/warmStart .*/d' ${new_script_path}
                                sed -i '/useActiveSetConvergence/d' ${new_script_path}
                                sed -i '/noActiveSetConvergence/d' ${new_script_path}
                                # Add lines for variables we are using
                                printf "constraintSolver ${constraintSolver}\n$(cat ${new_script_path})" > ${new_script_path}
                                printf "
                                    warmStart 1
                                    tol 1
                                    1e-3
                                    # Sweep parameters
                                    constraintOffset ${constraintOffset}
                                    constraintType ${constraintType}
                                    QPSolver ${QPSolver}
                                    ${activeSetConvergencePolicy}
                                    " | awk '{$1=$1;print}' >> ${new_script_path}

                                sbatch --output="${log_path}.out.txt" --error="${log_path}.err.txt" --time="$TIME_LIMIT" "$IPC_ROOT/tools/hpc-job.sh" ${IPC_BIN} ${new_script_path} "${example_name}"
                            done
                        done
                    done
                done
            else
                output_path="${RESULTS_DIR}/constraintSolver=${constraintSolver}/timestep=${timestep}/"
                mkdir -p ${output_path}
                cd ${output_path}

                # Create a directory for the logs
                mkdir -p "logs"
                log_path="$(realpath logs)/${example_name}__${constraintSolver}_${timestep}"

                new_script_path="${output_path}/input/collision/$1"
                mkdir -p $(dirname ${new_script_path})
                cp ${script_path} ${new_script_path}
                # Remove lines for variables we are using
                sed -i -E "s/time (.*) (.*)/time \1 ${timestep}/" ${new_script_path}
                sed -i '/constraintSolver .*/d' ${new_script_path}
                # Add lines for variables we are using
                printf "constraintSolver ${constraintSolver}\n$(cat ${new_script_path})" > ${new_script_path}

                sbatch --output="${log_path}.out.txt" --error="${log_path}.err.txt" --time="$TIME_LIMIT" "$IPC_ROOT/tools/hpc-job.sh" ${IPC_BIN} ${new_script_path} "${example_name}"
            fi
        done
    done
    cd ${calling_dir}
}

for scene in ${IPC_ROOT}/input/paperExamples/supplementB/*.txt; do
    run_example ${scene}
done
