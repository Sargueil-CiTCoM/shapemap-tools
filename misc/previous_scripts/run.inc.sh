#!/bin/bash

# DO NOT MODIFIED THIS SCRIPTS, use source run.inc.sh in another script.

if [ -n "$WAIT_PID" ]; then
    while [ -e /proc/$WAIT_PID ]
    do
        echo "Process: $WAIT_PID is still running"
        sleep 5
    done
    echo "Process $WAIT_PID has finished"
fi
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run simple ShapeMapper pipeline on a small subset of example data

#set -e # exit on first error (if any)

# Find the parent folder of this script,
# resolving (possibly nested) symlinks
SOURCE="${BASH_SOURàCE[0]}"
while [ -h "$SOURCE" ]; do
    BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$BASE_DIR/$SOURCE"
done
BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export PATH=${BASE_DIR}:${PATH}

cd ${BASE_DIR}

if $SPLIT_ARN then
    ARN_LIST=$(cat ${ARN_FILE}| cut -f1 -d$'\t')
fi

OLDIFS=$IFS;
for i in $LIBS; do
    IFS=',';
    set -- $i;
    IFS=$OLDIFS

    TITLE="$1";
    MODIFIED_ID=$2
    UNTREATED_ID=$3
    DENATURED_ID=$4
    SEQUENCE_PATH="$5"

    MODIFIED_R1="";
    MODIFIED_R2="";
    UNTREATED_R1="";
    UNTREATED_R2="";
    DENATURED_R1="";
    DENATURED_R2="";
    
    for dfidx in $(seq ${NDATA_FOLDER}); do
        DATA_FOLD_NAME="DATA_FOLDER$dfidx"
        echo ${!DATA_FOLD_NAME}
        if [ -d "${!DATA_FOLD_NAME}" ]; then

            MODIFIED_PATH="${!DATA_FOLD_NAME}/${MODIFIED_ID}/"
            UNTREATED_PATH="${!DATA_FOLD_NAME}/${UNTREATED_ID}/"
            DENATURED_PATH="${!DATA_FOLD_NAME}/${DENATURED_ID}/"

            #echo $MODIFIED_PATH $UNTREATED_PATH $DENATURED_PATH
            if [ -d "$MODIFIED_PATH" ]; then
                TMP="${MODIFIED_PATH}*_R1_*"
                MODIFIED_R1="${MODIFIED_R1} ${TMP}";
                TMP="${MODIFIED_PATH}*_R2_*"
                MODIFIED_R2="${MODIFIED_R2} ${TMP}";
            fi
            if [ -d "$UNTREATED_PATH" ]; then
                TMP="${UNTREATED_PATH}*_R1_*"
                UNTREATED_R1="${UNTREATED_R1} ${TMP}";
                TMP="${UNTREATED_PATH}*_R2_*"
                UNTREATED_R2="${UNTREATED_R2} ${TMP}";
            fi
            if [ -d "$DENATURED_PATH" ]; then
                TMP="${DENATURED_PATH}*_R1_*"
                DENATURED_R1="${DENATURED_R1} ${TMP}";
                TMP="${DENATURED_PATH}*_R2_*"
                DENATURED_R2="${DENATURED_R2} ${TMP}";
            fi
        fi
    done;
    TARGET=$SEQUENCE_PATH;
    OUTPUTDIR=${OUTPUT_PATH}/${SHAPEMAPPEROUT_NAME}/${TITLE};
    LOG=${OUTPUTDIR}/${TITLE}_shapemapper.log;

    echo "----------STARTING ${TITLE}-------------------"
    set -x;
    shapemapper  --name "$TITLE" \
        --target ${TARGET} \
        --overwrite \
        ${SHAPEMAPPER_ARGS} \
        --modified --R1 ${MODIFIED_R1} --R2 ${MODIFIED_R2} \
        --untreated --R1 ${UNTREATED_R1} --R2 ${UNTREATED_R2} \
        --denatured --R1 ${DENATURED_R1} --R2 ${DENATURED_R2} \
        --out "${OUTPUTDIR}" \
        --log "${LOG}" 2>&1 | tee shapemapper_script.log;

    EXIT_STATUS=$?
    set +x;

    if [ ${EXIT_STATUS} -ne 0 ] ; then
        mail -s "Shapemapper on ${TITLE} FAILED" francois-xavier.lyonnet@parisdescartes.fr <<< "Shapemapper on ${TITLE} FAILED"
    fi;
    echo "----------END ${TITLE}-------------------"
done;

