#!/bin/bash

echo "Running: $1"

CURRENTFILE=$(ls $1/output | wc -l)

TOTALFILES=$(ls $1 | wc -l)
TOTALFILES=$(($TOTALFILES - 2))

echo "Beginning at file: $CURRENTFILE out of $TOTALFILES"

while [ $CURRENTFILE -le $TOTALFILES ]
    do
    python OOPS.py $1/$CURRENTFILE.ini
    CURRENTFILE=$(ls $1/output | wc -l)
    done

echo "All files complete"