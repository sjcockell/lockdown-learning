#!/bin/bash

# This is a variable definition
VAR='hello world'
SLEEP=5

# This is variable use
echo ${0}
echo ${1}
sleep ${2}
echo "The number of arguments is: ${#}"
echo "The list of arguments is: ${@}"
echo ${VAR}
