#!/bin/bash

cat $1 | sed -e 's/43.\{8\}//g'
# need to change the first 2 numbers of the ID to match that of the file in question

# Usage: ./removeOID <file.nre>
