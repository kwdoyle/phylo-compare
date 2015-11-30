#!/bin/bash

cat $1 | sed -e 's/45.\{8\}//g' > $1_rm_oid.nre
# need to change the first 2 numbers of the ID to match that of the file in question
# this uses the sed command to search the file for "43" (or whatever the first 2 numbers
# of the object ID are) and then deletes the next 8 characters that follow it.

# Usage: ./removeOID <file.nre>
