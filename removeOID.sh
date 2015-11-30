#!/bin/bash

cat $1 | sed -e "s/$2.\{8\}//g" > tree_rm_oid.nre
# need to change the first 2 numbers of the ID to match that of the file in question
# this uses the sed command to search the file for "43" (or whatever the first 2 numbers
# of the object ID are) and then deletes the next 8 characters that follow it.

rm $1
mv tree_rm_oid.nre $1
# these two steps remove the old tree file and then rename the new file without object IDs
# to the old tree file's name.

# Usage: ./removeOID <file.nre> <first two numbers in object ID>
