#!/bin/bash

seq=$1

for (( i=0; i<${#seq}; i++ )); do
    case ${seq:$i:1} in

	A)
	    printf "T"
	    ;;

	T)
	    printf "A"
	    ;;

	G)
	    printf "C"
	    ;;

	C)
	    printf "G"
	    ;;

	*)
	    printf "N"
	    ;;
    esac
done
