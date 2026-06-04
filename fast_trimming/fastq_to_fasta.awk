#!/bin/awk -f

{
    # Header of the sequence
    if((NR-1)%4==0)
    {
	# Replace @ by the fasta format header character >
	sub("@",">",$0);
	# Replace all whitespaces by underscores
	gsub(" ","_",$0);
	print
    }
    # Sequence
    if((NR-2)%4==0)
    {
	print
    }
}
