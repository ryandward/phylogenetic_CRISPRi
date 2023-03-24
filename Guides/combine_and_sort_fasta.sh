#! /bin/bash

awk '
BEGIN {
    PROCINFO["sorted_in"] = "@ind_str_asc";
    FS = "\n";
    RS = ">";
}
NF {
    seq = $2;
    gsub("\n", "", seq);
    fname = FILENAME;
    gsub("\\..*$", "", fname);
    gsub("_.*$", "", fname);
    if (seq in seqs) {
        seqs[seq] = seqs[seq] "," fname;
    } else {
        seqs[seq] = fname;
    }
}
END {
    
    for (seq in seqs) {
        print ">" seqs[seq] "_" seq;
        print seq;
    }
}' *fasta
