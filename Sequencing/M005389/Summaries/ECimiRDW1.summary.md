``` logtalk
[14:30:30] Initializing heuristic barcode counting...                                                                                                                                                                      heuristicount.py:329
           Reading barcodes...                                                                                                                                                                                             heuristicount.py:337
           Sampling initial reads for orientation...                                                                                                                                                                       heuristicount.py:346
           Determining orientation of reads...                                                                                                                                                                             heuristicount.py:349
[14:30:32] Sampled 60,210 unique reads and found 27,476 forward and 32,528 reverse valid matches...                                                                                                                        heuristicount.py:353
           Swapping orientation...                                                                                                                                                                                         heuristicount.py:354
           Identifying forward flanking sequences...                                                                                                                                                                       heuristicount.py:359
           Identifying reverse flanking sequences...                                                                                                                                                                       heuristicount.py:367
[14:30:33] Executing high-throughput read analysis...                                                                                                                                                                      heuristicount.py:403
[14:31:31] Collating results...                                                                                                                                                                                            heuristicount.py:413
```

``` logtalk
                                                                                                                                                                                                                           heuristicount.py:502
                           heuristicount.py                               Summary                                                                                                                                                              
            ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━                                                                                                                                                             
                             Input & Config                                                                                                                                                                                                    
                                   Barcodes                        A_E_coli.fasta                                                                                                                                                              
                              Forward Reads   ECimiRDW1_S52_L005_R2_001.reads.zst                                                                                                                                                              
                              Reverse Reads   ECimiRDW1_S52_L005_R1_001.reads.zst                                                                                                                                                              
                                    Threads                                     6                                                                                                                                                              
                           Operating System                                 Linux                                                                                                                                                              
                                                                                                                                                                                                                                               
                                 Heuristics                                                                                                                                                                                                    
                             Barcode Length                                    20                                                                                                                                                              
                             Forward Offset                                   100                                                                                                                                                              
                             Reverse Offset                                    82                                                                                                                                                              
                             Forward Flanks                           TAGT...GTTT                                                                                                                                                              
                             Reverse Flanks                           AAAC...ACTA                                                                                                                                                              
                                                                                                                                                                                                                                               
                    Barcode Alignment Stats                                                                                                                                                                                                    
                          Declared Barcodes                                22,579                                                                                                                                                              
                   Seen Documented Barcodes                                22,574                                                                                                                                                              
                 Unseen Documented Barcodes                                     5                                                                                                                                                              
                Found Undocumented Barcodes                                61,471                                                                                                                                                              
                                Total Reads                            40,455,498                                                                                                                                                              
                   Documented Barcode Reads                            30,283,476                                                                                                                                                              
                 Undocumented Barcode Reads                               129,664                                                                                                                                                              
                        Documented Fraction                                 0.749                                                                                                                                                              
                      Undocumented Fraction                                 0.003                                                                                                                                                              
                                                                                                                                                                                                                                               
                  Top 5 Documented Barcodes                                                                                                                                                                                                    
                       CACCAATCTCAATGATCTTG                                 5,976                                                                                                                                                              
                       ATCGTCTGCTACCACACGCT                                 5,363                                                                                                                                                              
                       CGCCACTTTATGCATCGTAT                                 5,233                                                                                                                                                              
                       ATCGTCGACATATTGCTGCC                                 4,745                                                                                                                                                              
                       AGCCTTCGATTACGCCGCCC                                 4,489                                                                                                                                                              
                                                                                                                                                                                                                                               
                Top 5 Undocumented Barcodes                                                                                                                                                                                                    
                      TATCAATTCCGCCACGCCGA*                                   238                                                                                                                                                              
                      AGCGATCCCGGTGGGGCAGT*                                   223                                                                                                                                                              
                      GGCACAAACTCTGCCTGTTG*                                   207                                                                                                                                                              
                      TACGGGCCGGGAAATTGACA*                                   201                                                                                                                                                              
                      AATAGCTTAGTGCCCCATCT*                                   196                                                                                                                                                              
                                                                                                                                                                                                                                               
                            Finished at 2023-12-01 14:31:31.421461                                                                                                                                                                             
```