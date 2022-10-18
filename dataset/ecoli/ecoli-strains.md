# Ecoli Strains selection

genomic data was downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/167/). 

The idea is to select a pair of e.coli strains in order to serve as paternal and maternal data based on genome identity.

## E.coli strains 2003-3014 and F1 E4

Average identity ***97.55*** plus several SVs (image1.png).

- Escherichia coli strain 2003-3014 (CP101292.1) paternal
- Escherichia coli strain F1 E4 (CP040307) maternal

###  whole genome comparison alignment (dnadiff)

```
                               [CP101292.1]       [CP040307]
[Sequences]
TotalSeqs                          1                    1
AlignedSeqs               1(100.00%)           1(100.00%)
UnalignedSeqs               0(0.00%)             0(0.00%)

[Bases]
TotalBases                   5840137              5496151
AlignedBases         5012832(85.83%)      4871237(88.63%)
UnalignedBases        827305(14.17%)       624914(11.37%)

[Alignments]
1-to-1                           477                  477
TotalLength                  4724919              4724595
AvgLength                    9905.49              9904.81
AvgIdentity                    97.94                97.94

M-to-M                           938                  938
TotalLength                  5408040              5407723
AvgLength                    5765.50              5765.16
AvgIdentity                    97.55                97.55

[Feature Estimates]
Breakpoints                     1874                 1874
Relocations                      131                  143
Translocations                     0                    0
Inversions                        70                   68

Insertions                       749                  561
InsertionSum                 1252866               885647
InsertionAvg                 1672.72              1578.69

TandemIns                          3                    5
TandemInsSum                     718                  660
TandemInsAvg                  239.33               132.00

[SNPs]
TotalSNPs                      87246                87246
```

### genome plot (mummerplot)
![](plots/image.png)<!-- -->


## E.coli strains BW2952 and APEC O78

Average identity **98.70** no SVs.

- Escherichia coli BW2952 (CP001396.1) paternal 
- Escherichia coli APEC O78 (CP004009) maternal


###  whole genome comparison alignment (dnadiff)

```
                               [CP001396.1]      [CP004009]
[Sequences]
TotalSeqs                          1                    1
AlignedSeqs               1(100.00%)           1(100.00%)
UnalignedSeqs               0(0.00%)             0(0.00%)

[Bases]
TotalBases                   4578159              4798435
AlignedBases         4291796(93.75%)      4296259(89.53%)
UnalignedBases         286363(6.25%)       502176(10.47%)

[Alignments]
1-to-1                           250                  250
TotalLength                  4255961              4255679
AvgLength                   17023.84             17022.72
AvgIdentity                    98.70                98.70

M-to-M                           391                  391
TotalLength                  4386918              4386528
AvgLength                   11219.74             11218.74
AvgIdentity                    98.63                98.63

[Feature Estimates]
Breakpoints                      780                  780
Relocations                       56                   78
Translocations                     0                    0
Inversions                        18                   18

Insertions                       280                  265
InsertionSum                  342307               581148
InsertionAvg                 1222.53              2193.01

TandemIns                          6                    5
TandemInsSum                    1046                  752
TandemInsAvg                  174.33               150.40

[SNPs]
TotalSNPs                      51001                51001
```


### genome plot (mummerplot)
![](plots/image2.png)<!-- -->

## E.coli strains APEC O78 and ACN001 

Average identity **99.9** no SVs 

- Escherichia coli APEC O78 (CP004009) maternal
- Escherichia coli ACN001 (CP007442.1) paternal
 
 
###  whole genome comparison alignment (dnadiff)


```
                               [CP004009.1]       [CP004009]
[Sequences]
TotalSeqs                          1                    1
AlignedSeqs               1(100.00%)           1(100.00%)
UnalignedSeqs               0(0.00%)             0(0.00%)

[Bases]
TotalBases                   4798435              4936576
AlignedBases         4672551(97.38%)      4662853(94.46%)
UnalignedBases         125884(2.62%)        273723(5.54%)

[Alignments]
1-to-1                            60                   60
TotalLength                  4664831              4664928
AvgLength                   77747.18             77748.80
AvgIdentity                    99.95                99.95

M-to-M                           125                  125
TotalLength                  4770354              4770614
AvgLength                   38162.83             38164.91
AvgIdentity                    99.90                99.90

[Feature Estimates]
Breakpoints                      248                  248
Relocations                       24                   25
Translocations                     0                    0
Inversions                        10                   10

Insertions                        98                   72
InsertionSum                  184438               321974
InsertionAvg                 1882.02              4471.86

TandemIns                          2                    4
TandemInsSum                     138                  616
TandemInsAvg                   69.00               154.00

[SNPs]
TotalSNPs                       1755                 1755
```
 
 
### genome plot (mummerplot)
![](plots/image3.png)<!-- -->


# Conclusion

The pair of E.coli strains APEC 078 and ACN001 are the ones with the highest identity (99,9%), which is similar to the one expected for paternal/maternal data of human genomes.
