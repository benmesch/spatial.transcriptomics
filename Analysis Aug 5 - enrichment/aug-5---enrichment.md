For legacy 5330 project, compared group 1 to group 0 (cancer vs non
cancer sample log2( 10^gene means / 10^other gene means )

Then select only genes with absolute value &gt; 2, and sort before
passing into the cluster profiler function.

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 

    ## 

    ## clusterProfiler v3.16.0  For help: https://guangchuangyu.github.io/software/clusterProfiler
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

Load expression values (transformed and normalized, limited to 16k genes
from 31k originally.

    expression.vals <- readRDS('expression.vals.rds')
    g4.border.nodes <- readRDS('g4.border.nodes.rds')
    g8.border.nodes <- readRDS('g8.border.nodes.rds')
    gene.list <- readRDS('gene.list.rds')
    SCT <- readRDS('SCT.rds')

Cell types data (`expressoin.vals`) holds composition of 23 cell types
for all 1,072 spots. Transformed counts data (`SCT`) has normalized gene
counts of 16,017 genes for every spot.

    n <- dim(expression.vals)[1]
    dim(expression.vals)

    ## [1] 1072   23

    dim(SCT)

    ## [1] 16017  1072

Split the 1,072 samples into two groups “g0” and “g1”. Log10 transform
(this was done for 5330 project).

    print (c('border samples:', length(g4.border.nodes)))

    ## [1] "border samples:" "651"

    print (c('non border samples:', n-length(g4.border.nodes)))

    ## [1] "non border samples:" "421"

Take gene means for each of the 2 different groups of spots, and plot
the log-transformed ratio.

    g0.means <- rowMeans(SCT[,-g4.border.nodes])
    g1.means <- rowMeans(SCT[,g4.border.nodes])
    #create gene list
    gene.list.1.0 <- 10^g1.means / 10^g0.means #raw values are all log10 transormed. undo that.
    gene.list.1.0 <- log(gene.list.1.0, base=2)
    plot(sort(gene.list.1.0))
    abline(0,0)

![](aug-5---enrichment_files/figure-markdown_strict/unnamed-chunk-5-1.png)

Use bioconductor package to convert gene names into Entrez ids. Remove
duplicate Entrez located for 3 of the genes:

    annot <- AnnotationDbi::select(org.Hs.eg.db, keys=toupper(gene.list), columns=c("SYMBOL","ENTREZID")
                    , keytype="SYMBOL", multiVals="first")

    ## 'select()' returned 1:many mapping between keys and columns

    which(table((annot[,1]))>1)

    ## MEMO1  MMD2   TEC 
    ##  8609  8791 13869

    length(unique(annot[,1]))

    ## [1] 16017

    sum(annot[,1] %in% c("MEMO1","MMD2","TEC")) #6 entrez ids for these 3 gene symbols

    ## [1] 6

    #which(annot[,1] %in% c("MEMO1","MMD2","TEC"))
    annot[which(annot[,1] %in% c("MEMO1","MMD2","TEC")),]

    ##       SYMBOL  ENTREZID
    ## 4847     TEC      7006
    ## 4848     TEC 100124696
    ## 5442    MMD2    221938
    ## 5443    MMD2 100505381
    ## 14991  MEMO1      7795
    ## 14992  MEMO1     51072

Remove the duplicate rows, to stay at 16,017 gene rows.

    annot <- annot[-c(4848,5443,14992),]
    head(annot)

    ##   SYMBOL ENTREZID
    ## 1   XKR4   114786
    ## 2  SOX17    64321
    ## 3 MRPL15    29088
    ## 4 LYPLA1    10434
    ## 5  TCEA1     6917
    ## 6  RGS20     8601

Alternate plot of log-ratios for all 16k genes:

    hist((gene.list.1.0))

![](aug-5---enrichment_files/figure-markdown_strict/unnamed-chunk-8-1.png)
Save differentiall expressed genes. Limit to only genes that mapped to
an Entrez ID & had log-ratio &gt; 0.3

    differential.genes <- annot[!is.na(annot[,2]) & abs(gene.list.1.0) > 0.3,2]
    differential.genes

    ##  [1] "8476"   "149111" "3776"   "57512"  "51490"  "2934"   "25891"  "80025" 
    ##  [9] "4826"   "11255"  "7881"   "26289"  "8577"   "6456"   "8672"   "4681"  
    ## [17] "54997"  "84935"  "6863"   "644150" "6622"   "23233"  "2904"   "5582"  
    ## [25] "284339" "4793"   "27120"  "2558"   "2915"   "10368"  "1339"   "55117" 
    ## [33] "93058"  "729991" "51090"  "9467"   "2137"   "79873"  "2675"   "3609"  
    ## [41] "10518"  "2941"   "25907"  "22908"  "57659"  "3859"   "2670"   "1394"  
    ## [49] "2873"   "140767" "10814"  "6695"   "4208"   "10397"  "8214"   "257068"
    ## [57] "116448" "5121"   "10846"  "534"    "374868" "5977"   "119032" "51703"

    write.table(differential.genes, file="differential_entrez_logratio_abs_greater_point3.txt"
                , row.names = FALSE, col.names = FALSE, quote = FALSE)

save separate lists of differentially up and down expressed genes:

    differential.genes <- annot[!is.na(annot[,2]) & gene.list.1.0 > 0.3,2]
    length(differential.genes)

    ## [1] 41

    write.table(differential.genes, file="differential_entrez_logratio_greater_point3.txt"
                , row.names = FALSE, col.names = FALSE, quote = FALSE)

    differential.genes <- annot[!is.na(annot[,2]) & gene.list.1.0 < -0.3,2]
    length(differential.genes)

    ## [1] 23

    write.table(differential.genes, file="differential_entrez_logratio_less_point3.txt"
                , row.names = FALSE, col.names = FALSE, quote = FALSE)

### Differential genes (n=64)

    ggo <- groupGO(gene     = differential.genes,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   level    = 2,
                   readable = TRUE)
    ggo@result = ggo@result[order(ggo@result[,"Count"],decreasing=TRUE),]
    ggo@result

    ##                    ID                      Description Count GeneRatio
    ## GO:0005488 GO:0005488                          binding    21     21/23
    ## GO:0003824 GO:0003824               catalytic activity     6      6/23
    ## GO:0098772 GO:0098772     molecular function regulator     4      4/23
    ## GO:0005198 GO:0005198     structural molecule activity     2      2/23
    ## GO:0005215 GO:0005215             transporter activity     2      2/23
    ## GO:0060089 GO:0060089    molecular transducer activity     2      2/23
    ## GO:0140110 GO:0140110 transcription regulator activity     2      2/23
    ## GO:0045182 GO:0045182   translation regulator activity     1      1/23
    ## GO:0016209 GO:0016209             antioxidant activity     0      0/23
    ## GO:0031386 GO:0031386                      protein tag     0      0/23
    ## GO:0038024 GO:0038024          cargo receptor activity     0      0/23
    ## GO:0044183 GO:0044183        protein folding chaperone     0      0/23
    ## GO:0045735 GO:0045735      nutrient reservoir activity     0      0/23
    ## GO:0090729 GO:0090729                   toxin activity     0      0/23
    ## GO:0140104 GO:0140104       molecular carrier activity     0      0/23
    ## GO:0140299 GO:0140299   small molecule sensor activity     0      0/23
    ## GO:0140313 GO:0140313  molecular sequestering activity     0      0/23
    ##                                                                                                                                geneID
    ## GO:0005488 GSN/PANK2/EIF4G3/NBL1/TESC/TAC1/WIPF3/EXOC6B/GRIN2B/PRKCG/GRM5/CACNG3/PLLP/GSTA4/SACM1L/GFAP/CPLX2/NDRG1/OLIG1/PDE10A/DPF2
    ## GO:0003824                                                                                       PANK2/PRKCG/GSTA4/SACM1L/PDE10A/DPF2
    ## GO:0098772                                                                                                      NBL1/TESC/GRM5/CACNG3
    ## GO:0005198                                                                                                                  PLLP/GFAP
    ## GO:0005215                                                                                                              GRIN2B/CACNG3
    ## GO:0060089                                                                                                                GRIN2B/GRM5
    ## GO:0140110                                                                                                                 OLIG1/DPF2
    ## GO:0045182                                                                                                                     EIF4G3
    ## GO:0016209                                                                                                                           
    ## GO:0031386                                                                                                                           
    ## GO:0038024                                                                                                                           
    ## GO:0044183                                                                                                                           
    ## GO:0045735                                                                                                                           
    ## GO:0090729                                                                                                                           
    ## GO:0140104                                                                                                                           
    ## GO:0140299                                                                                                                           
    ## GO:0140313

    barplot(ggo,showCategory=8)

![](aug-5---enrichment_files/figure-markdown_strict/unnamed-chunk-11-1.png)

### Up expressed genes

    differential.genes <- annot[!is.na(annot[,2]) & gene.list.1.0 > 0.3,2]
    ggo <- groupGO(gene     = differential.genes,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   level    = 2,
                   readable = TRUE)
    ggo@result = ggo@result[order(ggo@result[,"Count"],decreasing=TRUE),]
    ggo@result

    ##                    ID                      Description Count GeneRatio
    ## GO:0005488 GO:0005488                          binding    28     28/41
    ## GO:0003824 GO:0003824               catalytic activity    10     10/41
    ## GO:0098772 GO:0098772     molecular function regulator     8      8/41
    ## GO:0005215 GO:0005215             transporter activity     6      6/41
    ## GO:0060089 GO:0060089    molecular transducer activity     5      5/41
    ## GO:0140110 GO:0140110 transcription regulator activity     3      3/41
    ## GO:0005198 GO:0005198     structural molecule activity     1      1/41
    ## GO:0016209 GO:0016209             antioxidant activity     0      0/41
    ## GO:0031386 GO:0031386                      protein tag     0      0/41
    ## GO:0038024 GO:0038024          cargo receptor activity     0      0/41
    ## GO:0044183 GO:0044183        protein folding chaperone     0      0/41
    ## GO:0045182 GO:0045182   translation regulator activity     0      0/41
    ## GO:0045735 GO:0045735      nutrient reservoir activity     0      0/41
    ## GO:0090729 GO:0090729                   toxin activity     0      0/41
    ## GO:0140104 GO:0140104       molecular carrier activity     0      0/41
    ## GO:0140299 GO:0140299   small molecule sensor activity     0      0/41
    ## GO:0140313 GO:0140313  molecular sequestering activity     0      0/41
    ##                                                                                                                                                                                      geneID
    ## GO:0005488 CDC42BPA/SPOUT1/PAMR1/KCNAB1/AK5/SH3GL2/SNCA/NFKBIB/DKKL1/GABRA5/COQ10A/BORCS8/SH3BP5/EXTL3/NUDT18/ILF3/CIB2/TMEM158/ZBTB4/CRHR1/GPS1/NRSN1/SPOCK1/MEF2C/PCP4/ATP9B/BORCS7/ACSL5
    ## GO:0003824                                                                                                                 CDC42BPA/SPOUT1/KCNAB1/AK5/SNCA/COX6A2/EXTL3/NUDT18/PLCXD2/ACSL5
    ## GO:0098772                                                                                                                                CNIH3/KCNAB1/SNCA/DKKL1/COX6A2/SH3BP5/GPS1/SPOCK1
    ## GO:0005215                                                                                                                                      KCNK2/KCNAB1/GABRA5/COX6A2/SLC6A15/ATP6V1G2
    ## GO:0060089                                                                                                                                                   GPR158/HRH3/GABRA5/GFRA2/CRHR1
    ## GO:0140110                                                                                                                                                               NFKBIB/ZBTB4/MEF2C
    ## GO:0005198                                                                                                                                                                            KRT12
    ## GO:0016209                                                                                                                                                                                 
    ## GO:0031386                                                                                                                                                                                 
    ## GO:0038024                                                                                                                                                                                 
    ## GO:0044183                                                                                                                                                                                 
    ## GO:0045182                                                                                                                                                                                 
    ## GO:0045735                                                                                                                                                                                 
    ## GO:0090729                                                                                                                                                                                 
    ## GO:0140104                                                                                                                                                                                 
    ## GO:0140299                                                                                                                                                                                 
    ## GO:0140313

### Down expressed genes

    differential.genes <- annot[!is.na(annot[,2]) & gene.list.1.0 < -0.3,2]
    ggo <- groupGO(gene     = differential.genes,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   level    = 3,
                   readable = TRUE)
    ggo@result = ggo@result[order(ggo@result[,"Count"],decreasing=TRUE),]
    ggo@result

    ##                    ID
    ## GO:0005515 GO:0005515
    ## GO:0043167 GO:0043167
    ## GO:0097159 GO:0097159
    ## GO:1901363 GO:1901363
    ## GO:0016740 GO:0016740
    ## GO:0036094 GO:0036094
    ## GO:0008144 GO:0008144
    ## GO:0097367 GO:0097367
    ## GO:0016787 GO:0016787
    ## GO:0022857 GO:0022857
    ## GO:0044877 GO:0044877
    ## GO:0038023 GO:0038023
    ## GO:0030234 GO:0030234
    ## GO:0003700 GO:0003700
    ## GO:0140096 GO:0140096
    ## GO:0005200 GO:0005200
    ## GO:0019911 GO:0019911
    ## GO:0033218 GO:0033218
    ## GO:0042165 GO:0042165
    ## GO:0090079 GO:0090079
    ## GO:0016247 GO:0016247
    ## GO:0030545 GO:0030545
    ## GO:0003712 GO:0003712
    ## GO:0004133 GO:0004133
    ## GO:0009975 GO:0009975
    ## GO:0016491 GO:0016491
    ## GO:0016829 GO:0016829
    ## GO:0016853 GO:0016853
    ## GO:0016874 GO:0016874
    ## GO:0032451 GO:0032451
    ## GO:0061783 GO:0061783
    ## GO:0140097 GO:0140097
    ## GO:0140098 GO:0140098
    ## GO:0003735 GO:0003735
    ## GO:0005199 GO:0005199
    ## GO:0005201 GO:0005201
    ## GO:0005212 GO:0005212
    ## GO:0005213 GO:0005213
    ## GO:0008147 GO:0008147
    ## GO:0008307 GO:0008307
    ## GO:0008316 GO:0008316
    ## GO:0016490 GO:0016490
    ## GO:0017056 GO:0017056
    ## GO:0030280 GO:0030280
    ## GO:0030281 GO:0030281
    ## GO:0030527 GO:0030527
    ## GO:0035804 GO:0035804
    ## GO:0039660 GO:0039660
    ## GO:0042302 GO:0042302
    ## GO:0043886 GO:0043886
    ## GO:0097099 GO:0097099
    ## GO:0097493 GO:0097493
    ## GO:0098918 GO:0098918
    ## GO:0005319 GO:0005319
    ## GO:0005344 GO:0005344
    ## GO:0140142 GO:0140142
    ## GO:0140318 GO:0140318
    ## GO:0000035 GO:0000035
    ## GO:0003682 GO:0003682
    ## GO:0003823 GO:0003823
    ## GO:0005549 GO:0005549
    ## GO:0008289 GO:0008289
    ## GO:0008430 GO:0008430
    ## GO:0015643 GO:0015643
    ## GO:0019808 GO:0019808
    ## GO:0030246 GO:0030246
    ## GO:0031409 GO:0031409
    ## GO:0033226 GO:0033226
    ## GO:0035731 GO:0035731
    ## GO:0042562 GO:0042562
    ## GO:0043176 GO:0043176
    ## GO:0043287 GO:0043287
    ## GO:0043515 GO:0043515
    ## GO:0046790 GO:0046790
    ## GO:0046812 GO:0046812
    ## GO:0046848 GO:0046848
    ## GO:0046904 GO:0046904
    ## GO:0048037 GO:0048037
    ## GO:0050436 GO:0050436
    ## GO:0050825 GO:0050825
    ## GO:0050840 GO:0050840
    ## GO:0050997 GO:0050997
    ## GO:0051540 GO:0051540
    ## GO:0060090 GO:0060090
    ## GO:0070026 GO:0070026
    ## GO:0072328 GO:0072328
    ## GO:0072341 GO:0072341
    ## GO:0097243 GO:0097243
    ## GO:0098631 GO:0098631
    ## GO:0140272 GO:0140272
    ## GO:0140355 GO:0140355
    ## GO:1901567 GO:1901567
    ## GO:1901681 GO:1901681
    ## GO:1904483 GO:1904483
    ## GO:1904493 GO:1904493
    ## GO:1904517 GO:1904517
    ## GO:1990300 GO:1990300
    ## GO:0004362 GO:0004362
    ## GO:0004601 GO:0004601
    ## GO:0004784 GO:0004784
    ## GO:0004791 GO:0004791
    ## GO:0032542 GO:0032542
    ## GO:0045174 GO:0045174
    ## GO:0050605 GO:0050605
    ## GO:0051920 GO:0051920
    ## GO:0004873 GO:0004873
    ## GO:0004998 GO:0004998
    ## GO:0005044 GO:0005044
    ## GO:0008196 GO:0008196
    ## GO:0016964 GO:0016964
    ## GO:0030226 GO:0030226
    ## GO:0030228 GO:0030228
    ## GO:0033568 GO:0033568
    ## GO:0061714 GO:0061714
    ## GO:0070287 GO:0070287
    ## GO:0004694 GO:0004694
    ## GO:0030371 GO:0030371
    ## GO:0045183 GO:0045183
    ## GO:0000156 GO:0000156
    ## GO:0009927 GO:0009927
    ## GO:0031992 GO:0031992
    ## GO:0005085 GO:0005085
    ## GO:0097690 GO:0097690
    ## GO:0140416 GO:0140416
    ## GO:0016530 GO:0016530
    ## GO:0032977 GO:0032977
    ## GO:0097163 GO:0097163
    ## GO:0140132 GO:0140132
    ## GO:0140414 GO:0140414
    ## GO:0001070 GO:0001070
    ## GO:0001072 GO:0001072
    ## GO:0001073 GO:0001073
    ## GO:0003711 GO:0003711
    ## GO:0016987 GO:0016987
    ## GO:0043856 GO:0043856
    ## GO:0140223 GO:0140223
    ## GO:0000155 GO:0000155
    ## GO:0005034 GO:0005034
    ## GO:0019826 GO:0019826
    ## GO:0035991 GO:0035991
    ## GO:0061891 GO:0061891
    ## GO:0070027 GO:0070027
    ## GO:0097063 GO:0097063
    ## GO:0097077 GO:0097077
    ## GO:0106219 GO:0106219
    ## GO:0106254 GO:0106254
    ## GO:0140442 GO:0140442
    ## GO:0140311 GO:0140311
    ## GO:0140314 GO:0140314
    ## GO:0140315 GO:0140315
    ## GO:0140319 GO:0140319
    ##                                                                Description
    ## GO:0005515                                                 protein binding
    ## GO:0043167                                                     ion binding
    ## GO:0097159                                 organic cyclic compound binding
    ## GO:1901363                                   heterocyclic compound binding
    ## GO:0016740                                            transferase activity
    ## GO:0036094                                          small molecule binding
    ## GO:0008144                                                    drug binding
    ## GO:0097367                                 carbohydrate derivative binding
    ## GO:0016787                                              hydrolase activity
    ## GO:0022857                              transmembrane transporter activity
    ## GO:0044877                              protein-containing complex binding
    ## GO:0038023                                     signaling receptor activity
    ## GO:0030234                                       enzyme regulator activity
    ## GO:0003700                       DNA-binding transcription factor activity
    ## GO:0140096                         catalytic activity, acting on a protein
    ## GO:0005200                          structural constituent of cytoskeleton
    ## GO:0019911                         structural constituent of myelin sheath
    ## GO:0033218                                                   amide binding
    ## GO:0042165                                        neurotransmitter binding
    ## GO:0090079            translation regulator activity, nucleic acid binding
    ## GO:0016247                                      channel regulator activity
    ## GO:0030545                                     receptor regulator activity
    ## GO:0003712                              transcription coregulator activity
    ## GO:0004133                            glycogen debranching enzyme activity
    ## GO:0009975                                                cyclase activity
    ## GO:0016491                                         oxidoreductase activity
    ## GO:0016829                                                  lyase activity
    ## GO:0016853                                              isomerase activity
    ## GO:0016874                                                 ligase activity
    ## GO:0032451                                            demethylase activity
    ## GO:0061783                                peptidoglycan muralytic activity
    ## GO:0140097                               catalytic activity, acting on DNA
    ## GO:0140098                               catalytic activity, acting on RNA
    ## GO:0003735                              structural constituent of ribosome
    ## GO:0005199                             structural constituent of cell wall
    ## GO:0005201                     extracellular matrix structural constituent
    ## GO:0005212                              structural constituent of eye lens
    ## GO:0005213                           structural constituent of egg chorion
    ## GO:0008147                                  structural constituent of bone
    ## GO:0008307                                structural constituent of muscle
    ## GO:0008316                    structural constituent of vitelline membrane
    ## GO:0016490                  structural constituent of peritrophic membrane
    ## GO:0017056                          structural constituent of nuclear pore
    ## GO:0030280                        structural constituent of skin epidermis
    ## GO:0030281                   structural constituent of cutaneous appendage
    ## GO:0030527                             structural constituent of chromatin
    ## GO:0035804                              structural constituent of egg coat
    ## GO:0039660                                structural constituent of virion
    ## GO:0042302                               structural constituent of cuticle
    ## GO:0043886                           structural constituent of carboxysome
    ## GO:0097099                               structural constituent of albumen
    ## GO:0097493              structural molecule activity conferring elasticity
    ## GO:0098918                               structural constituent of synapse
    ## GO:0005319                                      lipid transporter activity
    ## GO:0005344                                         oxygen carrier activity
    ## GO:0140142                              nucleocytoplasmic carrier activity
    ## GO:0140318                                    protein transporter activity
    ## GO:0000035                                                    acyl binding
    ## GO:0003682                                               chromatin binding
    ## GO:0003823                                                 antigen binding
    ## GO:0005549                                                 odorant binding
    ## GO:0008289                                                   lipid binding
    ## GO:0008430                                                selenium binding
    ## GO:0015643                                         toxic substance binding
    ## GO:0019808                                               polyamine binding
    ## GO:0030246                                            carbohydrate binding
    ## GO:0031409                                                 pigment binding
    ## GO:0033226                                 2-aminoethylphosphonate binding
    ## GO:0035731                                 dinitrosyl-iron complex binding
    ## GO:0042562                                                 hormone binding
    ## GO:0043176                                                   amine binding
    ## GO:0043287                                poly(3-hydroxyalkanoate) binding
    ## GO:0043515                                             kinetochore binding
    ## GO:0046790                                                  virion binding
    ## GO:0046812                                       host cell surface binding
    ## GO:0046848                                          hydroxyapatite binding
    ## GO:0046904                                         calcium oxalate binding
    ## GO:0048037                                                cofactor binding
    ## GO:0050436                                             microfibril binding
    ## GO:0050825                                                     ice binding
    ## GO:0050840                                    extracellular matrix binding
    ## GO:0050997                               quaternary ammonium group binding
    ## GO:0051540                                           metal cluster binding
    ## GO:0060090                                      molecular adaptor activity
    ## GO:0070026                                            nitric oxide binding
    ## GO:0072328                                                  alkene binding
    ## GO:0072341                                     modified amino acid binding
    ## GO:0097243                                               flavonoid binding
    ## GO:0098631                                 cell adhesion mediator activity
    ## GO:0140272                                       exogenous protein binding
    ## GO:0140355                                  cargo receptor ligand activity
    ## GO:1901567                                   fatty acid derivative binding
    ## GO:1901681                                         sulfur compound binding
    ## GO:1904483                                   synthetic cannabinoid binding
    ## GO:1904493                 tetrahydrofolyl-poly(glutamate) polymer binding
    ## GO:1904517                                               MgATP(2-) binding
    ## GO:1990300                                             cellulosome binding
    ## GO:0004362                        glutathione-disulfide reductase activity
    ## GO:0004601                                             peroxidase activity
    ## GO:0004784                                   superoxide dismutase activity
    ## GO:0004791                        thioredoxin-disulfide reductase activity
    ## GO:0032542                                           sulfiredoxin activity
    ## GO:0045174                  glutathione dehydrogenase (ascorbate) activity
    ## GO:0050605                                   superoxide reductase activity
    ## GO:0051920                                          peroxiredoxin activity
    ## GO:0004873                            asialoglycoprotein receptor activity
    ## GO:0004998                                   transferrin receptor activity
    ## GO:0005044                                     scavenger receptor activity
    ## GO:0008196                                  vitellogenin receptor activity
    ## GO:0016964                         alpha-2 macroglobulin receptor activity
    ## GO:0030226                                apolipoprotein receptor activity
    ## GO:0030228                          lipoprotein particle receptor activity
    ## GO:0033568                                   lactoferrin receptor activity
    ## GO:0061714                                    folic acid receptor activity
    ## GO:0070287                                      ferritin receptor activity
    ## GO:0004694 eukaryotic translation initiation factor 2alpha kinase activity
    ## GO:0030371                                  translation repressor activity
    ## GO:0045183           translation factor activity, non-nucleic acid binding
    ## GO:0000156                        phosphorelay response regulator activity
    ## GO:0009927                       histidine phosphotransfer kinase activity
    ## GO:0031992                                      energy transducer activity
    ## GO:0005085                      guanyl-nucleotide exchange factor activity
    ## GO:0097690           iron ion transmembrane transporter inhibitor activity
    ## GO:0140416             DNA-binding transcription factor inhibitor activity
    ## GO:0016530                                       metallochaperone activity
    ## GO:0032977                                     membrane insertase activity
    ## GO:0097163                                         sulfur carrier activity
    ## GO:0140132                            iron-sulfur cluster carrier activity
    ## GO:0140414                   phosphopantetheine-dependent carrier activity
    ## GO:0001070                    RNA-binding transcription regulator activity
    ## GO:0001072      transcription antitermination factor activity, RNA binding
    ## GO:0001073      transcription antitermination factor activity, DNA binding
    ## GO:0003711                     transcription elongation regulator activity
    ## GO:0016987                                           sigma factor activity
    ## GO:0043856                           anti-sigma factor antagonist activity
    ## GO:0140223                general transcription initiation factor activity
    ## GO:0000155                             phosphorelay sensor kinase activity
    ## GO:0005034                                             osmosensor activity
    ## GO:0019826                                          oxygen sensor activity
    ## GO:0035991                                    nitric oxide sensor activity
    ## GO:0061891                                     calcium ion sensor activity
    ## GO:0070027                                 carbon monoxide sensor activity
    ## GO:0097063                                     cadmium ion sensor activity
    ## GO:0097077                                      copper ion sensor activity
    ## GO:0106219                                        zinc ion sensor activity
    ## GO:0106254                                           lipid sensor activity
    ## GO:0140442                                        peroxide sensor activity
    ## GO:0140311                                   protein sequestering activity
    ## GO:0140314                               calcium ion sequestering activity
    ## GO:0140315                                  iron ion sequestering activity
    ## GO:0140319                                         receptor decoy activity
    ##            Count GeneRatio
    ## GO:0005515    17     17/23
    ## GO:0043167     7      7/23
    ## GO:0097159     6      6/23
    ## GO:1901363     6      6/23
    ## GO:0016740     4      4/23
    ## GO:0036094     4      4/23
    ## GO:0008144     3      3/23
    ## GO:0097367     3      3/23
    ## GO:0016787     2      2/23
    ## GO:0022857     2      2/23
    ## GO:0044877     2      2/23
    ## GO:0038023     2      2/23
    ## GO:0030234     2      2/23
    ## GO:0003700     2      2/23
    ## GO:0140096     1      1/23
    ## GO:0005200     1      1/23
    ## GO:0019911     1      1/23
    ## GO:0033218     1      1/23
    ## GO:0042165     1      1/23
    ## GO:0090079     1      1/23
    ## GO:0016247     1      1/23
    ## GO:0030545     1      1/23
    ## GO:0003712     1      1/23
    ## GO:0004133     0      0/23
    ## GO:0009975     0      0/23
    ## GO:0016491     0      0/23
    ## GO:0016829     0      0/23
    ## GO:0016853     0      0/23
    ## GO:0016874     0      0/23
    ## GO:0032451     0      0/23
    ## GO:0061783     0      0/23
    ## GO:0140097     0      0/23
    ## GO:0140098     0      0/23
    ## GO:0003735     0      0/23
    ## GO:0005199     0      0/23
    ## GO:0005201     0      0/23
    ## GO:0005212     0      0/23
    ## GO:0005213     0      0/23
    ## GO:0008147     0      0/23
    ## GO:0008307     0      0/23
    ## GO:0008316     0      0/23
    ## GO:0016490     0      0/23
    ## GO:0017056     0      0/23
    ## GO:0030280     0      0/23
    ## GO:0030281     0      0/23
    ## GO:0030527     0      0/23
    ## GO:0035804     0      0/23
    ## GO:0039660     0      0/23
    ## GO:0042302     0      0/23
    ## GO:0043886     0      0/23
    ## GO:0097099     0      0/23
    ## GO:0097493     0      0/23
    ## GO:0098918     0      0/23
    ## GO:0005319     0      0/23
    ## GO:0005344     0      0/23
    ## GO:0140142     0      0/23
    ## GO:0140318     0      0/23
    ## GO:0000035     0      0/23
    ## GO:0003682     0      0/23
    ## GO:0003823     0      0/23
    ## GO:0005549     0      0/23
    ## GO:0008289     0      0/23
    ## GO:0008430     0      0/23
    ## GO:0015643     0      0/23
    ## GO:0019808     0      0/23
    ## GO:0030246     0      0/23
    ## GO:0031409     0      0/23
    ## GO:0033226     0      0/23
    ## GO:0035731     0      0/23
    ## GO:0042562     0      0/23
    ## GO:0043176     0      0/23
    ## GO:0043287     0      0/23
    ## GO:0043515     0      0/23
    ## GO:0046790     0      0/23
    ## GO:0046812     0      0/23
    ## GO:0046848     0      0/23
    ## GO:0046904     0      0/23
    ## GO:0048037     0      0/23
    ## GO:0050436     0      0/23
    ## GO:0050825     0      0/23
    ## GO:0050840     0      0/23
    ## GO:0050997     0      0/23
    ## GO:0051540     0      0/23
    ## GO:0060090     0      0/23
    ## GO:0070026     0      0/23
    ## GO:0072328     0      0/23
    ## GO:0072341     0      0/23
    ## GO:0097243     0      0/23
    ## GO:0098631     0      0/23
    ## GO:0140272     0      0/23
    ## GO:0140355     0      0/23
    ## GO:1901567     0      0/23
    ## GO:1901681     0      0/23
    ## GO:1904483     0      0/23
    ## GO:1904493     0      0/23
    ## GO:1904517     0      0/23
    ## GO:1990300     0      0/23
    ## GO:0004362     0      0/23
    ## GO:0004601     0      0/23
    ## GO:0004784     0      0/23
    ## GO:0004791     0      0/23
    ## GO:0032542     0      0/23
    ## GO:0045174     0      0/23
    ## GO:0050605     0      0/23
    ## GO:0051920     0      0/23
    ## GO:0004873     0      0/23
    ## GO:0004998     0      0/23
    ## GO:0005044     0      0/23
    ## GO:0008196     0      0/23
    ## GO:0016964     0      0/23
    ## GO:0030226     0      0/23
    ## GO:0030228     0      0/23
    ## GO:0033568     0      0/23
    ## GO:0061714     0      0/23
    ## GO:0070287     0      0/23
    ## GO:0004694     0      0/23
    ## GO:0030371     0      0/23
    ## GO:0045183     0      0/23
    ## GO:0000156     0      0/23
    ## GO:0009927     0      0/23
    ## GO:0031992     0      0/23
    ## GO:0005085     0      0/23
    ## GO:0097690     0      0/23
    ## GO:0140416     0      0/23
    ## GO:0016530     0      0/23
    ## GO:0032977     0      0/23
    ## GO:0097163     0      0/23
    ## GO:0140132     0      0/23
    ## GO:0140414     0      0/23
    ## GO:0001070     0      0/23
    ## GO:0001072     0      0/23
    ## GO:0001073     0      0/23
    ## GO:0003711     0      0/23
    ## GO:0016987     0      0/23
    ## GO:0043856     0      0/23
    ## GO:0140223     0      0/23
    ## GO:0000155     0      0/23
    ## GO:0005034     0      0/23
    ## GO:0019826     0      0/23
    ## GO:0035991     0      0/23
    ## GO:0061891     0      0/23
    ## GO:0070027     0      0/23
    ## GO:0097063     0      0/23
    ## GO:0097077     0      0/23
    ## GO:0106219     0      0/23
    ## GO:0106254     0      0/23
    ## GO:0140442     0      0/23
    ## GO:0140311     0      0/23
    ## GO:0140314     0      0/23
    ## GO:0140315     0      0/23
    ## GO:0140319     0      0/23
    ##                                                                                                      geneID
    ## GO:0005515 GSN/NBL1/TESC/TAC1/WIPF3/EXOC6B/GRIN2B/GRM5/CACNG3/PLLP/GSTA4/SACM1L/GFAP/CPLX2/NDRG1/OLIG1/DPF2
    ## GO:0043167                                                          GSN/PANK2/TESC/GRIN2B/PRKCG/PDE10A/DPF2
    ## GO:0097159                                                             PANK2/EIF4G3/PRKCG/OLIG1/PDE10A/DPF2
    ## GO:1901363                                                             PANK2/EIF4G3/PRKCG/OLIG1/PDE10A/DPF2
    ## GO:0016740                                                                           PANK2/PRKCG/GSTA4/DPF2
    ## GO:0036094                                                                        PANK2/GRIN2B/PRKCG/PDE10A
    ## GO:0008144                                                                               PANK2/GRIN2B/PRKCG
    ## GO:0097367                                                                               PANK2/PRKCG/PDE10A
    ## GO:0016787                                                                                    SACM1L/PDE10A
    ## GO:0022857                                                                                    GRIN2B/CACNG3
    ## GO:0044877                                                                                        GSN/WIPF3
    ## GO:0038023                                                                                      GRIN2B/GRM5
    ## GO:0030234                                                                                        TESC/GRM5
    ## GO:0003700                                                                                       OLIG1/DPF2
    ## GO:0140096                                                                                            PRKCG
    ## GO:0005200                                                                                             GFAP
    ## GO:0019911                                                                                             PLLP
    ## GO:0033218                                                                                           GRIN2B
    ## GO:0042165                                                                                           GRIN2B
    ## GO:0090079                                                                                           EIF4G3
    ## GO:0016247                                                                                           CACNG3
    ## GO:0030545                                                                                             NBL1
    ## GO:0003712                                                                                             DPF2
    ## GO:0004133                                                                                                 
    ## GO:0009975                                                                                                 
    ## GO:0016491                                                                                                 
    ## GO:0016829                                                                                                 
    ## GO:0016853                                                                                                 
    ## GO:0016874                                                                                                 
    ## GO:0032451                                                                                                 
    ## GO:0061783                                                                                                 
    ## GO:0140097                                                                                                 
    ## GO:0140098                                                                                                 
    ## GO:0003735                                                                                                 
    ## GO:0005199                                                                                                 
    ## GO:0005201                                                                                                 
    ## GO:0005212                                                                                                 
    ## GO:0005213                                                                                                 
    ## GO:0008147                                                                                                 
    ## GO:0008307                                                                                                 
    ## GO:0008316                                                                                                 
    ## GO:0016490                                                                                                 
    ## GO:0017056                                                                                                 
    ## GO:0030280                                                                                                 
    ## GO:0030281                                                                                                 
    ## GO:0030527                                                                                                 
    ## GO:0035804                                                                                                 
    ## GO:0039660                                                                                                 
    ## GO:0042302                                                                                                 
    ## GO:0043886                                                                                                 
    ## GO:0097099                                                                                                 
    ## GO:0097493                                                                                                 
    ## GO:0098918                                                                                                 
    ## GO:0005319                                                                                                 
    ## GO:0005344                                                                                                 
    ## GO:0140142                                                                                                 
    ## GO:0140318                                                                                                 
    ## GO:0000035                                                                                                 
    ## GO:0003682                                                                                                 
    ## GO:0003823                                                                                                 
    ## GO:0005549                                                                                                 
    ## GO:0008289                                                                                                 
    ## GO:0008430                                                                                                 
    ## GO:0015643                                                                                                 
    ## GO:0019808                                                                                                 
    ## GO:0030246                                                                                                 
    ## GO:0031409                                                                                                 
    ## GO:0033226                                                                                                 
    ## GO:0035731                                                                                                 
    ## GO:0042562                                                                                                 
    ## GO:0043176                                                                                                 
    ## GO:0043287                                                                                                 
    ## GO:0043515                                                                                                 
    ## GO:0046790                                                                                                 
    ## GO:0046812                                                                                                 
    ## GO:0046848                                                                                                 
    ## GO:0046904                                                                                                 
    ## GO:0048037                                                                                                 
    ## GO:0050436                                                                                                 
    ## GO:0050825                                                                                                 
    ## GO:0050840                                                                                                 
    ## GO:0050997                                                                                                 
    ## GO:0051540                                                                                                 
    ## GO:0060090                                                                                                 
    ## GO:0070026                                                                                                 
    ## GO:0072328                                                                                                 
    ## GO:0072341                                                                                                 
    ## GO:0097243                                                                                                 
    ## GO:0098631                                                                                                 
    ## GO:0140272                                                                                                 
    ## GO:0140355                                                                                                 
    ## GO:1901567                                                                                                 
    ## GO:1901681                                                                                                 
    ## GO:1904483                                                                                                 
    ## GO:1904493                                                                                                 
    ## GO:1904517                                                                                                 
    ## GO:1990300                                                                                                 
    ## GO:0004362                                                                                                 
    ## GO:0004601                                                                                                 
    ## GO:0004784                                                                                                 
    ## GO:0004791                                                                                                 
    ## GO:0032542                                                                                                 
    ## GO:0045174                                                                                                 
    ## GO:0050605                                                                                                 
    ## GO:0051920                                                                                                 
    ## GO:0004873                                                                                                 
    ## GO:0004998                                                                                                 
    ## GO:0005044                                                                                                 
    ## GO:0008196                                                                                                 
    ## GO:0016964                                                                                                 
    ## GO:0030226                                                                                                 
    ## GO:0030228                                                                                                 
    ## GO:0033568                                                                                                 
    ## GO:0061714                                                                                                 
    ## GO:0070287                                                                                                 
    ## GO:0004694                                                                                                 
    ## GO:0030371                                                                                                 
    ## GO:0045183                                                                                                 
    ## GO:0000156                                                                                                 
    ## GO:0009927                                                                                                 
    ## GO:0031992                                                                                                 
    ## GO:0005085                                                                                                 
    ## GO:0097690                                                                                                 
    ## GO:0140416                                                                                                 
    ## GO:0016530                                                                                                 
    ## GO:0032977                                                                                                 
    ## GO:0097163                                                                                                 
    ## GO:0140132                                                                                                 
    ## GO:0140414                                                                                                 
    ## GO:0001070                                                                                                 
    ## GO:0001072                                                                                                 
    ## GO:0001073                                                                                                 
    ## GO:0003711                                                                                                 
    ## GO:0016987                                                                                                 
    ## GO:0043856                                                                                                 
    ## GO:0140223                                                                                                 
    ## GO:0000155                                                                                                 
    ## GO:0005034                                                                                                 
    ## GO:0019826                                                                                                 
    ## GO:0035991                                                                                                 
    ## GO:0061891                                                                                                 
    ## GO:0070027                                                                                                 
    ## GO:0097063                                                                                                 
    ## GO:0097077                                                                                                 
    ## GO:0106219                                                                                                 
    ## GO:0106254                                                                                                 
    ## GO:0140442                                                                                                 
    ## GO:0140311                                                                                                 
    ## GO:0140314                                                                                                 
    ## GO:0140315                                                                                                 
    ## GO:0140319

    differential.genes <- annot[!is.na(annot[,2]) & gene.list.1.0 > 0.3,2]
    ggo <- groupGO(gene     = differential.genes,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   level    = 3,
                   readable = TRUE)
    ggo@result = ggo@result[order(ggo@result[,"Count"],decreasing=TRUE),]
    ggo@result

    ##                    ID
    ## GO:0005515 GO:0005515
    ## GO:0043167 GO:0043167
    ## GO:0097159 GO:0097159
    ## GO:1901363 GO:1901363
    ## GO:0022857 GO:0022857
    ## GO:0036094 GO:0036094
    ## GO:0038023 GO:0038023
    ## GO:0030234 GO:0030234
    ## GO:0016740 GO:0016740
    ## GO:0008144 GO:0008144
    ## GO:0044877 GO:0044877
    ## GO:0097367 GO:0097367
    ## GO:0016491 GO:0016491
    ## GO:0016787 GO:0016787
    ## GO:0008289 GO:0008289
    ## GO:0033218 GO:0033218
    ## GO:0048037 GO:0048037
    ## GO:0016247 GO:0016247
    ## GO:0003700 GO:0003700
    ## GO:0016874 GO:0016874
    ## GO:0140096 GO:0140096
    ## GO:0003682 GO:0003682
    ## GO:0042165 GO:0042165
    ## GO:0042562 GO:0042562
    ## GO:0050840 GO:0050840
    ## GO:0005085 GO:0005085
    ## GO:0030545 GO:0030545
    ## GO:0003712 GO:0003712
    ## GO:0004133 GO:0004133
    ## GO:0009975 GO:0009975
    ## GO:0016829 GO:0016829
    ## GO:0016853 GO:0016853
    ## GO:0032451 GO:0032451
    ## GO:0061783 GO:0061783
    ## GO:0140097 GO:0140097
    ## GO:0140098 GO:0140098
    ## GO:0003735 GO:0003735
    ## GO:0005199 GO:0005199
    ## GO:0005200 GO:0005200
    ## GO:0005201 GO:0005201
    ## GO:0005212 GO:0005212
    ## GO:0005213 GO:0005213
    ## GO:0008147 GO:0008147
    ## GO:0008307 GO:0008307
    ## GO:0008316 GO:0008316
    ## GO:0016490 GO:0016490
    ## GO:0017056 GO:0017056
    ## GO:0019911 GO:0019911
    ## GO:0030280 GO:0030280
    ## GO:0030281 GO:0030281
    ## GO:0030527 GO:0030527
    ## GO:0035804 GO:0035804
    ## GO:0039660 GO:0039660
    ## GO:0042302 GO:0042302
    ## GO:0043886 GO:0043886
    ## GO:0097099 GO:0097099
    ## GO:0097493 GO:0097493
    ## GO:0098918 GO:0098918
    ## GO:0005319 GO:0005319
    ## GO:0005344 GO:0005344
    ## GO:0140142 GO:0140142
    ## GO:0140318 GO:0140318
    ## GO:0000035 GO:0000035
    ## GO:0003823 GO:0003823
    ## GO:0005549 GO:0005549
    ## GO:0008430 GO:0008430
    ## GO:0015643 GO:0015643
    ## GO:0019808 GO:0019808
    ## GO:0030246 GO:0030246
    ## GO:0031409 GO:0031409
    ## GO:0033226 GO:0033226
    ## GO:0035731 GO:0035731
    ## GO:0043176 GO:0043176
    ## GO:0043287 GO:0043287
    ## GO:0043515 GO:0043515
    ## GO:0046790 GO:0046790
    ## GO:0046812 GO:0046812
    ## GO:0046848 GO:0046848
    ## GO:0046904 GO:0046904
    ## GO:0050436 GO:0050436
    ## GO:0050825 GO:0050825
    ## GO:0050997 GO:0050997
    ## GO:0051540 GO:0051540
    ## GO:0060090 GO:0060090
    ## GO:0070026 GO:0070026
    ## GO:0072328 GO:0072328
    ## GO:0072341 GO:0072341
    ## GO:0097243 GO:0097243
    ## GO:0098631 GO:0098631
    ## GO:0140272 GO:0140272
    ## GO:0140355 GO:0140355
    ## GO:1901567 GO:1901567
    ## GO:1901681 GO:1901681
    ## GO:1904483 GO:1904483
    ## GO:1904493 GO:1904493
    ## GO:1904517 GO:1904517
    ## GO:1990300 GO:1990300
    ## GO:0004362 GO:0004362
    ## GO:0004601 GO:0004601
    ## GO:0004784 GO:0004784
    ## GO:0004791 GO:0004791
    ## GO:0032542 GO:0032542
    ## GO:0045174 GO:0045174
    ## GO:0050605 GO:0050605
    ## GO:0051920 GO:0051920
    ## GO:0004873 GO:0004873
    ## GO:0004998 GO:0004998
    ## GO:0005044 GO:0005044
    ## GO:0008196 GO:0008196
    ## GO:0016964 GO:0016964
    ## GO:0030226 GO:0030226
    ## GO:0030228 GO:0030228
    ## GO:0033568 GO:0033568
    ## GO:0061714 GO:0061714
    ## GO:0070287 GO:0070287
    ## GO:0004694 GO:0004694
    ## GO:0030371 GO:0030371
    ## GO:0045183 GO:0045183
    ## GO:0090079 GO:0090079
    ## GO:0000156 GO:0000156
    ## GO:0009927 GO:0009927
    ## GO:0031992 GO:0031992
    ## GO:0097690 GO:0097690
    ## GO:0140416 GO:0140416
    ## GO:0016530 GO:0016530
    ## GO:0032977 GO:0032977
    ## GO:0097163 GO:0097163
    ## GO:0140132 GO:0140132
    ## GO:0140414 GO:0140414
    ## GO:0001070 GO:0001070
    ## GO:0001072 GO:0001072
    ## GO:0001073 GO:0001073
    ## GO:0003711 GO:0003711
    ## GO:0016987 GO:0016987
    ## GO:0043856 GO:0043856
    ## GO:0140223 GO:0140223
    ## GO:0000155 GO:0000155
    ## GO:0005034 GO:0005034
    ## GO:0019826 GO:0019826
    ## GO:0035991 GO:0035991
    ## GO:0061891 GO:0061891
    ## GO:0070027 GO:0070027
    ## GO:0097063 GO:0097063
    ## GO:0097077 GO:0097077
    ## GO:0106219 GO:0106219
    ## GO:0106254 GO:0106254
    ## GO:0140442 GO:0140442
    ## GO:0140311 GO:0140311
    ## GO:0140314 GO:0140314
    ## GO:0140315 GO:0140315
    ## GO:0140319 GO:0140319
    ##                                                                Description
    ## GO:0005515                                                 protein binding
    ## GO:0043167                                                     ion binding
    ## GO:0097159                                 organic cyclic compound binding
    ## GO:1901363                                   heterocyclic compound binding
    ## GO:0022857                              transmembrane transporter activity
    ## GO:0036094                                          small molecule binding
    ## GO:0038023                                     signaling receptor activity
    ## GO:0030234                                       enzyme regulator activity
    ## GO:0016740                                            transferase activity
    ## GO:0008144                                                    drug binding
    ## GO:0044877                              protein-containing complex binding
    ## GO:0097367                                 carbohydrate derivative binding
    ## GO:0016491                                         oxidoreductase activity
    ## GO:0016787                                              hydrolase activity
    ## GO:0008289                                                   lipid binding
    ## GO:0033218                                                   amide binding
    ## GO:0048037                                                cofactor binding
    ## GO:0016247                                      channel regulator activity
    ## GO:0003700                       DNA-binding transcription factor activity
    ## GO:0016874                                                 ligase activity
    ## GO:0140096                         catalytic activity, acting on a protein
    ## GO:0003682                                               chromatin binding
    ## GO:0042165                                        neurotransmitter binding
    ## GO:0042562                                                 hormone binding
    ## GO:0050840                                    extracellular matrix binding
    ## GO:0005085                      guanyl-nucleotide exchange factor activity
    ## GO:0030545                                     receptor regulator activity
    ## GO:0003712                              transcription coregulator activity
    ## GO:0004133                            glycogen debranching enzyme activity
    ## GO:0009975                                                cyclase activity
    ## GO:0016829                                                  lyase activity
    ## GO:0016853                                              isomerase activity
    ## GO:0032451                                            demethylase activity
    ## GO:0061783                                peptidoglycan muralytic activity
    ## GO:0140097                               catalytic activity, acting on DNA
    ## GO:0140098                               catalytic activity, acting on RNA
    ## GO:0003735                              structural constituent of ribosome
    ## GO:0005199                             structural constituent of cell wall
    ## GO:0005200                          structural constituent of cytoskeleton
    ## GO:0005201                     extracellular matrix structural constituent
    ## GO:0005212                              structural constituent of eye lens
    ## GO:0005213                           structural constituent of egg chorion
    ## GO:0008147                                  structural constituent of bone
    ## GO:0008307                                structural constituent of muscle
    ## GO:0008316                    structural constituent of vitelline membrane
    ## GO:0016490                  structural constituent of peritrophic membrane
    ## GO:0017056                          structural constituent of nuclear pore
    ## GO:0019911                         structural constituent of myelin sheath
    ## GO:0030280                        structural constituent of skin epidermis
    ## GO:0030281                   structural constituent of cutaneous appendage
    ## GO:0030527                             structural constituent of chromatin
    ## GO:0035804                              structural constituent of egg coat
    ## GO:0039660                                structural constituent of virion
    ## GO:0042302                               structural constituent of cuticle
    ## GO:0043886                           structural constituent of carboxysome
    ## GO:0097099                               structural constituent of albumen
    ## GO:0097493              structural molecule activity conferring elasticity
    ## GO:0098918                               structural constituent of synapse
    ## GO:0005319                                      lipid transporter activity
    ## GO:0005344                                         oxygen carrier activity
    ## GO:0140142                              nucleocytoplasmic carrier activity
    ## GO:0140318                                    protein transporter activity
    ## GO:0000035                                                    acyl binding
    ## GO:0003823                                                 antigen binding
    ## GO:0005549                                                 odorant binding
    ## GO:0008430                                                selenium binding
    ## GO:0015643                                         toxic substance binding
    ## GO:0019808                                               polyamine binding
    ## GO:0030246                                            carbohydrate binding
    ## GO:0031409                                                 pigment binding
    ## GO:0033226                                 2-aminoethylphosphonate binding
    ## GO:0035731                                 dinitrosyl-iron complex binding
    ## GO:0043176                                                   amine binding
    ## GO:0043287                                poly(3-hydroxyalkanoate) binding
    ## GO:0043515                                             kinetochore binding
    ## GO:0046790                                                  virion binding
    ## GO:0046812                                       host cell surface binding
    ## GO:0046848                                          hydroxyapatite binding
    ## GO:0046904                                         calcium oxalate binding
    ## GO:0050436                                             microfibril binding
    ## GO:0050825                                                     ice binding
    ## GO:0050997                               quaternary ammonium group binding
    ## GO:0051540                                           metal cluster binding
    ## GO:0060090                                      molecular adaptor activity
    ## GO:0070026                                            nitric oxide binding
    ## GO:0072328                                                  alkene binding
    ## GO:0072341                                     modified amino acid binding
    ## GO:0097243                                               flavonoid binding
    ## GO:0098631                                 cell adhesion mediator activity
    ## GO:0140272                                       exogenous protein binding
    ## GO:0140355                                  cargo receptor ligand activity
    ## GO:1901567                                   fatty acid derivative binding
    ## GO:1901681                                         sulfur compound binding
    ## GO:1904483                                   synthetic cannabinoid binding
    ## GO:1904493                 tetrahydrofolyl-poly(glutamate) polymer binding
    ## GO:1904517                                               MgATP(2-) binding
    ## GO:1990300                                             cellulosome binding
    ## GO:0004362                        glutathione-disulfide reductase activity
    ## GO:0004601                                             peroxidase activity
    ## GO:0004784                                   superoxide dismutase activity
    ## GO:0004791                        thioredoxin-disulfide reductase activity
    ## GO:0032542                                           sulfiredoxin activity
    ## GO:0045174                  glutathione dehydrogenase (ascorbate) activity
    ## GO:0050605                                   superoxide reductase activity
    ## GO:0051920                                          peroxiredoxin activity
    ## GO:0004873                            asialoglycoprotein receptor activity
    ## GO:0004998                                   transferrin receptor activity
    ## GO:0005044                                     scavenger receptor activity
    ## GO:0008196                                  vitellogenin receptor activity
    ## GO:0016964                         alpha-2 macroglobulin receptor activity
    ## GO:0030226                                apolipoprotein receptor activity
    ## GO:0030228                          lipoprotein particle receptor activity
    ## GO:0033568                                   lactoferrin receptor activity
    ## GO:0061714                                    folic acid receptor activity
    ## GO:0070287                                      ferritin receptor activity
    ## GO:0004694 eukaryotic translation initiation factor 2alpha kinase activity
    ## GO:0030371                                  translation repressor activity
    ## GO:0045183           translation factor activity, non-nucleic acid binding
    ## GO:0090079            translation regulator activity, nucleic acid binding
    ## GO:0000156                        phosphorelay response regulator activity
    ## GO:0009927                       histidine phosphotransfer kinase activity
    ## GO:0031992                                      energy transducer activity
    ## GO:0097690           iron ion transmembrane transporter inhibitor activity
    ## GO:0140416             DNA-binding transcription factor inhibitor activity
    ## GO:0016530                                       metallochaperone activity
    ## GO:0032977                                     membrane insertase activity
    ## GO:0097163                                         sulfur carrier activity
    ## GO:0140132                            iron-sulfur cluster carrier activity
    ## GO:0140414                   phosphopantetheine-dependent carrier activity
    ## GO:0001070                    RNA-binding transcription regulator activity
    ## GO:0001072      transcription antitermination factor activity, RNA binding
    ## GO:0001073      transcription antitermination factor activity, DNA binding
    ## GO:0003711                     transcription elongation regulator activity
    ## GO:0016987                                           sigma factor activity
    ## GO:0043856                           anti-sigma factor antagonist activity
    ## GO:0140223                general transcription initiation factor activity
    ## GO:0000155                             phosphorelay sensor kinase activity
    ## GO:0005034                                             osmosensor activity
    ## GO:0019826                                          oxygen sensor activity
    ## GO:0035991                                    nitric oxide sensor activity
    ## GO:0061891                                     calcium ion sensor activity
    ## GO:0070027                                 carbon monoxide sensor activity
    ## GO:0097063                                     cadmium ion sensor activity
    ## GO:0097077                                      copper ion sensor activity
    ## GO:0106219                                        zinc ion sensor activity
    ## GO:0106254                                           lipid sensor activity
    ## GO:0140442                                        peroxide sensor activity
    ## GO:0140311                                   protein sequestering activity
    ## GO:0140314                               calcium ion sequestering activity
    ## GO:0140315                                  iron ion sequestering activity
    ## GO:0140319                                         receptor decoy activity
    ##            Count GeneRatio
    ## GO:0005515    20     20/41
    ## GO:0043167    13     13/41
    ## GO:0097159     9      9/41
    ## GO:1901363     9      9/41
    ## GO:0022857     6      6/41
    ## GO:0036094     6      6/41
    ## GO:0038023     5      5/41
    ## GO:0030234     5      5/41
    ## GO:0016740     4      4/41
    ## GO:0008144     4      4/41
    ## GO:0044877     4      4/41
    ## GO:0097367     4      4/41
    ## GO:0016491     3      3/41
    ## GO:0016787     2      2/41
    ## GO:0008289     2      2/41
    ## GO:0033218     2      2/41
    ## GO:0048037     2      2/41
    ## GO:0016247     2      2/41
    ## GO:0003700     2      2/41
    ## GO:0016874     1      1/41
    ## GO:0140096     1      1/41
    ## GO:0003682     1      1/41
    ## GO:0042165     1      1/41
    ## GO:0042562     1      1/41
    ## GO:0050840     1      1/41
    ## GO:0005085     1      1/41
    ## GO:0030545     1      1/41
    ## GO:0003712     1      1/41
    ## GO:0004133     0      0/41
    ## GO:0009975     0      0/41
    ## GO:0016829     0      0/41
    ## GO:0016853     0      0/41
    ## GO:0032451     0      0/41
    ## GO:0061783     0      0/41
    ## GO:0140097     0      0/41
    ## GO:0140098     0      0/41
    ## GO:0003735     0      0/41
    ## GO:0005199     0      0/41
    ## GO:0005200     0      0/41
    ## GO:0005201     0      0/41
    ## GO:0005212     0      0/41
    ## GO:0005213     0      0/41
    ## GO:0008147     0      0/41
    ## GO:0008307     0      0/41
    ## GO:0008316     0      0/41
    ## GO:0016490     0      0/41
    ## GO:0017056     0      0/41
    ## GO:0019911     0      0/41
    ## GO:0030280     0      0/41
    ## GO:0030281     0      0/41
    ## GO:0030527     0      0/41
    ## GO:0035804     0      0/41
    ## GO:0039660     0      0/41
    ## GO:0042302     0      0/41
    ## GO:0043886     0      0/41
    ## GO:0097099     0      0/41
    ## GO:0097493     0      0/41
    ## GO:0098918     0      0/41
    ## GO:0005319     0      0/41
    ## GO:0005344     0      0/41
    ## GO:0140142     0      0/41
    ## GO:0140318     0      0/41
    ## GO:0000035     0      0/41
    ## GO:0003823     0      0/41
    ## GO:0005549     0      0/41
    ## GO:0008430     0      0/41
    ## GO:0015643     0      0/41
    ## GO:0019808     0      0/41
    ## GO:0030246     0      0/41
    ## GO:0031409     0      0/41
    ## GO:0033226     0      0/41
    ## GO:0035731     0      0/41
    ## GO:0043176     0      0/41
    ## GO:0043287     0      0/41
    ## GO:0043515     0      0/41
    ## GO:0046790     0      0/41
    ## GO:0046812     0      0/41
    ## GO:0046848     0      0/41
    ## GO:0046904     0      0/41
    ## GO:0050436     0      0/41
    ## GO:0050825     0      0/41
    ## GO:0050997     0      0/41
    ## GO:0051540     0      0/41
    ## GO:0060090     0      0/41
    ## GO:0070026     0      0/41
    ## GO:0072328     0      0/41
    ## GO:0072341     0      0/41
    ## GO:0097243     0      0/41
    ## GO:0098631     0      0/41
    ## GO:0140272     0      0/41
    ## GO:0140355     0      0/41
    ## GO:1901567     0      0/41
    ## GO:1901681     0      0/41
    ## GO:1904483     0      0/41
    ## GO:1904493     0      0/41
    ## GO:1904517     0      0/41
    ## GO:1990300     0      0/41
    ## GO:0004362     0      0/41
    ## GO:0004601     0      0/41
    ## GO:0004784     0      0/41
    ## GO:0004791     0      0/41
    ## GO:0032542     0      0/41
    ## GO:0045174     0      0/41
    ## GO:0050605     0      0/41
    ## GO:0051920     0      0/41
    ## GO:0004873     0      0/41
    ## GO:0004998     0      0/41
    ## GO:0005044     0      0/41
    ## GO:0008196     0      0/41
    ## GO:0016964     0      0/41
    ## GO:0030226     0      0/41
    ## GO:0030228     0      0/41
    ## GO:0033568     0      0/41
    ## GO:0061714     0      0/41
    ## GO:0070287     0      0/41
    ## GO:0004694     0      0/41
    ## GO:0030371     0      0/41
    ## GO:0045183     0      0/41
    ## GO:0090079     0      0/41
    ## GO:0000156     0      0/41
    ## GO:0009927     0      0/41
    ## GO:0031992     0      0/41
    ## GO:0097690     0      0/41
    ## GO:0140416     0      0/41
    ## GO:0016530     0      0/41
    ## GO:0032977     0      0/41
    ## GO:0097163     0      0/41
    ## GO:0140132     0      0/41
    ## GO:0140414     0      0/41
    ## GO:0001070     0      0/41
    ## GO:0001072     0      0/41
    ## GO:0001073     0      0/41
    ## GO:0003711     0      0/41
    ## GO:0016987     0      0/41
    ## GO:0043856     0      0/41
    ## GO:0140223     0      0/41
    ## GO:0000155     0      0/41
    ## GO:0005034     0      0/41
    ## GO:0019826     0      0/41
    ## GO:0035991     0      0/41
    ## GO:0061891     0      0/41
    ## GO:0070027     0      0/41
    ## GO:0097063     0      0/41
    ## GO:0097077     0      0/41
    ## GO:0106219     0      0/41
    ## GO:0106254     0      0/41
    ## GO:0140442     0      0/41
    ## GO:0140311     0      0/41
    ## GO:0140314     0      0/41
    ## GO:0140315     0      0/41
    ## GO:0140319     0      0/41
    ##                                                                                                                                    geneID
    ## GO:0005515 CDC42BPA/SPOUT1/KCNAB1/SH3GL2/SNCA/NFKBIB/DKKL1/GABRA5/BORCS8/SH3BP5/NUDT18/ILF3/CIB2/ZBTB4/CRHR1/GPS1/NRSN1/MEF2C/PCP4/BORCS7
    ## GO:0043167                                                 CDC42BPA/PAMR1/KCNAB1/AK5/SNCA/EXTL3/NUDT18/CIB2/ZBTB4/SPOCK1/PCP4/ATP9B/ACSL5
    ## GO:0097159                                                                        CDC42BPA/SPOUT1/KCNAB1/AK5/ILF3/ZBTB4/MEF2C/ATP9B/ACSL5
    ## GO:1901363                                                                        CDC42BPA/SPOUT1/KCNAB1/AK5/ILF3/ZBTB4/MEF2C/ATP9B/ACSL5
    ## GO:0022857                                                                                    KCNK2/KCNAB1/GABRA5/COX6A2/SLC6A15/ATP6V1G2
    ## GO:0036094                                                                                          CDC42BPA/KCNAB1/AK5/ZBTB4/ATP9B/ACSL5
    ## GO:0038023                                                                                                 GPR158/HRH3/GABRA5/GFRA2/CRHR1
    ## GO:0030234                                                                                                 SNCA/COX6A2/SH3BP5/GPS1/SPOCK1
    ## GO:0016740                                                                                                      CDC42BPA/SPOUT1/AK5/EXTL3
    ## GO:0008144                                                                                                       CDC42BPA/AK5/ATP9B/ACSL5
    ## GO:0044877                                                                                                         SNCA/CIB2/CRHR1/SPOCK1
    ## GO:0097367                                                                                                       CDC42BPA/AK5/ATP9B/ACSL5
    ## GO:0016491                                                                                                             KCNAB1/SNCA/COX6A2
    ## GO:0016787                                                                                                                  NUDT18/PLCXD2
    ## GO:0008289                                                                                                                    SH3GL2/SNCA
    ## GO:0033218                                                                                                                  TMEM158/CRHR1
    ## GO:0048037                                                                                                                  KCNAB1/COQ10A
    ## GO:0016247                                                                                                                   CNIH3/KCNAB1
    ## GO:0003700                                                                                                                    ZBTB4/MEF2C
    ## GO:0016874                                                                                                                          ACSL5
    ## GO:0140096                                                                                                                       CDC42BPA
    ## GO:0003682                                                                                                                          MEF2C
    ## GO:0042165                                                                                                                          CRHR1
    ## GO:0042562                                                                                                                          CRHR1
    ## GO:0050840                                                                                                                         SPOCK1
    ## GO:0005085                                                                                                                         SH3BP5
    ## GO:0030545                                                                                                                          DKKL1
    ## GO:0003712                                                                                                                         NFKBIB
    ## GO:0004133                                                                                                                               
    ## GO:0009975                                                                                                                               
    ## GO:0016829                                                                                                                               
    ## GO:0016853                                                                                                                               
    ## GO:0032451                                                                                                                               
    ## GO:0061783                                                                                                                               
    ## GO:0140097                                                                                                                               
    ## GO:0140098                                                                                                                               
    ## GO:0003735                                                                                                                               
    ## GO:0005199                                                                                                                               
    ## GO:0005200                                                                                                                               
    ## GO:0005201                                                                                                                               
    ## GO:0005212                                                                                                                               
    ## GO:0005213                                                                                                                               
    ## GO:0008147                                                                                                                               
    ## GO:0008307                                                                                                                               
    ## GO:0008316                                                                                                                               
    ## GO:0016490                                                                                                                               
    ## GO:0017056                                                                                                                               
    ## GO:0019911                                                                                                                               
    ## GO:0030280                                                                                                                               
    ## GO:0030281                                                                                                                               
    ## GO:0030527                                                                                                                               
    ## GO:0035804                                                                                                                               
    ## GO:0039660                                                                                                                               
    ## GO:0042302                                                                                                                               
    ## GO:0043886                                                                                                                               
    ## GO:0097099                                                                                                                               
    ## GO:0097493                                                                                                                               
    ## GO:0098918                                                                                                                               
    ## GO:0005319                                                                                                                               
    ## GO:0005344                                                                                                                               
    ## GO:0140142                                                                                                                               
    ## GO:0140318                                                                                                                               
    ## GO:0000035                                                                                                                               
    ## GO:0003823                                                                                                                               
    ## GO:0005549                                                                                                                               
    ## GO:0008430                                                                                                                               
    ## GO:0015643                                                                                                                               
    ## GO:0019808                                                                                                                               
    ## GO:0030246                                                                                                                               
    ## GO:0031409                                                                                                                               
    ## GO:0033226                                                                                                                               
    ## GO:0035731                                                                                                                               
    ## GO:0043176                                                                                                                               
    ## GO:0043287                                                                                                                               
    ## GO:0043515                                                                                                                               
    ## GO:0046790                                                                                                                               
    ## GO:0046812                                                                                                                               
    ## GO:0046848                                                                                                                               
    ## GO:0046904                                                                                                                               
    ## GO:0050436                                                                                                                               
    ## GO:0050825                                                                                                                               
    ## GO:0050997                                                                                                                               
    ## GO:0051540                                                                                                                               
    ## GO:0060090                                                                                                                               
    ## GO:0070026                                                                                                                               
    ## GO:0072328                                                                                                                               
    ## GO:0072341                                                                                                                               
    ## GO:0097243                                                                                                                               
    ## GO:0098631                                                                                                                               
    ## GO:0140272                                                                                                                               
    ## GO:0140355                                                                                                                               
    ## GO:1901567                                                                                                                               
    ## GO:1901681                                                                                                                               
    ## GO:1904483                                                                                                                               
    ## GO:1904493                                                                                                                               
    ## GO:1904517                                                                                                                               
    ## GO:1990300                                                                                                                               
    ## GO:0004362                                                                                                                               
    ## GO:0004601                                                                                                                               
    ## GO:0004784                                                                                                                               
    ## GO:0004791                                                                                                                               
    ## GO:0032542                                                                                                                               
    ## GO:0045174                                                                                                                               
    ## GO:0050605                                                                                                                               
    ## GO:0051920                                                                                                                               
    ## GO:0004873                                                                                                                               
    ## GO:0004998                                                                                                                               
    ## GO:0005044                                                                                                                               
    ## GO:0008196                                                                                                                               
    ## GO:0016964                                                                                                                               
    ## GO:0030226                                                                                                                               
    ## GO:0030228                                                                                                                               
    ## GO:0033568                                                                                                                               
    ## GO:0061714                                                                                                                               
    ## GO:0070287                                                                                                                               
    ## GO:0004694                                                                                                                               
    ## GO:0030371                                                                                                                               
    ## GO:0045183                                                                                                                               
    ## GO:0090079                                                                                                                               
    ## GO:0000156                                                                                                                               
    ## GO:0009927                                                                                                                               
    ## GO:0031992                                                                                                                               
    ## GO:0097690                                                                                                                               
    ## GO:0140416                                                                                                                               
    ## GO:0016530                                                                                                                               
    ## GO:0032977                                                                                                                               
    ## GO:0097163                                                                                                                               
    ## GO:0140132                                                                                                                               
    ## GO:0140414                                                                                                                               
    ## GO:0001070                                                                                                                               
    ## GO:0001072                                                                                                                               
    ## GO:0001073                                                                                                                               
    ## GO:0003711                                                                                                                               
    ## GO:0016987                                                                                                                               
    ## GO:0043856                                                                                                                               
    ## GO:0140223                                                                                                                               
    ## GO:0000155                                                                                                                               
    ## GO:0005034                                                                                                                               
    ## GO:0019826                                                                                                                               
    ## GO:0035991                                                                                                                               
    ## GO:0061891                                                                                                                               
    ## GO:0070027                                                                                                                               
    ## GO:0097063                                                                                                                               
    ## GO:0097077                                                                                                                               
    ## GO:0106219                                                                                                                               
    ## GO:0106254                                                                                                                               
    ## GO:0140442                                                                                                                               
    ## GO:0140311                                                                                                                               
    ## GO:0140314                                                                                                                               
    ## GO:0140315                                                                                                                               
    ## GO:0140319

    differential.genes

    ##  [1] "8476"   "149111" "3776"   "57512"  "51490"  "25891"  "11255"  "7881"  
    ##  [9] "26289"  "6456"   "84935"  "6622"   "284339" "4793"   "27120"  "2558"  
    ## [17] "1339"   "55117"  "93058"  "729991" "9467"   "2137"   "79873"  "2675"  
    ## [25] "3609"   "10518"  "25907"  "57659"  "3859"   "1394"   "2873"   "140767"
    ## [33] "6695"   "4208"   "8214"   "257068" "5121"   "534"    "374868" "119032"
    ## [41] "51703"

    length(gene.list.1.0)

    ## [1] 16017

    length(gene.list)

    ## [1] 16017

    annot['ENTREZID']

    ##        ENTREZID
    ## 1        114786
    ## 2         64321
    ## 3         29088
    ## 4         10434
    ## 5          6917
    ## 6          8601
    ## 7         51606
    ## 8          4986
    ## 9          2831
    ## 10         9821
    ## 11         <NA>
    ## 12       389658
    ## 13         9705
    ## 14       115294
    ## 15         <NA>
    ## 16        54212
    ## 17        23212
    ## 18       137872
    ## 19         <NA>
    ## 20         4603
    ## 21        80124
    ## 22        23678
    ## 23       157777
    ## 24       641638
    ## 25    100129654
    ## 26       286187
    ## 27        10987
    ## 28        79848
    ## 29        10565
    ## 30        57094
    ## 31        80243
    ## 32         <NA>
    ## 33        23213
    ## 34        81796
    ## 35        10499
    ## 36        23471
    ## 37        51110
    ## 38         2138
    ## 39         <NA>
    ## 40         9242
    ## 41         9312
    ## 42         7013
    ## 43       157869
    ## 44         <NA>
    ## 45         6129
    ## 46       157506
    ## 47        27067
    ## 48        55284
    ## 49         6921
    ## 50        54968
    ## 51        23643
    ## 52        56704
    ## 53        54332
    ## 54         <NA>
    ## 55         <NA>
    ## 56        83690
    ## 57         <NA>
    ## 58         7021
    ## 59         4172
    ## 60         <NA>
    ## 61        85315
    ## 62       114327
    ## 63         9697
    ## 64        28978
    ## 65         2940
    ## 66         <NA>
    ## 67        56479
    ## 68        22999
    ## 69         <NA>
    ## 70        79627
    ## 71       135152
    ## 72        60682
    ## 73       135154
    ## 74        57579
    ## 75         1297
    ## 76         1310
    ## 77        55788
    ## 78          577
    ## 79        23469
    ## 80         7803
    ## 81       202559
    ## 82         5558
    ## 83        51715
    ## 84         9532
    ## 85         <NA>
    ## 86       221336
    ## 87          667
    ## 88        84317
    ## 89        92856
    ## 90        26469
    ## 91       205147
    ## 92        50649
    ## 93         <NA>
    ## 94       130074
    ## 95        55041
    ## 96         <NA>
    ## 97         9394
    ## 98        56886
    ## 99        93082
    ## 100       10865
    ## 101       55683
    ## 102       81562
    ## 103       26504
    ## 104       26505
    ## 105      200539
    ## 106       51239
    ## 107       54910
    ## 108       51252
    ## 109        1329
    ## 110       10120
    ## 111        <NA>
    ## 112        7535
    ## 113       23505
    ## 114      200403
    ## 115        3631
    ## 116      493753
    ## 117       25972
    ## 118       11320
    ## 119        <NA>
    ## 120       80705
    ## 121       51601
    ## 122      129531
    ## 123       51263
    ## 124       10190
    ## 125        9669
    ## 126       51455
    ## 127        3899
    ## 128      164832
    ## 129        9486
    ## 130       79031
    ## 131        4862
    ## 132        6160
    ## 133       11138
    ## 134       55571
    ## 135      284996
    ## 136      200407
    ## 137        <NA>
    ## 138      731220
    ## 139        9448
    ## 140        7850
    ## 141        3554
    ## 142        8808
    ## 143        8807
    ## 144      389015
    ## 145        6549
    ## 146       84804
    ## 147      130827
    ## 148        <NA>
    ## 149   100506421
    ## 150        5455
    ## 151   104940698
    ## 152       64965
    ## 153       11250
    ## 154        9392
    ## 155        <NA>
    ## 156        <NA>
    ## 157        <NA>
    ## 158        2274
    ## 159        8440
    ## 160        <NA>
    ## 161        <NA>
    ## 162       80146
    ## 163        7174
    ## 164       93081
    ## 165        <NA>
    ## 166       54841
    ## 167        2073
    ## 168       51454
    ## 169        1281
    ## 170        1290
    ## 171       84128
    ## 172       30061
    ## 173        <NA>
    ## 174        <NA>
    ## 175       57181
    ## 176        <NA>
    ## 177       23671
    ## 178        <NA>
    ## 179        8436
    ## 180       64859
    ## 181        4430
    ## 182        6775
    ## 183        <NA>
    ## 184        6772
    ## 185        2744
    ## 186        4664
    ## 187   100131211
    ## 188       54842
    ## 189        3628
    ## 190       26275
    ## 191        <NA>
    ## 192        2660
    ## 193        5378
    ## 194       94101
    ## 195       64172
    ## 196        <NA>
    ## 197       54529
    ## 198        <NA>
    ## 199        9262
    ## 200       57520
    ## 201        9330
    ## 202       80055
    ## 203       91526
    ## 204       23451
    ## 205       80219
    ## 206        3329
    ## 207        3336
    ## 208       25843
    ## 209      130132
    ## 210        <NA>
    ## 211       92935
    ## 212        5334
    ## 213        <NA>
    ## 214       23314
    ## 215        <NA>
    ## 216      129450
    ## 217       79568
    ## 218       26010
    ## 219      130535
    ## 220        <NA>
    ## 221         316
    ## 222        <NA>
    ## 223        <NA>
    ## 224        9689
    ## 225        1195
    ## 226       53938
    ## 227       60491
    ## 228        4999
    ## 229        <NA>
    ## 230      285172
    ## 231        4709
    ## 232        8837
    ## 233         841
    ## 234       66008
    ## 235       55437
    ## 236       65062
    ## 237        <NA>
    ## 238       58538
    ## 239       57679
    ## 240        8324
    ## 241        <NA>
    ## 242        7341
    ## 243       51602
    ## 244         659
    ## 245      150864
    ## 246      130026
    ## 247       55759
    ## 248       79800
    ## 249       65065
    ## 250       57404
    ## 251       10152
    ## 252       65059
    ## 253         940
    ## 254        1493
    ## 255      117583
    ## 256        8828
    ## 257       54891
    ## 258        <NA>
    ## 259        <NA>
    ## 260        4719
    ## 261        1933
    ## 262        2825
    ## 263       57683
    ## 264        8745
    ## 265      130752
    ## 266       22868
    ## 267        8609
    ## 268        <NA>
    ## 269        1385
    ## 270      151194
    ## 271      151195
    ## 272        7855
    ## 273      389072
    ## 274        1418
    ## 275        <NA>
    ## 276        3417
    ## 277      200576
    ## 278        4133
    ## 279      285175
    ## 280        6120
    ## 281      151050
    ## 282          33
    ## 283        4632
    ## 284       10314
    ## 285        1373
    ## 286        2066
    ## 287       22807
    ## 288      402117
    ## 289         580
    ## 290         471
    ## 291        2335
    ## 292       55686
    ## 293        <NA>
    ## 294       55825
    ## 295       92691
    ## 296        7520
    ## 297        <NA>
    ## 298       50485
    ## 299        <NA>
    ## 300        6168
    ## 301        3485
    ## 302        3488
    ## 303        <NA>
    ## 304        7145
    ## 305        3579
    ## 306       10109
    ## 307          14
    ## 308       25953
    ## 309       64114
    ## 310      375307
    ## 311        6556
    ## 312       58190
    ## 313        7429
    ## 314       57695
    ## 315        9125
    ## 316       84812
    ## 317        <NA>
    ## 318         617
    ## 319       64320
    ## 320       27148
    ## 321        9654
    ## 322        1593
    ## 323        7475
    ## 324       80326
    ## 325        8941
    ## 326        1412
    ## 327      255101
    ## 328       79840
    ## 329      151295
    ## 330       27013
    ## 331       79137
    ## 332      130617
    ## 333       10058
    ## 334       79065
    ## 335       55139
    ## 336       79411
    ## 337        8576
    ## 338        7277
    ## 339        3300
    ## 340        5798
    ## 341      389075
    ## 342       23549
    ## 343        1674
    ## 344       10290
    ## 345       29926
    ## 346       55515
    ## 347       79586
    ## 348       23363
    ## 349      130612
    ## 350        3623
    ## 351      114790
    ## 352        6508
    ## 353        <NA>
    ## 354        2043
    ## 355        <NA>
    ## 356      130367
    ## 357       10056
    ## 358        <NA>
    ## 359      116255
    ## 360        2181
    ## 361        <NA>
    ## 362       23704
    ## 363        7857
    ## 364      130340
    ## 365       57590
    ## 366       65080
    ## 367        5270
    ## 368        8452
    ## 369       55619
    ## 370       57624
    ## 371        3667
    ## 372       84236
    ## 373        1286
    ## 374       56947
    ## 375        3267
    ## 376        <NA>
    ## 377       80704
    ## 378       80309
    ## 379       55022
    ## 380       92737
    ## 381        9320
    ## 382      130888
    ## 383      151473
    ## 384        3431
    ## 385       11262
    ## 386        6672
    ## 387        <NA>
    ## 388       51719
    ## 389       81618
    ## 390        <NA>
    ## 391        <NA>
    ## 392      130560
    ## 393        <NA>
    ## 394        5707
    ## 395       80210
    ## 396       93010
    ## 397        4691
    ## 398        5757
    ## 399        5147
    ## 400       64708
    ## 401        4880
    ## 402      129563
    ## 403        9427
    ## 404      646960
    ## 405        9470
    ## 406       80303
    ## 407       26058
    ## 408        3769
    ## 409        <NA>
    ## 410       25791
    ## 411        4759
    ## 412        3635
    ## 413       55054
    ## 414        6295
    ## 415        8527
    ## 416       55230
    ## 417      414061
    ## 418       55355
    ## 419       10123
    ## 420       23677
    ## 421      116987
    ## 422      401036
    ## 423        <NA>
    ## 424       57007
    ## 425       10920
    ## 426        1293
    ## 427       79083
    ## 428        9208
    ## 429       10267
    ## 430      140739
    ## 431       51540
    ## 432      377007
    ## 433      151176
    ## 434       80895
    ## 435       55502
    ## 436        8864
    ## 437       26146
    ## 438       51665
    ## 439        <NA>
    ## 440      117581
    ## 441        9759
    ## 442        4705
    ## 443      150678
    ## 444      150677
    ## 445        <NA>
    ## 446        2817
    ## 447        <NA>
    ## 448       51281
    ## 449      285193
    ## 450        <NA>
    ## 451       57140
    ## 452        <NA>
    ## 453       11132
    ## 454         547
    ## 455       25992
    ## 456      130916
    ## 457        <NA>
    ## 458       23178
    ## 459        5510
    ## 460       50636
    ## 461        3069
    ## 462        <NA>
    ## 463        9855
    ## 464       10494
    ## 465         666
    ## 466       51078
    ## 467       23192
    ## 468        1841
    ## 469       84289
    ## 470      728294
    ## 471      129807
    ## 472      345757
    ## 473        7903
    ## 474      353189
    ## 475        <NA>
    ## 476        <NA>
    ## 477       23262
    ## 478       54826
    ## 479        5066
    ## 480       28316
    ## 481      220441
    ## 482        <NA>
    ## 483       23556
    ## 484        <NA>
    ## 485        8792
    ## 486        <NA>
    ## 487       54877
    ## 488       23239
    ## 489         596
    ## 490        2531
    ## 491        9525
    ## 492        5271
    ## 493        1005
    ## 494       28513
    ## 495        <NA>
    ## 496       92126
    ## 497        <NA>
    ## 498        7247
    ## 499       84365
    ## 500       23332
    ## 501       29842
    ## 502        2736
    ## 503        3625
    ## 504        5899
    ## 505        <NA>
    ## 506       79134
    ## 507       57669
    ## 508        5775
    ## 509       80775
    ## 510      200373
    ## 511        6344
    ## 512      140738
    ## 513        1622
    ## 514        <NA>
    ## 515       55240
    ## 516      165257
    ## 517        <NA>
    ## 518       51141
    ## 519       54520
    ## 520        <NA>
    ## 521        8886
    ## 522       57628
    ## 523       10096
    ## 524       80255
    ## 525      116372
    ## 526      344148
    ## 527        4249
    ## 528       81615
    ## 529        <NA>
    ## 530         905
    ## 531       22930
    ## 532       84083
    ## 533        <NA>
    ## 534       23518
    ## 535       23190
    ## 536        3938
    ## 537        4175
    ## 538        <NA>
    ## 539        7852
    ## 540       80731
    ## 541        <NA>
    ## 542        1604
    ## 543        <NA>
    ## 544        <NA>
    ## 545        5208
    ## 546       55432
    ## 547        <NA>
    ## 548        9214
    ## 549        9261
    ## 550        8444
    ## 551        1939
    ## 552       83593
    ## 553        9641
    ## 554        <NA>
    ## 555       23380
    ## 556      729533
    ## 557      338382
    ## 558      148811
    ## 559      254428
    ## 560        8934
    ## 561       64710
    ## 562       85414
    ## 563        2005
    ## 564      148808
    ## 565        5129
    ## 566        <NA>
    ## 567       93273
    ## 568       55220
    ## 569       81788
    ## 570        9911
    ## 571       25778
    ## 572        5929
    ## 573      388730
    ## 574        6900
    ## 575       23114
    ## 576       10446
    ## 577        4194
    ## 578        5287
    ## 579       84919
    ## 580       22874
    ## 581      127845
    ## 582       55224
    ## 583        9580
    ## 584        6635
    ## 585        <NA>
    ## 586        9877
    ## 587         493
    ## 588       26254
    ## 589        5549
    ## 590        2331
    ## 591        7832
    ## 592        <NA>
    ## 593        4608
    ## 594         134
    ## 595        8497
    ## 596       92703
    ## 597       51706
    ## 598       51094
    ## 599       59349
    ## 600        5877
    ## 601       10765
    ## 602        <NA>
    ## 603      127833
    ## 604        4660
    ## 605       29089
    ## 606       59352
    ## 607      127829
    ## 608        9283
    ## 609        6051
    ## 610       10440
    ## 611       25802
    ## 612      149345
    ## 613       55705
    ## 614       89796
    ## 615        <NA>
    ## 616        1465
    ## 617       23612
    ## 618        7135
    ## 619        3898
    ## 620        7139
    ## 621        5317
    ## 622       91156
    ## 623      252839
    ## 624         779
    ## 625       23046
    ## 626        <NA>
    ## 627        <NA>
    ## 628       23271
    ## 629       83479
    ## 630        <NA>
    ## 631        <NA>
    ## 632        <NA>
    ## 633        5788
    ## 634      140609
    ## 635       56956
    ## 636        <NA>
    ## 637      163486
    ## 638      360023
    ## 639      259266
    ## 640        3075
    ## 641      343450
    ## 642       79577
    ## 643        8707
    ## 644       51022
    ## 645        <NA>
    ## 646       51377
    ## 647        5997
    ## 648        6003
    ## 649        <NA>
    ## 650      339479
    ## 651        5321
    ## 652        5743
    ## 653        <NA>
    ## 654        7175
    ## 655       10216
    ## 656       83872
    ## 657        <NA>
    ## 658       10625
    ## 659       54823
    ## 660       81627
    ## 661        6045
    ## 662        <NA>
    ## 663        <NA>
    ## 664       80267
    ## 665        <NA>
    ## 666      116461
    ## 667       23127
    ## 668       23179
    ## 669       10092
    ## 670        4688
    ## 671        9887
    ## 672       23057
    ## 673        3918
    ## 674        3915
    ## 675       81626
    ## 676        1660
    ## 677       80896
    ## 678       85397
    ## 679        6004
    ## 680        6041
    ## 681      353299
    ## 682        <NA>
    ## 683        2752
    ## 684         777
    ## 685       51278
    ## 686        3140
    ## 687       10228
    ## 688        <NA>
    ## 689        9213
    ## 690       84320
    ## 691        5768
    ## 692        9857
    ## 693       26092
    ## 694      163590
    ## 695      148753
    ## 696      163589
    ## 697        <NA>
    ## 698        6646
    ## 699          27
    ## 700       64222
    ## 701        9917
    ## 702       55103
    ## 703        9068
    ## 704        9462
    ## 705        <NA>
    ## 706       89866
    ## 707       57795
    ## 708         460
    ## 709       60676
    ## 710       64326
    ## 711        7143
    ## 712        <NA>
    ## 713       63931
    ## 714       27101
    ## 715        9910
    ## 716      149041
    ## 717         462
    ## 718       84614
    ## 719       55157
    ## 720       91687
    ## 721       27252
    ## 722      339416
    ## 723        9588
    ## 724       51430
    ## 725        5279
    ## 726       26052
    ## 727        <NA>
    ## 728        8674
    ## 729        4653
    ## 730       23215
    ## 731        2326
    ## 732        2327
    ## 733        5396
    ## 734       92344
    ## 735      149281
    ## 736       22920
    ## 737       57147
    ## 738        <NA>
    ## 739       92342
    ## 740        2153
    ## 741       10560
    ## 742       57821
    ## 743        8548
    ## 744       29922
    ## 745         481
    ## 746      375035
    ## 747      261726
    ## 748       23432
    ## 749       55827
    ## 750       25874
    ## 751        9019
    ## 752       92241
    ## 753        8804
    ## 754        5451
    ## 755      387597
    ## 756        <NA>
    ## 757      117143
    ## 758       57645
    ## 759        <NA>
    ## 760        <NA>
    ## 761      149297
    ## 762        7371
    ## 763       54499
    ## 764         223
    ## 765        4259
    ## 766        6258
    ## 767        5087
    ## 768       83540
    ## 769        8490
    ## 770        5999
    ## 771      339512
    ## 772       51478
    ## 773        <NA>
    ## 774        4921
    ## 775        6675
    ## 776        <NA>
    ## 777      127933
    ## 778        <NA>
    ## 779        9722
    ## 780       25903
    ## 781       22926
    ## 782       11266
    ## 783        <NA>
    ## 784        2213
    ## 785        <NA>
    ## 786      257177
    ## 787        6391
    ## 788      654790
    ## 789       84134
    ## 790         336
    ## 791        2207
    ## 792        4720
    ## 793        9507
    ## 794        8703
    ## 795        5498
    ## 796       27005
    ## 797       51506
    ## 798        <NA>
    ## 799        9191
    ## 800        4817
    ## 801        5202
    ## 802      126823
    ## 803       81607
    ## 804      257106
    ## 805        7391
    ## 806   100131187
    ## 807        <NA>
    ## 808       50848
    ## 809        <NA>
    ## 810         962
    ## 811        6504
    ## 812        8832
    ## 813      114836
    ## 814       57216
    ## 815        4807
    ## 816       23385
    ## 817        1314
    ## 818        5824
    ## 819       50717
    ## 820        <NA>
    ## 821         844
    ## 822        <NA>
    ## 823       93185
    ## 824         477
    ## 825        3765
    ## 826        3766
    ## 827       93183
    ## 828       89886
    ## 829       57549
    ## 830        8407
    ## 831       25790
    ## 832      391123
    ## 833       54935
    ## 834        2532
    ## 835       57863
    ## 836        9447
    ## 837        <NA>
    ## 838        <NA>
    ## 839        <NA>
    ## 840        6708
    ## 841       56776
    ## 842       64388
    ## 843        6000
    ## 844        <NA>
    ## 845       23596
    ## 846        1122
    ## 847        9156
    ## 848      200150
    ## 849        <NA>
    ## 850        <NA>
    ## 851        9859
    ## 852       10806
    ## 853        <NA>
    ## 854       10000
    ## 855       10472
    ## 856        <NA>
    ## 857        <NA>
    ## 858       51029
    ## 859        <NA>
    ## 860      116228
    ## 861        3192
    ## 862       84288
    ## 863       55083
    ## 864       64754
    ## 865       64216
    ## 866      163882
    ## 867       51097
    ## 868       25909
    ## 869        8476
    ## 870       56997
    ## 871        <NA>
    ## 872        5664
    ## 873        3707
    ## 874      375057
    ## 875         142
    ## 876      286826
    ## 877       64746
    ## 878        <NA>
    ## 879        <NA>
    ## 880        <NA>
    ## 881      163859
    ## 882        7044
    ## 883       29920
    ## 884       10637
    ## 885        9725
    ## 886        2052
    ## 887        4931
    ## 888       29097
    ## 889       80232
    ## 890        <NA>
    ## 891      149111
    ## 892      127602
    ## 893        3930
    ## 894       55740
    ## 895        6726
    ## 896        8560
    ## 897        <NA>
    ## 898       23219
    ## 899        <NA>
    ## 900         824
    ## 901        <NA>
    ## 902       55061
    ## 903        <NA>
    ## 904       84976
    ## 905      148362
    ## 906       64853
    ## 907      375056
    ## 908        9015
    ## 909        <NA>
    ## 910       79802
    ## 911       11221
    ## 912        3142
    ## 913       54996
    ## 914        <NA>
    ## 915        4139
    ## 916       25782
    ## 917       55699
    ## 918       10380
    ## 919        <NA>
    ## 920        <NA>
    ## 921        <NA>
    ## 922       55532
    ## 923      127018
    ## 924        7042
    ## 925       51018
    ## 926        <NA>
    ## 927        <NA>
    ## 928       55105
    ## 929        2104
    ## 930       51133
    ## 931        3776
    ## 932        <NA>
    ## 933        1063
    ## 934        5784
    ## 935       56950
    ## 936        5629
    ## 937        <NA>
    ## 938       26750
    ## 939       90806
    ## 940       79805
    ## 941       28982
    ## 942      128387
    ## 943       25936
    ## 944       55509
    ## 945         467
    ## 946        <NA>
    ## 947       29937
    ## 948        <NA>
    ## 949        5525
    ## 950        <NA>
    ## 951       51514
    ## 952       25896
    ## 953        9926
    ## 954        4751
    ## 955        <NA>
    ## 956        7779
    ## 957        <NA>
    ## 958        7188
    ## 959        <NA>
    ## 960       55758
    ## 961        <NA>
    ## 962        3756
    ## 963       55733
    ## 964       56256
    ## 965        <NA>
    ## 966      255928
    ## 967        <NA>
    ## 968        <NA>
    ## 969        3664
    ## 970        <NA>
    ## 971       80342
    ## 972        3914
    ## 973        3290
    ## 974       50486
    ## 975       57172
    ## 976        5362
    ## 977         947
    ## 978        <NA>
    ## 979        4179
    ## 980        1379
    ## 981      221061
    ## 982        9397
    ## 983       10557
    ## 984      414149
    ## 985      644890
    ## 986       64421
    ## 987       79723
    ## 988       51182
    ## 989        <NA>
    ## 990      441549
    ## 991       83641
    ## 992       55691
    ## 993        8559
    ## 994      222389
    ## 995       22929
    ## 996        5264
    ## 997      221044
    ## 998       55388
    ## 999       10133
    ## 1000      83643
    ## 1001      57118
    ## 1002       8872
    ## 1003      11164
    ## 1004       <NA>
    ## 1005      55176
    ## 1006      55526
    ## 1007      26019
    ## 1008     254427
    ## 1009      79746
    ## 1010       9712
    ## 1011       <NA>
    ## 1012      10659
    ## 1013      83860
    ## 1014       <NA>
    ## 1015      22944
    ## 1016       3698
    ## 1017      80760
    ## 1018      57713
    ## 1019       5588
    ## 1020       <NA>
    ## 1021       5209
    ## 1022      84991
    ## 1023       3601
    ## 1024       <NA>
    ## 1025      54522
    ## 1026       8516
    ## 1027       <NA>
    ## 1028      80013
    ## 1029       9317
    ## 1030     389941
    ## 1031       6251
    ## 1032       1787
    ## 1033       7431
    ## 1034     338596
    ## 1035       9200
    ## 1036       8027
    ## 1037     653567
    ## 1038       4360
    ## 1039     221074
    ## 1040        783
    ## 1041     221078
    ## 1042     221079
    ## 1043      84898
    ## 1044      10529
    ## 1045       <NA>
    ## 1046     387640
    ## 1047       8028
    ## 1048      64215
    ## 1049      23412
    ## 1050        648
    ## 1051       9576
    ## 1052       <NA>
    ## 1053       5305
    ## 1054       <NA>
    ## 1055     219681
    ## 1056      22921
    ## 1057       <NA>
    ## 1058     220213
    ## 1059  100144434
    ## 1060       <NA>
    ## 1061      57584
    ## 1062       <NA>
    ## 1063     219670
    ## 1064      79896
    ## 1065      57512
    ## 1066       2572
    ## 1067      54518
    ## 1068      23590
    ## 1069      10006
    ## 1070      91452
    ## 1071      84930
    ## 1072      10730
    ## 1073       <NA>
    ## 1074      11249
    ## 1075     339745
    ## 1076       3176
    ## 1077       <NA>
    ## 1078      23550
    ## 1079        774
    ## 1080      79813
    ## 1081      92714
    ## 1082     116225
    ## 1083      92715
    ## 1084      64975
    ## 1085     375775
    ## 1086      26012
    ## 1087     441478
    ## 1088      54863
    ## 1089      25920
    ## 1090      10383
    ## 1091     142680
    ## 1092     727800
    ## 1093       <NA>
    ## 1094      27158
    ## 1095      94107
    ## 1096     286262
    ## 1097       8636
    ## 1098      29882
    ## 1099       2902
    ## 1100      11253
    ## 1101      29952
    ## 1102      91373
    ## 1103      89958
    ## 1104        954
    ## 1105      56654
    ## 1106       2529
    ## 1107         20
    ## 1108     286257
    ## 1109       5730
    ## 1110        733
    ## 1111      54461
    ## 1112       7186
    ## 1113       8721
    ## 1114     158056
    ## 1115      29085
    ## 1116       <NA>
    ## 1117      55684
    ## 1118       <NA>
    ## 1119      85014
    ## 1120       <NA>
    ## 1121       <NA>
    ## 1122       <NA>
    ## 1123      57582
    ## 1124     157922
    ## 1125      10422
    ## 1126     138151
    ## 1127       <NA>
    ## 1128     169714
    ## 1129     399693
    ## 1130      26086
    ## 1131     728489
    ## 1132      64170
    ## 1133       6621
    ## 1134       <NA>
    ## 1135      23203
    ## 1136      56623
    ## 1137       9919
    ## 1138       4851
    ## 1139      51162
    ## 1140      10555
    ## 1141       <NA>
    ## 1142       6838
    ## 1143       6837
    ## 1144       6130
    ## 1145       6834
    ## 1146       6835
    ## 1147       6836
    ## 1148      57109
    ## 1149      11093
    ## 1150      11094
    ## 1151      11182
    ## 1152     389827
    ## 1153       9719
    ## 1154     642968
    ## 1155       1621
    ## 1156       1757
    ## 1157       7410
    ## 1158     266655
    ## 1159       8019
    ## 1160      11091
    ## 1161       6256
    ## 1162       1289
    ## 1163      10439
    ## 1164       <NA>
    ## 1165       9858
    ## 1166       <NA>
    ## 1167      51116
    ## 1168       5900
    ## 1169       9328
    ## 1170       7248
    ## 1171     158067
    ## 1172       9329
    ## 1173      64794
    ## 1174       7270
    ## 1175      23064
    ## 1176      84628
    ## 1177       9442
    ## 1178       2889
    ## 1179      26995
    ## 1180      51117
    ## 1181      10999
    ## 1182      81605
    ## 1183       <NA>
    ## 1184      51148
    ## 1185       4957
    ## 1186       2733
    ## 1187       6709
    ## 1188      89891
    ## 1189       6418
    ## 1190      29941
    ## 1191      84885
    ## 1192      10444
    ## 1193      54662
    ## 1194       2021
    ## 1195      51490
    ## 1196        883
    ## 1197       <NA>
    ## 1198      56262
    ## 1199     254295
    ## 1200      22845
    ## 1201      23511
    ## 1202      56904
    ## 1203      84895
    ## 1204      57171
    ## 1205       1384
    ## 1206       5524
    ## 1207     389792
    ## 1208       <NA>
    ## 1209       <NA>
    ## 1210      28989
    ## 1211     140459
    ## 1212      51450
    ## 1213       9536
    ## 1214      27348
    ## 1215       1861
    ## 1216       <NA>
    ## 1217      10868
    ## 1218      23048
    ## 1219      57720
    ## 1220      23413
    ## 1221     256158
    ## 1222        445
    ## 1223       8939
    ## 1224      59335
    ## 1225      23404
    ## 1226         25
    ## 1227      84929
    ## 1228      10319
    ## 1229      83543
    ## 1230       8021
    ## 1231     286336
    ## 1232      84814
    ## 1233      84726
    ## 1234       <NA>
    ## 1235      10585
    ## 1236      83549
    ## 1237     375757
    ## 1238       2801
    ## 1239       1759
    ## 1240      25792
    ## 1241       <NA>
    ## 1242       3934
    ## 1243      80142
    ## 1244     114789
    ## 1245     203245
    ## 1246     399665
    ## 1247       8818
    ## 1248       <NA>
    ## 1249     138429
    ## 1250      27090
    ## 1251      30815
    ## 1252        203
    ## 1253       2022
    ## 1254       2356
    ## 1255       1025
    ## 1256      10044
    ## 1257       <NA>
    ## 1258     158248
    ## 1259      27433
    ## 1260     138428
    ## 1261     286207
    ## 1262       6812
    ## 1263       <NA>
    ## 1264       <NA>
    ## 1265      90678
    ## 1266       6136
    ## 1267      29988
    ## 1268      84253
    ## 1269       9649
    ## 1270      23452
    ## 1271     403341
    ## 1272      23099
    ## 1273       <NA>
    ## 1274      89853
    ## 1275       <NA>
    ## 1276       5090
    ## 1277      79109
    ## 1278      26130
    ## 1279       3309
    ## 1280      10244
    ## 1281      26190
    ## 1282       5711
    ## 1283      26147
    ## 1284       7185
    ## 1285       <NA>
    ## 1286      11064
    ## 1287      51552
    ## 1288       2934
    ## 1289       2040
    ## 1290       <NA>
    ## 1291     153090
    ## 1292     158135
    ## 1293       4702
    ## 1294     254956
    ## 1295      26468
    ## 1296      92400
    ## 1297      92399
    ## 1298       5742
    ## 1299       5082
    ## 1300      54542
    ## 1301      10773
    ## 1302      57684
    ## 1303      23637
    ## 1304      55342
    ## 1305       2844
    ## 1306     286204
    ## 1307      57706
    ## 1308       9355
    ## 1309       <NA>
    ## 1310       <NA>
    ## 1311      10783
    ## 1312       5695
    ## 1313       2649
    ## 1314     169611
    ## 1315       <NA>
    ## 1316     401551
    ## 1317      11224
    ## 1318      81873
    ## 1319       2800
    ## 1320     286205
    ## 1321       5537
    ## 1322      53353
    ## 1323      55843
    ## 1324      79712
    ## 1325       9839
    ## 1326       <NA>
    ## 1327       <NA>
    ## 1328         92
    ## 1329       5000
    ## 1330      55777
    ## 1331      26122
    ## 1332       3800
    ## 1333     130576
    ## 1334       <NA>
    ## 1335     130574
    ## 1336      27249
    ## 1337        390
    ## 1338     375287
    ## 1339       9111
    ## 1340       7130
    ## 1341      55183
    ## 1342       4703
    ## 1343      26225
    ## 1344        785
    ## 1345      10254
    ## 1346     114793
    ## 1347      55660
    ## 1348     151188
    ## 1349      56475
    ## 1350     114805
    ## 1351       <NA>
    ## 1352       3760
    ## 1353       4929
    ## 1354       2820
    ## 1355       <NA>
    ## 1356      57471
    ## 1357       9595
    ## 1358     130399
    ## 1359         90
    ## 1360     151531
    ## 1361     130940
    ## 1362       8502
    ## 1363      92196
    ## 1364      85461
    ## 1365     151525
    ## 1366      29994
    ## 1367       <NA>
    ## 1368       9936
    ## 1369       4065
    ## 1370      22925
    ## 1371       5937
    ## 1372      10010
    ## 1373      10213
    ## 1374      10716
    ## 1375      57282
    ## 1376       1803
    ## 1377       <NA>
    ## 1378       2191
    ## 1379      64135
    ## 1380      25801
    ## 1381      90134
    ## 1382      55137
    ## 1383       2888
    ## 1384      22837
    ## 1385       <NA>
    ## 1386       6328
    ## 1387       6326
    ## 1388      80034
    ## 1389       2591
    ## 1390      79809
    ## 1391       6323
    ## 1392       <NA>
    ## 1393       6335
    ## 1394       6332
    ## 1395       <NA>
    ## 1396       8708
    ## 1397      27347
    ## 1398     253782
    ## 1399     115677
    ## 1400      57405
    ## 1401     129880
    ## 1402      10324
    ## 1403      79675
    ## 1404       9360
    ## 1405     129881
    ## 1406     493911
    ## 1407     151230
    ## 1408       6741
    ## 1409      29081
    ## 1410       <NA>
    ## 1411     130507
    ## 1412     140469
    ## 1413     389058
    ## 1414     285141
    ## 1415       2571
    ## 1416       <NA>
    ## 1417      26003
    ## 1418       9874
    ## 1419      79828
    ## 1420      80067
    ## 1421      79901
    ## 1422       1781
    ## 1423       8604
    ## 1424       8520
    ## 1425     254042
    ## 1426       <NA>
    ## 1427       1745
    ## 1428       1746
    ## 1429       3655
    ## 1430       5163
    ## 1431      11069
    ## 1432      51776
    ## 1433      83879
    ## 1434       6670
    ## 1435       <NA>
    ## 1436      29789
    ## 1437  100131390
    ## 1438       9541
    ## 1439      79634
    ## 1440     151556
    ## 1441       7456
    ## 1442       1134
    ## 1443       <NA>
    ## 1444       <NA>
    ## 1445       1123
    ## 1446       1386
    ## 1447       <NA>
    ## 1448      80856
    ## 1449      10651
    ## 1450       <NA>
    ## 1451     220988
    ## 1452       4780
    ## 1453       8540
    ## 1454     150737
    ## 1455       <NA>
    ## 1456      50940
    ## 1457       <NA>
    ## 1458     129831
    ## 1459     114880
    ## 1460       8575
    ## 1461      51661
    ## 1462      65977
    ## 1463       7273
    ## 1464     285025
    ## 1465      91404
    ## 1466       <NA>
    ## 1467      57703
    ## 1468      10477
    ## 1469       3676
    ## 1470     375298
    ## 1471       4760
    ## 1472       <NA>
    ## 1473       5136
    ## 1474      54431
    ## 1475       2487
    ## 1476      10787
    ## 1477     142679
    ## 1478     129401
    ## 1479       <NA>
    ## 1480      55854
    ## 1481       3685
    ## 1482     165215
    ## 1483      10203
    ## 1484       7035
    ## 1485       <NA>
    ## 1486       1500
    ## 1487     280636
    ## 1488      51075
    ## 1489     219541
    ## 1490      25921
    ## 1491      10978
    ## 1492     219539
    ## 1493       <NA>
    ## 1494        710
    ## 1495       9246
    ## 1496      26519
    ## 1497      29015
    ## 1498       6749
    ## 1499      85456
    ## 1500        187
    ## 1501     219527
    ## 1502       <NA>
    ## 1503       <NA>
    ## 1504       5795
    ## 1505      23279
    ## 1506      23360
    ## 1507      79841
    ## 1508      23788
    ## 1509     114900
    ## 1510       4722
    ## 1511      55709
    ## 1512      10658
    ## 1513       5913
    ## 1514       5702
    ## 1515      91252
    ## 1516       6688
    ## 1517       4607
    ## 1518       8567
    ## 1519      10062
    ## 1520         53
    ## 1521       1643
    ## 1522       <NA>
    ## 1523       <NA>
    ## 1524      29763
    ## 1525      84364
    ## 1526       <NA>
    ## 1527       4038
    ## 1528       9793
    ## 1529       <NA>
    ## 1530        392
    ## 1531       9776
    ## 1532     283254
    ## 1533      55626
    ## 1534       4192
    ## 1535       8525
    ## 1536      90993
    ## 1537       <NA>
    ## 1538      51317
    ## 1539       9409
    ## 1540       <NA>
    ## 1541       9479
    ## 1542       1408
    ## 1543      55343
    ## 1544       8534
    ## 1545      57586
    ## 1546      56981
    ## 1547       <NA>
    ## 1548      90139
    ## 1549       3732
    ## 1550       2132
    ## 1551      84680
    ## 1552     390110
    ## 1553       <NA>
    ## 1554     221120
    ## 1555      51144
    ## 1556  100507261
    ## 1557      55761
    ## 1558       8539
    ## 1559      57689
    ## 1560       <NA>
    ## 1561       7189
    ## 1562      79899
    ## 1563      29099
    ## 1564     143458
    ## 1565      54765
    ## 1566      24147
    ## 1567      25891
    ## 1568       6506
    ## 1569        960
    ## 1570       8050
    ## 1571      51074
    ## 1572       2001
    ## 1573        847
    ## 1574      25841
    ## 1575      55226
    ## 1576       4076
    ## 1577       4005
    ## 1578      26273
    ## 1579       <NA>
    ## 1580       <NA>
    ## 1581       <NA>
    ## 1582       <NA>
    ## 1583       <NA>
    ## 1584      10114
    ## 1585       1479
    ## 1586      55346
    ## 1587      91614
    ## 1588      79832
    ## 1589     493860
    ## 1590      10480
    ## 1591       5954
    ## 1592       5080
    ## 1593      26610
    ## 1594     196294
    ## 1595     120526
    ## 1596       <NA>
    ## 1597        744
    ## 1598     120534
    ## 1599       3739
    ## 1600     196074
    ## 1601      81930
    ## 1602        627
    ## 1603      55327
    ## 1604      55366
    ## 1605      91057
    ## 1606       8424
    ## 1607     387758
    ## 1608      63982
    ## 1609     143662
    ## 1610     254531
    ## 1611      55505
    ## 1612       9990
    ## 1613      51234
    ## 1614      79768
    ## 1615      56851
    ## 1616      57099
    ## 1617       6263
    ## 1618     342184
    ## 1619      26585
    ## 1620       6447
    ## 1621       9824
    ## 1622      57369
    ## 1623       9716
    ## 1624       <NA>
    ## 1625      89978
    ## 1626       <NA>
    ## 1627       <NA>
    ## 1628       4212
    ## 1629       <NA>
    ## 1630       <NA>
    ## 1631       <NA>
    ## 1632       <NA>
    ## 1633       <NA>
    ## 1634     161742
    ## 1635     283742
    ## 1636      10125
    ## 1637       7057
    ## 1638     161835
    ## 1639      11245
    ## 1640     440275
    ## 1641       6727
    ## 1642      90427
    ## 1643        701
    ## 1644      56924
    ## 1645  100131244
    ## 1646       5330
    ## 1647  100505573
    ## 1648       <NA>
    ## 1649      85455
    ## 1650      90417
    ## 1651       3712
    ## 1652      22893
    ## 1653     113189
    ## 1654      90416
    ## 1655      27079
    ## 1656       5888
    ## 1657      55177
    ## 1658       2644
    ## 1659      55192
    ## 1660      84936
    ## 1661       6692
    ## 1662     171177
    ## 1663      57617
    ## 1664      54567
    ## 1665      79094
    ## 1666      54617
    ## 1667     161829
    ## 1668      11261
    ## 1669       <NA>
    ## 1670      11339
    ## 1671      51203
    ## 1672      51103
    ## 1673      23168
    ## 1674       3706
    ## 1675       4058
    ## 1676      26015
    ## 1677       7301
    ## 1678      23269
    ## 1679      23005
    ## 1680  100137047
    ## 1681      30844
    ## 1682     123745
    ## 1683      23339
    ## 1684      25963
    ## 1685       2595
    ## 1686        825
    ## 1687       <NA>
    ## 1688       8773
    ## 1689     255252
    ## 1690      55142
    ## 1691      57519
    ## 1692     146059
    ## 1693     146057
    ## 1694     197131
    ## 1695      80021
    ## 1696      23582
    ## 1697       2038
    ## 1698       9836
    ## 1699     161823
    ## 1700     146050
    ## 1701      27229
    ## 1702       <NA>
    ## 1703       4130
    ## 1704       9677
    ## 1705       <NA>
    ## 1706     117155
    ## 1707       2923
    ## 1708     619189
    ## 1709      10169
    ## 1710      25764
    ## 1711       <NA>
    ## 1712       <NA>
    ## 1713      79968
    ## 1714      84978
    ## 1715       <NA>
    ## 1716     113201
    ## 1717      51496
    ## 1718       <NA>
    ## 1719      80208
    ## 1720        567
    ## 1721       6652
    ## 1722      90525
    ## 1723       2628
    ## 1724       <NA>
    ## 1725       7782
    ## 1726      26258
    ## 1727      58472
    ## 1728      80031
    ## 1729     283652
    ## 1730      50804
    ## 1731     399697
    ## 1732       1854
    ## 1733       2200
    ## 1734      22995
    ## 1735     399694
    ## 1736      23741
    ## 1737       9728
    ## 1738       9318
    ## 1739       2585
    ## 1740       2252
    ## 1741      56986
    ## 1742      11001
    ## 1743       3067
    ## 1744       2553
    ## 1745       <NA>
    ## 1746       9101
    ## 1747     373509
    ## 1748      54822
    ## 1749      84888
    ## 1750      23431
    ## 1751        644
    ## 1752      23397
    ## 1753     150771
    ## 1754       <NA>
    ## 1755      23020
    ## 1756       9391
    ## 1757      55654
    ## 1758      56910
    ## 1759       1844
    ## 1760     150763
    ## 1761      51011
    ## 1762      30818
    ## 1763       <NA>
    ## 1764      64969
    ## 1765       4118
    ## 1766       7851
    ## 1767       4867
    ## 1768       <NA>
    ## 1769      10018
    ## 1770       <NA>
    ## 1771      64682
    ## 1772      10461
    ## 1773      84910
    ## 1774     129804
    ## 1775      84524
    ## 1776       <NA>
    ## 1777     376940
    ## 1778     150465
    ## 1779      84172
    ## 1780      84269
    ## 1781       <NA>
    ## 1782       6574
    ## 1783       <NA>
    ## 1784     150468
    ## 1785       3552
    ## 1786       3553
    ## 1787     140885
    ## 1788       5173
    ## 1789     140901
    ## 1790       7053
    ## 1791       6628
    ## 1792      10528
    ## 1793       3420
    ## 1794      57593
    ## 1795      56265
    ## 1796      64773
    ## 1797      64601
    ## 1798       5786
    ## 1799      64949
    ## 1800       5020
    ## 1801        551
    ## 1802      22888
    ## 1803      60493
    ## 1804       9762
    ## 1805      65992
    ## 1806       3704
    ## 1807      83959
    ## 1808       <NA>
    ## 1809       8455
    ## 1810      64096
    ## 1811      80332
    ## 1812     116835
    ## 1813       <NA>
    ## 1814      25876
    ## 1815       1059
    ## 1816        994
    ## 1817      55317
    ## 1818      57506
    ## 1819      80025
    ## 1820      11237
    ## 1821      54498
    ## 1822       <NA>
    ## 1823       <NA>
    ## 1824       <NA>
    ## 1825       5621
    ## 1826       9770
    ## 1827       9962
    ## 1828      29058
    ## 1829       5111
    ## 1830       <NA>
    ## 1831       8760
    ## 1832     128674
    ## 1833      56261
    ## 1834       <NA>
    ## 1835       1114
    ## 1836      51605
    ## 1837      84515
    ## 1838      54675
    ## 1839        650
    ## 1840      56255
    ## 1841      23236
    ## 1842       5332
    ## 1843      24141
    ## 1844       <NA>
    ## 1845      63926
    ## 1846       6616
    ## 1847       8195
    ## 1848     128710
    ## 1849        182
    ## 1850      22903
    ## 1851     140862
    ## 1852      55617
    ## 1853      51575
    ## 1854      79133
    ## 1855     140733
    ## 1856      23767
    ## 1857      55614
    ## 1858       6629
    ## 1859       <NA>
    ## 1860       5126
    ## 1861        631
    ## 1862      11034
    ## 1863       6238
    ## 1864      27131
    ## 1865      92667
    ## 1866      58495
    ## 1867       <NA>
    ## 1868      57325
    ## 1869      55184
    ## 1870      10621
    ## 1871      10741
    ## 1872      10483
    ## 1873     388789
    ## 1874      92675
    ## 1875      57419
    ## 1876      54453
    ## 1877      51126
    ## 1878      51340
    ## 1879      26074
    ## 1880       3642
    ## 1881      57186
    ## 1882      55857
    ## 1883      22803
    ## 1884       4821
    ## 1885       <NA>
    ## 1886       <NA>
    ## 1887       7056
    ## 1888      22918
    ## 1889      29107
    ## 1890      64412
    ## 1891      63908
    ## 1892       1471
    ## 1893      79953
    ## 1894       <NA>
    ## 1895       <NA>
    ## 1896       <NA>
    ## 1897       <NA>
    ## 1898       8530
    ## 1899      57136
    ## 1900      84532
    ## 1901        955
    ## 1902       5834
    ## 1903      26090
    ## 1904       9837
    ## 1905      22981
    ## 1906     140838
    ## 1907      55968
    ## 1908       2280
    ## 1909      27111
    ## 1910       9751
    ## 1911       <NA>
    ## 1912      55321
    ## 1913       9491
    ## 1914      83541
    ## 1915     113278
    ## 1916       <NA>
    ## 1917      85508
    ## 1918     140809
    ## 1919       6939
    ## 1920       1457
    ## 1921     128637
    ## 1922      10616
    ## 1923      57761
    ## 1924      80023
    ## 1925       6666
    ## 1926      85364
    ## 1927       <NA>
    ## 1928      28954
    ## 1929       <NA>
    ## 1930       <NA>
    ## 1931       <NA>
    ## 1932       <NA>
    ## 1933       3397
    ## 1934      84701
    ## 1935        598
    ## 1936      22974
    ## 1937       2307
    ## 1938     128853
    ## 1939     164395
    ## 1940      81572
    ## 1941     140706
    ## 1942       3055
    ## 1943       9777
    ## 1944       <NA>
    ## 1945       5326
    ## 1946      23509
    ## 1947       9371
    ## 1948     171023
    ## 1949     140688
    ## 1950     149951
    ## 1951       1789
    ## 1952      22919
    ## 1953     149954
    ## 1954      92747
    ## 1955      51654
    ## 1956       6640
    ## 1957       9139
    ## 1958      63941
    ## 1959       <NA>
    ## 1960       1869
    ## 1961      11264
    ## 1962       <NA>
    ## 1963     128866
    ## 1964      22913
    ## 1965       8894
    ## 1966        191
    ## 1967      83737
    ## 1968      83658
    ## 1969      84557
    ## 1970     128869
    ## 1971       <NA>
    ## 1972      23054
    ## 1973       2686
    ## 1974      55902
    ## 1975       2937
    ## 1976      57644
    ## 1977      26133
    ## 1978      55741
    ## 1979      10544
    ## 1980      10893
    ## 1981       <NA>
    ## 1982       3692
    ## 1983       <NA>
    ## 1984      55245
    ## 1985       8200
    ## 1986      11190
    ## 1987       <NA>
    ## 1988      51614
    ## 1989       6676
    ## 1990       8904
    ## 1991      10137
    ## 1992       9054
    ## 1993     140823
    ## 1994       9584
    ## 1995      51230
    ## 1996      51282
    ## 1997     140894
    ## 1998       <NA>
    ## 1999       2036
    ## 2000      25980
    ## 2001       <NA>
    ## 2002      22839
    ## 2003      10398
    ## 2004      60436
    ## 2005       <NA>
    ## 2006      84174
    ## 2007      57446
    ## 2008      79980
    ## 2009     140710
    ## 2010      25939
    ## 2011       5933
    ## 2012     140699
    ## 2013       6185
    ## 2014       2691
    ## 2015      63905
    ## 2016       6714
    ## 2017      10904
    ## 2018       4826
    ## 2019      56259
    ## 2020     128434
    ## 2021       9675
    ## 2022      58490
    ## 2023       7052
    ## 2024       <NA>
    ## 2025       3929
    ## 2026     128439
    ## 2027      57148
    ## 2028     140679
    ## 2029      79913
    ## 2030      26051
    ## 2031      81610
    ## 2032      60625
    ## 2033       9935
    ## 2034       7150
    ## 2035       5335
    ## 2036      23051
    ## 2037      64900
    ## 2038      84181
    ## 2039       <NA>
    ## 2040      11122
    ## 2041       6431
    ## 2042      26013
    ## 2043      10110
    ## 2044      51098
    ## 2045       4605
    ## 2046      84969
    ## 2047      57158
    ## 2048      51526
    ## 2049       <NA>
    ## 2050      78997
    ## 2051     128486
    ## 2052      79183
    ## 2053      10955
    ## 2054      11142
    ## 2055        100
    ## 2056      60598
    ## 2057     140730
    ## 2058       7529
    ## 2059      80336
    ## 2060      10953
    ## 2061       6789
    ## 2062       3787
    ## 2063       6590
    ## 2064       8785
    ## 2065      11317
    ## 2066       6385
    ## 2067      90196
    ## 2068       <NA>
    ## 2069      55861
    ## 2070      51604
    ## 2071       <NA>
    ## 2072      10406
    ## 2073     140686
    ## 2074     116092
    ## 2075      11065
    ## 2076       7125
    ## 2077      90203
    ## 2078      10005
    ## 2079     140831
    ## 2080      90204
    ## 2081     140825
    ## 2082       5476
    ## 2083       5360
    ## 2084      63935
    ## 2085       <NA>
    ## 2086       <NA>
    ## 2087       4318
    ## 2088      57468
    ## 2089      57727
    ## 2090        958
    ## 2091      64405
    ## 2092      51006
    ## 2093      63916
    ## 2094       <NA>
    ## 2095      64849
    ## 2096       <NA>
    ## 2097      81031
    ## 2098       2139
    ## 2099      23613
    ## 2100       8202
    ## 2101      55959
    ## 2102       <NA>
    ## 2103      57580
    ## 2104       <NA>
    ## 2105      10564
    ## 2106       1434
    ## 2107       6780
    ## 2108      55661
    ## 2109      57169
    ## 2110       3745
    ## 2111       5740
    ## 2112       9334
    ## 2113      23315
    ## 2114       9825
    ## 2115      55905
    ## 2116       6615
    ## 2117       7335
    ## 2118     387521
    ## 2119       1051
    ## 2120       5770
    ## 2121     140876
    ## 2122      84612
    ## 2123       8813
    ## 2124      27304
    ## 2125       3755
    ## 2126       4773
    ## 2127      10079
    ## 2128      55734
    ## 2129     128553
    ## 2130       <NA>
    ## 2131       8537
    ## 2132       5203
    ## 2133       <NA>
    ## 2134      55816
    ## 2135     140689
    ## 2136     116151
    ## 2137       6790
    ## 2138       1477
    ## 2139      57091
    ## 2140       <NA>
    ## 2141      51507
    ## 2142       7022
    ## 2143        655
    ## 2144       8480
    ## 2145      55544
    ## 2146      81030
    ## 2147      56937
    ## 2148      57403
    ## 2149       9217
    ## 2150       8675
    ## 2151      79716
    ## 2152       2778
    ## 2153       <NA>
    ## 2154      51497
    ## 2155       1522
    ## 2156       <NA>
    ## 2157      51012
    ## 2158       <NA>
    ## 2159       1908
    ## 2160       <NA>
    ## 2161       <NA>
    ## 2162       <NA>
    ## 2163       <NA>
    ## 2164       <NA>
    ## 2165       <NA>
    ## 2166       <NA>
    ## 2167       <NA>
    ## 2168       <NA>
    ## 2169       <NA>
    ## 2170       <NA>
    ## 2171       <NA>
    ## 2172       <NA>
    ## 2173       <NA>
    ## 2174       <NA>
    ## 2175       <NA>
    ## 2176       <NA>
    ## 2177       <NA>
    ## 2178       <NA>
    ## 2179       <NA>
    ## 2180       <NA>
    ## 2181     116154
    ## 2182      10388
    ## 2183       5509
    ## 2184      63939
    ## 2185       1002
    ## 2186       6874
    ## 2187     149986
    ## 2188       5688
    ## 2189      26039
    ## 2190      26164
    ## 2191      11255
    ## 2192       9885
    ## 2193      11047
    ## 2194       3911
    ## 2195       6227
    ## 2196      81928
    ## 2197       <NA>
    ## 2198      28231
    ## 2199       4923
    ## 2200      55257
    ## 2201      11054
    ## 2202       1299
    ## 2203      10732
    ## 2204      11083
    ## 2205      54994
    ## 2206      63910
    ## 2207       <NA>
    ## 2208       <NA>
    ## 2209      54915
    ## 2210     128414
    ## 2211      55738
    ## 2212      57642
    ## 2213       1137
    ## 2214       3785
    ## 2215       1917
    ## 2216      79144
    ## 2217      85441
    ## 2218      26205
    ## 2219      50861
    ## 2220      51750
    ## 2221      10139
    ## 2222      84619
    ## 2223      54923
    ## 2224     140685
    ## 2225       7165
    ## 2226      80331
    ## 2227      54963
    ## 2228       <NA>
    ## 2229      24148
    ## 2230     140700
    ## 2231      54345
    ## 2232       6919
    ## 2233      10287
    ## 2234       4987
    ## 2235       4661
    ## 2236      55251
    ## 2237      51728
    ## 2238      55190
    ## 2239     170685
    ## 2240      57477
    ## 2241       1184
    ## 2242     389856
    ## 2243      89801
    ## 2244  110437700
    ## 2245      28952
    ## 2246        778
    ## 2247       6855
    ## 2248       <NA>
    ## 2249       4007
    ## 2250       5355
    ## 2251      79917
    ## 2252      27238
    ## 2253      11152
    ## 2254      11230
    ## 2255      90060
    ## 2256       7030
    ## 2257      56850
    ## 2258       3750
    ## 2259      55593
    ## 2260      11040
    ## 2261       7355
    ## 2262      10084
    ## 2263      10245
    ## 2264      27344
    ## 2265      10013
    ## 2266     392465
    ## 2267       6839
    ## 2268      64743
    ## 2269       5935
    ## 2270       <NA>
    ## 2271       4943
    ## 2272      10682
    ## 2273      64840
    ## 2274      24140
    ## 2275      92745
    ## 2276       <NA>
    ## 2277       7504
    ## 2278       1536
    ## 2279       6990
    ## 2280       8406
    ## 2281       6103
    ## 2282       7102
    ## 2283      58526
    ## 2284      54880
    ## 2285      10159
    ## 2286       <NA>
    ## 2287       9282
    ## 2288       8239
    ## 2289       <NA>
    ## 2290       1654
    ## 2291       8573
    ## 2292       2857
    ## 2293       4128
    ## 2294       4129
    ## 2295       4693
    ## 2296      80258
    ## 2297     139341
    ## 2298       7403
    ## 2299      56548
    ## 2300      84679
    ## 2301       6102
    ## 2302       9767
    ## 2303      54539
    ## 2304       8241
    ## 2305       7317
    ## 2306       5127
    ## 2307       8237
    ## 2308        369
    ## 2309       6853
    ## 2310       7076
    ## 2311       5199
    ## 2312       2002
    ## 2313       8409
    ## 2314       <NA>
    ## 2315       <NA>
    ## 2316        186
    ## 2317      90293
    ## 2318      54521
    ## 2319     139818
    ## 2320       3597
    ## 2321     170261
    ## 2322      79836
    ## 2323      10857
    ## 2324       <NA>
    ## 2325        292
    ## 2326       <NA>
    ## 2327       <NA>
    ## 2328       7319
    ## 2329      55922
    ## 2330       <NA>
    ## 2331       6170
    ## 2332      65109
    ## 2333      79576
    ## 2334       4694
    ## 2335       <NA>
    ## 2336      10009
    ## 2337      55026
    ## 2338       3920
    ## 2339       8450
    ## 2340      28985
    ## 2341      29071
    ## 2342       2892
    ## 2343      57187
    ## 2344        331
    ## 2345      10735
    ## 2346       4068
    ## 2347      10178
    ## 2348     340578
    ## 2349     139170
    ## 2350       6594
    ## 2351       4952
    ## 2352       8862
    ## 2353      54440
    ## 2354      51114
    ## 2355      10813
    ## 2356      63035
    ## 2357       2000
    ## 2358       9131
    ## 2359       9363
    ## 2360       <NA>
    ## 2361       9016
    ## 2362      51634
    ## 2363      10495
    ## 2364     158763
    ## 2365       3547
    ## 2366      51765
    ## 2367      57826
    ## 2368      90161
    ## 2369       2239
    ## 2370       2719
    ## 2371     347475
    ## 2372      84295
    ## 2373       <NA>
    ## 2374     159090
    ## 2375      56180
    ## 2376       8933
    ## 2377      26071
    ## 2378     441518
    ## 2379       <NA>
    ## 2380       <NA>
    ## 2381       <NA>
    ## 2382       <NA>
    ## 2383       <NA>
    ## 2384     399668
    ## 2385     203522
    ## 2386      93380
    ## 2387      10479
    ## 2388       2273
    ## 2389      27336
    ## 2390       9459
    ## 2391      27316
    ## 2392      83550
    ## 2393       7547
    ## 2394       2258
    ## 2395       4168
    ## 2396     286410
    ## 2397       6658
    ## 2398       <NA>
    ## 2399      23641
    ## 2400     139065
    ## 2401      84631
    ## 2402       2332
    ## 2403       2334
    ## 2404       3423
    ## 2405       <NA>
    ## 2406      84548
    ## 2407      10046
    ## 2408       4534
    ## 2409       8776
    ## 2410      83692
    ## 2411       3149
    ## 2412     203547
    ## 2413       <NA>
    ## 2414      79057
    ## 2415       1260
    ## 2416       2556
    ## 2417      55879
    ## 2418       1069
    ## 2419      50814
    ## 2420       <NA>
    ## 2421     114824
    ## 2422      29944
    ## 2423       <NA>
    ## 2424       <NA>
    ## 2425       <NA>
    ## 2426       <NA>
    ## 2427     139735
    ## 2428      55559
    ## 2429        633
    ## 2430        492
    ## 2431       1852
    ## 2432     139728
    ## 2433       6535
    ## 2434      10134
    ## 2435        215
    ## 2436       5365
    ## 2437      26576
    ## 2438       3421
    ## 2439       6748
    ## 2440      57595
    ## 2441       3897
    ## 2442        393
    ## 2443       8260
    ## 2444       5973
    ## 2445       3054
    ## 2446       3654
    ## 2447       4204
    ## 2448       2652
    ## 2449       2316
    ## 2450       2010
    ## 2451       6134
    ## 2452       1774
    ## 2453       6901
    ## 2454        537
    ## 2455       2664
    ## 2456       9130
    ## 2457      55558
    ## 2458       8270
    ## 2459       8266
    ## 2460       8273
    ## 2461      60343
    ## 2462       8517
    ## 2463       <NA>
    ## 2464     139716
    ## 2465       1736
    ## 2466       4354
    ## 2467       2157
    ## 2468      65991
    ## 2469  100272147
    ## 2470       4515
    ## 2471      79184
    ## 2472       7411
    ## 2473     116442
    ## 2474       5358
    ## 2475       6907
    ## 2476       5613
    ## 2477       5638
    ## 2478       <NA>
    ## 2479      83604
    ## 2480       1756
    ## 2481     257397
    ## 2482       2710
    ## 2483       <NA>
    ## 2484      11141
    ## 2485     170302
    ## 2486       5422
    ## 2487       9468
    ## 2488       5165
    ## 2489       7543
    ## 2490       <NA>
    ## 2491      80311
    ## 2492      79135
    ## 2493       <NA>
    ## 2494       9500
    ## 2495      23708
    ## 2496     158586
    ## 2497     139886
    ## 2498      23229
    ## 2499     139285
    ## 2500      55906
    ## 2501     340554
    ## 2502      81887
    ## 2503       4478
    ## 2504       9843
    ## 2505        367
    ## 2506       4983
    ## 2507     286451
    ## 2508       9754
    ## 2509       1947
    ## 2510      64219
    ## 2511       <NA>
    ## 2512       1896
    ## 2513       3476
    ## 2514     347516
    ## 2515        407
    ## 2516      51248
    ## 2517       <NA>
    ## 2518      54857
    ## 2519       1741
    ## 2520      84889
    ## 2521      29934
    ## 2522       4303
    ## 2523       3561
    ## 2524       9968
    ## 2525      54413
    ## 2526       2705
    ## 2527       9203
    ## 2528       4841
    ## 2529      26548
    ## 2530       6872
    ## 2531       8473
    ## 2532       <NA>
    ## 2533     340527
    ## 2534     340526
    ## 2535       5303
    ## 2536       6191
    ## 2537       4435
    ## 2538      55869
    ## 2539       5255
    ## 2540       <NA>
    ## 2541       4674
    ## 2542      53344
    ## 2543     554203
    ## 2544  100302692
    ## 2545       6567
    ## 2546      51132
    ## 2547     340533
    ## 2548         22
    ## 2549     158866
    ## 2550     139599
    ## 2551      51260
    ## 2552      57692
    ## 2553       <NA>
    ## 2554       8823
    ## 2555        546
    ## 2556      84061
    ## 2557       1349
    ## 2558       <NA>
    ## 2559       5230
    ## 2560      51616
    ## 2561      10800
    ## 2562       <NA>
    ## 2563       2846
    ## 2564       <NA>
    ## 2565       9452
    ## 2566       <NA>
    ## 2567     254065
    ## 2568      79366
    ## 2569       6451
    ## 2570       <NA>
    ## 2571       5456
    ## 2572      27330
    ## 2573     139324
    ## 2574     139322
    ## 2575       <NA>
    ## 2576       1121
    ## 2577     117154
    ## 2578      56062
    ## 2579      27328
    ## 2580       4675
    ## 2581       <NA>
    ## 2582       1730
    ## 2583      57526
    ## 2584      64102
    ## 2585       7105
    ## 2586      27286
    ## 2587       1478
    ## 2588     402415
    ## 2589      79979
    ## 2590      59353
    ## 2591       2491
    ## 2592       1821
    ## 2593       <NA>
    ## 2594        695
    ## 2595       6173
    ## 2596       2717
    ## 2597       3188
    ## 2598  100131755
    ## 2599      51309
    ## 2600      54470
    ## 2601      51566
    ## 2602       9823
    ## 2603      84460
    ## 2604     158931
    ## 2605       <NA>
    ## 2606       9737
    ## 2607      64860
    ## 2608      80823
    ## 2609     114928
    ## 2610       <NA>
    ## 2611       <NA>
    ## 2612      84707
    ## 2613      56271
    ## 2614      90843
    ## 2615     340543
    ## 2616      55859
    ## 2617      56849
    ## 2618      51186
    ## 2619      27018
    ## 2620      85012
    ## 2621       9338
    ## 2622       9643
    ## 2623       5354
    ## 2624      51209
    ## 2625       <NA>
    ## 2626       <NA>
    ## 2627       <NA>
    ## 2628     401612
    ## 2629     644353
    ## 2630     139231
    ## 2631      26280
    ## 2632       <NA>
    ## 2633       <NA>
    ## 2634      79589
    ## 2635      54885
    ## 2636      79710
    ## 2637      84443
    ## 2638       5631
    ## 2639       1831
    ## 2640      11043
    ## 2641       5716
    ## 2642     115201
    ## 2643       1287
    ## 2644       8471
    ## 2645      55916
    ## 2646       <NA>
    ## 2647       2182
    ## 2648      84187
    ## 2649       9949
    ## 2650      91851
    ## 2651       5063
    ## 2652        827
    ## 2653       1641
    ## 2654      79868
    ## 2655       7224
    ## 2656     340595
    ## 2657     340596
    ## 2658     154796
    ## 2659       3358
    ## 2660      57631
    ## 2661       <NA>
    ## 2662       <NA>
    ## 2663      27301
    ## 2664        212
    ## 2665       5207
    ## 2666       7216
    ## 2667      10916
    ## 2668      54552
    ## 2669       2245
    ## 2670      90121
    ## 2671      65267
    ## 2672       <NA>
    ## 2673      54954
    ## 2674      23133
    ## 2675      10075
    ## 2676       3028
    ## 2677     158787
    ## 2678       8243
    ## 2679      23096
    ## 2680       8242
    ## 2681  102723508
    ## 2682      64061
    ## 2683      54328
    ## 2684        357
    ## 2685      28986
    ## 2686     139628
    ## 2687      10325
    ## 2688      11279
    ## 2689      29978
    ## 2690       <NA>
    ## 2691     550643
    ## 2692       <NA>
    ## 2693       6303
    ## 2694      23597
    ## 2695      10549
    ## 2696     139411
    ## 2697       <NA>
    ## 2698       5251
    ## 2699       6611
    ## 2700      51360
    ## 2701      22866
    ## 2702       6197
    ## 2703       1964
    ## 2704     256714
    ## 2705     256643
    ## 2706      30011
    ## 2707     389840
    ## 2708       5160
    ## 2709      10149
    ## 2710       5256
    ## 2711       5475
    ## 2712       6792
    ## 2713       <NA>
    ## 2714      10389
    ## 2715      10742
    ## 2716       4810
    ## 2717       9185
    ## 2718       5931
    ## 2719      55787
    ## 2720      94056
    ## 2721      56474
    ## 2722       <NA>
    ## 2723       8905
    ## 2724       8233
    ## 2725       <NA>
    ## 2726       <NA>
    ## 2727      59272
    ## 2728        660
    ## 2729       8544
    ## 2730       2277
    ## 2731       5277
    ## 2732     140456
    ## 2733     158747
    ## 2734       2187
    ## 2735       2742
    ## 2736      54960
    ## 2737       2824
    ## 2738       8481
    ## 2739       6399
    ## 2740       <NA>
    ## 2741     170082
    ## 2742      25975
    ## 2743       7114
    ## 2744      51284
    ## 2745       5634
    ## 2746       9758
    ## 2747      10943
    ## 2748        395
    ## 2749       3052
    ## 2750       4281
    ## 2751       <NA>
    ## 2752       <NA>
    ## 2753       <NA>
    ## 2754       <NA>
    ## 2755      79776
    ## 2756       5828
    ## 2757       <NA>
    ## 2758       5569
    ## 2759      51101
    ## 2760      11075
    ## 2761      23462
    ## 2762      28957
    ## 2763       7163
    ## 2764      65986
    ## 2765       <NA>
    ## 2766      55824
    ## 2767       2171
    ## 2768       3612
    ## 2769      79752
    ## 2770      64089
    ## 2771       <NA>
    ## 2772     138046
    ## 2773      85444
    ## 2774       1875
    ## 2775       <NA>
    ## 2776       <NA>
    ## 2777       <NA>
    ## 2778       <NA>
    ## 2779       <NA>
    ## 2780       <NA>
    ## 2781     253943
    ## 2782       <NA>
    ## 2783  100130155
    ## 2784      27319
    ## 2785       9420
    ## 2786      55156
    ## 2787       9650
    ## 2788       5150
    ## 2789       <NA>
    ## 2790       1392
    ## 2791       <NA>
    ## 2792       1356
    ## 2793      84343
    ## 2794       6596
    ## 2795       <NA>
    ## 2796       <NA>
    ## 2797       <NA>
    ## 2798      79718
    ## 2799     254827
    ## 2800      22871
    ## 2801       <NA>
    ## 2802       1894
    ## 2803      57552
    ## 2804       8743
    ## 2805      64778
    ## 2806       5337
    ## 2807      23043
    ## 2808       <NA>
    ## 2809      56648
    ## 2810     200916
    ## 2811       <NA>
    ## 2812       2122
    ## 2813      84517
    ## 2814      55892
    ## 2815       <NA>
    ## 2816       7095
    ## 2817      26996
    ## 2818      80012
    ## 2819       5584
    ## 2820       6498
    ## 2821       5010
    ## 2822      57709
    ## 2823       <NA>
    ## 2824       <NA>
    ## 2825      10242
    ## 2826      64393
    ## 2827       5290
    ## 2828       <NA>
    ## 2829      55669
    ## 2830      59345
    ## 2831         86
    ## 2832      57129
    ## 2833       4711
    ## 2834       8975
    ## 2835      51555
    ## 2836     151613
    ## 2837     339829
    ## 2838       <NA>
    ## 2839       8087
    ## 2840     131118
    ## 2841       <NA>
    ## 2842       6657
    ## 2843      23200
    ## 2844      54165
    ## 2845       <NA>
    ## 2846      56922
    ## 2847       <NA>
    ## 2848      28976
    ## 2849       <NA>
    ## 2850      84109
    ## 2851        308
    ## 2852       <NA>
    ## 2853       <NA>
    ## 2854       5393
    ## 2855        890
    ## 2856      55212
    ## 2857       7222
    ## 2858       <NA>
    ## 2859       3558
    ## 2860       <NA>
    ## 2861     166379
    ## 2862       2247
    ## 2863      11162
    ## 2864     166378
    ## 2865      10252
    ## 2866       <NA>
    ## 2867       <NA>
    ## 2868      57182
    ## 2869       <NA>
    ## 2870      79633
    ## 2871      27152
    ## 2872      22824
    ## 2873      10733
    ## 2874     256471
    ## 2875      80167
    ## 2876      55132
    ## 2877      10424
    ## 2878      79960
    ## 2879     132320
    ## 2880       <NA>
    ## 2881      57575
    ## 2882     132430
    ## 2883      54510
    ## 2884      23657
    ## 2885      25819
    ## 2886       1998
    ## 2887       4717
    ## 2888      80155
    ## 2889      83452
    ## 2890      80854
    ## 2891       4258
    ## 2892      55534
    ## 2893       2308
    ## 2894      57511
    ## 2895       <NA>
    ## 2896       <NA>
    ## 2897     387921
    ## 2898      80209
    ## 2899     341640
    ## 2900      51569
    ## 2901       7223
    ## 2902      10631
    ## 2903       <NA>
    ## 2904      11340
    ## 2905      29880
    ## 2906       4093
    ## 2907       5994
    ## 2908     400120
    ## 2909       8900
    ## 2910       <NA>
    ## 2911     728591
    ## 2912       9201
    ## 2913      26960
    ## 2914       4081
    ## 2915       4071
    ## 2916      25937
    ## 2917      51122
    ## 2918     389161
    ## 2919      11342
    ## 2920       5217
    ## 2921       <NA>
    ## 2922       <NA>
    ## 2923       9819
    ## 2924      27230
    ## 2925      83939
    ## 2926      51714
    ## 2927     131831
    ## 2928       6478
    ## 2929       <NA>
    ## 2930     116931
    ## 2931       9934
    ## 2932      53829
    ## 2933      64805
    ## 2934     285313
    ## 2935       <NA>
    ## 2936       4154
    ## 2937       5028
    ## 2938       <NA>
    ## 2939       5912
    ## 2940       <NA>
    ## 2941      26084
    ## 2942     170506
    ## 2943     344758
    ## 2944       4311
    ## 2945      23007
    ## 2946       <NA>
    ## 2947       9197
    ## 2948       8833
    ## 2949       7881
    ## 2950       6747
    ## 2951      25976
    ## 2952     389170
    ## 2953      57018
    ## 2954       <NA>
    ## 2955      51319
    ## 2956       4291
    ## 2957      85476
    ## 2958      56925
    ## 2959       5918
    ## 2960      64747
    ## 2961      29970
    ## 2962       3592
    ## 2963       <NA>
    ## 2964      57560
    ## 2965      10051
    ## 2966     286827
    ## 2967       3840
    ## 2968     151742
    ## 2969       8706
    ## 2970      51068
    ## 2971     165679
    ## 2972      22865
    ## 2973       <NA>
    ## 2974        590
    ## 2975      11235
    ## 2976       5274
    ## 2977      27333
    ## 2978      56884
    ## 2979       9693
    ## 2980       <NA>
    ## 2981      57600
    ## 2982       5481
    ## 2983       2110
    ## 2984       <NA>
    ## 2985      59350
    ## 2986      55314
    ## 2987       <NA>
    ## 2988       <NA>
    ## 2989       2891
    ## 2990       2743
    ## 2991      56034
    ## 2992       <NA>
    ## 2993       1519
    ## 2994       6999
    ## 2995       2983
    ## 2996       2982
    ## 2997      79884
    ## 2998       4887
    ## 2999       5356
    ## 3000       6423
    ## 3001       <NA>
    ## 3002       7097
    ## 3003      23240
    ## 3004      84057
    ## 3005      23321
    ## 3006      85462
    ## 3007      27236
    ## 3008     201799
    ## 3009      55294
    ## 3010       5188
    ## 3011     729830
    ## 3012       <NA>
    ## 3013     152503
    ## 3014       <NA>
    ## 3015        987
    ## 3016     166614
    ## 3017       <NA>
    ## 3018       <NA>
    ## 3019       <NA>
    ## 3020     115350
    ## 3021       2117
    ## 3022       9826
    ## 3023     375033
    ## 3024       4914
    ## 3025       5546
    ## 3026       3068
    ## 3027      79590
    ## 3028      51093
    ## 3029      81875
    ## 3030       1382
    ## 3031      10763
    ## 3032      63827
    ## 3033      60484
    ## 3034      54865
    ## 3035     128240
    ## 3036     128239
    ## 3037       4209
    ## 3038       <NA>
    ## 3039       <NA>
    ## 3040       <NA>
    ## 3041     128229
    ## 3042       7203
    ## 3043     112770
    ## 3044      84283
    ## 3045      23381
    ## 3046      79957
    ## 3047      11243
    ## 3048       9673
    ## 3049      64218
    ## 3050       4000
    ## 3051       <NA>
    ## 3052      92312
    ## 3053      28956
    ## 3054      56893
    ## 3055       6746
    ## 3056       9181
    ## 3057       <NA>
    ## 3058       <NA>
    ## 3059       6016
    ## 3060      23208
    ## 3061      54856
    ## 3062       <NA>
    ## 3063      55154
    ## 3064       7818
    ## 3065      55870
    ## 3066      23623
    ## 3067       2224
    ## 3068      57657
    ## 3069       1196
    ## 3070      10067
    ## 3071      10712
    ## 3072       2629
    ## 3073       4580
    ## 3074       7059
    ## 3075       4582
    ## 3076      80128
    ## 3077     200185
    ## 3078      54344
    ## 3079      55974
    ## 3080       1942
    ## 3081       1944
    ## 3082       1945
    ## 3083       8751
    ## 3084     127579
    ## 3085      51043
    ## 3086       <NA>
    ## 3087      80308
    ## 3088       1163
    ## 3089       6464
    ## 3090      90780
    ## 3091      57326
    ## 3092      10654
    ## 3093       3782
    ## 3094        103
    ## 3095       1141
    ## 3096       <NA>
    ## 3097      55585
    ## 3098     126669
    ## 3099       <NA>
    ## 3100      57198
    ## 3101      10456
    ## 3102       9898
    ## 3103       <NA>
    ## 3104       7170
    ## 3105       6232
    ## 3106       5872
    ## 3107      10899
    ## 3108      27173
    ## 3109     200186
    ## 3110       9909
    ## 3111      57459
    ## 3112      11000
    ## 3113      65123
    ## 3114       4881
    ## 3115       3608
    ## 3116      23557
    ## 3117      26097
    ## 3118       6271
    ## 3119       6284
    ## 3120     140576
    ## 3121       6274
    ## 3122       6275
    ## 3123       6276
    ## 3124       6277
    ## 3125     338324
    ## 3126       6279
    ## 3127       6280
    ## 3128       4014
    ## 3129       <NA>
    ## 3130       6698
    ## 3131       7062
    ## 3132       6282
    ## 3133       6281
    ## 3134     117145
    ## 3135  100191040
    ## 3136       6097
    ## 3137       <NA>
    ## 3138      11022
    ## 3139      65005
    ## 3140     284485
    ## 3141      11189
    ## 3142      81609
    ## 3143       7286
    ## 3144      57530
    ## 3145      23126
    ## 3146       5692
    ## 3147       8991
    ## 3148       5993
    ## 3149       <NA>
    ## 3150       5298
    ## 3151       <NA>
    ## 3152       <NA>
    ## 3153       5710
    ## 3154       8394
    ## 3155       6944
    ## 3156      29765
    ## 3157      79005
    ## 3158     388695
    ## 3159      79626
    ## 3160      10500
    ## 3161     126626
    ## 3162       <NA>
    ## 3163      10962
    ## 3164      56882
    ## 3165       <NA>
    ## 3166      58497
    ## 3167      55793
    ## 3168       8416
    ## 3169       <NA>
    ## 3170      29956
    ## 3171       9869
    ## 3172       <NA>
    ## 3173        405
    ## 3174       1513
    ## 3175       1520
    ## 3176      55204
    ## 3177       2029
    ## 3178       4170
    ## 3179      54507
    ## 3180       1893
    ## 3181      80222
    ## 3182      23248
    ## 3183       9129
    ## 3184      54460
    ## 3185       <NA>
    ## 3186     148523
    ## 3187       <NA>
    ## 3188      51107
    ## 3189       <NA>
    ## 3190      81611
    ## 3191      51177
    ## 3192      11311
    ## 3193      56957
    ## 3194      10903
    ## 3195      10262
    ## 3196       9900
    ## 3197      51027
    ## 3198       <NA>
    ## 3199       <NA>
    ## 3200       <NA>
    ## 3201       <NA>
    ## 3202       <NA>
    ## 3203       <NA>
    ## 3204       <NA>
    ## 3205       <NA>
    ## 3206       <NA>
    ## 3207      10628
    ## 3208      84265
    ## 3209     284615
    ## 3210     128077
    ## 3211       9939
    ## 3212       8799
    ## 3213       8515
    ## 3214     148741
    ## 3215      10401
    ## 3216     200035
    ## 3217      10623
    ## 3218      27246
    ## 3219       <NA>
    ## 3220      51205
    ## 3221        607
    ## 3222       9557
    ## 3223       2330
    ## 3224       5565
    ## 3225       9659
    ## 3226       9554
    ## 3227       <NA>
    ## 3228       4853
    ## 3229       3158
    ## 3230      26227
    ## 3231       <NA>
    ## 3232      10352
    ## 3233       6913
    ## 3234      10885
    ## 3235      54834
    ## 3236       <NA>
    ## 3237      10905
    ## 3238      80263
    ## 3239       8458
    ## 3240       5738
    ## 3241       3321
    ## 3242        476
    ## 3243      55356
    ## 3244       4808
    ## 3245        845
    ## 3246      81839
    ## 3247       <NA>
    ## 3248       4803
    ## 3249       <NA>
    ## 3250      10100
    ## 3251      80143
    ## 3252       7812
    ## 3253       4893
    ## 3254       <NA>
    ## 3255     163259
    ## 3256      10286
    ## 3257      51592
    ## 3258     148281
    ## 3259      56944
    ## 3260     204851
    ## 3261      64858
    ## 3262      10717
    ## 3263     440603
    ## 3264      26191
    ## 3265      54665
    ## 3266       <NA>
    ## 3267      10745
    ## 3268     260425
    ## 3269       9860
    ## 3270       6566
    ## 3271     333926
    ## 3272        389
    ## 3273       4343
    ## 3274        829
    ## 3275      54879
    ## 3276       7482
    ## 3277      55917
    ## 3278       <NA>
    ## 3279       3752
    ## 3280      11218
    ## 3281       <NA>
    ## 3282       <NA>
    ## 3283       5906
    ## 3284        140
    ## 3285       <NA>
    ## 3286       <NA>
    ## 3287      79084
    ## 3288       5016
    ## 3289     128344
    ## 3290       <NA>
    ## 3291       <NA>
    ## 3292      10390
    ## 3293     128338
    ## 3294      55791
    ## 3295        963
    ## 3296       <NA>
    ## 3297       3737
    ## 3298       3744
    ## 3299      10542
    ## 3300       9122
    ## 3301      64783
    ## 3302       3749
    ## 3303     388662
    ## 3304        257
    ## 3305      85369
    ## 3306      10768
    ## 3307       1435
    ## 3308       <NA>
    ## 3309       2949
    ## 3310       <NA>
    ## 3311       <NA>
    ## 3312       2946
    ## 3313       2944
    ## 3314       2948
    ## 3315        271
    ## 3316       <NA>
    ## 3317       2773
    ## 3318      83873
    ## 3319      57463
    ## 3320     284613
    ## 3321     127002
    ## 3322     284612
    ## 3323       5686
    ## 3324       6272
    ## 3325       <NA>
    ## 3326      84722
    ## 3327       1952
    ## 3328       <NA>
    ## 3329       <NA>
    ## 3330       <NA>
    ## 3331      56900
    ## 3332       6884
    ## 3333      22911
    ## 3334      23155
    ## 3335      29899
    ## 3336       6814
    ## 3337     163479
    ## 3338      55119
    ## 3339     113802
    ## 3340     284611
    ## 3341      29957
    ## 3342      10451
    ## 3343      22854
    ## 3344      55170
    ## 3345       <NA>
    ## 3346       <NA>
    ## 3347      55599
    ## 3348       1301
    ## 3349     118427
    ## 3350       1901
    ## 3351       <NA>
    ## 3352      51611
    ## 3353     148867
    ## 3354       2135
    ## 3355       7412
    ## 3356      54112
    ## 3357       8556
    ## 3358       8634
    ## 3359       1629
    ## 3360      54482
    ## 3361     163786
    ## 3362       <NA>
    ## 3363      64645
    ## 3364      23443
    ## 3365        178
    ## 3366     391059
    ## 3367      54873
    ## 3368       9890
    ## 3369     163404
    ## 3370      51375
    ## 3371       <NA>
    ## 3372       1806
    ## 3373      58155
    ## 3374      25950
    ## 3375       <NA>
    ## 3376     199857
    ## 3377       1266
    ## 3378       <NA>
    ## 3379       2152
    ## 3380       5825
    ## 3381       9411
    ## 3382         24
    ## 3383       2730
    ## 3384      30836
    ## 3385       8412
    ## 3386      54874
    ## 3387       8654
    ## 3388       <NA>
    ## 3389      54532
    ## 3390     171024
    ## 3391       9871
    ## 3392      57721
    ## 3393       <NA>
    ## 3394       8492
    ## 3395       9348
    ## 3396     133022
    ## 3397      64579
    ## 3398       <NA>
    ## 3399      79642
    ## 3400        817
    ## 3401        287
    ## 3402       <NA>
    ## 3403      51574
    ## 3404      55345
    ## 3405      63973
    ## 3406      80216
    ## 3407      92610
    ## 3408      55435
    ## 3409     132720
    ## 3410       2028
    ## 3411      79071
    ## 3412       1950
    ## 3413      54433
    ## 3414      81579
    ## 3415        839
    ## 3416      55013
    ## 3417      10427
    ## 3418      84570
    ## 3419      64850
    ## 3420      58505
    ## 3421       6164
    ## 3422      51176
    ## 3423       3033
    ## 3424     113612
    ## 3425     166929
    ## 3426       9061
    ## 3427      27123
    ## 3428       9255
    ## 3429      93627
    ## 3430     255743
    ## 3431      79807
    ## 3432      57117
    ## 3433       <NA>
    ## 3434      27068
    ## 3435      54790
    ## 3436       <NA>
    ## 3437       <NA>
    ## 3438       <NA>
    ## 3439      80319
    ## 3440       6870
    ## 3441       1062
    ## 3442      56898
    ## 3443     133308
    ## 3444     150159
    ## 3445     493856
    ## 3446       <NA>
    ## 3447       7323
    ## 3448       4126
    ## 3449       4790
    ## 3450      64116
    ## 3451      55024
    ## 3452       5530
    ## 3453      51705
    ## 3454     115265
    ## 3455       <NA>
    ## 3456       <NA>
    ## 3457      79982
    ## 3458       <NA>
    ## 3459       8649
    ## 3460      27071
    ## 3461       <NA>
    ## 3462       <NA>
    ## 3463       4547
    ## 3464      93587
    ## 3465       <NA>
    ## 3466        128
    ## 3467      23173
    ## 3468       1977
    ## 3469      10098
    ## 3470       5910
    ## 3471       <NA>
    ## 3472       8633
    ## 3473        658
    ## 3474      10611
    ## 3475     115362
    ## 3476     388646
    ## 3477       2635
    ## 3478       2634
    ## 3479      56267
    ## 3480       2959
    ## 3481       5586
    ## 3482       <NA>
    ## 3483       <NA>
    ## 3484       <NA>
    ## 3485       8543
    ## 3486       9653
    ## 3487       9403
    ## 3488      51100
    ## 3489       <NA>
    ## 3490      57489
    ## 3491     255631
    ## 3492      54680
    ## 3493       <NA>
    ## 3494       <NA>
    ## 3495      23576
    ## 3496       8915
    ## 3497       <NA>
    ## 3498      84144
    ## 3499     117178
    ## 3500       1486
    ## 3501  100505741
    ## 3502       2787
    ## 3503      80135
    ## 3504     391051
    ## 3505       5567
    ## 3506      79739
    ## 3507       <NA>
    ## 3508      23266
    ## 3509       <NA>
    ## 3510      64123
    ## 3511      10561
    ## 3512       5737
    ## 3513      54810
    ## 3514       <NA>
    ## 3515      11080
    ## 3516       8880
    ## 3517      91624
    ## 3518     374986
    ## 3519      23032
    ## 3520      26009
    ## 3521      26289
    ## 3522       <NA>
    ## 3523      10026
    ## 3524      81849
    ## 3525     256435
    ## 3526       <NA>
    ## 3527       5876
    ## 3528         34
    ## 3529     204962
    ## 3530     431707
    ## 3531     127253
    ## 3532       1429
    ## 3533     127254
    ## 3534      51086
    ## 3535       8790
    ## 3536     127255
    ## 3537       <NA>
    ## 3538     257194
    ## 3539       <NA>
    ## 3540       9406
    ## 3541       5733
    ## 3542       1491
    ## 3543      81573
    ## 3544       9295
    ## 3545      55631
    ## 3546       <NA>
    ## 3547      57554
    ## 3548       6121
    ## 3549      79971
    ## 3550     137695
    ## 3551      96764
    ## 3552       4067
    ## 3553       6224
    ## 3554      79145
    ## 3555       <NA>
    ## 3556       5179
    ## 3557      54928
    ## 3558      90362
    ## 3559     137886
    ## 3560       6386
    ## 3561       8439
    ## 3562       9760
    ## 3563       <NA>
    ## 3564       5862
    ## 3565      55636
    ## 3566     157807
    ## 3567        444
    ## 3568       <NA>
    ## 3569      79666
    ## 3570     137682
    ## 3571       <NA>
    ## 3572       9134
    ## 3573      55656
    ## 3574     286148
    ## 3575      25962
    ## 3576      25788
    ## 3577  100861412
    ## 3578       2669
    ## 3579      54704
    ## 3580      91147
    ## 3581       <NA>
    ## 3582       <NA>
    ## 3583     137392
    ## 3584     286144
    ## 3585        862
    ## 3586     115111
    ## 3587      51633
    ## 3588       <NA>
    ## 3589       <NA>
    ## 3590      64168
    ## 3591     169200
    ## 3592        793
    ## 3593       1666
    ## 3594       4683
    ## 3595        734
    ## 3596       8767
    ## 3597       4325
    ## 3598       8895
    ## 3599      51115
    ## 3600      11059
    ## 3601       7274
    ## 3602       8836
    ## 3603     286183
    ## 3604        892
    ## 3605  100130890
    ## 3606      85015
    ## 3607      25957
    ## 3608      51805
    ## 3609      84553
    ## 3610      26235
    ## 3611       5454
    ## 3612       <NA>
    ## 3613     253714
    ## 3614     114792
    ## 3615       <NA>
    ## 3616      29078
    ## 3617      23376
    ## 3618      10690
    ## 3619      79694
    ## 3620       2045
    ## 3621       6885
    ## 3622      60468
    ## 3623       9994
    ## 3624      23195
    ## 3625      57226
    ## 3626      22881
    ## 3627      58528
    ## 3628      51465
    ## 3629       2570
    ## 3630     135295
    ## 3631      10957
    ## 3632       8732
    ## 3633       1268
    ## 3634      81833
    ## 3635      55122
    ## 3636      23595
    ## 3637      57038
    ## 3638      10559
    ## 3639     154313
    ## 3640      57150
    ## 3641       <NA>
    ## 3642      79817
    ## 3643       <NA>
    ## 3644     158038
    ## 3645       <NA>
    ## 3646         48
    ## 3647      23586
    ## 3648      10210
    ## 3649  100129250
    ## 3650       4712
    ## 3651     401498
    ## 3652      54840
    ## 3653       <NA>
    ## 3654       3301
    ## 3655      55234
    ## 3656       2683
    ## 3657        573
    ## 3658      51510
    ## 3659       4799
    ## 3660      65083
    ## 3661      54926
    ## 3662      55833
    ## 3663      25853
    ## 3664      51271
    ## 3665        318
    ## 3666       <NA>
    ## 3667       <NA>
    ## 3668     203259
    ## 3669       <NA>
    ## 3670     375704
    ## 3671       1271
    ## 3672       <NA>
    ## 3673     138716
    ## 3674      11258
    ## 3675      10280
    ## 3676       2592
    ## 3677       <NA>
    ## 3678       <NA>
    ## 3679       <NA>
    ## 3680       <NA>
    ## 3681       <NA>
    ## 3682       <NA>
    ## 3683       6363
    ## 3684       <NA>
    ## 3685      23349
    ## 3686      25822
    ## 3687       7415
    ## 3688       2189
    ## 3689      84720
    ## 3690      30968
    ## 3691      80256
    ## 3692       <NA>
    ## 3693      10497
    ## 3694       9853
    ## 3695     730112
    ## 3696       7016
    ## 3697        971
    ## 3698     203260
    ## 3699      84904
    ## 3700       7169
    ## 3701       7094
    ## 3702      10488
    ## 3703      57704
    ## 3704       9827
    ## 3705       4882
    ## 3706      26206
    ## 3707      84681
    ## 3708     392307
    ## 3709      51754
    ## 3710     646962
    ## 3711     158376
    ## 3712       8434
    ## 3713     152007
    ## 3714       1211
    ## 3715      10020
    ## 3716     152006
    ## 3717       9833
    ## 3718      84186
    ## 3719       9380
    ## 3720       9925
    ## 3721       <NA>
    ## 3722      64425
    ## 3723      26267
    ## 3724     401505
    ## 3725      22844
    ## 3726     158234
    ## 3727      51010
    ## 3728      79269
    ## 3729      92014
    ## 3730       6461
    ## 3731        219
    ## 3732     347252
    ## 3733      23424
    ## 3734       7111
    ## 3735     158427
    ## 3736       4686
    ## 3737       7507
    ## 3738      51531
    ## 3739      10541
    ## 3740      54187
    ## 3741       7464
    ## 3742      55357
    ## 3743       9568
    ## 3744     203286
    ## 3745      79695
    ## 3746       1306
    ## 3747       7046
    ## 3748      85365
    ## 3749      10952
    ## 3750       8013
    ## 3751      55014
    ## 3752      23071
    ## 3753      27130
    ## 3754      54881
    ## 3755      91283
    ## 3756       8577
    ## 3757     347273
    ## 3758      54886
    ## 3759      54534
    ## 3760       <NA>
    ## 3761        229
    ## 3762       <NA>
    ## 3763      56254
    ## 3764     116443
    ## 3765       5535
    ## 3766       <NA>
    ## 3767      10592
    ## 3768      55335
    ## 3769         19
    ## 3770      23446
    ## 3771      83856
    ## 3772       2218
    ## 3773      55151
    ## 3774       <NA>
    ## 3775       5887
    ## 3776       9314
    ## 3777       <NA>
    ## 3778       <NA>
    ## 3779       8727
    ## 3780      23731
    ## 3781      23732
    ## 3782      54566
    ## 3783       5774
    ## 3784       <NA>
    ## 3785       <NA>
    ## 3786       <NA>
    ## 3787      79987
    ## 3788       1902
    ## 3789       <NA>
    ## 3790      22949
    ## 3791     548645
    ## 3792       2790
    ## 3793       7357
    ## 3794      64420
    ## 3795       9991
    ## 3796      84263
    ## 3797       <NA>
    ## 3798      58493
    ## 3799     401548
    ## 3800       <NA>
    ## 3801       7539
    ## 3802       1318
    ## 3803      23307
    ## 3804       1317
    ## 3805     246184
    ## 3806       9128
    ## 3807     114987
    ## 3808      54836
    ## 3809      81932
    ## 3810        210
    ## 3811      54107
    ## 3812       5998
    ## 3813       <NA>
    ## 3814      85301
    ## 3815      80709
    ## 3816      25861
    ## 3817       9550
    ## 3818     203197
    ## 3819       3371
    ## 3820      23245
    ## 3821      22954
    ## 3822       7099
    ## 3823       1620
    ## 3824      55755
    ## 3825       1955
    ## 3826       7088
    ## 3827     257019
    ## 3828      23081
    ## 3829       <NA>
    ## 3830      90871
    ## 3831       5789
    ## 3832     286343
    ## 3833       8777
    ## 3834       4781
    ## 3835       <NA>
    ## 3836     340481
    ## 3837     158326
    ## 3838     158219
    ## 3839       6619
    ## 3840      11168
    ## 3841     203238
    ## 3842      54875
    ## 3843       6456
    ## 3844      92949
    ## 3845      10670
    ## 3846      54801
    ## 3847        123
    ## 3848      55667
    ## 3849       6194
    ## 3850     340485
    ## 3851      25769
    ## 3852       4300
    ## 3853      54914
    ## 3854     401494
    ## 3855      55958
    ## 3856       <NA>
    ## 3857       4507
    ## 3858       1030
    ## 3859       1993
    ## 3860     286319
    ## 3861      79886
    ## 3862       <NA>
    ## 3863       9373
    ## 3864      80173
    ## 3865       7010
    ## 3866      54586
    ## 3867     114803
    ## 3868       3725
    ## 3869       <NA>
    ## 3870      55277
    ## 3871      51361
    ## 3872       <NA>
    ## 3873       <NA>
    ## 3874       <NA>
    ## 3875       <NA>
    ## 3876       4774
    ## 3877      83941
    ## 3878      10207
    ## 3879     163782
    ## 3880       7398
    ## 3881      85440
    ## 3882      27329
    ## 3883      84938
    ## 3884       <NA>
    ## 3885      27022
    ## 3886      29929
    ## 3887      23421
    ## 3888      84455
    ## 3889      55276
    ## 3890       4919
    ## 3891      57685
    ## 3892      55225
    ## 3893       3716
    ## 3894       <NA>
    ## 3895        205
    ## 3896       9829
    ## 3897      54741
    ## 3898       3953
    ## 3899       5142
    ## 3900      84251
    ## 3901       <NA>
    ## 3902      10022
    ## 3903      79819
    ## 3904      57708
    ## 3905      23169
    ## 3906       <NA>
    ## 3907     115209
    ## 3908       1600
    ## 3909     199920
    ## 3910       5563
    ## 3911       8613
    ## 3912      23358
    ## 3913     255738
    ## 3914       <NA>
    ## 3915       1718
    ## 3916       <NA>
    ## 3917      55001
    ## 3918      25973
    ## 3919       7268
    ## 3920     374977
    ## 3921      26027
    ## 3922      23648
    ## 3923      51253
    ## 3924     606495
    ## 3925       <NA>
    ## 3926     127428
    ## 3927       9528
    ## 3928     115353
    ## 3929      51668
    ## 3930      54432
    ## 3931      55706
    ## 3932     148979
    ## 3933      63948
    ## 3934       <NA>
    ## 3935       7804
    ## 3936       <NA>
    ## 3937       4116
    ## 3938       <NA>
    ## 3939       1376
    ## 3940     127435
    ## 3941       6342
    ## 3942      55268
    ## 3943      79699
    ## 3944      65260
    ## 3945       2882
    ## 3946       <NA>
    ## 3947      84950
    ## 3948     200014
    ## 3949       <NA>
    ## 3950       9372
    ## 3951      91408
    ## 3952      51060
    ## 3953     112970
    ## 3954       5865
    ## 3955       <NA>
    ## 3956     114883
    ## 3957       2060
    ## 3958      22996
    ## 3959      26994
    ## 3960       1031
    ## 3961      11124
    ## 3962       <NA>
    ## 3963      63950
    ## 3964       1996
    ## 3965      84871
    ## 3966      79656
    ## 3967      54558
    ## 3968     388630
    ## 3969       2306
    ## 3970      51727
    ## 3971       6491
    ## 3972       6886
    ## 3973     260293
    ## 3974       <NA>
    ## 3975       1580
    ## 3976       9813
    ## 3977      64756
    ## 3978     148932
    ## 3979       8569
    ## 3980       2166
    ## 3981     387338
    ## 3982       7388
    ## 3983      10489
    ## 3984       8438
    ## 3985      55624
    ## 3986     541468
    ## 3987      10103
    ## 3988       8503
    ## 3989      23139
    ## 3990       3652
    ## 3991      51249
    ## 3992      60313
    ## 3993     149483
    ## 3994       4678
    ## 3995      10327
    ## 3996       5052
    ## 3997      25974
    ## 3998     126661
    ## 3999      10420
    ## 4000     114034
    ## 4001       4595
    ## 4002      84842
    ## 4003      57643
    ## 4004       7389
    ## 4005      79654
    ## 4006       8891
    ## 4007       8643
    ## 4008     149478
    ## 4009       1263
    ## 4010       6202
    ## 4011      11004
    ## 4012      79639
    ## 4013      55182
    ## 4014       <NA>
    ## 4015      79033
    ## 4016      55929
    ## 4017       6536
    ## 4018     149473
    ## 4019       8704
    ## 4020        533
    ## 4021       1802
    ## 4022       9670
    ## 4023       9048
    ## 4024       6487
    ## 4025       9682
    ## 4026       5792
    ## 4027      23334
    ## 4028     112950
    ## 4029      64834
    ## 4030        991
    ## 4031       4352
    ## 4032       7075
    ## 4033     128218
    ## 4034     149465
    ## 4035      10969
    ## 4036       6513
    ## 4037       <NA>
    ## 4038     114625
    ## 4039     374969
    ## 4040       <NA>
    ## 4041      64175
    ## 4042     149461
    ## 4043       4904
    ## 4044      10465
    ## 4045     728621
    ## 4046      79717
    ## 4047      84217
    ## 4048     284716
    ## 4049       <NA>
    ## 4050      22887
    ## 4051      59269
    ## 4052  100132074
    ## 4053      22955
    ## 4054       <NA>
    ## 4055     163732
    ## 4056       9132
    ## 4057       4802
    ## 4058       9783
    ## 4059      64789
    ## 4060     339559
    ## 4061      64744
    ## 4062       1298
    ## 4063      10269
    ## 4064     127391
    ## 4065       6018
    ## 4066       5538
    ## 4067      10487
    ## 4068      84879
    ## 4069       4610
    ## 4070      54802
    ## 4071      10450
    ## 4072      51440
    ## 4073      26508
    ## 4074       8761
    ## 4075      23499
    ## 4076       4725
    ## 4077      79647
    ## 4078       <NA>
    ## 4079      26292
    ## 4080      64121
    ## 4081       <NA>
    ## 4082       5453
    ## 4083      51118
    ## 4084       2275
    ## 4085      10946
    ## 4086       3633
    ## 4087       4520
    ## 4088       <NA>
    ## 4089      79693
    ## 4090     149175
    ## 4091     284656
    ## 4092      55143
    ## 4093       <NA>
    ## 4094     284654
    ## 4095      29889
    ## 4096      79753
    ## 4097      64769
    ## 4098      80149
    ## 4099       2899
    ## 4100       <NA>
    ## 4101       1441
    ## 4102      64960
    ## 4103     127700
    ## 4104      84967
    ## 4105      83931
    ## 4106      55194
    ## 4107      79729
    ## 4108       9967
    ## 4109      55700
    ## 4110      27095
    ## 4111       1296
    ## 4112      54936
    ## 4113      27285
    ## 4114     192669
    ## 4115      26523
    ## 4116     192670
    ## 4117       <NA>
    ## 4118       5690
    ## 4119      23154
    ## 4120       <NA>
    ## 4121       9202
    ## 4122       6421
    ## 4123      79830
    ## 4124       9204
    ## 4125  100506144
    ## 4126      58512
    ## 4127     113444
    ## 4128       2701
    ## 4129       <NA>
    ## 4130       <NA>
    ## 4131     114784
    ## 4132       7579
    ## 4133       <NA>
    ## 4134       1912
    ## 4135       <NA>
    ## 4136      55223
    ## 4137     113451
    ## 4138        204
    ## 4139       <NA>
    ## 4140     127544
    ## 4141       3208
    ## 4142     252995
    ## 4143      64766
    ## 4144       <NA>
    ## 4145       <NA>
    ## 4146       <NA>
    ## 4147       5928
    ## 4148       <NA>
    ## 4149     339487
    ## 4150     653121
    ## 4151     728116
    ## 4152      55108
    ## 4153  100128071
    ## 4154      65108
    ## 4155       3065
    ## 4156       3932
    ## 4157      84734
    ## 4158       8668
    ## 4159      56063
    ## 4160     149069
    ## 4161      55721
    ## 4162      79140
    ## 4163     200081
    ## 4164      23633
    ## 4165      55116
    ## 4166      10657
    ## 4167       8073
    ## 4168       <NA>
    ## 4169      90853
    ## 4170       <NA>
    ## 4171        576
    ## 4172       1307
    ## 4173     553115
    ## 4174       3061
    ## 4175      64129
    ## 4176     347735
    ## 4177       2170
    ## 4178      51538
    ## 4179       9410
    ## 4180      79570
    ## 4181       9698
    ## 4182       9672
    ## 4183       7805
    ## 4184       4146
    ## 4185      10076
    ## 4186      51102
    ## 4187       6429
    ## 4188       2035
    ## 4189       4985
    ## 4190      51441
    ## 4191      10691
    ## 4192       6883
    ## 4193     115273
    ## 4194      85028
    ## 4195      54952
    ## 4196       1104
    ## 4197       <NA>
    ## 4198      65979
    ## 4199      54797
    ## 4200      83667
    ## 4201       <NA>
    ## 4202      22826
    ## 4203       5724
    ## 4204       2140
    ## 4205      55113
    ## 4206      27293
    ## 4207       6118
    ## 4208       9473
    ## 4209       5511
    ## 4210      23673
    ## 4211     199870
    ## 4212      27245
    ## 4213      10163
    ## 4214       2827
    ## 4215     388611
    ## 4216       9064
    ## 4217      84958
    ## 4218      84065
    ## 4219      23038
    ## 4220       6548
    ## 4221     388610
    ## 4222      10726
    ## 4223      63906
    ## 4224      54707
    ## 4225       2810
    ## 4226      84243
    ## 4227      55650
    ## 4228       8289
    ## 4229       6195
    ## 4230       3151
    ## 4231      79947
    ## 4232      55057
    ## 4233       1043
    ## 4234      91544
    ## 4235      83442
    ## 4236      64793
    ## 4237       <NA>
    ## 4238       <NA>
    ## 4239     149420
    ## 4240       2134
    ## 4241       5051
    ## 4242       <NA>
    ## 4243       3925
    ## 4244     164091
    ## 4245      79000
    ## 4246      56181
    ## 4247      57190
    ## 4248      57134
    ## 4249      26119
    ## 4250       <NA>
    ## 4251       6007
    ## 4252      23585
    ## 4253      57035
    ## 4254      25949
    ## 4255      25932
    ## 4256      10250
    ## 4257     400746
    ## 4258      11123
    ## 4259      57185
    ## 4260      90529
    ## 4261      57822
    ## 4262     127294
    ## 4263      10772
    ## 4264      55629
    ## 4265       1269
    ## 4266       2517
    ## 4267       3155
    ## 4268       2582
    ## 4269      11313
    ## 4270      57095
    ## 4271       6924
    ## 4272       6135
    ## 4273       3399
    ## 4274       1870
    ## 4275      55616
    ## 4276       6920
    ## 4277       <NA>
    ## 4278      10236
    ## 4279       3352
    ## 4280       7798
    ## 4281      23028
    ## 4282       2048
    ## 4283        713
    ## 4284        714
    ## 4285        712
    ## 4286       2046
    ## 4287       9923
    ## 4288      54361
    ## 4289        998
    ## 4290       <NA>
    ## 4291       3339
    ## 4292      84196
    ## 4293       5909
    ## 4294        249
    ## 4295       1889
    ## 4296       8672
    ## 4297       <NA>
    ## 4298      50809
    ## 4299     400745
    ## 4300      57576
    ## 4301       1650
    ## 4302      65018
    ## 4303        978
    ## 4304     163933
    ## 4305       <NA>
    ## 4306      79594
    ## 4307      55450
    ## 4308     127731
    ## 4309     127733
    ## 4310       5322
    ## 4311      23252
    ## 4312     255104
    ## 4313       3362
    ## 4314       <NA>
    ## 4315       4681
    ## 4316       <NA>
    ## 4317       <NA>
    ## 4318        832
    ## 4319       <NA>
    ## 4320       <NA>
    ## 4321      51154
    ## 4322      23065
    ## 4323      23352
    ## 4324       <NA>
    ## 4325     126917
    ## 4326       8659
    ## 4327     127707
    ## 4328      84966
    ## 4329       <NA>
    ## 4330      55160
    ## 4331       <NA>
    ## 4332      55920
    ## 4333     353238
    ## 4334      11240
    ## 4335       <NA>
    ## 4336       6390
    ## 4337      23400
    ## 4338       4237
    ## 4339       9696
    ## 4340      55707
    ## 4341     374955
    ## 4342      26099
    ## 4343      54455
    ## 4344     128272
    ## 4345       <NA>
    ## 4346       1969
    ## 4347     348487
    ## 4348       <NA>
    ## 4349       1188
    ## 4350      27129
    ## 4351     149563
    ## 4352       7709
    ## 4353      23013
    ## 4354      54751
    ## 4355     388595
    ## 4356     284723
    ## 4357      23207
    ## 4358      84301
    ## 4359      79814
    ## 4360      23341
    ## 4361        842
    ## 4362      79180
    ## 4363       <NA>
    ## 4364     114827
    ## 4365      55092
    ## 4366      23254
    ## 4367       7799
    ## 4368      10630
    ## 4369     126755
    ## 4370     391002
    ## 4371       9249
    ## 4372      55187
    ## 4373       7133
    ## 4374        943
    ## 4375       <NA>
    ## 4376       <NA>
    ## 4377       <NA>
    ## 4378       <NA>
    ## 4379       <NA>
    ## 4380       <NA>
    ## 4381       <NA>
    ## 4382      60672
    ## 4383       <NA>
    ## 4384       9927
    ## 4385       5351
    ## 4386       <NA>
    ## 4387       4878
    ## 4388       <NA>
    ## 4389       1185
    ## 4390       4524
    ## 4391      57085
    ## 4392     374946
    ## 4393      10459
    ## 4394      26270
    ## 4395      93611
    ## 4396      26232
    ## 4397      57540
    ## 4398      29914
    ## 4399       2475
    ## 4400       5394
    ## 4401       6723
    ## 4402      10747
    ## 4403      23435
    ## 4404      54897
    ## 4405       <NA>
    ## 4406       5195
    ## 4407       1676
    ## 4408       1325
    ## 4409     378708
    ## 4410       5226
    ## 4411      23095
    ## 4412      10277
    ## 4413     116362
    ## 4414      64802
    ## 4415      84328
    ## 4416      56998
    ## 4417      22883
    ## 4418       5293
    ## 4419     199953
    ## 4420      84275
    ## 4421      80176
    ## 4422       9563
    ## 4423      80045
    ## 4424       6518
    ## 4425       2023
    ## 4426        473
    ## 4427      50651
    ## 4428      54206
    ## 4429      11315
    ## 4430       3604
    ## 4431      23261
    ## 4432       8863
    ## 4433       9341
    ## 4434      55735
    ## 4435      90326
    ## 4436     148479
    ## 4437       9903
    ## 4438       3104
    ## 4439      79707
    ## 4440      57449
    ## 4441       8718
    ## 4442      83715
    ## 4443      11332
    ## 4444     387509
    ## 4445     390992
    ## 4446      23463
    ## 4447     388591
    ## 4448       6146
    ## 4449      26038
    ## 4450       8514
    ## 4451     261734
    ## 4452      55966
    ## 4453       <NA>
    ## 4454       <NA>
    ## 4455       1677
    ## 4456       9731
    ## 4457      57470
    ## 4458     388588
    ## 4459       <NA>
    ## 4460      49856
    ## 4461       <NA>
    ## 4462       1953
    ## 4463      63976
    ## 4464       <NA>
    ## 4465  100287898
    ## 4466      79258
    ## 4467       <NA>
    ## 4468     388585
    ## 4469      55229
    ## 4470       9651
    ## 4471       5192
    ## 4472      11079
    ## 4473      79906
    ## 4474       6497
    ## 4475     199990
    ## 4476       5590
    ## 4477       2563
    ## 4478      85452
    ## 4479     339456
    ## 4480       2782
    ## 4481      65220
    ## 4482       <NA>
    ## 4483        984
    ## 4484       <NA>
    ## 4485     142678
    ## 4486     643988
    ## 4487      29101
    ## 4488       <NA>
    ## 4489     339453
    ## 4490      55210
    ## 4491       <NA>
    ## 4492      64856
    ## 4493     643965
    ## 4494      55052
    ## 4495      81669
    ## 4496      54998
    ## 4497      54587
    ## 4498       1855
    ## 4499      83756
    ## 4500      80772
    ## 4501      54973
    ## 4502     126789
    ## 4503     116983
    ## 4504     118424
    ## 4505     388581
    ## 4506     126792
    ## 4507      51150
    ## 4508       7293
    ## 4509       <NA>
    ## 4510       <NA>
    ## 4511     401934
    ## 4512     375790
    ## 4513       9636
    ## 4514       <NA>
    ## 4515      84808
    ## 4516      84069
    ## 4517      26155
    ## 4518     148398
    ## 4519       1021
    ## 4520     257415
    ## 4521       <NA>
    ## 4522      84060
    ## 4523       5189
    ## 4524       <NA>
    ## 4525      57798
    ## 4526      54467
    ## 4527        889
    ## 4528     401387
    ## 4529       <NA>
    ## 4530      10142
    ## 4531       <NA>
    ## 4532       <NA>
    ## 4533       8321
    ## 4534       5218
    ## 4535       <NA>
    ## 4536       9069
    ## 4537      85865
    ## 4538       <NA>
    ## 4539      79846
    ## 4540     261729
    ## 4541      26872
    ## 4542       <NA>
    ## 4543      79689
    ## 4544       6717
    ## 4545      53616
    ## 4546      10926
    ## 4547      55972
    ## 4548     154661
    ## 4549       <NA>
    ## 4550       <NA>
    ## 4551       5244
    ## 4552      54677
    ## 4553      79161
    ## 4554       9988
    ## 4555       <NA>
    ## 4556       2913
    ## 4557     223117
    ## 4558      10371
    ## 4559       9723
    ## 4560      27445
    ## 4561       <NA>
    ## 4562        781
    ## 4563       3082
    ## 4564      10512
    ## 4565        948
    ## 4566       2770
    ## 4567       <NA>
    ## 4568       9863
    ## 4569      57157
    ## 4570      85025
    ## 4571     222194
    ## 4572       <NA>
    ## 4573       5782
    ## 4574      54103
    ## 4575      57639
    ## 4576      10875
    ## 4577     222234
    ## 4578     222235
    ## 4579      10234
    ## 4580      83787
    ## 4581     222236
    ## 4582       9512
    ## 4583      27000
    ## 4584       5701
    ## 4585       5649
    ## 4586       5001
    ## 4587       <NA>
    ## 4588       <NA>
    ## 4589     375612
    ## 4590       <NA>
    ## 4591      55904
    ## 4592       6733
    ## 4593      54517
    ## 4594      60561
    ## 4595      54543
    ## 4596       <NA>
    ## 4597      84668
    ## 4598      55975
    ## 4599       <NA>
    ## 4600       3757
    ## 4601       4846
    ## 4602       <NA>
    ## 4603     285973
    ## 4604      11194
    ## 4605       9311
    ## 4606       1020
    ## 4607       6522
    ## 4608      10922
    ## 4609      83590
    ## 4610       <NA>
    ## 4611     116988
    ## 4612      10061
    ## 4613      54480
    ## 4614       6604
    ## 4615      51667
    ## 4616     349136
    ## 4617     155051
    ## 4618       6009
    ## 4619      51422
    ## 4620       <NA>
    ## 4621      63917
    ## 4622      58508
    ## 4623       <NA>
    ## 4624       7516
    ## 4625      57180
    ## 4626       <NA>
    ## 4627       1804
    ## 4628      22976
    ## 4629       3361
    ## 4630       3638
    ## 4631       2020
    ## 4632     155435
    ## 4633       6469
    ## 4634       <NA>
    ## 4635       <NA>
    ## 4636     140545
    ## 4637      64327
    ## 4638       <NA>
    ## 4639      64434
    ## 4640       9690
    ## 4641      10049
    ## 4642       7298
    ## 4643       <NA>
    ## 4644     150946
    ## 4645       3030
    ## 4646       3032
    ## 4647      85465
    ## 4648      92749
    ## 4649       9381
    ## 4650       3777
    ## 4651      54978
    ## 4652       1058
    ## 4653      56896
    ## 4654      22924
    ## 4655      54867
    ## 4656      60509
    ## 4657  100128731
    ## 4658      11117
    ## 4659       3795
    ## 4660      10669
    ## 4661      10113
    ## 4662       8884
    ## 4663      51374
    ## 4664        790
    ## 4665       7781
    ## 4666      57159
    ## 4667       7349
    ## 4668       4358
    ## 4669       <NA>
    ## 4670       2976
    ## 4671       8890
    ## 4672       9784
    ## 4673       <NA>
    ## 4674       5496
    ## 4675      29959
    ## 4676     200634
    ## 4677      26160
    ## 4678      64838
    ## 4679       <NA>
    ## 4680      11321
    ## 4681       9913
    ## 4682      22950
    ## 4683       <NA>
    ## 4684       9553
    ## 4685      64080
    ## 4686       9577
    ## 4687       2355
    ## 4688     151056
    ## 4689       5500
    ## 4690       7525
    ## 4691      23761
    ## 4692     253143
    ## 4693       9681
    ## 4694       7533
    ## 4695       1487
    ## 4696      10296
    ## 4697      57654
    ## 4698     152877
    ## 4699       7884
    ## 4700      92305
    ## 4701      10460
    ## 4702       2261
    ## 4703       3954
    ## 4704       7468
    ## 4705       7469
    ## 4706       <NA>
    ## 4707       <NA>
    ## 4708     339983
    ## 4709     353497
    ## 4710      79441
    ## 4711      10608
    ## 4712      57732
    ## 4713       6047
    ## 4714       8603
    ## 4715      79155
    ## 4716       6452
    ## 4717        118
    ## 4718      10227
    ## 4719       8602
    ## 4720       2868
    ## 4721       3064
    ## 4722     345222
    ## 4723       6002
    ## 4724     285489
    ## 4725       4043
    ## 4726        152
    ## 4727       8532
    ## 4728     152992
    ## 4729       8310
    ## 4730      94031
    ## 4731      54436
    ## 4732      84448
    ## 4733      60312
    ## 4734      57537
    ## 4735       <NA>
    ## 4736      80273
    ## 4737      93624
    ## 4738     257236
    ## 4739      57533
    ## 4740       <NA>
    ## 4741      55330
    ## 4742      93621
    ## 4743      23324
    ## 4744       5522
    ## 4745       7466
    ## 4746     152789
    ## 4747       1400
    ## 4748       2121
    ## 4749     132884
    ## 4750      55351
    ## 4751       4487
    ## 4752      53407
    ## 4753      27065
    ## 4754     166793
    ## 4755      55646
    ## 4756      85013
    ## 4757       1816
    ## 4758      56606
    ## 4759       9948
    ## 4760       <NA>
    ## 4761       9957
    ## 4762       9364
    ## 4763       <NA>
    ## 4764       <NA>
    ## 4765     132864
    ## 4766     114905
    ## 4767      57545
    ## 4768      26234
    ## 4769        683
    ## 4770        952
    ## 4771       9982
    ## 4772       8842
    ## 4773     202018
    ## 4774       <NA>
    ## 4775       9079
    ## 4776       5860
    ## 4777      51056
    ## 4778      80306
    ## 4779      27146
    ## 4780       <NA>
    ## 4781     254251
    ## 4782       9353
    ## 4783       <NA>
    ## 4784      80333
    ## 4785     166647
    ## 4786      10891
    ## 4787       1665
    ## 4788       6649
    ## 4789      91050
    ## 4790      55203
    ## 4791      51091
    ## 4792      55300
    ## 4793      29063
    ## 4794      29945
    ## 4795      23231
    ## 4796     389203
    ## 4797       <NA>
    ## 4798       3516
    ## 4799      55296
    ## 4800      57620
    ## 4801       5099
    ## 4802     116984
    ## 4803      57495
    ## 4804       <NA>
    ## 4805     768211
    ## 4806       5236
    ## 4807      23216
    ## 4808      51274
    ## 4809       7096
    ## 4810      10333
    ## 4811      92689
    ## 4812      51088
    ## 4813      57728
    ## 4814       5981
    ## 4815       6133
    ## 4816      11019
    ## 4817       7358
    ## 4818     201895
    ## 4819       3093
    ## 4820      23244
    ## 4821      55728
    ## 4822        399
    ## 4823      54502
    ## 4824      79730
    ## 4825        323
    ## 4826       7345
    ## 4827      22998
    ## 4828      55161
    ## 4829      10463
    ## 4830     389206
    ## 4831       <NA>
    ## 4832     152573
    ## 4833      10396
    ## 4834       <NA>
    ## 4835     386617
    ## 4836      60558
    ## 4837     132789
    ## 4838       2565
    ## 4839       2555
    ## 4840     170712
    ## 4841       2557
    ## 4842       2560
    ## 4843       <NA>
    ## 4844      54951
    ## 4845     152518
    ## 4846     326340
    ## 4847       7006
    ## 4849      57606
    ## 4850     201780
    ## 4851     285527
    ## 4852      54940
    ## 4853     132299
    ## 4854      80157
    ## 4855      23142
    ## 4856       6443
    ## 4857      64854
    ## 4858      65997
    ## 4859     152579
    ## 4860      81608
    ## 4861      84708
    ## 4862      26511
    ## 4863       5156
    ## 4864       3815
    ## 4865       3791
    ## 4866      79644
    ## 4867      55858
    ## 4868       9575
    ## 4869       <NA>
    ## 4870      55763
    ## 4871       9662
    ## 4872       <NA>
    ## 4873     132949
    ## 4874       5471
    ## 4875      10606
    ## 4876       6731
    ## 4877      84525
    ## 4878       5978
    ## 4879      84273
    ## 4880       5431
    ## 4881       3490
    ## 4882      23284
    ## 4883       2044
    ## 4884       <NA>
    ## 4885      55236
    ## 4886      91746
    ## 4887       3512
    ## 4888      57050
    ## 4889      22902
    ## 4890       2926
    ## 4891      92597
    ## 4892       1633
    ## 4893       8671
    ## 4894       9508
    ## 4895     285521
    ## 4896      26057
    ## 4897       <NA>
    ## 4898        174
    ## 4899     166824
    ## 4900       6374
    ## 4901       5473
    ## 4902       5196
    ## 4903     441024
    ## 4904        374
    ## 4905      25849
    ## 4906      25898
    ## 4907       8999
    ## 4908       9908
    ## 4909       8615
    ## 4910       5470
    ## 4911      27163
    ## 4912      55153
    ## 4913        419
    ## 4914      53371
    ## 4915        950
    ## 4916       8987
    ## 4917     339965
    ## 4918      57619
    ## 4919     345079
    ## 4920       <NA>
    ## 4921      10983
    ## 4922        901
    ## 4923     246175
    ## 4924      65008
    ## 4925      80144
    ## 4926        306
    ## 4927       <NA>
    ## 4928      55589
    ## 4929     152559
    ## 4930     118429
    ## 4931      56978
    ## 4932       <NA>
    ## 4933       2250
    ## 4934       <NA>
    ## 4935        651
    ## 4936       5593
    ## 4937     153020
    ## 4938       <NA>
    ## 4939       <NA>
    ## 4940       <NA>
    ## 4941       3184
    ## 4942       <NA>
    ## 4943       <NA>
    ## 4944       9987
    ## 4945      58478
    ## 4946     441027
    ## 4947       <NA>
    ## 4948      22872
    ## 4949       <NA>
    ## 4950     132660
    ## 4951      51138
    ## 4952      51316
    ## 4953      27235
    ## 4954      10855
    ## 4955     113510
    ## 4956      51023
    ## 4957      84142
    ## 4958      84803
    ## 4959       <NA>
    ## 4960       1040
    ## 4961      23001
    ## 4962       <NA>
    ## 4963      83478
    ## 4964       5602
    ## 4965       5783
    ## 4966       4299
    ## 4967      57563
    ## 4968     345275
    ## 4969      51170
    ## 4970      53343
    ## 4971       <NA>
    ## 4972       8404
    ## 4973       <NA>
    ## 4974       1758
    ## 4975       3381
    ## 4976       6696
    ## 4977       5311
    ## 4978       <NA>
    ## 4979       <NA>
    ## 4980       <NA>
    ## 4981       <NA>
    ## 4982       <NA>
    ## 4983     115361
    ## 4984      23507
    ## 4985      84230
    ## 4986      55144
    ## 4987       <NA>
    ## 4988     343472
    ## 4989       <NA>
    ## 4990       <NA>
    ## 4991     164045
    ## 4992       8317
    ## 4993       7049
    ## 4994        676
    ## 4995       <NA>
    ## 4996     253152
    ## 4997       <NA>
    ## 4998       <NA>
    ## 4999       <NA>
    ## 5000      11146
    ## 5001      79871
    ## 5002       7813
    ## 5003       6125
    ## 5004       <NA>
    ## 5005      22823
    ## 5006      50999
    ## 5007     343099
    ## 5008       1810
    ## 5009      54872
    ## 5010       <NA>
    ## 5011       <NA>
    ## 5012       <NA>
    ## 5013       <NA>
    ## 5014       <NA>
    ## 5015       <NA>
    ## 5016      10336
    ## 5017      10815
    ## 5018       2580
    ## 5019      84286
    ## 5020       1609
    ## 5021       3425
    ## 5022      10861
    ## 5023      53834
    ## 5024      64109
    ## 5025       <NA>
    ## 5026       <NA>
    ## 5027       <NA>
    ## 5028       <NA>
    ## 5029       <NA>
    ## 5030       <NA>
    ## 5031       <NA>
    ## 5032       <NA>
    ## 5033      55344
    ## 5034       8225
    ## 5035       <NA>
    ## 5036      55743
    ## 5037       2802
    ## 5038       <NA>
    ## 5039      23141
    ## 5040     192111
    ## 5041       5827
    ## 5042       5426
    ## 5043      57666
    ## 5044       <NA>
    ## 5045       <NA>
    ## 5046      50614
    ## 5047      79050
    ## 5048     317781
    ## 5049       <NA>
    ## 5050      57634
    ## 5051       <NA>
    ## 5052      80324
    ## 5053       <NA>
    ## 5054       8408
    ## 5055     150274
    ## 5056      11200
    ## 5057      23331
    ## 5058      23760
    ## 5059       <NA>
    ## 5060       4330
    ## 5061       <NA>
    ## 5062       <NA>
    ## 5063     440823
    ## 5064       <NA>
    ## 5065       1413
    ## 5066       1414
    ## 5067       8459
    ## 5068      24144
    ## 5069     402055
    ## 5070      89781
    ## 5071      57168
    ## 5072      23544
    ## 5073       <NA>
    ## 5074       <NA>
    ## 5075       <NA>
    ## 5076      84700
    ## 5077        157
    ## 5078       1417
    ## 5079       <NA>
    ## 5080     129049
    ## 5081       9671
    ## 5082       1240
    ## 5083      11153
    ## 5084       9733
    ## 5085      23479
    ## 5086     338773
    ## 5087       6404
    ## 5088      23603
    ## 5089      54434
    ## 5090      55530
    ## 5091      84749
    ## 5092     121642
    ## 5093       7374
    ## 5094         32
    ## 5095     283446
    ## 5096      83892
    ## 5097      89910
    ## 5098     326625
    ## 5099       4598
    ## 5100      84915
    ## 5101      59341
    ## 5102      51228
    ## 5103      84260
    ## 5104       9815
    ## 5105      88455
    ## 5106       <NA>
    ## 5107       <NA>
    ## 5108       <NA>
    ## 5109       <NA>
    ## 5110       <NA>
    ## 5111     121665
    ## 5112       <NA>
    ## 5113         35
    ## 5114      84747
    ## 5115       9761
    ## 5116       9478
    ## 5117      51367
    ## 5118       9921
    ## 5119      84274
    ## 5120       8655
    ## 5121       <NA>
    ## 5122       8683
    ## 5123     283459
    ## 5124      51499
    ## 5125       1337
    ## 5126       4440
    ## 5127      23409
    ## 5128       5829
    ## 5129       6175
    ## 5130       <NA>
    ## 5131       <NA>
    ## 5132      11021
    ## 5133      92558
    ## 5134      11113
    ## 5135       5564
    ## 5136       <NA>
    ## 5137     160777
    ## 5138      26353
    ## 5139      84530
    ## 5140      64426
    ## 5141      51347
    ## 5142       5037
    ## 5143      54621
    ## 5144      55884
    ## 5145       5985
    ## 5146       <NA>
    ## 5147     283455
    ## 5148       4842
    ## 5149      23014
    ## 5150      54997
    ## 5151      26259
    ## 5152       8739
    ## 5153      84900
    ## 5154       <NA>
    ## 5155      23389
    ## 5156       <NA>
    ## 5157       6926
    ## 5158       9904
    ## 5159       <NA>
    ## 5160     113675
    ## 5161       <NA>
    ## 5162     196463
    ## 5163      80024
    ## 5164      53373
    ## 5165     115811
    ## 5166      84934
    ## 5167      79039
    ## 5168       8437
    ## 5169       1840
    ## 5170       4939
    ## 5171       <NA>
    ## 5172       <NA>
    ## 5173      22895
    ## 5174       5781
    ## 5175       6128
    ## 5176     283450
    ## 5177      10906
    ## 5178       <NA>
    ## 5179      80018
    ## 5180      10961
    ## 5181      89894
    ## 5182       8759
    ## 5183       8550
    ## 5184        217
    ## 5185       <NA>
    ## 5186      80724
    ## 5187       8315
    ## 5188       6311
    ## 5189      10019
    ## 5190       <NA>
    ## 5191      23316
    ## 5192       4633
    ## 5193     160762
    ## 5194       5501
    ## 5195      84329
    ## 5196      79600
    ## 5197     160760
    ## 5198     144715
    ## 5199      51699
    ## 5200      29902
    ## 5201      51184
    ## 5202      10094
    ## 5203      51434
    ## 5204        488
    ## 5205      28981
    ## 5206       5027
    ## 5207       5025
    ## 5208      10645
    ## 5209      51433
    ## 5210      80196
    ## 5211      84678
    ## 5212       <NA>
    ## 5213      84876
    ## 5214     144404
    ## 5215      54509
    ## 5216       <NA>
    ## 5217       <NA>
    ## 5218      23067
    ## 5219       5715
    ## 5220     144406
    ## 5221       <NA>
    ## 5222        605
    ## 5223      22877
    ## 5224      56616
    ## 5225      79369
    ## 5226      65082
    ## 5227       6249
    ## 5228      55596
    ## 5229      65117
    ## 5230       9735
    ## 5231      27198
    ## 5232       8562
    ## 5233      84660
    ## 5234       9026
    ## 5235      79720
    ## 5236      23457
    ## 5237      79676
    ## 5238      51329
    ## 5239      57605
    ## 5240      10198
    ## 5241       <NA>
    ## 5242       8099
    ## 5243      55206
    ## 5244     387893
    ## 5245     196383
    ## 5246      11066
    ## 5247     353116
    ## 5248      10959
    ## 5249      57696
    ## 5250       1967
    ## 5251       2967
    ## 5252      23545
    ## 5253     196385
    ## 5254      80212
    ## 5255       <NA>
    ## 5256     144347
    ## 5257       9612
    ## 5258        949
    ## 5259       7316
    ## 5260       <NA>
    ## 5261      57647
    ## 5262     140707
    ## 5263      65985
    ## 5264     114795
    ## 5265       <NA>
    ## 5266      92293
    ## 5267     121260
    ## 5268     144423
    ## 5269     121256
    ## 5270       <NA>
    ## 5271      11211
    ## 5272       9271
    ## 5273      23504
    ## 5274       2054
    ## 5275       5901
    ## 5276     283383
    ## 5277       6433
    ## 5278       4326
    ## 5279       <NA>
    ## 5280      51373
    ## 5281       2631
    ## 5282       5723
    ## 5283      25870
    ## 5284       5260
    ## 5285      51142
    ## 5286      58486
    ## 5287       <NA>
    ## 5288     154807
    ## 5289       2990
    ## 5290        435
    ## 5291      27297
    ## 5292       8460
    ## 5293     154881
    ## 5294      27342
    ## 5295      55069
    ## 5296       <NA>
    ## 5297      51119
    ## 5298      55253
    ## 5299      83698
    ## 5300       <NA>
    ## 5301      64409
    ## 5302      26053
    ## 5303       <NA>
    ## 5304       <NA>
    ## 5305      81554
    ## 5306      84163
    ## 5307     653361
    ## 5308       2969
    ## 5309       9569
    ## 5310       7461
    ## 5311       5982
    ## 5312       7462
    ## 5313       7458
    ## 5314       3984
    ## 5315       2006
    ## 5316       <NA>
    ## 5317     155368
    ## 5318       <NA>
    ## 5319      83451
    ## 5320       6804
    ## 5321     114049
    ## 5322      84277
    ## 5323     155382
    ## 5324      51085
    ## 5325      26608
    ## 5326       9275
    ## 5327       9031
    ## 5328       8326
    ## 5329      55695
    ## 5330       9883
    ## 5331       3092
    ## 5332       6369
    ## 5333      57414
    ## 5334       5447
    ## 5335      83862
    ## 5336       4191
    ## 5337     222183
    ## 5338       3315
    ## 5339       7532
    ## 5340     136853
    ## 5341     113878
    ## 5342       <NA>
    ## 5343      10156
    ## 5344       5439
    ## 5345       <NA>
    ## 5346     222229
    ## 5347      54784
    ## 5348      80228
    ## 5349      79706
    ## 5350      10603
    ## 5351       1523
    ## 5352     136227
    ## 5353      64792
    ## 5354       <NA>
    ## 5355      51024
    ## 5356      10467
    ## 5357       8985
    ## 5358       7425
    ## 5359       1174
    ## 5360       5054
    ## 5361      81844
    ## 5362       4584
    ## 5363         43
    ## 5364     402682
    ## 5365      51593
    ## 5366       7205
    ## 5367      56996
    ## 5368       2050
    ## 5369      10248
    ## 5370       <NA>
    ## 5371      64599
    ## 5372       2783
    ## 5373      51412
    ## 5374       7036
    ## 5375      64598
    ## 5376       5118
    ## 5377       <NA>
    ## 5378       4034
    ## 5379       <NA>
    ## 5380  100316904
    ## 5381       <NA>
    ## 5382       3268
    ## 5383     222950
    ## 5384      81628
    ## 5385     221908
    ## 5386      56257
    ## 5387      55063
    ## 5388      29992
    ## 5389       <NA>
    ## 5390     349149
    ## 5391       7586
    ## 5392       7589
    ## 5393       <NA>
    ## 5394      10980
    ## 5395       4176
    ## 5396       9179
    ## 5397       6878
    ## 5398     245812
    ## 5399     255374
    ## 5400     389541
    ## 5401       <NA>
    ## 5402       <NA>
    ## 5403      79690
    ## 5404      10734
    ## 5405       <NA>
    ## 5406       <NA>
    ## 5407       <NA>
    ## 5408      56975
    ## 5409       5154
    ## 5410       <NA>
    ## 5411       5575
    ## 5412      54919
    ## 5413      23353
    ## 5414      51608
    ## 5415      11033
    ## 5416      90639
    ## 5417       <NA>
    ## 5418     115330
    ## 5419       2852
    ## 5420      90637
    ## 5421      79778
    ## 5422      26173
    ## 5423       7975
    ## 5424      84262
    ## 5425     392617
    ## 5426       8379
    ## 5427       4521
    ## 5428      29960
    ## 5429      29886
    ## 5430       8662
    ## 5431      55501
    ## 5432       3955
    ## 5433      80727
    ## 5434      23288
    ## 5435     221927
    ## 5436     155185
    ## 5437       2768
    ## 5438     221935
    ## 5439     221937
    ## 5440       9907
    ## 5441      55698
    ## 5442     221938
    ## 5444      26100
    ## 5445     222962
    ## 5446      84629
    ## 5447      80028
    ## 5448         60
    ## 5449       6624
    ## 5450      54476
    ## 5451      57786
    ## 5452       <NA>
    ## 5453       <NA>
    ## 5454       <NA>
    ## 5455       <NA>
    ## 5456       <NA>
    ## 5457      55146
    ## 5458       <NA>
    ## 5459     392862
    ## 5460       <NA>
    ## 5461      11014
    ## 5462     221955
    ## 5463       5879
    ## 5464      84792
    ## 5465       9265
    ## 5466      84132
    ## 5467       <NA>
    ## 5468      27102
    ## 5469       7965
    ## 5470       5395
    ## 5471     222967
    ## 5472      51622
    ## 5473     654231
    ## 5474      22853
    ## 5475      25851
    ## 5476       <NA>
    ## 5477      25798
    ## 5478      55971
    ## 5479       4885
    ## 5480     222865
    ## 5481       8295
    ## 5482      57154
    ## 5483      10552
    ## 5484      10095
    ## 5485      11333
    ## 5486       8896
    ## 5487      26024
    ## 5488      10898
    ## 5489       <NA>
    ## 5490       <NA>
    ## 5491      23660
    ## 5492       <NA>
    ## 5493       <NA>
    ## 5494       <NA>
    ## 5495       6049
    ## 5496       1024
    ## 5497      10810
    ## 5498       2835
    ## 5499     219333
    ## 5500       6144
    ## 5501     387496
    ## 5502       2971
    ## 5503     219402
    ## 5504     222484
    ## 5505      51082
    ## 5506     219409
    ## 5507       2322
    ## 5508     255967
    ## 5509       2321
    ## 5510      51371
    ## 5511     283537
    ## 5512       <NA>
    ## 5513      23281
    ## 5514       6541
    ## 5515       5412
    ## 5516      84056
    ## 5517       <NA>
    ## 5518       <NA>
    ## 5519       3146
    ## 5520       <NA>
    ## 5521       <NA>
    ## 5522      10208
    ## 5523       <NA>
    ## 5524        241
    ## 5525       <NA>
    ## 5526      84935
    ## 5527     122046
    ## 5528       <NA>
    ## 5529       <NA>
    ## 5530      10808
    ## 5531     145173
    ## 5532     122042
    ## 5533      10129
    ## 5534       <NA>
    ## 5535     646799
    ## 5536        675
    ## 5537      90634
    ## 5538       <NA>
    ## 5539       <NA>
    ## 5540      10443
    ## 5541      23047
    ## 5542       <NA>
    ## 5543       9365
    ## 5544      90627
    ## 5545       <NA>
    ## 5546       5983
    ## 5547     219285
    ## 5548      55610
    ## 5549       2791
    ## 5550      10282
    ## 5551       1278
    ## 5552      64921
    ## 5553       8910
    ## 5554      23089
    ## 5555      55607
    ## 5556       5446
    ## 5557       5445
    ## 5558      51666
    ## 5559       5166
    ## 5560       1780
    ## 5561      10165
    ## 5562       7979
    ## 5563       <NA>
    ## 5564       1750
    ## 5565       1749
    ## 5566      57001
    ## 5567       6863
    ## 5568        440
    ## 5569      56913
    ## 5570      54468
    ## 5571       6119
    ## 5572     729852
    ## 5573     113263
    ## 5574       3382
    ## 5575      30010
    ## 5576       4697
    ## 5577       9678
    ## 5578     221981
    ## 5579      54664
    ## 5580      64418
    ## 5581     154743
    ## 5582      54329
    ## 5583       <NA>
    ## 5584       <NA>
    ## 5585      93986
    ## 5586      29969
    ## 5587       <NA>
    ## 5588      26136
    ## 5589        858
    ## 5590        857
    ## 5591       4233
    ## 5592        830
    ## 5593       7982
    ## 5594       7472
    ## 5595       1080
    ## 5596      83992
    ## 5597      51691
    ## 5598       3751
    ## 5599      23554
    ## 5600      54556
    ## 5601      79974
    ## 5602      10447
    ## 5603       5803
    ## 5604      10157
    ## 5605      93664
    ## 5606       6561
    ## 5607     154865
    ## 5608       4698
    ## 5609       8976
    ## 5610     730130
    ## 5611       2861
    ## 5612       <NA>
    ## 5613       2918
    ## 5614       <NA>
    ## 5615       <NA>
    ## 5616      79571
    ## 5617        381
    ## 5618      27044
    ## 5619      64101
    ## 5620       <NA>
    ## 5621      55131
    ## 5622     401399
    ## 5623       3614
    ## 5624      29923
    ## 5625     346653
    ## 5626        813
    ## 5627      64753
    ## 5628       2318
    ## 5629       9296
    ## 5630     375616
    ## 5631       3663
    ## 5632      23534
    ## 5633     340348
    ## 5634       6608
    ## 5635      23382
    ## 5636      57464
    ## 5637       4899
    ## 5638       7328
    ## 5639      51530
    ## 5640      23008
    ## 5641       <NA>
    ## 5642      84928
    ## 5643       1358
    ## 5644      95681
    ## 5645       4232
    ## 5646      26958
    ## 5647     136259
    ## 5648       <NA>
    ## 5649       4289
    ## 5650       5420
    ## 5651       <NA>
    ## 5652       <NA>
    ## 5653       <NA>
    ## 5654      91584
    ## 5655       <NA>
    ## 5656      54927
    ## 5657      60412
    ## 5658     136332
    ## 5659      84912
    ## 5660       <NA>
    ## 5661       <NA>
    ## 5662      57016
    ## 5663        669
    ## 5664        800
    ## 5665       <NA>
    ## 5666     340351
    ## 5667      55281
    ## 5668       <NA>
    ## 5669      29062
    ## 5670       <NA>
    ## 5671       4850
    ## 5672      23165
    ## 5673       <NA>
    ## 5674      26266
    ## 5675     389558
    ## 5676     136319
    ## 5677       <NA>
    ## 5678       1129
    ## 5679       5764
    ## 5680       9162
    ## 5681      64764
    ## 5682       8805
    ## 5683      50617
    ## 5684       <NA>
    ## 5685      92092
    ## 5686      56829
    ## 5687      79989
    ## 5688     254048
    ## 5689     154791
    ## 5690      51631
    ## 5691     154790
    ## 5692      28996
    ## 5693       6916
    ## 5694      64761
    ## 5695      80853
    ## 5696      84255
    ## 5697      23608
    ## 5698       <NA>
    ## 5699      27147
    ## 5700      90956
    ## 5701       4708
    ## 5702        673
    ## 5703      51650
    ## 5704  100507421
    ## 5705      55750
    ## 5706       <NA>
    ## 5707       6742
    ## 5708      23601
    ## 5709       <NA>
    ## 5710       <NA>
    ## 5711       <NA>
    ## 5712      28638
    ## 5713       2051
    ## 5714      55503
    ## 5715     373156
    ## 5716        835
    ## 5717       1180
    ## 5718       9715
    ## 5719       <NA>
    ## 5720       7791
    ## 5721     285966
    ## 5722       9747
    ## 5723       7984
    ## 5724      27010
    ## 5725      26047
    ## 5726       8454
    ## 5727       2146
    ## 5728       9601
    ## 5729       <NA>
    ## 5730       <NA>
    ## 5731       <NA>
    ## 5732       <NA>
    ## 5733       <NA>
    ## 5734       <NA>
    ## 5735       <NA>
    ## 5736       <NA>
    ## 5737      84626
    ## 5738       <NA>
    ## 5739      23145
    ## 5740       <NA>
    ## 5741     155066
    ## 5742      65999
    ## 5743       5919
    ## 5744       <NA>
    ## 5745      29803
    ## 5746       <NA>
    ## 5747       <NA>
    ## 5748     155038
    ## 5749       <NA>
    ## 5750      55303
    ## 5751     474344
    ## 5752     170575
    ## 5753      55340
    ## 5754      28959
    ## 5755      55365
    ## 5756       <NA>
    ## 5757      10457
    ## 5758     115416
    ## 5759      29896
    ## 5760      90693
    ## 5761     340277
    ## 5762       4852
    ## 5763      51678
    ## 5764       1687
    ## 5765      26031
    ## 5766      54205
    ## 5767       <NA>
    ## 5768       <NA>
    ## 5769       9603
    ## 5770       3181
    ## 5771      11335
    ## 5772      29887
    ## 5773       8935
    ## 5774      11112
    ## 5775       8887
    ## 5776     221895
    ## 5777       9586
    ## 5778       9865
    ## 5779       <NA>
    ## 5780       <NA>
    ## 5781       1124
    ## 5782     222171
    ## 5783     644150
    ## 5784       9805
    ## 5785      55033
    ## 5786      84725
    ## 5787     222166
    ## 5788     223082
    ## 5789      10392
    ## 5790      79017
    ## 5791       <NA>
    ## 5792       1395
    ## 5793      11185
    ## 5794      84182
    ## 5795        358
    ## 5796       <NA>
    ## 5797        117
    ## 5798      63974
    ## 5799       <NA>
    ## 5800       <NA>
    ## 5801      10842
    ## 5802       5137
    ## 5803      23658
    ## 5804      23080
    ## 5805      25948
    ## 5806      11328
    ## 5807       <NA>
    ## 5808     152926
    ## 5809      55008
    ## 5810  100996939
    ## 5811      55915
    ## 5812      81552
    ## 5813       9429
    ## 5814       8916
    ## 5815     266812
    ## 5816      10144
    ## 5817     166815
    ## 5818     285513
    ## 5819       6622
    ## 5820     401145
    ## 5821       2895
    ## 5822      56916
    ## 5823      27306
    ## 5824      79931
    ## 5825      79625
    ## 5826      11107
    ## 5827       4085
    ## 5828      55970
    ## 5829       1647
    ## 5830      26135
    ## 5831       4070
    ## 5832       3514
    ## 5833      22934
    ## 5834       9451
    ## 5835      55258
    ## 5836     150572
    ## 5837      51315
    ## 5838       <NA>
    ## 5839      64795
    ## 5840       7844
    ## 5841      51652
    ## 5842      55818
    ## 5843      65055
    ## 5844      51318
    ## 5845      10989
    ## 5846      55037
    ## 5847      25885
    ## 5848       8869
    ## 5849      84913
    ## 5850      10713
    ## 5851       <NA>
    ## 5852       <NA>
    ## 5853     129303
    ## 5854      51255
    ## 5855      10791
    ## 5856       8673
    ## 5857       2677
    ## 5858       4144
    ## 5859     284948
    ## 5860        822
    ## 5861      84173
    ## 5862      54884
    ## 5863       <NA>
    ## 5864      83439
    ## 5865      56888
    ## 5866       9168
    ## 5867       1768
    ## 5868       8802
    ## 5869       <NA>
    ## 5870       1496
    ## 5871     347730
    ## 5872      80059
    ## 5873       6936
    ## 5874       9801
    ## 5875       <NA>
    ## 5876      84141
    ## 5877       6869
    ## 5878      56655
    ## 5879       3099
    ## 5880       <NA>
    ## 5881      10505
    ## 5882       1796
    ## 5883      84695
    ## 5884      27429
    ## 5885        550
    ## 5886     165545
    ## 5887      84759
    ## 5888      85474
    ## 5889       <NA>
    ## 5890       <NA>
    ## 5891     116540
    ## 5892       7841
    ## 5893      23559
    ## 5894      83444
    ## 5895       6242
    ## 5896      84058
    ## 5897       <NA>
    ## 5898       1639
    ## 5899       <NA>
    ## 5900      10797
    ## 5901      55233
    ## 5902     388962
    ## 5903     200424
    ## 5904       1716
    ## 5905      10617
    ## 5906       <NA>
    ## 5907     113419
    ## 5908      55577
    ## 5909     400961
    ## 5910       <NA>
    ## 5911       8291
    ## 5912      56603
    ## 5913      23233
    ## 5914       6697
    ## 5915       2016
    ## 5916      94097
    ## 5917      26056
    ## 5918      10322
    ## 5919      84279
    ## 5920      10574
    ## 5921     150726
    ## 5922       1961
    ## 5923       7840
    ## 5924       <NA>
    ## 5925       9027
    ## 5926       <NA>
    ## 5927       <NA>
    ## 5928       <NA>
    ## 5929      51002
    ## 5930       8446
    ## 5931       <NA>
    ## 5932        119
    ## 5933       7039
    ## 5934      84908
    ## 5935       6637
    ## 5936      51449
    ## 5937       7072
    ## 5938       <NA>
    ## 5939       <NA>
    ## 5940       <NA>
    ## 5941       5093
    ## 5942       <NA>
    ## 5943     151516
    ## 5944       4084
    ## 5945      11017
    ## 5946       <NA>
    ## 5947      64395
    ## 5948        307
    ## 5949      22848
    ## 5950      27247
    ## 5951       2673
    ## 5952      84168
    ## 5953       <NA>
    ## 5954       9938
    ## 5955     200558
    ## 5956       <NA>
    ## 5957       <NA>
    ## 5958      79825
    ## 5959       2815
    ## 5960     339122
    ## 5961      57461
    ## 5962       7555
    ## 5963      22820
    ## 5964      56941
    ## 5965       <NA>
    ## 5966       <NA>
    ## 5967       <NA>
    ## 5968       6184
    ## 5969       2624
    ## 5970      60678
    ## 5971       8607
    ## 5972      29927
    ## 5973     166348
    ## 5974       <NA>
    ## 5975      11343
    ## 5976      80325
    ## 5977      50512
    ## 5978       <NA>
    ## 5979       4171
    ## 5980     131601
    ## 5981       5361
    ## 5982      84303
    ## 5983     114112
    ## 5984      79364
    ## 5985     348807
    ## 5986      28999
    ## 5987       <NA>
    ## 5988       <NA>
    ## 5989      10840
    ## 5990      54946
    ## 5991       9922
    ## 5992      23225
    ## 5993      79885
    ## 5994       2199
    ## 5995       7476
    ## 5996       <NA>
    ## 5997     131474
    ## 5998      79188
    ## 5999       7508
    ## 6000      27258
    ## 6001       6533
    ## 6002      80852
    ## 6003      51244
    ## 6004     152273
    ## 6005       7182
    ## 6006      64432
    ## 6007      64145
    ## 6008       7200
    ## 6009     166336
    ## 6010      56999
    ## 6011       9223
    ## 6012     115286
    ## 6013      26018
    ## 6014      84541
    ## 6015       8801
    ## 6016       <NA>
    ## 6017       <NA>
    ## 6018     285203
    ## 6019       7110
    ## 6020       9039
    ## 6021      10550
    ## 6022      23150
    ## 6023       4286
    ## 6024       <NA>
    ## 6025      27086
    ## 6026     317649
    ## 6027       2850
    ## 6028      23429
    ## 6029      55164
    ## 6030     727936
    ## 6031     151987
    ## 6032      23024
    ## 6033       5067
    ## 6034      10752
    ## 6035      27255
    ## 6036       <NA>
    ## 6037     152330
    ## 6038      51095
    ## 6039      51185
    ## 6040      57633
    ## 6041       6419
    ## 6042     285362
    ## 6043       3708
    ## 6044       <NA>
    ## 6045       8553
    ## 6046      55207
    ## 6047       9695
    ## 6048       <NA>
    ## 6049       2917
    ## 6050      29995
    ## 6051       <NA>
    ## 6052        859
    ## 6053       5021
    ## 6054      56852
    ## 6055       <NA>
    ## 6056       9901
    ## 6057      25917
    ## 6058       <NA>
    ## 6059      55209
    ## 6060     375323
    ## 6061      64419
    ## 6062     151835
    ## 6063       7862
    ## 6064       <NA>
    ## 6065       4968
    ## 6066       8536
    ## 6067      10474
    ## 6068      10093
    ## 6069      26140
    ## 6070     285367
    ## 6071      84522
    ## 6072      84818
    ## 6073      78987
    ## 6074     285368
    ## 6075      55831
    ## 6076       2177
    ## 6077      55845
    ## 6078       7428
    ## 6079       3656
    ## 6080       9797
    ## 6081       <NA>
    ## 6082       <NA>
    ## 6083       6396
    ## 6084        491
    ## 6085       <NA>
    ## 6086       6538
    ## 6087       6529
    ## 6088       3269
    ## 6089      10533
    ## 6090       9686
    ## 6091     132001
    ## 6092       6854
    ## 6093       7079
    ## 6094       5468
    ## 6095      80746
    ## 6096      23609
    ## 6097       5894
    ## 6098      55287
    ## 6099      23066
    ## 6100       6161
    ## 6101      90288
    ## 6102       8930
    ## 6103      55764
    ## 6104      23129
    ## 6105      23023
    ## 6106       <NA>
    ## 6107      93550
    ## 6108       <NA>
    ## 6109        240
    ## 6110       <NA>
    ## 6111       <NA>
    ## 6112      83937
    ## 6113       <NA>
    ## 6114       6387
    ## 6115       <NA>
    ## 6116       <NA>
    ## 6117       <NA>
    ## 6118       3185
    ## 6119     221002
    ## 6120      55454
    ## 6121       5979
    ## 6122       9790
    ## 6123       <NA>
    ## 6124       <NA>
    ## 6125      22852
    ## 6126        775
    ## 6127     196513
    ## 6128      93589
    ## 6129     654429
    ## 6130      79602
    ## 6131      81029
    ## 6132     144699
    ## 6133      23085
    ## 6134       5893
    ## 6135      65125
    ## 6136       4815
    ## 6137      84318
    ## 6138       5927
    ## 6139      23765
    ## 6140      27439
    ## 6141      27440
    ## 6142       <NA>
    ## 6143      27443
    ## 6144      83733
    ## 6145        529
    ## 6146      23786
    ## 6147        637
    ## 6148      57553
    ## 6149       <NA>
    ## 6150      55670
    ## 6151      51807
    ## 6152      11274
    ## 6153       6540
    ## 6154       6539
    ## 6155     440073
    ## 6156          2
    ## 6157      10219
    ## 6158       4074
    ## 6159       1911
    ## 6160      57494
    ## 6161        339
    ## 6162       6515
    ## 6163      55810
    ## 6164        719
    ## 6165      25977
    ## 6166       <NA>
    ## 6167       <NA>
    ## 6168       <NA>
    ## 6169       9332
    ## 6170       5830
    ## 6171       9746
    ## 6172      51279
    ## 6173       <NA>
    ## 6174       <NA>
    ## 6175      10162
    ## 6176      10436
    ## 6177      11331
    ## 6178       5777
    ## 6179       <NA>
    ## 6180       1822
    ## 6181       2026
    ## 6182      10233
    ## 6183      84727
    ## 6184       7167
    ## 6185       8078
    ## 6186      83461
    ## 6187      10536
    ## 6188      27239
    ## 6189        920
    ## 6190       3902
    ## 6191       5763
    ## 6192       <NA>
    ## 6193       8079
    ## 6194      50813
    ## 6195     196500
    ## 6196       <NA>
    ## 6197      51147
    ## 6198      84519
    ## 6199       1108
    ## 6200       4839
    ## 6201      25900
    ## 6202       2597
    ## 6203       9918
    ## 6204      51258
    ## 6205       6843
    ## 6206      55080
    ## 6207       4055
    ## 6208       6337
    ## 6209       7132
    ## 6210        928
    ## 6211       7450
    ## 6212      57101
    ## 6213       4908
    ## 6214       3741
    ## 6215       3736
    ## 6216       3742
    ## 6217       <NA>
    ## 6218       4704
    ## 6219      10566
    ## 6220      10635
    ## 6221       <NA>
    ## 6222      57103
    ## 6223        894
    ## 6224       <NA>
    ## 6225      57097
    ## 6226      84766
    ## 6227      56341
    ## 6228       <NA>
    ## 6229     441631
    ## 6230      10867
    ## 6231       <NA>
    ## 6232       <NA>
    ## 6233       7289
    ## 6234      83695
    ## 6235       2305
    ## 6236  101929469
    ## 6237      83714
    ## 6238      55846
    ## 6239       2288
    ## 6240       <NA>
    ## 6241     144568
    ## 6242       <NA>
    ## 6243      29121
    ## 6244       <NA>
    ## 6245     160364
    ## 6246      23710
    ## 6247      55110
    ## 6248       8531
    ## 6249  100129361
    ## 6250       2120
    ## 6251       4040
    ## 6252       <NA>
    ## 6253      54682
    ## 6254     118426
    ## 6255      80824
    ## 6256       1389
    ## 6257       2842
    ## 6258       <NA>
    ## 6259       1027
    ## 6260       <NA>
    ## 6261      81575
    ## 6262      51202
    ## 6263      55507
    ## 6264      50865
    ## 6265      57613
    ## 6266      83445
    ## 6267       <NA>
    ## 6268       2012
    ## 6269       2904
    ## 6270       <NA>
    ## 6271      55729
    ## 6272       2984
    ## 6273       <NA>
    ## 6274      51729
    ## 6275     440087
    ## 6276       4256
    ## 6277     121506
    ## 6278        397
    ## 6279       5149
    ## 6280      85004
    ## 6281       5800
    ## 6282       2059
    ## 6283      11171
    ## 6284      51071
    ## 6285       4257
    ## 6286      55885
    ## 6287       <NA>
    ## 6288       <NA>
    ## 6289      89869
    ## 6290      54477
    ## 6291     121536
    ## 6292       <NA>
    ## 6293       5139
    ## 6294      53919
    ## 6295       <NA>
    ## 6296       <NA>
    ## 6297      79912
    ## 6298       5965
    ## 6299      51026
    ## 6300      80763
    ## 6301       3945
    ## 6302       3764
    ## 6303      10060
    ## 6304      55907
    ## 6305       <NA>
    ## 6306       6489
    ## 6307       9847
    ## 6308      55500
    ## 6309       <NA>
    ## 6310       6660
    ## 6311        586
    ## 6312       <NA>
    ## 6313       <NA>
    ## 6314       4033
    ## 6315      55259
    ## 6316     144363
    ## 6317       3845
    ## 6318       <NA>
    ## 6319       <NA>
    ## 6320      11228
    ## 6321      79365
    ## 6322       8082
    ## 6323       3709
    ## 6324      55726
    ## 6325      26127
    ## 6326      51768
    ## 6327       9412
    ## 6328       <NA>
    ## 6329      23012
    ## 6330      56938
    ## 6331     341346
    ## 6332       8496
    ## 6333     387849
    ## 6334      60488
    ## 6335      57542
    ## 6336       5744
    ## 6337      55297
    ## 6338      55711
    ## 6339       <NA>
    ## 6340      51290
    ## 6341      83857
    ## 6342      10526
    ## 6343      65981
    ## 6344       <NA>
    ## 6345     160518
    ## 6346     254013
    ## 6347     196394
    ## 6348       <NA>
    ## 6349        636
    ## 6350       <NA>
    ## 6351      91663
    ## 6352       5582
    ## 6353      59284
    ## 6354      59283
    ## 6355       <NA>
    ## 6356      59285
    ## 6357     126014
    ## 6358       4696
    ## 6359      29844
    ## 6360      26121
    ## 6361       4849
    ## 6362      79165
    ## 6363     147798
    ## 6364      79143
    ## 6365      79042
    ## 6366       6203
    ## 6367       <NA>
    ## 6368       3903
    ## 6369       <NA>
    ## 6370      57348
    ## 6371     114823
    ## 6372      94059
    ## 6373     148170
    ## 6374     353514
    ## 6375     112724
    ## 6376      54869
    ## 6377      54776
    ## 6378       7138
    ## 6379       7137
    ## 6380     352909
    ## 6381       6861
    ## 6382       5794
    ## 6383     255043
    ## 6384      22870
    ## 6385      23640
    ## 6386      84446
    ## 6387     284417
    ## 6388      84787
    ## 6389     125965
    ## 6390       3589
    ## 6391     388564
    ## 6392       6158
    ## 6393      27338
    ## 6394     729956
    ## 6395       <NA>
    ## 6396       <NA>
    ## 6397       <NA>
    ## 6398      57106
    ## 6399     284297
    ## 6400  100130827
    ## 6401       <NA>
    ## 6402      84922
    ## 6403       <NA>
    ## 6404       <NA>
    ## 6405      29903
    ## 6406       <NA>
    ## 6407       <NA>
    ## 6408       <NA>
    ## 6409      11338
    ## 6410      29924
    ## 6411       <NA>
    ## 6412       <NA>
    ## 6413       <NA>
    ## 6414       <NA>
    ## 6415       <NA>
    ## 6416       <NA>
    ## 6417       <NA>
    ## 6418       <NA>
    ## 6419       <NA>
    ## 6420     140612
    ## 6421     147670
    ## 6422       <NA>
    ## 6423       5178
    ## 6424      57663
    ## 6425       <NA>
    ## 6426       <NA>
    ## 6427       <NA>
    ## 6428       <NA>
    ## 6429       <NA>
    ## 6430       1183
    ## 6431       <NA>
    ## 6432     284307
    ## 6433       <NA>
    ## 6434      65982
    ## 6435       <NA>
    ## 6436       <NA>
    ## 6437       <NA>
    ## 6438     342945
    ## 6439       6193
    ## 6440     646862
    ## 6441       <NA>
    ## 6442       <NA>
    ## 6443      10998
    ## 6444      10155
    ## 6445      84878
    ## 6446      27243
    ## 6447       9040
    ## 6448       7593
    ## 6449     374920
    ## 6450       3978
    ## 6451      56344
    ## 6452       <NA>
    ## 6453       6415
    ## 6454      29997
    ## 6455      30846
    ## 6456      29998
    ## 6457       <NA>
    ## 6458       8775
    ## 6459      11133
    ## 6460       6543
    ## 6461      56917
    ## 6462       9704
    ## 6463      27202
    ## 6464        728
    ## 6465     255783
    ## 6466      26093
    ## 6467      27113
    ## 6468      10055
    ## 6469      23211
    ## 6470       <NA>
    ## 6471      54958
    ## 6472       4861
    ## 6473       2909
    ## 6474       <NA>
    ## 6475       1175
    ## 6476       6510
    ## 6477      79147
    ## 6478      29888
    ## 6479      25865
    ## 6480       <NA>
    ## 6481     147906
    ## 6482      94235
    ## 6483        808
    ## 6484       <NA>
    ## 6485       <NA>
    ## 6486       <NA>
    ## 6487      83987
    ## 6488       5536
    ## 6489      64344
    ## 6490       <NA>
    ## 6491       <NA>
    ## 6492       8993
    ## 6493     729440
    ## 6494       4858
    ## 6495     339345
    ## 6496     339344
    ## 6497      26145
    ## 6498       8189
    ## 6499      81492
    ## 6500       1762
    ## 6501       1760
    ## 6502     147912
    ## 6503      23403
    ## 6504      54814
    ## 6505       6633
    ## 6506       2696
    ## 6507      24139
    ## 6508       2828
    ## 6509      80207
    ## 6510       7408
    ## 6511       6253
    ## 6512       2354
    ## 6513       2067
    ## 6514      10849
    ## 6515      10848
    ## 6516       2068
    ## 6517     147700
    ## 6518       1158
    ## 6519       <NA>
    ## 6520      57787
    ## 6521      90332
    ## 6522     388552
    ## 6523      79090
    ## 6524     284353
    ## 6525     284352
    ## 6526      79760
    ## 6527       <NA>
    ## 6528      11129
    ## 6529       5971
    ## 6530       1209
    ## 6531        341
    ## 6532        348
    ## 6533      10452
    ## 6534       <NA>
    ## 6535       5819
    ## 6536       4059
    ## 6537        602
    ## 6538       <NA>
    ## 6539       5817
    ## 6540       <NA>
    ## 6541       <NA>
    ## 6542       <NA>
    ## 6543       <NA>
    ## 6544       <NA>
    ## 6545       <NA>
    ## 6546       <NA>
    ## 6547       <NA>
    ## 6548       <NA>
    ## 6549       <NA>
    ## 6550       3783
    ## 6551      56006
    ## 6552       <NA>
    ## 6553       5329
    ## 6554     199731
    ## 6555       <NA>
    ## 6556     126298
    ## 6557       7515
    ## 6558       <NA>
    ## 6559      23474
    ## 6560       <NA>
    ## 6561       6223
    ## 6562        973
    ## 6563       9138
    ## 6564      10567
    ## 6565        478
    ## 6566       2901
    ## 6567       <NA>
    ## 6568       5452
    ## 6569     162989
    ## 6570       <NA>
    ## 6571       2931
    ## 6572       <NA>
    ## 6573       2077
    ## 6574      23152
    ## 6575       5050
    ## 6576     284339
    ## 6577       1954
    ## 6578       <NA>
    ## 6579       3991
    ## 6580        634
    ## 6581      55101
    ## 6582     374907
    ## 6583        593
    ## 6584      56915
    ## 6585     641649
    ## 6586      80776
    ## 6587       7040
    ## 6588      90324
    ## 6589      11100
    ## 6590        558
    ## 6591      29785
    ## 6592       <NA>
    ## 6593     112398
    ## 6594      53916
    ## 6595       8190
    ## 6596       6626
    ## 6597       <NA>
    ## 6598      80271
    ## 6599      79934
    ## 6600       <NA>
    ## 6601       9253
    ## 6602       8425
    ## 6603      92799
    ## 6604      57731
    ## 6605        645
    ## 6606      29946
    ## 6607      29950
    ## 6608     147746
    ## 6609      23646
    ## 6610       <NA>
    ## 6611        208
    ## 6612     148014
    ## 6613       4294
    ## 6614       <NA>
    ## 6615       <NA>
    ## 6616       <NA>
    ## 6617       <NA>
    ## 6618       <NA>
    ## 6619       <NA>
    ## 6620       <NA>
    ## 6621       <NA>
    ## 6622       5704
    ## 6623       2091
    ## 6624       9149
    ## 6625     163126
    ## 6626     126272
    ## 6627     348303
    ## 6628      10683
    ## 6629      92609
    ## 6630       <NA>
    ## 6631       6217
    ## 6632      64857
    ## 6633       7538
    ## 6634      55588
    ## 6635      54623
    ## 6636      55095
    ## 6637       9535
    ## 6638      57622
    ## 6639      10298
    ## 6640       <NA>
    ## 6641     126433
    ## 6642     115290
    ## 6643       6183
    ## 6644      54938
    ## 6645     643669
    ## 6646       4793
    ## 6647      22933
    ## 6648     126432
    ## 6649       3191
    ## 6650       <NA>
    ## 6651       1891
    ## 6652       3960
    ## 6653       3963
    ## 6654     147968
    ## 6655       <NA>
    ## 6656         81
    ## 6657      27335
    ## 6658      11184
    ## 6659       6261
    ## 6660       <NA>
    ## 6661     115727
    ## 6662     147965
    ## 6663       <NA>
    ## 6664     399473
    ## 6665     199720
    ## 6666       5714
    ## 6667       <NA>
    ## 6668       9424
    ## 6669      90522
    ## 6670       <NA>
    ## 6671      10653
    ## 6672      94274
    ## 6673       8193
    ## 6674      23094
    ## 6675       <NA>
    ## 6676      22835
    ## 6677       <NA>
    ## 6678       <NA>
    ## 6679       <NA>
    ## 6680       <NA>
    ## 6681       <NA>
    ## 6682       <NA>
    ## 6683       <NA>
    ## 6684      57677
    ## 6685       <NA>
    ## 6686       <NA>
    ## 6687       <NA>
    ## 6688       <NA>
    ## 6689       <NA>
    ## 6690       1346
    ## 6691        826
    ## 6692       1155
    ## 6693       5438
    ## 6694     284403
    ## 6695      25999
    ## 6696      84964
    ## 6697     163183
    ## 6698     644096
    ## 6699      79414
    ## 6700       7305
    ## 6701      10870
    ## 6702      84807
    ## 6703        333
    ## 6704      84063
    ## 6705       4868
    ## 6706     115703
    ## 6707     148137
    ## 6708     126393
    ## 6709      55957
    ## 6710      55851
    ## 6711     199746
    ## 6712      79713
    ## 6713       9757
    ## 6714      11045
    ## 6715       1340
    ## 6716      79171
    ## 6717       <NA>
    ## 6718      23354
    ## 6719        495
    ## 6720      10430
    ## 6721       <NA>
    ## 6722       <NA>
    ## 6723     374897
    ## 6724      93099
    ## 6725       4099
    ## 6726       7392
    ## 6727       <NA>
    ## 6728      51599
    ## 6729      53827
    ## 6730      53822
    ## 6731       5348
    ## 6732     163175
    ## 6733       5349
    ## 6734       3249
    ## 6735       6324
    ## 6736      57655
    ## 6737     126374
    ## 6738      10054
    ## 6739      84306
    ## 6740       <NA>
    ## 6741       <NA>
    ## 6742       <NA>
    ## 6743      26065
    ## 6744       <NA>
    ## 6745      79047
    ## 6746      64377
    ## 6747       5184
    ## 6748       1054
    ## 6749       1050
    ## 6750      56301
    ## 6751       4037
    ## 6752      55094
    ## 6753      85415
    ## 6754      91442
    ## 6755      84902
    ## 6756      91646
    ## 6757     390916
    ## 6758       <NA>
    ## 6759     388531
    ## 6760      84079
    ## 6761       9141
    ## 6762     147991
    ## 6763       <NA>
    ## 6764       <NA>
    ## 6765      57616
    ## 6766       <NA>
    ## 6767       <NA>
    ## 6768       8725
    ## 6769        898
    ## 6770       <NA>
    ## 6771      79156
    ## 6772      10775
    ## 6773       <NA>
    ## 6774       <NA>
    ## 6775       <NA>
    ## 6776     342865
    ## 6777       <NA>
    ## 6778       <NA>
    ## 6779       <NA>
    ## 6780       <NA>
    ## 6781       <NA>
    ## 6782       <NA>
    ## 6783       <NA>
    ## 6784       <NA>
    ## 6785       <NA>
    ## 6786       <NA>
    ## 6787       <NA>
    ## 6788       2109
    ## 6789     147645
    ## 6790     402665
    ## 6791        945
    ## 6792       <NA>
    ## 6793       <NA>
    ## 6794       <NA>
    ## 6795      90353
    ## 6796      11202
    ## 6797       5653
    ## 6798       <NA>
    ## 6799       <NA>
    ## 6800       6320
    ## 6801      50944
    ## 6802       <NA>
    ## 6803      84258
    ## 6804      94030
    ## 6805     554235
    ## 6806     126119
    ## 6807     284361
    ## 6808     112703
    ## 6809       5424
    ## 6810       7376
    ## 6811       9476
    ## 6812       <NA>
    ## 6813       3748
    ## 6814      79784
    ## 6815       <NA>
    ## 6816      51231
    ## 6817      22809
    ## 6818      23636
    ## 6819     259307
    ## 6820      79735
    ## 6821      84335
    ## 6822      11284
    ## 6823      53635
    ## 6824      81857
    ## 6825      80199
    ## 6826        160
    ## 6827      60385
    ## 6828     126129
    ## 6829       <NA>
    ## 6830       3276
    ## 6831       <NA>
    ## 6832      83596
    ## 6833       3661
    ## 6834      58506
    ## 6835       6237
    ## 6836      57479
    ## 6837       5639
    ## 6838      51070
    ## 6839      57333
    ## 6840       2217
    ## 6841       6205
    ## 6842      23521
    ## 6843       <NA>
    ## 6844     126133
    ## 6845      55011
    ## 6846      57030
    ## 6847  100507003
    ## 6848     147872
    ## 6849      27120
    ## 6850       8463
    ## 6851        951
    ## 6852      28968
    ## 6853      54795
    ## 6854       8541
    ## 6855       <NA>
    ## 6856      64130
    ## 6857       6625
    ## 6858      10856
    ## 6859       2997
    ## 6860       <NA>
    ## 6861        581
    ## 6862      27294
    ## 6863       4924
    ## 6864      23645
    ## 6865      57664
    ## 6866       <NA>
    ## 6867        587
    ## 6868       <NA>
    ## 6869      54922
    ## 6870       <NA>
    ## 6871     284358
    ## 6872       <NA>
    ## 6873     126147
    ## 6874       <NA>
    ## 6875       1628
    ## 6876      56848
    ## 6877       6141
    ## 6878       6820
    ## 6879       <NA>
    ## 6880     114783
    ## 6881       9266
    ## 6882       3770
    ## 6883      83743
    ## 6884       <NA>
    ## 6885       2906
    ## 6886      10945
    ## 6887      55260
    ## 6888       2014
    ## 6889      93233
    ## 6890        368
    ## 6891      23420
    ## 6892       3767
    ## 6893       6833
    ## 6894       3746
    ## 6895      26297
    ## 6896       7166
    ## 6897     113174
    ## 6898       6288
    ## 6899      11234
    ## 6900       2965
    ## 6901       3939
    ## 6902       7251
    ## 6903      55293
    ## 6904     144108
    ## 6905     144110
    ## 6906      84867
    ## 6907      54503
    ## 6908      79733
    ## 6909      89797
    ## 6910      10553
    ## 6911      10196
    ## 6912       4745
    ## 6913      57084
    ## 6914       2188
    ## 6915       2620
    ## 6916     258010
    ## 6917     338645
    ## 6918       <NA>
    ## 6919     114791
    ## 6920      23191
    ## 6921      81614
    ## 6922     123606
    ## 6923       8924
    ## 6924       4948
    ## 6925       2567
    ## 6926       <NA>
    ## 6927       2558
    ## 6928       2562
    ## 6929      57194
    ## 6930       7337
    ## 6931       <NA>
    ## 6932       6638
    ## 6933       <NA>
    ## 6934       <NA>
    ## 6935       <NA>
    ## 6936       <NA>
    ## 6937       4692
    ## 6938      54551
    ## 6939       7681
    ## 6940       1139
    ## 6941     161725
    ## 6942       <NA>
    ## 6943      51621
    ## 6944      54893
    ## 6945      10199
    ## 6946      84693
    ## 6947        321
    ## 6948      23359
    ## 6949      56160
    ## 6950       7082
    ## 6951       <NA>
    ## 6952      80213
    ## 6953       5046
    ## 6954       6627
    ## 6955      55829
    ## 6956      22856
    ## 6957      79705
    ## 6958        220
    ## 6959     140460
    ## 6960      55180
    ## 6961     204219
    ## 6962     170691
    ## 6963       <NA>
    ## 6964     145748
    ## 6965       4205
    ## 6966     123355
    ## 6967      64927
    ## 6968      23336
    ## 6969       <NA>
    ## 6970       3480
    ## 6971     283777
    ## 6972      91947
    ## 6973       7026
    ## 6974      55784
    ## 6975       <NA>
    ## 6976       <NA>
    ## 6977      56963
    ## 6978       1106
    ## 6979     400451
    ## 6980       8128
    ## 6981      28232
    ## 6982       9899
    ## 6983      11214
    ## 6984       <NA>
    ## 6985      64410
    ## 6986       4916
    ## 6987       <NA>
    ## 6988      26589
    ## 6989      64963
    ## 6990      55070
    ## 6991       <NA>
    ## 6992      64782
    ## 6993       3669
    ## 6994        176
    ## 6995     145864
    ## 6996       4240
    ## 6997      11057
    ## 6998       6017
    ## 6999      55215
    ## 7000       5428
    ## 7001     254559
    ## 7002       <NA>
    ## 7003       <NA>
    ## 7004      51458
    ## 7005     374654
    ## 7006       <NA>
    ## 7007       8800
    ## 7008      56964
    ## 7009     145873
    ## 7010        290
    ## 7011      10239
    ## 7012     348110
    ## 7013       <NA>
    ## 7014       3418
    ## 7015      10509
    ## 7016      10519
    ## 7017     390637
    ## 7018      51335
    ## 7019      26276
    ## 7020      91433
    ## 7021       9055
    ## 7022      55898
    ## 7023     374659
    ## 7024       4122
    ## 7025       2242
    ## 7026       5045
    ## 7027        641
    ## 7028      64784
    ## 7029       8826
    ## 7030      54993
    ## 7031      84942
    ## 7032       4828
    ## 7033      23478
    ## 7034       <NA>
    ## 7035      57538
    ## 7036       9154
    ## 7037       5151
    ## 7038       6218
    ## 7039      64506
    ## 7040       8120
    ## 7041       <NA>
    ## 7042     123720
    ## 7043       9455
    ## 7044       <NA>
    ## 7045       <NA>
    ## 7046      53339
    ## 7047      53346
    ## 7048      50810
    ## 7049       6457
    ## 7050     283726
    ## 7051      79631
    ## 7052      84206
    ## 7053       <NA>
    ## 7054       <NA>
    ## 7055      80765
    ## 7056       3603
    ## 7057     161502
    ## 7058       <NA>
    ## 7059      59274
    ## 7060      23184
    ## 7061      57214
    ## 7062      58489
    ## 7063       9915
    ## 7064       <NA>
    ## 7065       2184
    ## 7066      54469
    ## 7067       2346
    ## 7068       2915
    ## 7069       1075
    ## 7070      65084
    ## 7071       8322
    ## 7072      11098
    ## 7073       <NA>
    ## 7074      10873
    ## 7075      60494
    ## 7076      51501
    ## 7077       8726
    ## 7078       8301
    ## 7079      54843
    ## 7080     220388
    ## 7081      58487
    ## 7082      84233
    ## 7083      55863
    ## 7084       1740
    ## 7085      60492
    ## 7086     338699
    ## 7087      51585
    ## 7088      27314
    ## 7089     220042
    ## 7090       5547
    ## 7091     220382
    ## 7092      26011
    ## 7093      79731
    ## 7094       9846
    ## 7095      57558
    ## 7096     283219
    ## 7097      79053
    ## 7098       4718
    ## 7099       7069
    ## 7100      65987
    ## 7101      92105
    ## 7102      28971
    ## 7103      51773
    ## 7104       1207
    ## 7105     282679
    ## 7106       5058
    ## 7107       4647
    ## 7108        726
    ## 7109       4975
    ## 7110     192134
    ## 7111      55331
    ## 7112      25987
    ## 7113       2615
    ## 7114      56946
    ## 7115       5612
    ## 7116       7481
    ## 7117       7405
    ## 7118       <NA>
    ## 7119       <NA>
    ## 7120      84649
    ## 7121       4135
    ## 7122        871
    ## 7123      81544
    ## 7124     283212
    ## 7125       6188
    ## 7126        408
    ## 7127  100507050
    ## 7128       <NA>
    ## 7129      11309
    ## 7130      10825
    ## 7131       9789
    ## 7132     143570
    ## 7133     254225
    ## 7134      10714
    ## 7135     387787
    ## 7136      10008
    ## 7137     283209
    ## 7138     283208
    ## 7139      51400
    ## 7140      26005
    ## 7141       7352
    ## 7142       7351
    ## 7143     374407
    ## 7144      51287
    ## 7145      51642
    ## 7146       5870
    ## 7147      58473
    ## 7148      23201
    ## 7149      84957
    ## 7150       9828
    ## 7151       5031
    ## 7152       9873
    ## 7153      89849
    ## 7154      10809
    ## 7155     116985
    ## 7156       5138
    ## 7157      81570
    ## 7158       3636
    ## 7159       <NA>
    ## 7160       2348
    ## 7161      25906
    ## 7162      55004
    ## 7163       <NA>
    ## 7164       4926
    ## 7165      10068
    ## 7166      55298
    ## 7167       <NA>
    ## 7168       <NA>
    ## 7169      57053
    ## 7170       4928
    ## 7171      27315
    ## 7172        391
    ## 7173       6786
    ## 7174       6240
    ## 7175       6737
    ## 7176      55128
    ## 7177       <NA>
    ## 7178       <NA>
    ## 7179       <NA>
    ## 7180       <NA>
    ## 7181       <NA>
    ## 7182       <NA>
    ## 7183      85363
    ## 7184       <NA>
    ## 7185       <NA>
    ## 7186       <NA>
    ## 7187       <NA>
    ## 7188      84067
    ## 7189       1262
    ## 7190        887
    ## 7191     112464
    ## 7192       6609
    ## 7193        322
    ## 7194       3263
    ## 7195      10612
    ## 7196      23647
    ## 7197      26515
    ## 7198       <NA>
    ## 7199     144132
    ## 7200      23378
    ## 7201       3611
    ## 7202       6881
    ## 7203       1200
    ## 7204       8642
    ## 7205      63875
    ## 7206     143425
    ## 7207     283298
    ## 7208       8495
    ## 7209       <NA>
    ## 7210       8665
    ## 7211       7275
    ## 7212      79608
    ## 7213       4004
    ## 7214      65975
    ## 7215       <NA>
    ## 7216       9866
    ## 7217       6157
    ## 7218       <NA>
    ## 7219      56672
    ## 7220      56674
    ## 7221      56675
    ## 7222      23258
    ## 7223     440026
    ## 7224       <NA>
    ## 7225      10527
    ## 7226       <NA>
    ## 7227       7465
    ## 7228      23075
    ## 7229      81846
    ## 7230        133
    ## 7231        272
    ## 7232      50862
    ## 7233      10894
    ## 7234      10335
    ## 7235       9646
    ## 7236       1982
    ## 7237     374378
    ## 7238      55031
    ## 7239      27122
    ## 7240       9645
    ## 7241      84953
    ## 7242      55742
    ## 7243       7003
    ## 7244     644943
    ## 7245        406
    ## 7246      84280
    ## 7247      84188
    ## 7248      10418
    ## 7249      22800
    ## 7250       1315
    ## 7251       5682
    ## 7252       5140
    ## 7253        796
    ## 7254     387755
    ## 7255      55553
    ## 7256       <NA>
    ## 7257       <NA>
    ## 7258     144100
    ## 7259       6207
    ## 7260       <NA>
    ## 7261       5286
    ## 7262       4925
    ## 7263       <NA>
    ## 7264      64131
    ## 7265       6210
    ## 7266      23204
    ## 7267      23049
    ## 7268      51760
    ## 7269     162073
    ## 7270      10229
    ## 7271      79905
    ## 7272       <NA>
    ## 7273      51573
    ## 7274       9738
    ## 7275       <NA>
    ## 7276     400506
    ## 7277     124152
    ## 7278       <NA>
    ## 7279      51704
    ## 7280     124274
    ## 7281       7369
    ## 7282      55623
    ## 7283     112479
    ## 7284      81691
    ## 7285     123879
    ## 7286      57149
    ## 7287      57146
    ## 7288       1428
    ## 7289       7385
    ## 7290     255762
    ## 7291     730094
    ## 7292     146177
    ## 7293      29904
    ## 7294      55718
    ## 7295       1039
    ## 7296  105371130
    ## 7297      51108
    ## 7298      10261
    ## 7299       9956
    ## 7300      57478
    ## 7301      91949
    ## 7302      23062
    ## 7303     124454
    ## 7304      56061
    ## 7305       4706
    ## 7306      79728
    ## 7307      84516
    ## 7308       5347
    ## 7309       5579
    ## 7310      10368
    ## 7311       5930
    ## 7312       <NA>
    ## 7313      27327
    ## 7314     115584
    ## 7315      55114
    ## 7316      51451
    ## 7317     342357
    ## 7318       <NA>
    ## 7319       9951
    ## 7320      79831
    ## 7321     197370
    ## 7322       <NA>
    ## 7323      50615
    ## 7324       2975
    ## 7325       <NA>
    ## 7326     146395
    ## 7327      23214
    ## 7328     388228
    ## 7329      27040
    ## 7330      83985
    ## 7331      84901
    ## 7332      79874
    ## 7333        487
    ## 7334      25970
    ## 7335       7284
    ## 7336      11273
    ## 7337       8663
    ## 7338       1201
    ## 7339      55911
    ## 7340      26471
    ## 7341       <NA>
    ## 7342     112869
    ## 7343       6817
    ## 7344      79008
    ## 7345     552900
    ## 7346      11151
    ## 7347       5595
    ## 7348      79153
    ## 7349      83719
    ## 7350       6911
    ## 7351       5531
    ## 7352        226
    ## 7353       <NA>
    ## 7354       <NA>
    ## 7355       8448
    ## 7356     283899
    ## 7357       8479
    ## 7358       9344
    ## 7359     124446
    ## 7360     253980
    ## 7361     253982
    ## 7362      26470
    ## 7363       <NA>
    ## 7364      10423
    ## 7365       9961
    ## 7366       <NA>
    ## 7367     112476
    ## 7368       4150
    ## 7369       3835
    ## 7370       <NA>
    ## 7371      23475
    ## 7372      10421
    ## 7373      26000
    ## 7374      29895
    ## 7375       <NA>
    ## 7376       <NA>
    ## 7377       <NA>
    ## 7378      79077
    ## 7379      22928
    ## 7380       3683
    ## 7381       <NA>
    ## 7382       <NA>
    ## 7383       <NA>
    ## 7384       <NA>
    ## 7385       <NA>
    ## 7386       <NA>
    ## 7387       <NA>
    ## 7388       <NA>
    ## 7389      78994
    ## 7390      64319
    ## 7391       <NA>
    ## 7392      10847
    ## 7393       <NA>
    ## 7394       5261
    ## 7395      90835
    ## 7396       9810
    ## 7397       <NA>
    ## 7398       9274
    ## 7399       1489
    ## 7400      54620
    ## 7401      93129
    ## 7402       9739
    ## 7403     112755
    ## 7404       <NA>
    ## 7405       <NA>
    ## 7406       <NA>
    ## 7407     339105
    ## 7408      79001
    ## 7409      10295
    ## 7410      84148
    ## 7411     146547
    ## 7412       2521
    ## 7413      29108
    ## 7414       3684
    ## 7415       3687
    ## 7416       3681
    ## 7417       1339
    ## 7418       <NA>
    ## 7419      79798
    ## 7420       7041
    ## 7421       <NA>
    ## 7422       6001
    ## 7423       7073
    ## 7424       9531
    ## 7425      22876
    ## 7426      79892
    ## 7427      11196
    ## 7428     196051
    ## 7429       <NA>
    ## 7430      55717
    ## 7431       2263
    ## 7432      11101
    ## 7433       <NA>
    ## 7434      54780
    ## 7435      10579
    ## 7436     118663
    ## 7437       <NA>
    ## 7438      59338
    ## 7439       5654
    ## 7440      50624
    ## 7441     196792
    ## 7442       <NA>
    ## 7443     118672
    ## 7444       <NA>
    ## 7445      64376
    ## 7446         36
    ## 7447       9184
    ## 7448       2849
    ## 7449     119587
    ## 7450      51363
    ## 7451       4942
    ## 7452      64077
    ## 7453       9679
    ## 7454     399818
    ## 7455      23172
    ## 7456      54764
    ## 7457       1488
    ## 7458      26098
    ## 7459       7390
    ## 7460      56647
    ## 7461      55760
    ## 7462      92565
    ## 7463       8038
    ## 7464       <NA>
    ## 7465       1793
    ## 7466       <NA>
    ## 7467       5791
    ## 7468       4288
    ## 7469       <NA>
    ## 7470       4255
    ## 7471     253738
    ## 7472       <NA>
    ## 7473      10539
    ## 7474     256536
    ## 7475       <NA>
    ## 7476      55844
    ## 7477        664
    ## 7478     282973
    ## 7479      10570
    ## 7480     282974
    ## 7481      80313
    ## 7482     170394
    ## 7483       3632
    ## 7484      84504
    ## 7485      54777
    ## 7486      84435
    ## 7487      85442
    ## 7488       <NA>
    ## 7489        101
    ## 7490      10844
    ## 7491       <NA>
    ## 7492      50632
    ## 7493     282969
    ## 7494       1892
    ## 7495     196743
    ## 7496      92170
    ## 7497     503542
    ## 7498      93426
    ## 7499       <NA>
    ## 7500       <NA>
    ## 7501      51272
    ## 7502      60626
    ## 7503      23410
    ## 7504       5719
    ## 7505       <NA>
    ## 7506     171389
    ## 7507      80162
    ## 7508      10581
    ## 7509       8519
    ## 7510      10410
    ## 7511       <NA>
    ## 7512     338707
    ## 7513      11187
    ## 7514      59307
    ## 7515      81490
    ## 7516       6050
    ## 7517       3265
    ## 7518     115399
    ## 7519       8045
    ## 7520      57661
    ## 7521       3665
    ## 7522       1815
    ## 7523      10522
    ## 7524     283232
    ## 7525      64787
    ## 7526       6888
    ## 7527     347862
    ## 7528      51286
    ## 7529      79751
    ## 7530      55367
    ## 7531       6181
    ## 7532      57104
    ## 7533     283229
    ## 7534        977
    ## 7535       5441
    ## 7536       7106
    ## 7537      66005
    ## 7538        161
    ## 7539       4583
    ## 7540      54472
    ## 7541       9024
    ## 7542      81532
    ## 7543       1850
    ## 7544     402778
    ## 7545       1509
    ## 7546       7136
    ## 7547       4046
    ## 7548       6150
    ## 7549     283120
    ## 7550       3481
    ## 7551       7054
    ## 7552        430
    ## 7553        975
    ## 7554      10078
    ## 7555       3784
    ## 7556      10984
    ## 7557       <NA>
    ## 7558       1028
    ## 7559       5002
    ## 7560       4676
    ## 7561       <NA>
    ## 7562       <NA>
    ## 7563       <NA>
    ## 7564     114879
    ## 7565     116534
    ## 7566      55191
    ## 7567       1717
    ## 7568      22941
    ## 7569       2017
    ## 7570       8500
    ## 7571       8772
    ## 7572      55107
    ## 7573       2248
    ## 7574       <NA>
    ## 7575       <NA>
    ## 7576        595
    ## 7577     219931
    ## 7578       <NA>
    ## 7579      81706
    ## 7580      57480
    ## 7581      25902
    ## 7582       9590
    ## 7583      57621
    ## 7584       <NA>
    ## 7585      55005
    ## 7586      79624
    ## 7587       2099
    ## 7588      23345
    ## 7589      80177
    ## 7590       7432
    ## 7591      26271
    ## 7592      54516
    ## 7593      26575
    ## 7594      26034
    ## 7595     154043
    ## 7596      84918
    ## 7597       5110
    ## 7598     348995
    ## 7599       9113
    ## 7600      11104
    ## 7601     116254
    ## 7602      85313
    ## 7603      23118
    ## 7604       <NA>
    ## 7605      10090
    ## 7606      23328
    ## 7607     389432
    ## 7608     134957
    ## 7609      79747
    ## 7610      10981
    ## 7611       2911
    ## 7612     257218
    ## 7613      84085
    ## 7614       7957
    ## 7615       7402
    ## 7616       <NA>
    ## 7617      83443
    ## 7618       5325
    ## 7619      84946
    ## 7620       9749
    ## 7621       2519
    ## 7622       8504
    ## 7623     134637
    ## 7624      51390
    ## 7625       3097
    ## 7626       <NA>
    ## 7627       <NA>
    ## 7628      57211
    ## 7629      51534
    ## 7630      10370
    ## 7631       <NA>
    ## 7632     167838
    ## 7633      51696
    ## 7634      58527
    ## 7635      85021
    ## 7636      25901
    ## 7637       <NA>
    ## 7638      57224
    ## 7639      23593
    ## 7640      57221
    ## 7641       <NA>
    ## 7642      64065
    ## 7643       7128
    ## 7644       3459
    ## 7645      53832
    ## 7646     340146
    ## 7647       5191
    ## 7648       4217
    ## 7649       9053
    ## 7650       <NA>
    ## 7651       9774
    ## 7652      27115
    ## 7653      54806
    ## 7654       4602
    ## 7655      10767
    ## 7656       6446
    ## 7657       <NA>
    ## 7658     154091
    ## 7659       9519
    ## 7660       2070
    ## 7661       6206
    ## 7662     116843
    ## 7663       8417
    ## 7664      26002
    ## 7665       <NA>
    ## 7666       5167
    ## 7667       9439
    ## 7668       <NA>
    ## 7669       9465
    ## 7670       2037
    ## 7671       <NA>
    ## 7672     114801
    ## 7673       <NA>
    ## 7674     154075
    ## 7675      84456
    ## 7676      93663
    ## 7677       3908
    ## 7678       5796
    ## 7679     387104
    ## 7680       <NA>
    ## 7681       <NA>
    ## 7682      55862
    ## 7683      81847
    ## 7684      84870
    ## 7685     387103
    ## 7686      60487
    ## 7687       <NA>
    ## 7688     135114
    ## 7689     135112
    ## 7690      23493
    ## 7691      51020
    ## 7692       7164
    ## 7693     154214
    ## 7694     154215
    ## 7695      10345
    ## 7696     134829
    ## 7697     345895
    ## 7698       <NA>
    ## 7699      51389
    ## 7700       <NA>
    ## 7701       <NA>
    ## 7702       <NA>
    ## 7703      29940
    ## 7704       <NA>
    ## 7705       7259
    ## 7706     221294
    ## 7707       <NA>
    ## 7708      23270
    ## 7709       2444
    ## 7710       <NA>
    ## 7711     222537
    ## 7712       3066
    ## 7713       4082
    ## 7714       3910
    ## 7715     619208
    ## 7716      51175
    ## 7717       2534
    ## 7718       <NA>
    ## 7719      10758
    ## 7720       <NA>
    ## 7721       5980
    ## 7722       <NA>
    ## 7723       <NA>
    ## 7724     117247
    ## 7725      84154
    ## 7726     112495
    ## 7727        262
    ## 7728      23097
    ## 7729       8528
    ## 7730     728464
    ## 7731      51362
    ## 7732       8936
    ## 7733       2830
    ## 7734       9896
    ## 7735       9841
    ## 7736      64780
    ## 7737       6610
    ## 7738     285755
    ## 7739       8763
    ## 7740       <NA>
    ## 7741     285753
    ## 7742      27244
    ## 7743      84071
    ## 7744       2309
    ## 7745     246269
    ## 7746       8724
    ## 7747       7101
    ## 7748      28962
    ## 7749      11231
    ## 7750     256380
    ## 7751       <NA>
    ## 7752      55084
    ## 7753       <NA>
    ## 7754      57107
    ## 7755      57673
    ## 7756       <NA>
    ## 7757       <NA>
    ## 7758       <NA>
    ## 7759      55278
    ## 7760      84816
    ## 7761        202
    ## 7762       9474
    ## 7763        639
    ## 7764       5550
    ## 7765      64208
    ## 7766      11149
    ## 7767     389421
    ## 7768      57531
    ## 7769       <NA>
    ## 7770       2898
    ## 7771      10973
    ## 7772       <NA>
    ## 7773     285761
    ## 7774      57120
    ## 7775     116150
    ## 7776     222553
    ## 7777     387119
    ## 7778       5350
    ## 7779     254394
    ## 7780      25842
    ## 7781      79632
    ## 7782       <NA>
    ## 7783     221322
    ## 7784       <NA>
    ## 7785       2697
    ## 7786       3298
    ## 7787      57515
    ## 7788       5570
    ## 7789       2173
    ## 7790      10924
    ## 7791       9648
    ## 7792       3987
    ## 7793       5903
    ## 7794     165055
    ## 7795     344558
    ## 7796       <NA>
    ## 7797      65124
    ## 7798       5033
    ## 7799       <NA>
    ## 7800      90550
    ## 7801      10367
    ## 7802      54788
    ## 7803      54541
    ## 7804     119504
    ## 7805      51008
    ## 7806       9806
    ## 7807       9469
    ## 7808       5660
    ## 7809      64072
    ## 7810      64115
    ## 7811       <NA>
    ## 7812      55315
    ## 7813     219699
    ## 7814       5092
    ## 7815       8879
    ## 7816     219793
    ## 7817     140766
    ## 7818      27143
    ## 7819       4838
    ## 7820       1979
    ## 7821       <NA>
    ## 7822      55222
    ## 7823       5464
    ## 7824      56681
    ## 7825     219743
    ## 7826      84883
    ## 7827       <NA>
    ## 7828       1305
    ## 7829     219738
    ## 7830      23555
    ## 7831       3098
    ## 7832       <NA>
    ## 7833      80201
    ## 7834       6832
    ## 7835       <NA>
    ## 7836       9559
    ## 7837       5552
    ## 7838       <NA>
    ## 7839       <NA>
    ## 7840       9188
    ## 7841      79009
    ## 7842     219736
    ## 7843      55749
    ## 7844      80312
    ## 7845       8034
    ## 7846       1763
    ## 7847      55680
    ## 7848       3189
    ## 7849       <NA>
    ## 7850       <NA>
    ## 7851     220202
    ## 7852      84665
    ## 7853      26091
    ## 7854      23411
    ## 7855      56521
    ## 7856      29119
    ## 7857       <NA>
    ## 7858     347731
    ## 7859       <NA>
    ## 7860     221035
    ## 7861     221037
    ## 7862      29982
    ## 7863       1959
    ## 7864      84890
    ## 7865       <NA>
    ## 7866     219790
    ## 7867      84159
    ## 7868     219621
    ## 7869       9886
    ## 7870       <NA>
    ## 7871        983
    ## 7872        288
    ## 7873       8030
    ## 7874  100507027
    ## 7875     220963
    ## 7876     220965
    ## 7877      84457
    ## 7878      80114
    ## 7879       7019
    ## 7880       7321
    ## 7881      55847
    ## 7882     253430
    ## 7883      11130
    ## 7884      65217
    ## 7885       2781
    ## 7886       9609
    ## 7887        613
    ## 7888      23384
    ## 7889        135
    ## 7890      83606
    ## 7891      51733
    ## 7892       6634
    ## 7893     388886
    ## 7894       2678
    ## 7895       2687
    ## 7896      56241
    ## 7897      23523
    ## 7898       1652
    ## 7899       <NA>
    ## 7900       2952
    ## 7901       2953
    ## 7902       4282
    ## 7903      91319
    ## 7904       6598
    ## 7905       4320
    ## 7906     400916
    ## 7907       <NA>
    ## 7908       3275
    ## 7909       6285
    ## 7910      23181
    ## 7911       5116
    ## 7912      54059
    ## 7913       8888
    ## 7914       4047
    ## 7915       1292
    ## 7916       1291
    ## 7917      54039
    ## 7918       6573
    ## 7919      80781
    ## 7920      23275
    ## 7921        104
    ## 7922      85395
    ## 7923       3689
    ## 7924        754
    ## 7925       6612
    ## 7926       7327
    ## 7927      81543
    ## 7928       7226
    ## 7929       <NA>
    ## 7930       5211
    ## 7931      29947
    ## 7932       <NA>
    ## 7933       <NA>
    ## 7934       5822
    ## 7935       7109
    ## 7936       <NA>
    ## 7937      56894
    ## 7938       8568
    ## 7939       1476
    ## 7940       <NA>
    ## 7941       8566
    ## 7942      10994
    ## 7943      85360
    ## 7944     126402
    ## 7945       6511
    ## 7946       <NA>
    ## 7947       8612
    ## 7948      54531
    ## 7949      51298
    ## 7950     126567
    ## 7951      25759
    ## 7952     284451
    ## 7953      91978
    ## 7954        997
    ## 7955       3004
    ## 7956        682
    ## 7957        610
    ## 7958       5442
    ## 7959      27006
    ## 7960      55658
    ## 7961      10272
    ## 7962     400668
    ## 7963       5064
    ## 7964     126353
    ## 7965       5725
    ## 7966      79948
    ## 7967       5657
    ## 7968      10025
    ## 7969      91300
    ## 7970      84634
    ## 7971       1820
    ## 7972      57418
    ## 7973     116444
    ## 7974      91304
    ## 7975       1265
    ## 7976      10347
    ## 7977      23526
    ## 7978       5434
    ## 7979       2879
    ## 7980      22904
    ## 7981       6794
    ## 7982     255057
    ## 7983       <NA>
    ## 7984      90007
    ## 7985       1153
    ## 7986       <NA>
    ## 7987       1943
    ## 7988       <NA>
    ## 7989     374291
    ## 7990       2593
    ## 7991      26528
    ## 7992       6209
    ## 7993      10297
    ## 7994       <NA>
    ## 7995      54760
    ## 7996      92840
    ## 7997     339366
    ## 7998     126520
    ## 7999     399664
    ## 8000      53615
    ## 8001      10975
    ## 8002       6929
    ## 8003      57455
    ## 8004      83855
    ## 8005       <NA>
    ## 8006      81926
    ## 8007       <NA>
    ## 8008       1455
    ## 8009      55643
    ## 8010       2872
    ## 8011     126308
    ## 8012     113177
    ## 8013       8943
    ## 8014      84444
    ## 8015      55111
    ## 8016       8175
    ## 8017     126306
    ## 8018       4946
    ## 8019     645191
    ## 8020      51690
    ## 8021      56928
    ## 8022     360200
    ## 8023      26517
    ## 8024      84823
    ## 8025       4616
    ## 8026       2788
    ## 8027     148252
    ## 8028       <NA>
    ## 8029      29985
    ## 8030       6449
    ## 8031       7064
    ## 8032      84699
    ## 8033       5605
    ## 8034      51341
    ## 8035      51588
    ## 8036       1938
    ## 8037       1613
    ## 8038      85300
    ## 8039      23217
    ## 8040       4145
    ## 8041     116541
    ## 8042       9546
    ## 8043      27134
    ## 8044      23396
    ## 8045       <NA>
    ## 8046      58509
    ## 8047       6915
    ## 8048      10362
    ## 8049     126321
    ## 8050       <NA>
    ## 8051      51343
    ## 8052      83475
    ## 8053       <NA>
    ## 8054     284422
    ## 8055       4782
    ## 8056      60680
    ## 8057       <NA>
    ## 8058      56926
    ## 8059       8698
    ## 8060       2769
    ## 8061       2767
    ## 8062       <NA>
    ## 8063       7089
    ## 8064      79816
    ## 8065       <NA>
    ## 8066      51548
    ## 8067     170961
    ## 8068       <NA>
    ## 8069       <NA>
    ## 8070       <NA>
    ## 8071       <NA>
    ## 8072       <NA>
    ## 8073       <NA>
    ## 8074       <NA>
    ## 8075       <NA>
    ## 8076       <NA>
    ## 8077       <NA>
    ## 8078       6996
    ## 8079      83468
    ## 8080      29915
    ## 8081       4801
    ## 8082       7296
    ## 8083      50515
    ## 8084      84102
    ## 8085       <NA>
    ## 8086     160428
    ## 8087      23325
    ## 8088      55198
    ## 8089       <NA>
    ## 8090       <NA>
    ## 8091       9891
    ## 8092      10970
    ## 8093     255394
    ## 8094      55703
    ## 8095       5992
    ## 8096      55188
    ## 8097       <NA>
    ## 8098      90488
    ## 8099      80298
    ## 8100       1407
    ## 8101     121551
    ## 8102      11137
    ## 8103      11108
    ## 8104      51493
    ## 8105       <NA>
    ## 8106      25793
    ## 8107       8224
    ## 8108       7078
    ## 8109       7184
    ## 8110       <NA>
    ## 8111      51559
    ## 8112        429
    ## 8113       3479
    ## 8114       5367
    ## 8115      79023
    ## 8116      51019
    ## 8117      55332
    ## 8118      79158
    ## 8119      56994
    ## 8120      50511
    ## 8121       4604
    ## 8122        400
    ## 8123      27340
    ## 8124     121601
    ## 8125     283431
    ## 8126     246213
    ## 8127      55681
    ## 8128      64431
    ## 8129      23074
    ## 8130      56899
    ## 8131        317
    ## 8132     121457
    ## 8133       5250
    ## 8134       7112
    ## 8135     196475
    ## 8136     121441
    ## 8137     144535
    ## 8138       5128
    ## 8139       2004
    ## 8140       <NA>
    ## 8141       4048
    ## 8142     144193
    ## 8143       6636
    ## 8144      59277
    ## 8145      10988
    ## 8146      55591
    ## 8147      55785
    ## 8148       7181
    ## 8149      55967
    ## 8150      57458
    ## 8151      51134
    ## 8152      10154
    ## 8153       8738
    ## 8154       8835
    ## 8155      28977
    ## 8156       7334
    ## 8157      11163
    ## 8158       <NA>
    ## 8159       8411
    ## 8160        694
    ## 8161       1634
    ## 8162       4060
    ## 8163       <NA>
    ## 8164        490
    ## 8165     282809
    ## 8166       8693
    ## 8167       <NA>
    ## 8168       1848
    ## 8169       <NA>
    ## 8170       <NA>
    ## 8171     160418
    ## 8172      80184
    ## 8173       <NA>
    ## 8174      25834
    ## 8175       <NA>
    ## 8176       <NA>
    ## 8177       4922
    ## 8178       9182
    ## 8179      55117
    ## 8180     160335
    ## 8181      84190
    ## 8182      29080
    ## 8183       <NA>
    ## 8184       8499
    ## 8185       <NA>
    ## 8186       8825
    ## 8187     283310
    ## 8188       4659
    ## 8189       <NA>
    ## 8190       6857
    ## 8191      89795
    ## 8192       1466
    ## 8193      23390
    ## 8194     114882
    ## 8195      79738
    ## 8196       4673
    ## 8197      22822
    ## 8198      11103
    ## 8199      11010
    ## 8200      84698
    ## 8201       3747
    ## 8202     552889
    ## 8203      29953
    ## 8204     121278
    ## 8205      64786
    ## 8206      23011
    ## 8207      55266
    ## 8208      83591
    ## 8209     196441
    ## 8210       8549
    ## 8211       <NA>
    ## 8212       7103
    ## 8213       5801
    ## 8214       5787
    ## 8215      27345
    ## 8216       <NA>
    ## 8217       4848
    ## 8218     117177
    ## 8219      10576
    ## 8220      10818
    ## 8221       8089
    ## 8222       <NA>
    ## 8223      11052
    ## 8224       1368
    ## 8225       4193
    ## 8226      55508
    ## 8227      57122
    ## 8228       5908
    ## 8229      56890
    ## 8230       8445
    ## 8231      55832
    ## 8232      23426
    ## 8233      92797
    ## 8234      11213
    ## 8235      51643
    ## 8236      84298
    ## 8237     253827
    ## 8238      23592
    ## 8239      11197
    ## 8240      23329
    ## 8241       2799
    ## 8242     283349
    ## 8243      29110
    ## 8244      11260
    ## 8245       <NA>
    ## 8246      57522
    ## 8247       <NA>
    ## 8248        552
    ## 8249      57460
    ## 8250      23041
    ## 8251       9958
    ## 8252       <NA>
    ## 8253       9194
    ## 8254     121227
    ## 8255      91419
    ## 8256      10106
    ## 8257      10102
    ## 8258       4234
    ## 8259       <NA>
    ## 8260       1019
    ## 8261       6302
    ## 8262       <NA>
    ## 8263     116986
    ## 8264      10956
    ## 8265       2583
    ## 8266      65012
    ## 8267     115557
    ## 8268     196403
    ## 8269       <NA>
    ## 8270      79837
    ## 8271       3798
    ## 8272      10540
    ## 8273     114785
    ## 8274       1649
    ## 8275       <NA>
    ## 8276      64333
    ## 8277       2735
    ## 8278      83729
    ## 8279      22864
    ## 8280     246329
    ## 8281      56901
    ## 8282       6472
    ## 8283      11247
    ## 8284       4035
    ## 8285       6778
    ## 8286       4665
    ## 8287      23306
    ## 8288       4640
    ## 8289       <NA>
    ## 8290       9880
    ## 8291      11318
    ## 8292       <NA>
    ## 8293     121214
    ## 8294       5557
    ## 8295       4666
    ## 8296      10728
    ## 8297       <NA>
    ## 8298       <NA>
    ## 8299      11176
    ## 8300       5939
    ## 8301      27165
    ## 8302     283377
    ## 8303       8914
    ## 8304       6773
    ## 8305      51561
    ## 8306       9924
    ## 8307      10330
    ## 8308       1431
    ## 8309      93058
    ## 8310     283373
    ## 8311     283375
    ## 8312      79035
    ## 8313      10193
    ## 8314       6601
    ## 8315       4637
    ## 8316     140465
    ## 8317       <NA>
    ## 8318      23344
    ## 8319      84872
    ## 8320       6171
    ## 8321       5036
    ## 8322       2065
    ## 8323       6231
    ## 8324      64375
    ## 8325       <NA>
    ## 8326       6821
    ## 8327       5869
    ## 8328       1017
    ## 8329       6490
    ## 8330       1606
    ## 8331      84305
    ## 8332       4327
    ## 8333     440104
    ## 8334      85406
    ## 8335      84324
    ## 8336      29095
    ## 8337      10220
    ## 8338        967
    ## 8339       5959
    ## 8340       2647
    ## 8341       3679
    ## 8342      58158
    ## 8343       <NA>
    ## 8344       <NA>
    ## 8345       3643
    ## 8346      23370
    ## 8347      92960
    ## 8348       <NA>
    ## 8349      57192
    ## 8350      10908
    ## 8351      57662
    ## 8352      56949
    ## 8353  100131801
    ## 8354       6813
    ## 8355      56729
    ## 8356     199675
    ## 8357       <NA>
    ## 8358       <NA>
    ## 8359     126003
    ## 8360     339390
    ## 8361       <NA>
    ## 8362     115704
    ## 8363      80164
    ## 8364       5609
    ## 8365  100507588
    ## 8366       6618
    ## 8367     404217
    ## 8368      10469
    ## 8369       1994
    ## 8370       6370
    ## 8371      79603
    ## 8372       <NA>
    ## 8373      79801
    ## 8374       1948
    ## 8375      55082
    ## 8376     728215
    ## 8377       3981
    ## 8378      84945
    ## 8379      10673
    ## 8380      23026
    ## 8381       <NA>
    ## 8382       8660
    ## 8383       <NA>
    ## 8384       1282
    ## 8385       1284
    ## 8386      55647
    ## 8387      55739
    ## 8388      79587
    ## 8389       3621
    ## 8390      55608
    ## 8391       8874
    ## 8392       6656
    ## 8393      10426
    ## 8394      23250
    ## 8395      23263
    ## 8396      55795
    ## 8397       8451
    ## 8398       3916
    ## 8399      79774
    ## 8400      55208
    ## 8401      55002
    ## 8402       7027
    ## 8403     348013
    ## 8404       2621
    ## 8405      22821
    ## 8406       <NA>
    ## 8407       8881
    ## 8408      65110
    ## 8409     283489
    ## 8410      55352
    ## 8411      26260
    ## 8412     157695
    ## 8413     157697
    ## 8414       9228
    ## 8415       2055
    ## 8416       9639
    ## 8417       9920
    ## 8418      64478
    ## 8419      79648
    ## 8420        285
    ## 8421      55326
    ## 8422       <NA>
    ## 8423       1672
    ## 8424        540
    ## 8425     440138
    ## 8426       4752
    ## 8427      26586
    ## 8428      51028
    ## 8429      55901
    ## 8430      10166
    ## 8431      10240
    ## 8432     114926
    ## 8433       6575
    ## 8434       7419
    ## 8435       5423
    ## 8436       3551
    ## 8437       5327
    ## 8438      10947
    ## 8439       <NA>
    ## 8440       7994
    ## 8441        286
    ## 8442       <NA>
    ## 8443     137964
    ## 8444      84296
    ## 8445      51125
    ## 8446       6422
    ## 8447      79698
    ## 8448      56892
    ## 8449     169355
    ## 8450       3620
    ## 8451       8749
    ## 8452     255926
    ## 8453     203102
    ## 8454       8754
    ## 8455      83877
    ## 8456     203100
    ## 8457       <NA>
    ## 8458      59339
    ## 8459       6867
    ## 8460       2260
    ## 8461     137994
    ## 8462       <NA>
    ## 8463      54904
    ## 8464      84513
    ## 8465      23259
    ## 8466       <NA>
    ## 8467       9530
    ## 8468      27257
    ## 8469       6770
    ## 8470       9070
    ## 8471     157855
    ## 8472     138050
    ## 8473      84197
    ## 8474       2339
    ## 8475      84376
    ## 8476      81790
    ## 8477       <NA>
    ## 8478      55145
    ## 8479       <NA>
    ## 8480      11160
    ## 8481       <NA>
    ## 8482      11212
    ## 8483      25960
    ## 8484      55290
    ## 8485      80223
    ## 8486     137362
    ## 8487        155
    ## 8488       <NA>
    ## 8489       1978
    ## 8490       1142
    ## 8491       8973
    ## 8492     137970
    ## 8493      78986
    ## 8494      79845
    ## 8495      80185
    ## 8496      84549
    ## 8497      84750
    ## 8498       3084
    ## 8499       7486
    ## 8500      29942
    ## 8501       <NA>
    ## 8502      56154
    ## 8503       5516
    ## 8504       7993
    ## 8505       2936
    ## 8506       2961
    ## 8507  100507341
    ## 8508      11030
    ## 8509      10671
    ## 8510     619373
    ## 8511       <NA>
    ## 8512      23484
    ## 8513      51669
    ## 8514       1846
    ## 8515       8658
    ## 8516      79660
    ## 8517      90459
    ## 8518       9258
    ## 8519       <NA>
    ## 8520     137075
    ## 8521     157285
    ## 8522      91694
    ## 8523       <NA>
    ## 8524      10395
    ## 8525       <NA>
    ## 8526     137868
    ## 8527       <NA>
    ## 8528       7991
    ## 8529     286097
    ## 8530      51201
    ## 8531      29883
    ## 8532     137492
    ## 8533       9108
    ## 8534       6542
    ## 8535       5157
    ## 8536      57509
    ## 8537       5108
    ## 8538        427
    ## 8539       2483
    ## 8540       <NA>
    ## 8541       2195
    ## 8542       <NA>
    ## 8543      25854
    ## 8544       7098
    ## 8545       8470
    ## 8546      27295
    ## 8547     256309
    ## 8548       <NA>
    ## 8549      55325
    ## 8550     353322
    ## 8551      83891
    ## 8552      57587
    ## 8553        291
    ## 8554     391723
    ## 8555       2180
    ## 8556      79682
    ## 8557     201973
    ## 8558        836
    ## 8559       3660
    ## 8560     133121
    ## 8561      56977
    ## 8562      60684
    ## 8563       <NA>
    ## 8564       3622
    ## 8565      55602
    ## 8566       <NA>
    ## 8567      80014
    ## 8568      53842
    ## 8569       1635
    ## 8570      55714
    ## 8571        175
    ## 8572      55247
    ## 8573       7424
    ## 8574      60559
    ## 8575     140458
    ## 8576       <NA>
    ## 8577     116966
    ## 8578       2823
    ## 8579       8001
    ## 8580       3248
    ## 8581      80817
    ## 8582      26269
    ## 8583      11341
    ## 8584       8819
    ## 8585       3148
    ## 8586      51809
    ## 8587     442117
    ## 8588       <NA>
    ## 8589       9848
    ## 8590      54969
    ## 8591       1182
    ## 8592       <NA>
    ## 8593       4750
    ## 8594      57630
    ## 8595      84869
    ## 8596      23022
    ## 8597      55601
    ## 8598      50859
    ## 8599       1363
    ## 8600       6307
    ## 8601      11275
    ## 8602     201931
    ## 8603       <NA>
    ## 8604      55319
    ## 8605       4889
    ## 8606       4886
    ## 8607      92345
    ## 8608          9
    ## 8609         10
    ## 8610      23362
    ## 8611      55790
    ## 8612      55174
    ## 8613       4023
    ## 8614        526
    ## 8615      11178
    ## 8616       <NA>
    ## 8617       <NA>
    ## 8618       <NA>
    ## 8619       <NA>
    ## 8620       <NA>
    ## 8621       <NA>
    ## 8622       <NA>
    ## 8623      57130
    ## 8624      51291
    ## 8625       9170
    ## 8626      80714
    ## 8627     148113
    ## 8628     374887
    ## 8629      51079
    ## 8630      83983
    ## 8631      54815
    ## 8632      23383
    ## 8633      57794
    ## 8634      53345
    ## 8635     404037
    ## 8636       1463
    ## 8637       8625
    ## 8638     126382
    ## 8639     729991
    ## 8640      54929
    ## 8641     284439
    ## 8642      93436
    ## 8643      10147
    ## 8644       9454
    ## 8645      54555
    ## 8646      11316
    ## 8647      10715
    ## 8648       5976
    ## 8649      23373
    ## 8650      55295
    ## 8651      25789
    ## 8652       9244
    ## 8653      55049
    ## 8654       7311
    ## 8655      79036
    ## 8656      23770
    ## 8657       8178
    ## 8658      51477
    ## 8659     170463
    ## 8660     126364
    ## 8661      54858
    ## 8662      25804
    ## 8663       3727
    ## 8664       5143
    ## 8665       5864
    ## 8666      84769
    ## 8667       5296
    ## 8668       <NA>
    ## 8669      23031
    ## 8670      27106
    ## 8671       3780
    ## 8672     115098
    ## 8673       6528
    ## 8674       6142
    ## 8675      55201
    ## 8676      93323
    ## 8677       4650
    ## 8678      55850
    ## 8679      79629
    ## 8680       2063
    ## 8681      83878
    ## 8682      29086
    ## 8683      79575
    ## 8684      64981
    ## 8685      79016
    ## 8686      57719
    ## 8687      84705
    ## 8688      83483
    ## 8689        684
    ## 8690      93343
    ## 8691  100130519
    ## 8692     376497
    ## 8693      25796
    ## 8694       <NA>
    ## 8695      79709
    ## 8696      23025
    ## 8697       3718
    ## 8698      10331
    ## 8699      23149
    ## 8700       <NA>
    ## 8701       <NA>
    ## 8702       <NA>
    ## 8703       <NA>
    ## 8704       <NA>
    ## 8705       7171
    ## 8706       <NA>
    ## 8707       4218
    ## 8708      26017
    ## 8709       8907
    ## 8710      10365
    ## 8711      58513
    ## 8712     125972
    ## 8713       <NA>
    ## 8714      10523
    ## 8715      79939
    ## 8716       <NA>
    ## 8717       9441
    ## 8718      79086
    ## 8719      79041
    ## 8720     284434
    ## 8721      23309
    ## 8722       9215
    ## 8723      10042
    ## 8724      10043
    ## 8725       3162
    ## 8726       4174
    ## 8727      23551
    ## 8728       4306
    ## 8729      79658
    ## 8730       <NA>
    ## 8731      90826
    ## 8732      55751
    ## 8733       1909
    ## 8734     494115
    ## 8735      84068
    ## 8736      11157
    ## 8737       <NA>
    ## 8738     166785
    ## 8739       4086
    ## 8740      54726
    ## 8741       6059
    ## 8742      10393
    ## 8743      64399
    ## 8744       2993
    ## 8745       8467
    ## 8746       2549
    ## 8747       <NA>
    ## 8748      84640
    ## 8749       8821
    ## 8750       3600
    ## 8751       <NA>
    ## 8752      57484
    ## 8753       <NA>
    ## 8754      23158
    ## 8755     255520
    ## 8756       1047
    ## 8757      60592
    ## 8758       4713
    ## 8759       9524
    ## 8760       3337
    ## 8761      10755
    ## 8762       5585
    ## 8763       5731
    ## 8764       <NA>
    ## 8765        976
    ## 8766       <NA>
    ## 8767      22859
    ## 8768       <NA>
    ## 8769      55723
    ## 8770       5566
    ## 8771      90378
    ## 8772       <NA>
    ## 8773     342979
    ## 8774       9466
    ## 8775       <NA>
    ## 8776       <NA>
    ## 8777       5989
    ## 8778      90379
    ## 8779      54862
    ## 8780       <NA>
    ## 8781     342977
    ## 8782       <NA>
    ## 8783      65249
    ## 8784       <NA>
    ## 8785      84245
    ## 8786      81576
    ## 8787        773
    ## 8788       9592
    ## 8789       <NA>
    ## 8790     112939
    ## 8791      55621
    ## 8792       4784
    ## 8793       4066
    ## 8794     199699
    ## 8795      90480
    ## 8796       5886
    ## 8797        811
    ## 8798       2193
    ## 8799     256126
    ## 8800       2639
    ## 8801      10661
    ## 8802       <NA>
    ## 8803      22983
    ## 8804      83546
    ## 8805      10535
    ## 8806       7001
    ## 8807       3726
    ## 8808      29911
    ## 8809       <NA>
    ## 8810      79002
    ## 8811      30000
    ## 8812      84261
    ## 8813       1725
    ## 8814      84292
    ## 8815      51398
    ## 8816       4125
    ## 8817       <NA>
    ## 8818      55737
    ## 8819      23594
    ## 8820      91807
    ## 8821       <NA>
    ## 8822      84706
    ## 8823      10294
    ## 8824      81831
    ## 8825      81533
    ## 8826       5257
    ## 8827      94160
    ## 8828      83752
    ## 8829       <NA>
    ## 8830       <NA>
    ## 8831       9683
    ## 8832       <NA>
    ## 8833        869
    ## 8834       <NA>
    ## 8835     255919
    ## 8836      55027
    ## 8837       <NA>
    ## 8838       <NA>
    ## 8839        113
    ## 8840      29117
    ## 8841      85407
    ## 8842     124460
    ## 8843       1540
    ## 8844       6299
    ## 8845      27324
    ## 8846       <NA>
    ## 8847      80205
    ## 8848       5934
    ## 8849      64400
    ## 8850      23322
    ## 8851      79068
    ## 8852       4313
    ## 8853      54947
    ## 8854       <NA>
    ## 8855     221223
    ## 8856       <NA>
    ## 8857       2775
    ## 8858        267
    ## 8859      11051
    ## 8860      55239
    ## 8861        583
    ## 8862       4504
    ## 8863       <NA>
    ## 8864       <NA>
    ## 8865       <NA>
    ## 8866       9688
    ## 8867       <NA>
    ## 8868       9709
    ## 8869      84166
    ## 8870     221184
    ## 8871       <NA>
    ## 8872      89970
    ## 8873      23568
    ## 8874      51090
    ## 8875       <NA>
    ## 8876       6367
    ## 8877       6376
    ## 8878       6361
    ## 8879      57019
    ## 8880      57017
    ## 8881       5432
    ## 8882      55715
    ## 8883      92922
    ## 8884       9289
    ## 8885     222487
    ## 8886      10300
    ## 8887       3801
    ## 8888       1258
    ## 8889       <NA>
    ## 8890      79650
    ## 8891       4324
    ## 8892      29105
    ## 8893       1459
    ## 8894      29070
    ## 8895       <NA>
    ## 8896      64785
    ## 8897      65009
    ## 8898      79918
    ## 8899      23019
    ## 8900       <NA>
    ## 8901       <NA>
    ## 8902       <NA>
    ## 8903      55238
    ## 8904       2806
    ## 8905       1006
    ## 8906       <NA>
    ## 8907       1009
    ## 8908       1003
    ## 8909     146227
    ## 8910       7084
    ## 8911      51192
    ## 8912     123920
    ## 8913     146223
    ## 8914       1783
    ## 8915       <NA>
    ## 8916     283847
    ## 8917       8883
    ## 8918       <NA>
    ## 8919      57546
    ## 8920       6236
    ## 8921       <NA>
    ## 8922        865
    ## 8923       <NA>
    ## 8924      84752
    ## 8925       8717
    ## 8926      55336
    ## 8927       3299
    ## 8928       8996
    ## 8929       <NA>
    ## 8930       <NA>
    ## 8931       1874
    ## 8932      79767
    ## 8933      26231
    ## 8934      29100
    ## 8935      29109
    ## 8936       6553
    ## 8937      25894
    ## 8938      55282
    ## 8939      51673
    ## 8940      29800
    ## 8941       9114
    ## 8942        181
    ## 8943      79567
    ## 8944       <NA>
    ## 8945      10664
    ## 8946     146206
    ## 8947      65057
    ## 8948      50855
    ## 8949      84080
    ## 8950       <NA>
    ## 8951      81577
    ## 8952      57610
    ## 8953      80152
    ## 8954      57215
    ## 8955      10204
    ## 8956      23644
    ## 8957     123904
    ## 8958       5681
    ## 8959       1506
    ## 8960       5699
    ## 8961       3931
    ## 8962       6560
    ## 8963      54920
    ## 8964      55794
    ## 8965       4775
    ## 8966      23659
    ## 8967       9057
    ## 8968      84138
    ## 8969      54496
    ## 8970      55512
    ## 8971     146198
    ## 8972       <NA>
    ## 8973        999
    ## 8974      79613
    ## 8975       3038
    ## 8976      54921
    ## 8977      84916
    ## 8978       6645
    ## 8979      27183
    ## 8980      84342
    ## 8981      51388
    ## 8982       7014
    ## 8983      80777
    ## 8984      10725
    ## 8985       1728
    ## 8986      28987
    ## 8987      11060
    ## 8988       5713
    ## 8989        463
    ## 8990       <NA>
    ## 8991       <NA>
    ## 8992      83449
    ## 8993       9785
    ## 8994      54957
    ## 8995       3240
    ## 8996       1723
    ## 8997       9798
    ## 8998       <NA>
    ## 8999     342371
    ## 9000        164
    ## 9001      23035
    ## 9002      91862
    ## 9003       <NA>
    ## 9004        794
    ## 9005      55783
    ## 9006      55697
    ## 9007       <NA>
    ## 9008     146433
    ## 9009       <NA>
    ## 9010      23450
    ## 9011      25839
    ## 9012       <NA>
    ## 9013       6483
    ## 9014      55308
    ## 9015      11269
    ## 9016       <NA>
    ## 9017     348174
    ## 9018      55066
    ## 9019       2734
    ## 9020      55159
    ## 9021      79152
    ## 9022      79726
    ## 9023      84937
    ## 9024     197257
    ## 9025     162239
    ## 9026       9564
    ## 9027      10428
    ## 9028       <NA>
    ## 9029      23563
    ## 9030      79583
    ## 9031      11345
    ## 9032      23536
    ## 9033       <NA>
    ## 9034      54386
    ## 9035       <NA>
    ## 9036       <NA>
    ## 9037      85445
    ## 9038      22879
    ## 9039     170692
    ## 9040     283927
    ## 9041      57687
    ## 9042      51741
    ## 9043       <NA>
    ## 9044       4094
    ## 9045       <NA>
    ## 9046      83657
    ## 9047     124359
    ## 9048      56942
    ## 9049      55839
    ## 9050      23300
    ## 9051       <NA>
    ## 9052       2653
    ## 9053       8139
    ## 9054      80790
    ## 9055       5336
    ## 9056      93517
    ## 9057       3294
    ## 9058      10200
    ## 9059       1012
    ## 9060       3281
    ## 9061      23417
    ## 9062      29948
    ## 9063      54550
    ## 9064       8720
    ## 9065      83693
    ## 9066       9013
    ## 9067      93107
    ## 9068      58189
    ## 9069       9914
    ## 9070       <NA>
    ## 9071      23406
    ## 9072      79786
    ## 9073       9100
    ## 9074      83716
    ## 9075      55625
    ## 9076       <NA>
    ## 9077     339145
    ## 9078       <NA>
    ## 9079      23199
    ## 9080      51659
    ## 9081       <NA>
    ## 9082      10328
    ## 9083       1327
    ## 9084       3394
    ## 9085       2294
    ## 9086      64779
    ## 9087       2303
    ## 9088       <NA>
    ## 9089      79791
    ## 9090      81631
    ## 9091      23174
    ## 9092       <NA>
    ## 9093       <NA>
    ## 9094      57338
    ## 9095      54758
    ## 9096       <NA>
    ## 9097       8140
    ## 9098       <NA>
    ## 9099      54971
    ## 9100       <NA>
    ## 9101     161882
    ## 9102       <NA>
    ## 9103     124245
    ## 9104       1535
    ## 9105       4597
    ## 9106     333929
    ## 9107     115992
    ## 9108     348180
    ## 9109       9780
    ## 9110      81620
    ## 9111        353
    ## 9112       2588
    ## 9113      51693
    ## 9114        863
    ## 9115       <NA>
    ## 9116     197322
    ## 9117       1013
    ## 9118      29123
    ## 9119       6687
    ## 9120       6137
    ## 9121      27132
    ## 9122       5119
    ## 9123     124045
    ## 9124       8558
    ## 9125     124044
    ## 9126       <NA>
    ## 9127       9605
    ## 9128       <NA>
    ## 9129       2175
    ## 9130      84501
    ## 9131      22980
    ## 9132      10381
    ## 9133      54849
    ## 9134       <NA>
    ## 9135      79007
    ## 9136       2622
    ## 9137      58480
    ## 9138       5867
    ## 9139       <NA>
    ## 9140     126731
    ## 9141         58
    ## 9142       <NA>
    ## 9143      55746
    ## 9144       <NA>
    ## 9145      23456
    ## 9146      27097
    ## 9147       9816
    ## 9148      79605
    ## 9149      22796
    ## 9150        183
    ## 9151      10753
    ## 9152       <NA>
    ## 9153      79573
    ## 9154      64801
    ## 9155     375061
    ## 9156     440730
    ## 9157       <NA>
    ## 9158       8443
    ## 9159     149371
    ## 9160      83932
    ## 9161      54583
    ## 9162       7257
    ## 9163      57568
    ## 9164      54627
    ## 9165      84284
    ## 9166      80003
    ## 9167       <NA>
    ## 9168      84451
    ## 9169       3775
    ## 9170     148641
    ## 9171     388753
    ## 9172       6894
    ## 9173       <NA>
    ## 9174     359948
    ## 9175       9804
    ## 9176      23029
    ## 9177      56288
    ## 9178       8829
    ## 9179       3688
    ## 9180       <NA>
    ## 9181       <NA>
    ## 9182       <NA>
    ## 9183       <NA>
    ## 9184       <NA>
    ## 9185       <NA>
    ## 9186       <NA>
    ## 9187       <NA>
    ## 9188       <NA>
    ## 9189       <NA>
    ## 9190       <NA>
    ## 9191       <NA>
    ## 9192       <NA>
    ## 9193       <NA>
    ## 9194       <NA>
    ## 9195       <NA>
    ## 9196       2317
    ## 9197      57406
    ## 9198      11102
    ## 9199      54899
    ## 9200       5162
    ## 9201     200845
    ## 9202      11170
    ## 9203       <NA>
    ## 9204       2272
    ## 9205       <NA>
    ## 9206       5793
    ## 9207       <NA>
    ## 9208      55079
    ## 9209       8618
    ## 9210     132204
    ## 9211      80145
    ## 9212       6314
    ## 9213       9861
    ## 9214       3563
    ## 9215       9497
    ## 9216     152110
    ## 9217     116135
    ## 9218      54995
    ## 9219       <NA>
    ## 9220      55768
    ## 9221       7155
    ## 9222       5915
    ## 9223       7068
    ## 9224       9975
    ## 9225       6138
    ## 9226      28512
    ## 9227       7324
    ## 9228       7325
    ## 9229       <NA>
    ## 9230      22795
    ## 9231      51637
    ## 9232       <NA>
    ## 9233      54331
    ## 9234      55776
    ## 9235       8645
    ## 9236      25961
    ## 9237      11319
    ## 9238       <NA>
    ## 9239      23234
    ## 9240      51021
    ## 9241     118491
    ## 9242        310
    ## 9243     118490
    ## 9244       5532
    ## 9245       <NA>
    ## 9246     159195
    ## 9247       <NA>
    ## 9248      79933
    ## 9249       9632
    ## 9250       <NA>
    ## 9251     170384
    ## 9252     118487
    ## 9253      23053
    ## 9254       8509
    ## 9255        818
    ## 9256       <NA>
    ## 9257       5328
    ## 9258       7414
    ## 9259      26985
    ## 9260        132
    ## 9261      23522
    ## 9262      51207
    ## 9263     142891
    ## 9264       7417
    ## 9265     118881
    ## 9266       <NA>
    ## 9267       3778
    ## 9268       9231
    ## 9269      11128
    ## 9270       6229
    ## 9271       <NA>
    ## 9272       <NA>
    ## 9273      57178
    ## 9274      10105
    ## 9275     219654
    ## 9276        311
    ## 9277       <NA>
    ## 9278       7871
    ## 9279     201627
    ## 9280       <NA>
    ## 9281        378
    ## 9282     201626
    ## 9283       <NA>
    ## 9284      26060
    ## 9285      54756
    ## 9286      50650
    ## 9287       <NA>
    ## 9288     285331
    ## 9289      26059
    ## 9290       7474
    ## 9291      55799
    ## 9292      58515
    ## 9293      93973
    ## 9294      55540
    ## 9295      55349
    ## 9296        776
    ## 9297      55802
    ## 9298       7086
    ## 9299       5580
    ## 9300      91869
    ## 9301      51460
    ## 9302       <NA>
    ## 9303     389125
    ## 9304       3699
    ## 9305       6787
    ## 9306      28972
    ## 9307      55830
    ## 9308      26354
    ## 9309      55193
    ## 9310     440957
    ## 9311      64943
    ## 9312      23166
    ## 9313      11188
    ## 9314       7134
    ## 9315      56920
    ## 9316      51533
    ## 9317       8314
    ## 9318      25981
    ## 9319       <NA>
    ## 9320      23473
    ## 9321       9467
    ## 9322     131965
    ## 9323      85403
    ## 9324       8292
    ## 9325      26061
    ## 9326        686
    ## 9327      23243
    ## 9328     117248
    ## 9329     285381
    ## 9330      92106
    ## 9331       8031
    ## 9332  100287932
    ## 9333       8505
    ## 9334       <NA>
    ## 9335      55753
    ## 9336       6572
    ## 9337       2074
    ## 9338       <NA>
    ## 9339     196740
    ## 9340      57705
    ## 9341      58504
    ## 9342       5599
    ## 9343     143162
    ## 9344      26095
    ## 9345       2662
    ## 9346       <NA>
    ## 9347       9721
    ## 9348      83849
    ## 9349       <NA>
    ## 9350       <NA>
    ## 9351       2746
    ## 9352       6623
    ## 9353      79812
    ## 9354        657
    ## 9355      11155
    ## 9356       <NA>
    ## 9357       <NA>
    ## 9358      23063
    ## 9359       <NA>
    ## 9360       2894
    ## 9361      54462
    ## 9362      92211
    ## 9363       <NA>
    ## 9364      27069
    ## 9365      10718
    ## 9366     387694
    ## 9367      81619
    ## 9368       <NA>
    ## 9369      84332
    ## 9370       <NA>
    ## 9371       <NA>
    ## 9372       <NA>
    ## 9373       <NA>
    ## 9374       5729
    ## 9375       5732
    ## 9376      57544
    ## 9377     283554
    ## 9378       <NA>
    ## 9379       5706
    ## 9380       6815
    ## 9381      64841
    ## 9382      10979
    ## 9383       <NA>
    ## 9384      80821
    ## 9385       <NA>
    ## 9386        652
    ## 9387       1033
    ## 9388      10175
    ## 9389       <NA>
    ## 9390       2764
    ## 9391      10668
    ## 9392       <NA>
    ## 9393       <NA>
    ## 9394       <NA>
    ## 9395       2643
    ## 9396      11169
    ## 9397     122809
    ## 9398      93487
    ## 9399       3958
    ## 9400       9787
    ## 9401      55030
    ## 9402      22863
    ## 9403       3895
    ## 9404      57161
    ## 9405      54916
    ## 9406      10640
    ## 9407      55745
    ## 9408     122830
    ## 9409     341880
    ## 9410       <NA>
    ## 9411      91875
    ## 9412      57820
    ## 9413      10038
    ## 9414       7011
    ## 9415       <NA>
    ## 9416     123103
    ## 9417      55644
    ## 9418        328
    ## 9419       <NA>
    ## 9420       4860
    ## 9421       <NA>
    ## 9422       6038
    ## 9423        283
    ## 9424       6039
    ## 9425       6035
    ## 9426      64745
    ## 9427      29986
    ## 9428      57447
    ## 9429      55701
    ## 9430       <NA>
    ## 9431       <NA>
    ## 9432       3183
    ## 9433      57096
    ## 9434       <NA>
    ## 9435      57680
    ## 9436      84932
    ## 9437       9878
    ## 9438      56339
    ## 9439       6297
    ## 9440      28755
    ## 9441       1603
    ## 9442      63874
    ## 9443       5018
    ## 9444       9056
    ## 9445     122704
    ## 9446       4323
    ## 9447      26020
    ## 9448     161253
    ## 9449      10419
    ## 9450      54930
    ## 9451      84962
    ## 9452       <NA>
    ## 9453       5693
    ## 9454      64403
    ## 9455      22985
    ## 9456       <NA>
    ## 9457      23428
    ## 9458      57594
    ## 9459      90673
    ## 9460        599
    ## 9461       8106
    ## 9462      51310
    ## 9463      10278
    ## 9464     116173
    ## 9465       4625
    ## 9466      25983
    ## 9467       <NA>
    ## 9468      85446
    ## 9469       8906
    ## 9470      84502
    ## 9471      10901
    ## 9472      90668
    ## 9473       9362
    ## 9474       4901
    ## 9475       5106
    ## 9476      80344
    ## 9477       5720
    ## 9478      51016
    ## 9479       5721
    ## 9480      55072
    ## 9481      10379
    ## 9482       9985
    ## 9483      79711
    ## 9484      10548
    ## 9485       <NA>
    ## 9486     283629
    ## 9487     145553
    ## 9488       4738
    ## 9489      51292
    ## 9490      26277
    ## 9491       5875
    ## 9492     115817
    ## 9493     161424
    ## 9494      27141
    ## 9495       <NA>
    ## 9496     196883
    ## 9497       4776
    ## 9498      57523
    ## 9499     643866
    ## 9500      23351
    ## 9501      56948
    ## 9502      56163
    ## 9503      55835
    ## 9504        143
    ## 9505      54737
    ## 9506       <NA>
    ## 9507      55269
    ## 9508       9205
    ## 9509       7750
    ## 9510       2706
    ## 9511       <NA>
    ## 9512      10804
    ## 9513      51084
    ## 9514       8100
    ## 9515      53342
    ## 9516     221143
    ## 9517      64328
    ## 9518      26524
    ## 9519      10284
    ## 9520     221150
    ## 9521      78988
    ## 9522     253832
    ## 9523     221154
    ## 9524       2254
    ## 9525      55213
    ## 9526       <NA>
    ## 9527       <NA>
    ## 9528       <NA>
    ## 9529      83852
    ## 9530      81617
    ## 9531      81602
    ## 9532     387914
    ## 9533      51761
    ## 9534       <NA>
    ## 9535       9107
    ## 9536     219287
    ## 9537     221178
    ## 9538       4285
    ## 9539      55504
    ## 9540      26278
    ## 9541       6445
    ## 9542     115761
    ## 9543      84650
    ## 9544       3839
    ## 9545      57213
    ## 9546      10206
    ## 9547     220107
    ## 9548      79621
    ## 9549     220108
    ## 9550      26512
    ## 9551       <NA>
    ## 9552       <NA>
    ## 9553     115825
    ## 9554       1508
    ## 9555       2222
    ## 9556     252969
    ## 9557      83648
    ## 9558      66036
    ## 9559     286046
    ## 9560      54984
    ## 9561      83595
    ## 9562       <NA>
    ## 9563       <NA>
    ## 9564       4482
    ## 9565       <NA>
    ## 9566       <NA>
    ## 9567      23303
    ## 9568      79618
    ## 9569      55756
    ## 9570       2137
    ## 9571       7976
    ## 9572     157574
    ## 9573       <NA>
    ## 9574       5368
    ## 9575      55140
    ## 9576      55872
    ## 9577     157570
    ## 9578      55246
    ## 9579      51435
    ## 9580       1191
    ## 9581       2053
    ## 9582       1135
    ## 9583       2185
    ## 9584      23087
    ## 9585      81551
    ## 9586        148
    ## 9587       1808
    ## 9588      10687
    ## 9589        665
    ## 9590       5520
    ## 9591      64641
    ## 9592      54793
    ## 9593       2796
    ## 9594      80005
    ## 9595       <NA>
    ## 9596       4741
    ## 9597       4747
    ## 9598       6781
    ## 9599       4824
    ## 9600      51312
    ## 9601       4017
    ## 9602     203069
    ## 9603      91782
    ## 9604       <NA>
    ## 9605       8795
    ## 9606      23221
    ## 9607       1960
    ## 9608      55909
    ## 9609      57805
    ## 9610       <NA>
    ## 9611      64236
    ## 9612      10174
    ## 9613       5533
    ## 9614      23516
    ## 9615      55124
    ## 9616        661
    ## 9617       9796
    ## 9618        649
    ## 9619       6440
    ## 9620     203190
    ## 9621      80346
    ## 9622      55806
    ## 9623      79873
    ## 9624      64760
    ## 9625       2039
    ## 9626       8822
    ## 9627      10361
    ## 9628      23039
    ## 9629       9046
    ## 9630       2675
    ## 9631      22862
    ## 9632      57105
    ## 9633       1102
    ## 9634       5925
    ## 9635      10161
    ## 9636       9445
    ## 9637      29079
    ## 9638      55270
    ## 9639       8803
    ## 9640       3356
    ## 9641       2098
    ## 9642      23143
    ## 9643      80183
    ## 9644       3936
    ## 9645      23091
    ## 9646     283514
    ## 9647      83548
    ## 9648     253512
    ## 9649       <NA>
    ## 9650       7178
    ## 9651       2963
    ## 9652     386618
    ## 9653      55425
    ## 9654      26747
    ## 9655       <NA>
    ## 9656       8848
    ## 9657     387923
    ## 9658       <NA>
    ## 9659     144811
    ## 9660     160857
    ## 9661      55068
    ## 9662      29103
    ## 9663      94240
    ## 9664      11215
    ## 9665     160851
    ## 9666      23078
    ## 9667       <NA>
    ## 9668      28984
    ## 9669      79612
    ## 9670       9617
    ## 9671      84078
    ## 9672      89890
    ## 9673      11193
    ## 9674       1997
    ## 9675      10910
    ## 9676      11061
    ## 9677       <NA>
    ## 9678       5100
    ## 9679       <NA>
    ## 9680      10562
    ## 9681      27253
    ## 9682       <NA>
    ## 9683       <NA>
    ## 9684      81624
    ## 9685      81550
    ## 9686      64881
    ## 9687       5101
    ## 9688      57626
    ## 9689       1602
    ## 9690       <NA>
    ## 9691     440145
    ## 9692      79866
    ## 9693      22894
    ## 9694      10464
    ## 9695        688
    ## 9696      11278
    ## 9697       9882
    ## 9698     170622
    ## 9699       7347
    ## 9700       4008
    ## 9701     115207
    ## 9702       1203
    ## 9703      26224
    ## 9704      23077
    ## 9705       <NA>
    ## 9706     122060
    ## 9707       1910
    ## 9708       <NA>
    ## 9709      64062
    ## 9710      54602
    ## 9711       <NA>
    ## 9712      10253
    ## 9713     114798
    ## 9714      84189
    ## 9715       <NA>
    ## 9716      26050
    ## 9717       <NA>
    ## 9718       2262
    ## 9719      10082
    ## 9720       1638
    ## 9721      23483
    ## 9722     160897
    ## 9723      11166
    ## 9724       <NA>
    ## 9725      10257
    ## 9726       9071
    ## 9727      22873
    ## 9728       <NA>
    ## 9729       5611
    ## 9730      55757
    ## 9731      27199
    ## 9732      10150
    ## 9733       5911
    ## 9734       3843
    ## 9735      10160
    ## 9736       8428
    ## 9737      23348
    ## 9738       <NA>
    ## 9739     337867
    ## 9740       1880
    ## 9741       9375
    ## 9742     171425
    ## 9743       <NA>
    ## 9744      85416
    ## 9745       7546
    ## 9746       5095
    ## 9747      87769
    ## 9748      84899
    ## 9749     259232
    ## 9750       9358
    ## 9751       2259
    ## 9752      91801
    ## 9753     143884
    ## 9754       2977
    ## 9755      60496
    ## 9756     143879
    ## 9757      84437
    ## 9758       2893
    ## 9759        834
    ## 9760        837
    ## 9761  100506742
    ## 9762      80310
    ## 9763      79659
    ## 9764      84259
    ## 9765     114908
    ## 9766        329
    ## 9767        330
    ## 9768      10413
    ## 9769       <NA>
    ## 9770      57562
    ## 9771       7225
    ## 9772       <NA>
    ## 9773       5241
    ## 9774       <NA>
    ## 9775     143872
    ## 9776      53942
    ## 9777       8690
    ## 9778      79780
    ## 9779      84441
    ## 9780       8898
    ## 9781       9702
    ## 9782     143684
    ## 9783     143686
    ## 9784      23052
    ## 9785      55693
    ## 9786      51503
    ## 9787     154810
    ## 9788      54851
    ## 9789       <NA>
    ## 9790      10888
    ## 9791      24145
    ## 9792     387804
    ## 9793       9440
    ## 9794       <NA>
    ## 9795      79101
    ## 9796      85459
    ## 9797      56935
    ## 9798     120103
    ## 9799     120114
    ## 9800      26973
    ## 9801      10003
    ## 9802       <NA>
    ## 9803       <NA>
    ## 9804       <NA>
    ## 9805       <NA>
    ## 9806       <NA>
    ## 9807       <NA>
    ## 9808       <NA>
    ## 9809       <NA>
    ## 9810       <NA>
    ## 9811       <NA>
    ## 9812       <NA>
    ## 9813      54850
    ## 9814      59286
    ## 9815       5300
    ## 9816      93145
    ## 9817      50509
    ## 9818       <NA>
    ## 9819      83854
    ## 9820      56342
    ## 9821       8666
    ## 9822       1786
    ## 9823       9294
    ## 9824      51073
    ## 9825       3383
    ## 9826       <NA>
    ## 9827       3386
    ## 9828       7087
    ## 9829  100125288
    ## 9830       <NA>
    ## 9831       <NA>
    ## 9832       <NA>
    ## 9833       7297
    ## 9834      11140
    ## 9835       5141
    ## 9836       <NA>
    ## 9837       9817
    ## 9838      53637
    ## 9839      84971
    ## 9840      65095
    ## 9841       1032
    ## 9842      10053
    ## 9843       <NA>
    ## 9844      57153
    ## 9845       3609
    ## 9846      81890
    ## 9847       1785
    ## 9848      11018
    ## 9849       <NA>
    ## 9850      10498
    ## 9851      78992
    ## 9852      90580
    ## 9853       6597
    ## 9854       3949
    ## 9855     147841
    ## 9856      25959
    ## 9857      57572
    ## 9858       9545
    ## 9859     374882
    ## 9860     126075
    ## 9861      64748
    ## 9862     126074
    ## 9863       2057
    ## 9864      57139
    ## 9865     115948
    ## 9866       5589
    ## 9867       1995
    ## 9868       <NA>
    ## 9869       <NA>
    ## 9870      51295
    ## 9871       1264
    ## 9872      84337
    ## 9873         54
    ## 9874       <NA>
    ## 9875       <NA>
    ## 9876       <NA>
    ## 9877       <NA>
    ## 9878      54443
    ## 9879       6100
    ## 9880       <NA>
    ## 9881      27241
    ## 9882     168667
    ## 9883     387129
    ## 9884      23333
    ## 9885     283417
    ## 9886      64224
    ## 9887       <NA>
    ## 9888      80820
    ## 9889       <NA>
    ## 9890       <NA>
    ## 9891      27087
    ## 9892      27034
    ## 9893      29087
    ## 9894     112936
    ## 9895      23310
    ## 9896      83700
    ## 9897      22997
    ## 9898     219938
    ## 9899       4978
    ## 9900      50863
    ## 9901     399979
    ## 9902     170689
    ## 9903      11095
    ## 9904      29068
    ## 9905       6768
    ## 9906        334
    ## 9907       <NA>
    ## 9908      56980
    ## 9909       4798
    ## 9910       <NA>
    ## 9911       8538
    ## 9912       9743
    ## 9913       3762
    ## 9914       2313
    ## 9915       2113
    ## 9916      84623
    ## 9917       6484
    ## 9918      28960
    ## 9919     114609
    ## 9920       <NA>
    ## 9921      55572
    ## 9922      79607
    ## 9923      84881
    ## 9924      50937
    ## 9925      29118
    ## 9926      83480
    ## 9927     219844
    ## 9928     399967
    ## 9929       1111
    ## 9930       3703
    ## 9931       9538
    ## 9932       9638
    ## 9933      63876
    ## 9934     219854
    ## 9935     219855
    ## 9936      80071
    ## 9937     220296
    ## 9938      54538
    ## 9939      64221
    ## 9940      79684
    ## 9941      90952
    ## 9942      23584
    ## 9943       4900
    ## 9944      54414
    ## 9945      53340
    ## 9946      84897
    ## 9947       <NA>
    ## 9948       4013
    ## 9949       <NA>
    ## 9950       <NA>
    ## 9951       <NA>
    ## 9952      55800
    ## 9953      57476
    ## 9954       <NA>
    ## 9955      79827
    ## 9956       3312
    ## 9957      84959
    ## 9958       <NA>
    ## 9959       <NA>
    ## 9960       <NA>
    ## 9961       <NA>
    ## 9962       6653
    ## 9963       6309
    ## 9964     219899
    ## 9965       2900
    ## 9966      23365
    ## 9967       <NA>
    ## 9968      25833
    ## 9969     220323
    ## 9970       <NA>
    ## 9971       5818
    ## 9972       7070
    ## 9973       <NA>
    ## 9974       9099
    ## 9975      79102
    ## 9976       <NA>
    ## 9977       <NA>
    ## 9978       4162
    ## 9979        867
    ## 9980     283152
    ## 9981      79849
    ## 9982      79671
    ## 9983      64137
    ## 9984       <NA>
    ## 9985      25988
    ## 9986       9854
    ## 9987       1798
    ## 9988       <NA>
    ## 9989       3145
    ## 9990      55823
    ## 9991      10525
    ## 9992       2542
    ## 9993      51399
    ## 9994       6230
    ## 9995     338657
    ## 9996       <NA>
    ## 9997       <NA>
    ## 9998     283149
    ## 9999       <NA>
    ## 10000      1656
    ## 10001     23187
    ## 10002       372
    ## 10003     56912
    ## 10004     84866
    ## 10005    143941
    ## 10006      4297
    ## 10007      <NA>
    ## 10008      <NA>
    ## 10009      9354
    ## 10010       915
    ## 10011     10205
    ## 10012    196264
    ## 10013    120425
    ## 10014      6327
    ## 10015      6330
    ## 10016     56649
    ## 10017      <NA>
    ## 10018      3587
    ## 10019     53826
    ## 10020       486
    ## 10021     57453
    ## 10022     22897
    ## 10023     23621
    ## 10024      <NA>
    ## 10025    257160
    ## 10026      9159
    ## 10027      6876
    ## 10028      <NA>
    ## 10029     51092
    ## 10030      <NA>
    ## 10031      5049
    ## 10032     23387
    ## 10033       335
    ## 10034       345
    ## 10035      8882
    ## 10036     84811
    ## 10037      <NA>
    ## 10038      <NA>
    ## 10039     23705
    ## 10040      <NA>
    ## 10041     54827
    ## 10042    120406
    ## 10043     25996
    ## 10044     10179
    ## 10045      <NA>
    ## 10046      4837
    ## 10047      7704
    ## 10048      3359
    ## 10049     57646
    ## 10050      9183
    ## 10051     80975
    ## 10052      1813
    ## 10053     54970
    ## 10054      4684
    ## 10055      <NA>
    ## 10056      5805
    ## 10057     83875
    ## 10058      3606
    ## 10059      6392
    ## 10060     26521
    ## 10061      <NA>
    ## 10062    120379
    ## 10063      1737
    ## 10064     85458
    ## 10065      <NA>
    ## 10066      3316
    ## 10067      1410
    ## 10068      <NA>
    ## 10069     91893
    ## 10070     79796
    ## 10071      5519
    ## 10072     23235
    ## 10073    143903
    ## 10074    120376
    ## 10075      <NA>
    ## 10076     57569
    ## 10077      <NA>
    ## 10078      2230
    ## 10079      5962
    ## 10080     85463
    ## 10081      <NA>
    ## 10082      1662
    ## 10083     23086
    ## 10084      <NA>
    ## 10085      <NA>
    ## 10086       472
    ## 10087      4863
    ## 10088        38
    ## 10089      8065
    ## 10090      <NA>
    ## 10091      <NA>
    ## 10092     54733
    ## 10093      6588
    ## 10094     55531
    ## 10095    388121
    ## 10096    342035
    ## 10097     23312
    ## 10098     10518
    ## 10099      3419
    ## 10100     23205
    ## 10101     55466
    ## 10102     80349
    ## 10103      1381
    ## 10104      3658
    ## 10105    123688
    ## 10106      <NA>
    ## 10107      5685
    ## 10108      1138
    ## 10109      1136
    ## 10110      1143
    ## 10111     92912
    ## 10112     26263
    ## 10113    145957
    ## 10114    123591
    ## 10115      2108
    ## 10116     49855
    ## 10117      5955
    ## 10118      9051
    ## 10119     10099
    ## 10120     79834
    ## 10121      <NA>
    ## 10122     10363
    ## 10123     84894
    ## 10124      1464
    ## 10125    257364
    ## 10126     55272
    ## 10127     10073
    ## 10128      5780
    ## 10129      <NA>
    ## 10130     25942
    ## 10131      <NA>
    ## 10132      <NA>
    ## 10133      4123
    ## 10134     79661
    ## 10135     54939
    ## 10136      <NA>
    ## 10137     60490
    ## 10138    192683
    ## 10139     54913
    ## 10140      9377
    ## 10141     57184
    ## 10142      4351
    ## 10143     10066
    ## 10144     25989
    ## 10145    594855
    ## 10146      1445
    ## 10147     80153
    ## 10148      1198
    ## 10149     10620
    ## 10150     84993
    ## 10151      8482
    ## 10152      <NA>
    ## 10153      1583
    ## 10154     80125
    ## 10155     64220
    ## 10156      3671
    ## 10157     57611
    ## 10158      5371
    ## 10159      9399
    ## 10160      4016
    ## 10161      <NA>
    ## 10162     80381
    ## 10163     27020
    ## 10164    283677
    ## 10165     10021
    ## 10166      4756
    ## 10167     83440
    ## 10168       585
    ## 10169     25820
    ## 10170    338949
    ## 10171      3073
    ## 10172     60677
    ## 10173     56965
    ## 10174      5315
    ## 10175      <NA>
    ## 10176    123228
    ## 10177      4649
    ## 10178     79875
    ## 10179      <NA>
    ## 10180     54839
    ## 10181      <NA>
    ## 10182     55323
    ## 10183     55075
    ## 10184      7090
    ## 10185      6176
    ## 10186      9493
    ## 10187     54852
    ## 10188     26035
    ## 10189      <NA>
    ## 10190      8125
    ## 10191     10391
    ## 10192     22801
    ## 10193     10116
    ## 10194     54982
    ## 10195     91860
    ## 10196      8554
    ## 10197      5607
    ## 10198      <NA>
    ## 10199     79719
    ## 10200      4088
    ## 10201      4091
    ## 10202     55055
    ## 10203      6124
    ## 10204     10302
    ## 10205      5604
    ## 10206      <NA>
    ## 10207     54962
    ## 10208    115752
    ## 10209     84465
    ## 10210      8766
    ## 10211     10260
    ## 10212     81556
    ## 10213     51495
    ## 10214     54878
    ## 10215     57722
    ## 10216      9543
    ## 10217     54956
    ## 10218     10845
    ## 10219     10081
    ## 10220     51285
    ## 10221    123263
    ## 10222      <NA>
    ## 10223     51324
    ## 10224    348094
    ## 10225     80301
    ## 10226    348093
    ## 10227      4947
    ## 10228      <NA>
    ## 10229      9325
    ## 10230      9768
    ## 10231     53944
    ## 10232      5479
    ## 10233     79856
    ## 10234      6642
    ## 10235      <NA>
    ## 10236     23604
    ## 10237      8925
    ## 10238    283807
    ## 10239      9960
    ## 10240      <NA>
    ## 10241     83464
    ## 10242      <NA>
    ## 10243     51762
    ## 10244     51065
    ## 10245    114294
    ## 10246      7168
    ## 10247     83660
    ## 10248      <NA>
    ## 10249    388125
    ## 10250    145741
    ## 10251     54832
    ## 10252      6095
    ## 10253     79664
    ## 10254       302
    ## 10255       663
    ## 10256      2958
    ## 10257    145773
    ## 10258      4643
    ## 10259      9133
    ## 10260     54778
    ## 10261     79811
    ## 10262     54629
    ## 10263      <NA>
    ## 10264       102
    ## 10265      3990
    ## 10266       366
    ## 10267      8854
    ## 10268     81488
    ## 10269 100820829
    ## 10270     84952
    ## 10271      6938
    ## 10272      <NA>
    ## 10273      <NA>
    ## 10274     55329
    ## 10275    374618
    ## 10276     64864
    ## 10277      4734
    ## 10278    283659
    ## 10279     26108
    ## 10280      <NA>
    ## 10281      <NA>
    ## 10282      9236
    ## 10283      9488
    ## 10284      5873
    ## 10285     51187
    ## 10286    440279
    ## 10287     56204
    ## 10288     10776
    ## 10289      4644
    ## 10290     55930
    ## 10291     10681
    ## 10292      <NA>
    ## 10293      5597
    ## 10294    123169
    ## 10295     29766
    ## 10296     29767
    ## 10297    256586
    ## 10298     29106
    ## 10299       653
    ## 10300     54511
    ## 10301      3062
    ## 10302     90523
    ## 10303     55227
    ## 10304      2729
    ## 10305     60481
    ## 10306     26268
    ## 10307     22858
    ## 10308      <NA>
    ## 10309      2941
    ## 10310      <NA>
    ## 10311    441161
    ## 10312      <NA>
    ## 10313      <NA>
    ## 10314     25821
    ## 10315      1915
    ## 10316     26503
    ## 10317    135228
    ## 10318      <NA>
    ## 10319      <NA>
    ## 10320      1303
    ## 10321      1347
    ## 10322     55754
    ## 10323     27145
    ## 10324     26054
    ## 10325      4646
    ## 10326      3351
    ## 10327      <NA>
    ## 10328 101928601
    ## 10329    134728
    ## 10330     55023
    ## 10331      9324
    ## 10332      <NA>
    ## 10333    167691
    ## 10334     83699
    ## 10335      6785
    ## 10336      7272
    ## 10337       594
    ## 10338      <NA>
    ## 10339     25998
    ## 10340      7162
    ## 10341      <NA>
    ## 10342      <NA>
    ## 10343      5238
    ## 10344    112611
    ## 10345      4199
    ## 10346      <NA>
    ## 10347    167681
    ## 10348      9892
    ## 10349    134701
    ## 10350     51167
    ## 10351    112609
    ## 10352     22832
    ## 10353      9096
    ## 10354      4907
    ## 10355     57231
    ## 10356     10492
    ## 10357      <NA>
    ## 10358      <NA>
    ## 10359      <NA>
    ## 10360      <NA>
    ## 10361      <NA>
    ## 10362      <NA>
    ## 10363      <NA>
    ## 10364     10588
    ## 10365      <NA>
    ## 10366     23423
    ## 10367    390616
    ## 10368      5923
    ## 10369      1512
    ## 10370     10933
    ## 10371     23102
    ## 10372      <NA>
    ## 10373      7545
    ## 10374     84107
    ## 10375      <NA>
    ## 10376      5359
    ## 10377     57047
    ## 10378     57088
    ## 10379      5352
    ## 10380      <NA>
    ## 10381    285195
    ## 10382      9435
    ## 10383     23350
    ## 10384      <NA>
    ## 10385    344838
    ## 10386     26577
    ## 10387      7220
    ## 10388      5357
    ## 10389      <NA>
    ## 10390       545
    ## 10391      <NA>
    ## 10392     54464
    ## 10393    256356
    ## 10394      7029
    ## 10395       483
    ## 10396      9616
    ## 10397      5922
    ## 10398    253461
    ## 10399     92370
    ## 10400     92369
    ## 10401     55186
    ## 10402     64084
    ## 10403    349565
    ## 10404      5947
    ## 10405      9276
    ## 10406     56945
    ## 10407      <NA>
    ## 10408       668
    ## 10409     55179
    ## 10410      5291
    ## 10411     80321
    ## 10412     83850
    ## 10413     22808
    ## 10414     25852
    ## 10415     51163
    ## 10416    199221
    ## 10417     53833
    ## 10418      4690
    ## 10419     80723
    ## 10420     10274
    ## 10421      5096
    ## 10422     55167
    ## 10423      5523
    ## 10424      2047
    ## 10425    339855
    ## 10426     80254
    ## 10427     25847
    ## 10428     51421
    ## 10429      6259
    ## 10430      6578
    ## 10431     51560
    ## 10432     58477
    ## 10433      <NA>
    ## 10434     11073
    ## 10435     55573
    ## 10436      8419
    ## 10437     66000
    ## 10438     27031
    ## 10439     79876
    ## 10440     84129
    ## 10441     23317
    ## 10442    131034
    ## 10443     11222
    ## 10444    131870
    ## 10445     79858
    ## 10446     28990
    ## 10447     27032
    ## 10448     30849
    ## 10449      <NA>
    ## 10450      <NA>
    ## 10451    132158
    ## 10452     80335
    ## 10453    132160
    ## 10454     11344
    ## 10455       211
    ## 10456      <NA>
    ## 10457     25886
    ## 10458      1849
    ## 10459      6159
    ## 10460        95
    ## 10461     25864
    ## 10462     84836
    ## 10463     57060
    ## 10464    118442
    ## 10465     10039
    ## 10466      9136
    ## 10467      2912
    ## 10468     51368
    ## 10469     23132
    ## 10470      9730
    ## 10471      7873
    ## 10472     29890
    ## 10473      <NA>
    ## 10474      1795
    ## 10475      7867
    ## 10476      1154
    ## 10477     51409
    ## 10478      <NA>
    ## 10479      9254
    ## 10480     11070
    ## 10481     11068
    ## 10482     10641
    ## 10483     51364
    ## 10484     11186
    ## 10485      <NA>
    ## 10486     11334
    ## 10487      8692
    ## 10488      3373
    ## 10489      <NA>
    ## 10490      8372
    ## 10491      7866
    ## 10492      <NA>
    ## 10493      7869
    ## 10494      2771
    ## 10495     10991
    ## 10496      6405
    ## 10497     10181
    ## 10498     10180
    ## 10499     84315
    ## 10500      4486
    ## 10501     79012
    ## 10502     10293
    ## 10503      7318
    ## 10504      <NA>
    ## 10505    389118
    ## 10506      9807
    ## 10507     29925
    ## 10508     63891
    ## 10509    386724
    ## 10510       327
    ## 10511      8927
    ## 10512      1605
    ## 10513     84276
    ## 10514       275
    ## 10515      6988
    ## 10516       387
    ## 10517      2876
    ## 10518      7375
    ## 10519      <NA>
    ## 10520    339834
    ## 10521    200942
    ## 10522     64925
    ## 10523      3913
    ## 10524     10869
    ## 10525      <NA>
    ## 10526     54870
    ## 10527      3615
    ## 10528     25915
    ## 10529      <NA>
    ## 10530     55152
    ## 10531     11180
    ## 10532     54681
    ## 10533     10425
    ## 10534       788
    ## 10535      5576
    ## 10536     51447
    ## 10537     51517
    ## 10538      1951
    ## 10539     65010
    ## 10540      7384
    ## 10541      1294
    ## 10542      5210
    ## 10543     51246
    ## 10544     84126
    ## 10545     51372
    ## 10546     79714
    ## 10547      <NA>
    ## 10548      5364
    ## 10549    646424
    ## 10550     10201
    ## 10551       820
    ## 10552       993
    ## 10553      4134
    ## 10554     22907
    ## 10555      6599
    ## 10556     10675
    ## 10557     54859
    ## 10558     22937
    ## 10559     25930
    ## 10560      <NA>
    ## 10561     23276
    ## 10562     64147
    ## 10563     29072
    ## 10564      <NA>
    ## 10565     23218
    ## 10566    151903
    ## 10567      5745
    ## 10568      <NA>
    ## 10569     29122
    ## 10570    259236
    ## 10571    259173
    ## 10572      6997
    ## 10573     79442
    ## 10574      9034
    ## 10575      9209
    ## 10576      4292
    ## 10577      9852
    ## 10578      9881
    ## 10579     85443
    ## 10580      6769
    ## 10581     10777
    ## 10582      <NA>
    ## 10583      <NA>
    ## 10584      <NA>
    ## 10585     10015
    ## 10586     23122
    ## 10587      7342
    ## 10588     25827
    ## 10589      <NA>
    ## 10590     26032
    ## 10591     10491
    ## 10592      2720
    ## 10593     25904
    ## 10594     51143
    ## 10595     54918
    ## 10596    112616
    ## 10597      <NA>
    ## 10598    152189
    ## 10599     23171
    ## 10600    114884
    ## 10601      <NA>
    ## 10602    201595
    ## 10603    339896
    ## 10604      7048
    ## 10605     27303
    ## 10606     64343
    ## 10607    152100
    ## 10608      <NA>
    ## 10609      8320
    ## 10610      2803
    ## 10611      3680
    ## 10612      <NA>
    ## 10613     10217
    ## 10614     50853
    ## 10615      5333
    ## 10616      9940
    ## 10617      <NA>
    ## 10618      9943
    ## 10619      4615
    ## 10620      <NA>
    ## 10621      9942
    ## 10622        93
    ## 10623      9941
    ## 10624      6331
    ## 10625      6336
    ## 10626     57599
    ## 10627      1524
    ## 10628     64689
    ## 10629     64651
    ## 10630     54977
    ## 10631      3921
    ## 10632      4336
    ## 10633      <NA>
    ## 10634     25924
    ## 10635      <NA>
    ## 10636     10289
    ## 10637       956
    ## 10638      9045
    ## 10639      <NA>
    ## 10640      <NA>
    ## 10641      <NA>
    ## 10642      1499
    ## 10643     54986
    ## 10644     22906
    ## 10645       885
    ## 10646      <NA>
    ## 10647    131375
    ## 10648      7433
    ## 10649      9117
    ## 10650     51188
    ## 10651      4820
    ## 10652      <NA>
    ## 10653      <NA>
    ## 10654    131377
    ## 10655     57467
    ## 10656    152206
    ## 10657     25994
    ## 10658      <NA>
    ## 10659      <NA>
    ## 10660      1238
    ## 10661      <NA>
    ## 10662      <NA>
    ## 10663     84892
    ## 10664     54861
    ## 10665     55129
    ## 10666      <NA>
    ## 10667     51099
    ## 10668      <NA>
    ## 10669    285343
    ## 10670      <NA>
    ## 10671     55888
    ## 10672      <NA>
    ## 10673      <NA>
    ## 10674      <NA>
    ## 10675    131616
    ## 10676     51304
    ## 10677     23016
    ## 10678      7123
    ## 10679     64866
    ## 10680     25907
    ## 10681     23395
    ## 10682      8994
    ## 10683     22908
    ## 10684      <NA>
    ## 10685     54585
    ## 10686     79443
    ## 10687      1230
    ## 10688    729230
    ## 10689      1234
    ## 10690      <NA>
    ## 10691      <NA>
    ## 10692      <NA>
    ## 10693      9814
    ## 10694      4733
    ## 10695     56478
    ## 10696     23598
    ## 10697    113791
    ## 10698      3985
    ## 10699     91445
    ## 10700     50487
    ## 10701     27124
    ## 10702    140606
    ## 10703      6525
    ## 10704     55000
    ## 10705      <NA>
    ## 10706     23762
    ## 10707    150290
    ## 10708    339665
    ## 10709      6948
    ## 10710     23481
    ## 10711      9514
    ## 10712     51537
    ## 10713     23541
    ## 10714    200312
    ## 10715    550631
    ## 10716     10291
    ## 10717     83874
    ## 10718      <NA>
    ## 10719      8897
    ## 10720     84164
    ## 10721     29796
    ## 10722     55954
    ## 10723    164633
    ## 10724      4771
    ## 10725      8508
    ## 10726      8563
    ## 10727      4744
    ## 10728       162
    ## 10729     10634
    ## 10730     10633
    ## 10731      2130
    ## 10732     25807
    ## 10733    129080
    ## 10734     83999
    ## 10735     84133
    ## 10736      7494
    ## 10737    150275
    ## 10738     64951
    ## 10739     55665
    ## 10740     28988
    ## 10741      5224
    ## 10742     27434
    ## 10743       165
    ## 10744      5425
    ## 10745      2645
    ## 10746     10652
    ## 10747       816
    ## 10748     23386
    ## 10749     29881
    ## 10750     54606
    ## 10751    222068
    ## 10752      4967
    ## 10753     83637
    ## 10754      5478
    ## 10755      <NA>
    ## 10756      5814
    ## 10757      <NA>
    ## 10758     64005
    ## 10759     83605
    ## 10760     23148
    ## 10761      9238
    ## 10762     10268
    ## 10763       107
    ## 10764      3486
    ## 10765     64759
    ## 10766      3364
    ## 10767      <NA>
    ## 10768      7378
    ## 10769    375567
    ## 10770     11055
    ## 10771     10320
    ## 10772     63979
    ## 10773      1644
    ## 10774      2887
    ## 10775     23242
    ## 10776    222008
    ## 10777     23480
    ## 10778      1956
    ## 10779      5341
    ## 10780     25927
    ## 10781      5534
    ## 10782    116143
    ## 10783     56902
    ## 10784     10438
    ## 10785     54465
    ## 10786      <NA>
    ## 10787      4211
    ## 10788    200734
    ## 10789     10097
    ## 10790      5861
    ## 10791     23177
    ## 10792      6509
    ## 10793      9792
    ## 10794     54812
    ## 10795     29094
    ## 10796     57162
    ## 10797     51542
    ## 10798      7360
    ## 10799      4190
    ## 10800     51057
    ## 10801      5013
    ## 10802     23301
    ## 10803    200728
    ## 10804      <NA>
    ## 10805    150684
    ## 10806      <NA>
    ## 10807     10575
    ## 10808     84140
    ## 10809      7514
    ## 10810      9736
    ## 10811      <NA>
    ## 10812      <NA>
    ## 10813      5194
    ## 10814    150962
    ## 10815      5966
    ## 10816     64895
    ## 10817      <NA>
    ## 10818     53335
    ## 10819      <NA>
    ## 10820     55120
    ## 10821      7444
    ## 10822    114800
    ## 10823      2202
    ## 10824     87178
    ## 10825     57223
    ## 10826    112942
    ## 10827     55704
    ## 10828      <NA>
    ## 10829      4528
    ## 10830      6233
    ## 10831     57142
    ## 10832    400954
    ## 10833      <NA>
    ## 10834      6711
    ## 10835        98
    ## 10836     23198
    ## 10837     10936
    ## 10838     51130
    ## 10839     27248
    ## 10840    494143
    ## 10841      8614
    ## 10842     91272
    ## 10843      <NA>
    ## 10844     80315
    ## 10845     51617
    ## 10846     79622
    ## 10847     64285
    ## 10848      8131
    ## 10849      4350
    ## 10850      <NA>
    ## 10851      <NA>
    ## 10852      <NA>
    ## 10853      <NA>
    ## 10854    285590
    ## 10855     92181
    ## 10856      6793
    ## 10857     23291
    ## 10858      8817
    ## 10859      4869
    ## 10860     64901
    ## 10861     30820
    ## 10862      3779
    ## 10863      3937
    ## 10864      1794
    ## 10865      <NA>
    ## 10866     54908
    ## 10867      6586
    ## 10868     79646
    ## 10869    345630
    ## 10870      <NA>
    ## 10871     23286
    ## 10872     57451
    ## 10873     27430
    ## 10874      3161
    ## 10875    134492
    ## 10876       900
    ## 10877      2566
    ## 10878      2554
    ## 10879      2561
    ## 10880     23120
    ## 10881      9232
    ## 10882     10569
    ## 10883    114898
    ## 10884     79616
    ## 10885      <NA>
    ## 10886    114825
    ## 10887      7265
    ## 10888       147
    ## 10889    134510
    ## 10890    153830
    ## 10891      1879
    ## 10892      9685
    ## 10893    134353
    ## 10894     54974
    ## 10895     11063
    ## 10896      <NA>
    ## 10897      8728
    ## 10898    348938
    ## 10899     26999
    ## 10900    408263
    ## 10901      3702
    ## 10902      9443
    ## 10903     84868
    ## 10904      <NA>
    ## 10905      6444
    ## 10906     10399
    ## 10907     90933
    ## 10908      <NA>
    ## 10909     81786
    ## 10910      <NA>
    ## 10911      <NA>
    ## 10912      <NA>
    ## 10913      <NA>
    ## 10914    643836
    ## 10915      4245
    ## 10916      2324
    ## 10917     92304
    ## 10918     57472
    ## 10919      9945
    ## 10920      5601
    ## 10921    255426
    ## 10922     55819
    ## 10923     23061
    ## 10924     51149
    ## 10925      8878
    ## 10926     11282
    ## 10927      4056
    ## 10928      <NA>
    ## 10929      9794
    ## 10930       821
    ## 10931    646019
    ## 10932      3187
    ## 10933     80230
    ## 10934      9509
    ## 10935      <NA>
    ## 10936      <NA>
    ## 10937      <NA>
    ## 10938     80108
    ## 10939      <NA>
    ## 10940      <NA>
    ## 10941     57396
    ## 10942     91522
    ## 10943     85007
    ## 10944      3182
    ## 10945     55651
    ## 10946     64777
    ## 10947     23138
    ## 10948      <NA>
    ## 10949      <NA>
    ## 10950     10802
    ## 10951     51128
    ## 10952     23338
    ## 10953      <NA>
    ## 10954     91368
    ## 10955      7320
    ## 10956     51265
    ## 10957      5515
    ## 10958      <NA>
    ## 10959      6932
    ## 10960      <NA>
    ## 10961      7416
    ## 10962      <NA>
    ## 10963     23105
    ## 10964      3308
    ## 10965     54819
    ## 10966     27125
    ## 10967     27089
    ## 10968      2661
    ## 10969    134549
    ## 10970    134548
    ## 10971      <NA>
    ## 10972      <NA>
    ## 10973     11127
    ## 10974     10111
    ## 10975      3659
    ## 10976      <NA>
    ## 10977      6584
    ## 10978      <NA>
    ## 10979      6583
    ## 10980      8572
    ## 10981      8974
    ## 10982     23305
    ## 10983     96459
    ## 10984     51735
    ## 10985     56990
    ## 10986     90624
    ## 10987      3094
    ## 10988      2878
    ## 10989     10318
    ## 10990       309
    ## 10991     26112
    ## 10992      2760
    ## 10993    206358
    ## 10994      6678
    ## 10995       475
    ## 10996     10146
    ## 10997      2890
    ## 10998     10827
    ## 10999      4238
    ## 11000     55568
    ## 11001     79685
    ## 11002     23367
    ## 11003      9337
    ## 11004     25929
    ## 11005     29093
    ## 11006      <NA>
    ## 11007      <NA>
    ## 11008      <NA>
    ## 11009      <NA>
    ## 11010     80851
    ## 11011      <NA>
    ## 11012      <NA>
    ## 11013      <NA>
    ## 11014      <NA>
    ## 11015      <NA>
    ## 11016    149603
    ## 11017      <NA>
    ## 11018      <NA>
    ## 11019     51127
    ## 11020     81559
    ## 11021     84033
    ## 11022    200205
    ## 11023     57165
    ## 11024      2987
    ## 11025    128308
    ## 11026      <NA>
    ## 11027       375
    ## 11028      7483
    ## 11029    116841
    ## 11030     65094
    ## 11031      <NA>
    ## 11032      <NA>
    ## 11033    114548
    ## 11034     23164
    ## 11035    201163
    ## 11036      <NA>
    ## 11037      8533
    ## 11038     56953
    ## 11039      <NA>
    ## 11040     55090
    ## 11041     51655
    ## 11042     10400
    ## 11043     10743
    ## 11044      6720
    ## 11045    146691
    ## 11046     83450
    ## 11047     91647
    ## 11048     79018
    ## 11049      1819
    ## 11050      <NA>
    ## 11051     54890
    ## 11052      3996
    ## 11053      2314
    ## 11054    125170
    ## 11055      7156
    ## 11056    140775
    ## 11057      6470
    ## 11058     25979
    ## 11059      8834
    ## 11060      <NA>
    ## 11061    256302
    ## 11062      <NA>
    ## 11063      5606
    ## 11064      3768
    ## 11065     23495
    ## 11066     23326
    ## 11067       218
    ## 11068       224
    ## 11069     55244
    ## 11070      <NA>
    ## 11071      7732
    ## 11072      4239
    ## 11073      5598
    ## 11074      <NA>
    ## 11075     27077
    ## 11076     22905
    ## 11077     10750
    ## 11078      5636
    ## 11079      9706
    ## 11080      <NA>
    ## 11081     11216
    ## 11082     92521
    ## 11083       136
    ## 11084    125150
    ## 11085     54902
    ## 11086      9611
    ## 11087      9487
    ## 11088    201161
    ## 11089      7314
    ## 11090     51393
    ## 11091    388341
    ## 11092      <NA>
    ## 11093      <NA>
    ## 11094      <NA>
    ## 11095     10626
    ## 11096     10517
    ## 11097     51030
    ## 11098      5376
    ## 11099      9953
    ## 11100      1352
    ## 11101      <NA>
    ## 11102      9955
    ## 11103     60528
    ## 11104      9912
    ## 11105      6416
    ## 11106      <NA>
    ## 11107      1770
    ## 11108    388336
    ## 11109    644139
    ## 11110    388335
    ## 11111     56985
    ## 11112      6341
    ## 11113      4621
    ## 11114      4620
    ## 11115      4619
    ## 11116 100128560
    ## 11117      4622
    ## 11118      4626
    ## 11119      8735
    ## 11120      8522
    ## 11121      9340
    ## 11122      5957
    ## 11123    124739
    ## 11124    146845
    ## 11125      9482
    ## 11126      9423
    ## 11127     23533
    ## 11128    146850
    ## 11129    146849
    ## 11130      4628
    ## 11131     81565
    ## 11132      6154
    ## 11133     22899
    ## 11134    399512
    ## 11135     29098
    ## 11136      5198
    ## 11137     80169
    ## 11138      <NA>
    ## 11139     54785
    ## 11140     84314
    ## 11141      6844
    ## 11142      5187
    ## 11143      <NA>
    ## 11144     84667
    ## 11145     59344
    ## 11146       242
    ## 11147      <NA>
    ## 11148      <NA>
    ## 11149    116840
    ## 11150     58485
    ## 11151      9196
    ## 11152      <NA>
    ## 11153      1107
    ## 11154    124637
    ## 11155     84316
    ## 11156     92162
    ## 11157     23135
    ## 11158    146754
    ## 11159      <NA>
    ## 11160      1949
    ## 11161     55135
    ## 11162      <NA>
    ## 11163       482
    ## 11164    112483
    ## 11165      9513
    ## 11166      6665
    ## 11167      9526
    ## 11168       968
    ## 11169      1973
    ## 11170     26168
    ## 11171      <NA>
    ## 11172      <NA>
    ## 11173      8742
    ## 11174      5430
    ## 11175     57659
    ## 11176      1140
    ## 11177      2256
    ## 11178     57555
    ## 11179    254863
    ## 11180     57048
    ## 11181    147040
    ## 11182      <NA>
    ## 11183     84461
    ## 11184      2874
    ## 11185      1984
    ## 11186     51087
    ## 11187      6517
    ## 11188     23587
    ## 11189     23399
    ## 11190     11337
    ## 11191     79142
    ## 11192      1856
    ## 11193        37
    ## 11194      1742
    ## 11195       432
    ## 11196      <NA>
    ## 11197     10462
    ## 11198    162515
    ## 11199    201232
    ## 11200    255877
    ## 11201      <NA>
    ## 11202      <NA>
    ## 11203    440400
    ## 11204       239
    ## 11205     27043
    ## 11206       409
    ## 11207    400569
    ## 11208     58191
    ## 11209     84225
    ## 11210      9032
    ## 11211      5694
    ## 11212      5338
    ## 11213     50488
    ## 11214      8402
    ## 11215     26001
    ## 11216      5216
    ## 11217      2027
    ## 11218      9552
    ## 11219     23125
    ## 11220    388324
    ## 11221     10749
    ## 11222    124961
    ## 11223      9135
    ## 11224      4927
    ## 11225     84268
    ## 11226       708
    ## 11227     56919
    ## 11228     51009
    ## 11229     79003
    ## 11230      <NA>
    ## 11231      <NA>
    ## 11232     23302
    ## 11233     54478
    ## 11234     83394
    ## 11235      <NA>
    ## 11236     84817
    ## 11237     51003
    ## 11238    284111
    ## 11239     54739
    ## 11240     83659
    ## 11241    342527
    ## 11242     10514
    ## 11243    124976
    ## 11244    201305
    ## 11245      7326
    ## 11246     51479
    ## 11247    124936
    ## 11248     23140
    ## 11249       489
    ## 11250      5023
    ## 11251      <NA>
    ## 11252     84254
    ## 11253     55421
    ## 11254      3682
    ## 11255     83903
    ## 11256      5026
    ## 11257     83460
    ## 11258     30851
    ## 11259      1497
    ## 11260     23729
    ## 11261    162514
    ## 11262       443
    ## 11263     23108
    ## 11264 101928991
    ## 11265     23277
    ## 11266      5048
    ## 11267     79066
    ## 11268      4335
    ## 11269      9905
    ## 11270     55720
    ## 11271     63826
    ## 11272     23293
    ## 11273      3090
    ## 11274    124641
    ## 11275    146760
    ## 11276      6117
    ## 11277    114826
    ## 11278      5176
    ## 11279    124997
    ## 11280    727910
    ## 11281     10594
    ## 11282     83547
    ## 11283      8578
    ## 11284      <NA>
    ## 11285    124935
    ## 11286      5306
    ## 11287     51763
    ## 11288      4641
    ## 11289      1398
    ## 11290      7531
    ## 11291      8447
    ## 11292      9501
    ## 11293      <NA>
    ## 11294    359845
    ## 11295     55275
    ## 11296     51031
    ## 11297      <NA>
    ## 11298     50628
    ## 11299      <NA>
    ## 11300     55178
    ## 11301     64359
    ## 11302     29928
    ## 11303        29
    ## 11304    727857
    ## 11305      <NA>
    ## 11306      9527
    ## 11307      1362
    ## 11308       642
    ## 11309     84081
    ## 11310    374786
    ## 11311     85464
    ## 11312     84940
    ## 11313    124930
    ## 11314      <NA>
    ## 11315     28964
    ## 11316      <NA>
    ## 11317    116236
    ## 11318     57551
    ## 11319     57532
    ## 11320    399687
    ## 11321     51268
    ## 11322    124925
    ## 11323     57649
    ## 11324    147015
    ## 11325      2319
    ## 11326     26284
    ## 11327     55731
    ## 11328      9618
    ## 11329    284086
    ## 11330    116238
    ## 11331      6147
    ## 11332     83871
    ## 11333    147011
    ## 11334      <NA>
    ## 11335      6388
    ## 11336      <NA>
    ## 11337      <NA>
    ## 11338     10615
    ## 11339       230
    ## 11340     94005
    ## 11341      9094
    ## 11342    113235
    ## 11343     23098
    ## 11344      7448
    ## 11345    645832
    ## 11346    147007
    ## 11347     26073
    ## 11348      7126
    ## 11349     90410
    ## 11350     27346
    ## 11351     51701
    ## 11352      <NA>
    ## 11353    201229
    ## 11354      3965
    ## 11355      8844
    ## 11356     26118
    ## 11357      4763
    ## 11358      <NA>
    ## 11359      4974
    ## 11360      <NA>
    ## 11361      2123
    ## 11362     84440
    ## 11363     55813
    ## 11364     23512
    ## 11365     51379
    ## 11366     79915
    ## 11367     79736
    ## 11368     55803
    ## 11369     84282
    ## 11370     55288
    ## 11371    162494
    ## 11372      <NA>
    ## 11373      <NA>
    ## 11374      5717
    ## 11375      8851
    ## 11376      4642
    ## 11377     26022
    ## 11378        40
    ## 11379      <NA>
    ## 11380    124842
    ## 11381      <NA>
    ## 11382      <NA>
    ## 11383      3980
    ## 11384    117584
    ## 11385      5892
    ## 11386     54475
    ## 11387    146862
    ## 11388    162394
    ## 11389      <NA>
    ## 11390      <NA>
    ## 11391      <NA>
    ## 11392      5193
    ## 11393       163
    ## 11394     91608
    ## 11395      <NA>
    ## 11396     79148
    ## 11397      8148
    ## 11398      6352
    ## 11399      <NA>
    ## 11400      <NA>
    ## 11401      6348
    ## 11402      6351
    ## 11403      <NA>
    ## 11404      <NA>
    ## 11405     63897
    ## 11406     11056
    ## 11407     11276
    ## 11408     11072
    ## 11409      6871
    ## 11410        31
    ## 11411     26574
    ## 11412     79922
    ## 11413     79154
    ## 11414     79893
    ## 11415    284098
    ## 11416     80179
    ## 11417      9326
    ## 11418      <NA>
    ## 11419     84669
    ## 11420      <NA>
    ## 11421     10513
    ## 11422      8493
    ## 11423     54828
    ## 11424      6909
    ## 11425      <NA>
    ## 11426     57508
    ## 11427      9969
    ## 11428     51136
    ## 11429      6198
    ## 11430     51174
    ## 11431     81671
    ## 11432     51651
    ## 11433      1213
    ## 11434     79665
    ## 11435    388403
    ## 11436    284161
    ## 11437     55181
    ## 11438    348235
    ## 11439      4591
    ## 11440     22843
    ## 11441      5889
    ## 11442      <NA>
    ## 11443      9110
    ## 11444    124535
    ## 11445     54894
    ## 11446      <NA>
    ## 11447      9256
    ## 11448      4025
    ## 11449     54903
    ## 11450      8288
    ## 11451      <NA>
    ## 11452    140735
    ## 11453      <NA>
    ## 11454      6426
    ## 11455      <NA>
    ## 11456      7716
    ## 11457    404093
    ## 11458     51649
    ## 11459    124540
    ## 11460      <NA>
    ## 11461      8165
    ## 11462     59342
    ## 11463      8161
    ## 11464      7706
    ## 11465      8526
    ## 11466      <NA>
    ## 11467      <NA>
    ## 11468      9241
    ## 11469     58488
    ## 11470     55273
    ## 11471      <NA>
    ## 11472     23531
    ## 11473      3131
    ## 11474    252983
    ## 11475      1353
    ## 11476     10040
    ## 11477      <NA>
    ## 11478     51096
    ## 11479     54799
    ## 11480      4831
    ## 11481      4830
    ## 11482      9043
    ## 11483     10140
    ## 11484    124857
    ## 11485     51747
    ## 11486     91369
    ## 11487      8714
    ## 11488      8913
    ## 11489     55040
    ## 11490     84073
    ## 11491     55316
    ## 11492     80221
    ## 11493     55379
    ## 11494     51264
    ## 11495     64132
    ## 11496      1277
    ## 11497      6442
    ## 11498     84687
    ## 11499    201191
    ## 11500      5164
    ## 11501      3675
    ## 11502     11143
    ## 11503     81558
    ## 11504     10237
    ## 11505      8405
    ## 11506     11248
    ## 11507      4804
    ## 11508      5245
    ## 11509      <NA>
    ## 11510      <NA>
    ## 11511      <NA>
    ## 11512     51225
    ## 11513      2793
    ## 11514     11267
    ## 11515     65264
    ## 11516      <NA>
    ## 11517      8631
    ## 11518     29916
    ## 11519     10951
    ## 11520      4779
    ## 11521     51226
    ## 11522     80279
    ## 11523      <NA>
    ## 11524     55163
    ## 11525      6668
    ## 11526     90507
    ## 11527     90506
    ## 11528    124995
    ## 11529    114881
    ## 11530     30009
    ## 11531      9755
    ## 11532      3837
    ## 11533      9520
    ## 11534     84311
    ## 11535     30837
    ## 11536     57636
    ## 11537     80725
    ## 11538 100170841
    ## 11539      4302
    ## 11540    284106
    ## 11541      7703
    ## 11542      5691
    ## 11543      8396
    ## 11544     54883
    ## 11545      9349
    ## 11546      3927
    ## 11547      <NA>
    ## 11548     57125
    ## 11549    390790
    ## 11550       782
    ## 11551      6143
    ## 11552    342667
    ## 11553     84961
    ## 11554      5469
    ## 11555     51755
    ## 11556      4761
    ## 11557     84152
    ## 11558     10948
    ## 11559      8557
    ## 11560      5409
    ## 11561     93210
    ## 11562      2064
    ## 11563     84299
    ## 11564      2886
    ## 11565     94103
    ## 11566      <NA>
    ## 11567      5709
    ## 11568      9862
    ## 11569      7067
    ## 11570      9572
    ## 11571    339287
    ## 11572      <NA>
    ## 11573     22794
    ## 11574     51195
    ## 11575    147179
    ## 11576       990
    ## 11577      5914
    ## 11578      7153
    ## 11579      3487
    ## 11580     84951
    ## 11581      1236
    ## 11582      6605
    ## 11583    125113
    ## 11584    162605
    ## 11585      3858
    ## 11586      3859
    ## 11587      <NA>
    ## 11588      3880
    ## 11589      3857
    ## 11590      3872
    ## 11591     10209
    ## 11592      9001
    ## 11593      3728
    ## 11594     10609
    ## 11595     60681
    ## 11596    115024
    ## 11597     55175
    ## 11598        47
    ## 11599     83538
    ## 11600      1267
    ## 11601      7266
    ## 11602     28511
    ## 11603      <NA>
    ## 11604     79132
    ## 11605      2648
    ## 11606      5878
    ## 11607     23415
    ## 11608      3060
    ## 11609     84514
    ## 11610      6777
    ## 11611      6776
    ## 11612      6774
    ## 11613    284119
    ## 11614       535
    ## 11615      4669
    ## 11616      3292
    ## 11617     80347
    ## 11618      6945
    ## 11619     29893
    ## 11620    162427
    ## 11621      7283
    ## 11622     27175
    ## 11623     79990
    ## 11624      8506
    ## 11625      2826
    ## 11626      2145
    ## 11627     10266
    ## 11628     84313
    ## 11629     65266
    ## 11630     28958
    ## 11631    124817
    ## 11632      8678
    ## 11633     10197
    ## 11634       314
    ## 11635      8639
    ## 11636      <NA>
    ## 11637     80755
    ## 11638 100885848
    ## 11639    146923
    ## 11640      6155
    ## 11641      3430
    ## 11642     10493
    ## 11643      8153
    ## 11644       672
    ## 11645      4077
    ## 11646    113277
    ## 11647    201299
    ## 11648       379
    ## 11649      1659
    ## 11650      2118
    ## 11651      4222
    ## 11652      1845
    ## 11653      4356
    ## 11654      4355
    ## 11655     84336
    ## 11656    124801
    ## 11657     92579
    ## 11658     10014
    ## 11659      <NA>
    ## 11660     92591
    ## 11661     79089
    ## 11662     56970
    ## 11663      7343
    ## 11664      6521
    ## 11665     10900
    ## 11666     51629
    ## 11667      2896
    ## 11668    284069
    ## 11669      3674
    ## 11670     23131
    ## 11671      <NA>
    ## 11672      2535
    ## 11673      <NA>
    ## 11674    124808
    ## 11675      4185
    ## 11676     10052
    ## 11677     51751
    ## 11678      9343
    ## 11679    388389
    ## 11680      2670
    ## 11681     10882
    ## 11682      <NA>
    ## 11683     79877
    ## 11684      4836
    ## 11685    113026
    ## 11686     79777
    ## 11687     10614
    ## 11688    124790
    ## 11689       752
    ## 11690      <NA>
    ## 11691    124783
    ## 11692      9020
    ## 11693    201176
    ## 11694      <NA>
    ## 11695      <NA>
    ## 11696      9842
    ## 11697      <NA>
    ## 11698    388394
    ## 11699      9570
    ## 11700      7473
    ## 11701      4905
    ## 11702      <NA>
    ## 11703      1394
    ## 11704      4137
    ## 11705    284058
    ## 11706       996
    ## 11707      4635
    ## 11708      3690
    ## 11709    146779
    ## 11710      <NA>
    ## 11711      <NA>
    ## 11712     11011
    ## 11713      9902
    ## 11714      <NA>
    ## 11715     26115
    ## 11716      1534
    ## 11717      1636
    ## 11718     81033
    ## 11719     10238
    ## 11720      <NA>
    ## 11721     51204
    ## 11722      4215
    ## 11723     80774
    ## 11724     92335
    ## 11725     57003
    ## 11726     11325
    ## 11727    117246
    ## 11728      5705
    ## 11729      6603
    ## 11730      <NA>
    ## 11731       974
    ## 11732      3384
    ## 11733      2081
    ## 11734     55852
    ## 11735      5175
    ## 11736     11232
    ## 11737      1655
    ## 11738     90799
    ## 11739     64750
    ## 11740      3838
    ## 11741      <NA>
    ## 11742      2186
    ## 11743     25926
    ## 11744     26207
    ## 11745      5718
    ## 11746      9931
    ## 11747     27092
    ## 11748     27091
    ## 11749      5578
    ## 11750       350
    ## 11751    201134
    ## 11752      8313
    ## 11753      8787
    ## 11754      <NA>
    ## 11755     10672
    ## 11756     51321
    ## 11757      9120
    ## 11758     22901
    ## 11759     55062
    ## 11760      5573
    ## 11761     54757
    ## 11762      <NA>
    ## 11763      <NA>
    ## 11764     10350
    ## 11765     23460
    ## 11766     23461
    ## 11767      5608
    ## 11768      3773
    ## 11769      3759
    ## 11770      6662
    ## 11771      <NA>
    ## 11772    201266
    ## 11773      6752
    ## 11774      9382
    ## 11775     84923
    ## 11776      <NA>
    ## 11777     23580
    ## 11778     54549
    ## 11779      6169
    ## 11780     94015
    ## 11781      <NA>
    ## 11782      <NA>
    ## 11783    388419
    ## 11784     55890
    ## 11785     11314
    ## 11786      <NA>
    ## 11787    326624
    ## 11788      9368
    ## 11789     26151
    ## 11790     54868
    ## 11791      2905
    ## 11792      2232
    ## 11793    283985
    ## 11794     92736
    ## 11795    124590
    ## 11796    283987
    ## 11797     30850
    ## 11798      3396
    ## 11799      <NA>
    ## 11800     23510
    ## 11801     79637
    ## 11802     30833
    ## 11803     51155
    ## 11804      6613
    ## 11805     79902
    ## 11806     23163
    ## 11807     51081
    ## 11808     57409
    ## 11809     60386
    ## 11810      2885
    ## 11811      9772
    ## 11812     57513
    ## 11813    283989
    ## 11814      3993
    ## 11815     80022
    ## 11816      9400
    ## 11817     29115
    ## 11818      3691
    ## 11819      2584
    ## 11820      <NA>
    ## 11821     85451
    ## 11822    201294
    ## 11823     23558
    ## 11824     91107
    ## 11825    201292
    ## 11826     64978
    ## 11827     85302
    ## 11828        51
    ## 11829 100134934
    ## 11830      2125
    ## 11831      6730
    ## 11832     23265
    ## 11833      2302
    ## 11834    114804
    ## 11835    283991
    ## 11836      5635
    ## 11837      8877
    ## 11838     63893
    ## 11839        15
    ## 11840     79651
    ## 11841    114757
    ## 11842    768206
    ## 11843      <NA>
    ## 11844     10610
    ## 11845    439921
    ## 11846     23210
    ## 11847    124512
    ## 11848      6427
    ## 11849     79157
    ## 11850    146664
    ## 11851    654434
    ## 11852      <NA>
    ## 11853      6397
    ## 11854      <NA>
    ## 11855      <NA>
    ## 11856     57690
    ## 11857     11322
    ## 11858    147138
    ## 11859      9144
    ## 11860      7083
    ## 11861    125061
    ## 11862       332
    ## 11863      <NA>
    ## 11864      9021
    ## 11865      9489
    ## 11866      <NA>
    ## 11867      9267
    ## 11868     57602
    ## 11869      7077
    ## 11870      3959
    ## 11871    124583
    ## 11872    114897
    ## 11873     64772
    ## 11874    146713
    ## 11875     84733
    ## 11876     57332
    ## 11877      8535
    ## 11878    125058
    ## 11879     55036
    ## 11880      2548
    ## 11881      9775
    ## 11882     79092
    ## 11883      6448
    ## 11884    284129
    ## 11885     57674
    ## 11886    284131
    ## 11887      4884
    ## 11888     57521
    ## 11889     79643
    ## 11890      <NA>
    ## 11891     10458
    ## 11892      9625
    ## 11893     22994
    ## 11894    146705
    ## 11895    284184
    ## 11896    124565
    ## 11897     57597
    ## 11898      <NA>
    ## 11899        71
    ## 11900     25794
    ## 11901     80233
    ## 11902     55666
    ## 11903      5148
    ## 11904    339229
    ## 11905    339230
    ## 11906    339231
    ## 11907      9146
    ## 11908      6182
    ## 11909      1468
    ## 11910    348262
    ## 11911      5034
    ## 11912       396
    ## 11913     10189
    ## 11914     51529
    ## 11915    256933
    ## 11916      5833
    ## 11917     51547
    ## 11918      4097
    ## 11919      <NA>
    ## 11920      5831
    ## 11921    255275
    ## 11922    147111
    ## 11923      <NA>
    ## 11924     79058
    ## 11925    201254
    ## 11926    201255
    ## 11927      5881
    ## 11928     51181
    ## 11929      <NA>
    ## 11930      <NA>
    ## 11931      5986
    ## 11932      2873
    ## 11933     64118
    ## 11934      2194
    ## 11935    284001
    ## 11936      9123
    ## 11937      1453
    ## 11938       924
    ## 11939      2837
    ## 11940     79701
    ## 11941      <NA>
    ## 11942      <NA>
    ## 11943     26502
    ## 11944      3607
    ## 11945     56270
    ## 11946     10966
    ## 11947     79672
    ## 11948     64122
    ## 11949      6904
    ## 11950    146712
    ## 11951    284207
    ## 11952      <NA>
    ## 11953      <NA>
    ## 11954      2665
    ## 11955      <NA>
    ## 11956     79754
    ## 11957     10276
    ## 11958      <NA>
    ## 11959      <NA>
    ## 11960      <NA>
    ## 11961      1316
    ## 11962     10531
    ## 11963      5214
    ## 11964       105
    ## 11965     22884
    ## 11966      3422
    ## 11967     23560
    ## 11968     23185
    ## 11969     22982
    ## 11970      <NA>
    ## 11971     10771
    ## 11972      1131
    ## 11973      6262
    ## 11974      4548
    ## 11975        88
    ## 11976     55127
    ## 11977      3964
    ## 11978      <NA>
    ## 11979      7107
    ## 11980      4811
    ## 11981      1130
    ## 11982      2786
    ## 11983    148789
    ## 11984      6905
    ## 11985      9453
    ## 11986     51742
    ## 11987     23072
    ## 11988      <NA>
    ## 11989     64983
    ## 11990      5683
    ## 11991      <NA>
    ## 11992      <NA>
    ## 11993      2737
    ## 11994      <NA>
    ## 11995      3624
    ## 11996     79783
    ## 11997    136647
    ## 11998      8621
    ## 11999      5898
    ## 12000      <NA>
    ## 12001     11281
    ## 12002     27072
    ## 12003       273
    ## 12004     83930
    ## 12005     54749
    ## 12006      6424
    ## 12007      9844
    ## 12008      <NA>
    ## 12009      5987
    ## 12010    257202
    ## 12011      9753
    ## 12012     80317
    ## 12013     84547
    ## 12014      7741
    ## 12015    222698
    ## 12016    387032
    ## 12017      <NA>
    ## 12018      <NA>
    ## 12019      7745
    ## 12020      <NA>
    ## 12021      <NA>
    ## 12022      <NA>
    ## 12023     10279
    ## 12024      <NA>
    ## 12025      <NA>
    ## 12026     29777
    ## 12027      <NA>
    ## 12028      <NA>
    ## 12029      <NA>
    ## 12030      <NA>
    ## 12031      <NA>
    ## 12032      <NA>
    ## 12033      <NA>
    ## 12034      <NA>
    ## 12035      3077
    ## 12036      <NA>
    ## 12037     10590
    ## 12038     55604
    ## 12039      9750
    ## 12040      <NA>
    ## 12041     51053
    ## 12042      <NA>
    ## 12043     55856
    ## 12044     51567
    ## 12045      <NA>
    ## 12046      7915
    ## 12047      2822
    ## 12048     57380
    ## 12049      <NA>
    ## 12050    140767
    ## 12051      <NA>
    ## 12052      <NA>
    ## 12053      <NA>
    ## 12054      6659
    ## 12055      <NA>
    ## 12056     54901
    ## 12057      1871
    ## 12058    154141
    ## 12059      <NA>
    ## 12060      7386
    ## 12061     56940
    ## 12062      3662
    ## 12063      <NA>
    ## 12064     55770
    ## 12065     94234
    ## 12066      2295
    ## 12067      2296
    ## 12068      2762
    ## 12069     56897
    ## 12070      <NA>
    ## 12071      <NA>
    ## 12072      5272
    ## 12073      <NA>
    ## 12074      <NA>
    ## 12075      <NA>
    ## 12076      <NA>
    ## 12077      4835
    ## 12078      8737
    ## 12079       670
    ## 12080      7280
    ## 12081    347733
    ## 12082    389362
    ## 12083     63027
    ## 12084    221749
    ## 12085      <NA>
    ## 12086      8899
    ## 12087    222826
    ## 12088     10455
    ## 12089      9425
    ## 12090     10799
    ## 12091    648791
    ## 12092     57128
    ## 12093     10667
    ## 12094     51299
    ## 12095      2162
    ## 12096      9450
    ## 12097      <NA>
    ## 12098      6239
    ## 12099      6745
    ## 12100    285782
    ## 12101     83732
    ## 12102    154007
    ## 12103       654
    ## 12104     81567
    ## 12105     63915
    ## 12106      9521
    ## 12107      <NA>
    ## 12108     51000
    ## 12109      <NA>
    ## 12110      2651
    ## 12111      <NA>
    ## 12112      <NA>
    ## 12113     55003
    ## 12114     51522
    ## 12115      4117
    ## 12116     54898
    ## 12117    221710
    ## 12118      4739
    ## 12119 100113407
    ## 12120     84830
    ## 12121      3096
    ## 12122      1906
    ## 12123    221692
    ## 12124     51256
    ## 12125     54438
    ## 12126     23408
    ## 12127     51406
    ## 12128     10048
    ## 12129      <NA>
    ## 12130     63933
    ## 12131    221687
    ## 12132      <NA>
    ## 12133      9308
    ## 12134      <NA>
    ## 12135      <NA>
    ## 12136      <NA>
    ## 12137      3720
    ## 12138     84062
    ## 12139      <NA>
    ## 12140      <NA>
    ## 12141     29116
    ## 12142      2766
    ## 12143      6310
    ## 12144      <NA>
    ## 12145      <NA>
    ## 12146    221662
    ## 12147     10486
    ## 12148     51439
    ## 12149      9972
    ## 12150     63971
    ## 12151    378884
    ## 12152      7172
    ## 12153    221656
    ## 12154      7913
    ## 12155    255488
    ## 12156      <NA>
    ## 12157      <NA>
    ## 12158      3400
    ## 12159      <NA>
    ## 12160      <NA>
    ## 12161      <NA>
    ## 12162      <NA>
    ## 12163    138639
    ## 12164      5253
    ## 12165      <NA>
    ## 12166     23196
    ## 12167    158293
    ## 12168     65268
    ## 12169      4814
    ## 12170     84270
    ## 12171    203328
    ## 12172     89846
    ## 12173     23299
    ## 12174     64768
    ## 12175    401541
    ## 12176      <NA>
    ## 12177      1842
    ## 12178     54829
    ## 12179      4958
    ## 12180      4969
    ## 12181     55035
    ## 12182      <NA>
    ## 12183      <NA>
    ## 12184     10927
    ## 12185    158046
    ## 12186      <NA>
    ## 12187      1903
    ## 12188     53358
    ## 12189      1164
    ## 12190     79048
    ## 12191     10507
    ## 12192     10912
    ## 12193     54769
    ## 12194      6850
    ## 12195       549
    ## 12196      4783
    ## 12197      4920
    ## 12198     10558
    ## 12199      <NA>
    ## 12200      4488
    ## 12201      1812
    ## 12202     94081
    ## 12203      3274
    ## 12204      <NA>
    ## 12205     10814
    ## 12206     84321
    ## 12207    375484
    ## 12208      <NA>
    ## 12209    285598
    ## 12210     51491
    ## 12211    192286
    ## 12212      1212
    ## 12213     23197
    ## 12214     22838
    ## 12215    114787
    ## 12216      6620
    ## 12217     26262
    ## 12218     90249
    ## 12219      3101
    ## 12220     51720
    ## 12221      <NA>
    ## 12222     64324
    ## 12223     53917
    ## 12224     27166
    ## 12225     83463
    ## 12226     10960
    ## 12227     10636
    ## 12228    345456
    ## 12229      2161
    ## 12230      2870
    ## 12231     80758
    ## 12232      1627
    ## 12233      9260
    ## 12234     79930
    ## 12235     51428
    ## 12236     54540
    ## 12237     54732
    ## 12238     11285
    ## 12239      <NA>
    ## 12240      9879
    ## 12241      <NA>
    ## 12242     79770
    ## 12243     84105
    ## 12244      <NA>
    ## 12245    497189
    ## 12246      9547
    ## 12247      <NA>
    ## 12248      7045
    ## 12249      4090
    ## 12250      <NA>
    ## 12251     57113
    ## 12252      6695
    ## 12253     26249
    ## 12254     10949
    ## 12255    414328
    ## 12256     29979
    ## 12257     80318
    ## 12258     55582
    ## 12259      <NA>
    ## 12260      3190
    ## 12261     80010
    ## 12262      4915
    ## 12263      <NA>
    ## 12264      <NA>
    ## 12265      <NA>
    ## 12266      <NA>
    ## 12267     23287
    ## 12268     60560
    ## 12269     51280
    ## 12270     81689
    ## 12271      <NA>
    ## 12272      <NA>
    ## 12273      <NA>
    ## 12274      2619
    ## 12275      1612
    ## 12276      <NA>
    ## 12277      <NA>
    ## 12278      <NA>
    ## 12279      <NA>
    ## 12280      <NA>
    ## 12281      <NA>
    ## 12282      <NA>
    ## 12283      2203
    ## 12284      <NA>
    ## 12285      2176
    ## 12286      5727
    ## 12287    375748
    ## 12288     11046
    ## 12289      <NA>
    ## 12290     22927
    ## 12291      8555
    ## 12292      <NA>
    ## 12293      1514
    ## 12294     23552
    ## 12295      <NA>
    ## 12296     79937
    ## 12297     84641
    ## 12298      <NA>
    ## 12299      <NA>
    ## 12300      <NA>
    ## 12301      7381
    ## 12302      <NA>
    ## 12303     51001
    ## 12304      9791
    ## 12305      <NA>
    ## 12306      <NA>
    ## 12307      <NA>
    ## 12308      <NA>
    ## 12309      <NA>
    ## 12310      <NA>
    ## 12311      <NA>
    ## 12312      <NA>
    ## 12313      <NA>
    ## 12314      <NA>
    ## 12315      <NA>
    ## 12316      <NA>
    ## 12317      <NA>
    ## 12318      <NA>
    ## 12319      <NA>
    ## 12320      <NA>
    ## 12321      <NA>
    ## 12322      <NA>
    ## 12323      <NA>
    ## 12324      <NA>
    ## 12325      <NA>
    ## 12326      <NA>
    ## 12327      <NA>
    ## 12328      <NA>
    ## 12329      <NA>
    ## 12330      <NA>
    ## 12331      <NA>
    ## 12332      <NA>
    ## 12333      4552
    ## 12334     79072
    ## 12335      <NA>
    ## 12336       108
    ## 12337      <NA>
    ## 12338      <NA>
    ## 12339     54888
    ## 12340      6715
    ## 12341      <NA>
    ## 12342    134111
    ## 12343     84246
    ## 12344     23379
    ## 12345    170690
    ## 12346      <NA>
    ## 12347      4726
    ## 12348     64979
    ## 12349     79888
    ## 12350     81037
    ## 12351      7015
    ## 12352     10723
    ## 12353     85409
    ## 12354      9319
    ## 12355     65980
    ## 12356     11076
    ## 12357     55722
    ## 12358      6550
    ## 12359     11336
    ## 12360     57491
    ## 12361     10016
    ## 12362      6389
    ## 12363    133957
    ## 12364    389257
    ## 12365      <NA>
    ## 12366      <NA>
    ## 12367      <NA>
    ## 12368      <NA>
    ## 12369      <NA>
    ## 12370     51752
    ## 12371       831
    ## 12372      5122
    ## 12373     22936
    ## 12374      2745
    ## 12375     22836
    ## 12376    317671
    ## 12377    285601
    ## 12378    153642
    ## 12379      9652
    ## 12380     79772
    ## 12381     84250
    ## 12382      <NA>
    ## 12383     83989
    ## 12384      7025
    ## 12385      <NA>
    ## 12386      <NA>
    ## 12387     57561
    ## 12388      <NA>
    ## 12389     84059
    ## 12390    116068
    ## 12391     10622
    ## 12392    153364
    ## 12393      <NA>
    ## 12394      1070
    ## 12395      4208
    ## 12396      <NA>
    ## 12397      <NA>
    ## 12398      <NA>
    ## 12399      <NA>
    ## 12400    153396
    ## 12401       902
    ## 12402      5921
    ## 12403      1350
    ## 12404     10085
    ## 12405      1404
    ## 12406      1462
    ## 12407      7518
    ## 12408      <NA>
    ## 12409     92270
    ## 12410      6228
    ## 12411     83734
    ## 12412      <NA>
    ## 12413     23635
    ## 12414      <NA>
    ## 12415     84240
    ## 12416      5924
    ## 12417      4437
    ## 12418      1719
    ## 12419      <NA>
    ## 12420    340120
    ## 12421    167555
    ## 12422      9765
    ## 12423    256987
    ## 12424      7060
    ## 12425    345778
    ## 12426    202333
    ## 12427      <NA>
    ## 12428      <NA>
    ## 12429      <NA>
    ## 12430      9456
    ## 12431    133746
    ## 12432     23743
    ## 12433     29958
    ## 12434      <NA>
    ## 12435       411
    ## 12436     10184
    ## 12437      9522
    ## 12438      8546
    ## 12439      6902
    ## 12440     55255
    ## 12441      <NA>
    ## 12442      8622
    ## 12443     84327
    ## 12444     55109
    ## 12445      1393
    ## 12446      2149
    ## 12447     10788
    ## 12448      2151
    ## 12449     22987
    ## 12450      <NA>
    ## 12451      <NA>
    ## 12452    134359
    ## 12453     51426
    ## 12454      <NA>
    ## 12455      3156
    ## 12456    256006
    ## 12457     51301
    ## 12458     26049
    ## 12459      <NA>
    ## 12460     10412
    ## 12461     84340
    ## 12462      3074
    ## 12463      8507
    ## 12464      <NA>
    ## 12465      <NA>
    ## 12466      <NA>
    ## 12467     64283
    ## 12468     84135
    ## 12469     57763
    ## 12470       689
    ## 12471      2297
    ## 12472    134288
    ## 12473    134285
    ## 12474    115548
    ## 12475      <NA>
    ## 12476      3842
    ## 12477      <NA>
    ## 12478     79810
    ## 12479     23107
    ## 12480      4131
    ## 12481      9607
    ## 12482     64087
    ## 12483     55814
    ## 12484      <NA>
    ## 12485      6606
    ## 12486      <NA>
    ## 12487      <NA>
    ## 12488      2966
    ## 12489 100506658
    ## 12490    153562
    ## 12491      5884
    ## 12492 102157402
    ## 12493      6880
    ## 12494    202243
    ## 12495      1022
    ## 12496     92259
    ## 12497     64946
    ## 12498       891
    ## 12499     64924
    ## 12500      <NA>
    ## 12501      5295
    ## 12502      4064
    ## 12503    375449
    ## 12504    140890
    ## 12505     55914
    ## 12506     57486
    ## 12507     54557
    ## 12508     80006
    ## 12509       373
    ## 12510     23398
    ## 12511     64105
    ## 12512     11174
    ## 12513     10283
    ## 12514    285672
    ## 12515      <NA>
    ## 12516    401190
    ## 12517    285671
    ## 12518      3350
    ## 12519     51194
    ## 12520     27292
    ## 12521      3796
    ## 12522      <NA>
    ## 12523     91942
    ## 12524    643155
    ## 12525      1161
    ## 12526     79993
    ## 12527     55789
    ## 12528      5144
    ## 12529    115827
    ## 12530     10769
    ## 12531      <NA>
    ## 12532     65056
    ## 12533    166968
    ## 12534      4214
    ## 12535     79722
    ## 12536      3572
    ## 12537     54514
    ## 12538    153129
    ## 12539      8611
    ## 12540      <NA>
    ## 12541     54505
    ## 12542     10309
    ## 12543    493869
    ## 12544      3001
    ## 12545     11082
    ## 12546      <NA>
    ## 12547    112574
    ## 12548      8988
    ## 12549     54622
    ## 12550      4724
    ## 12551     10468
    ## 12552      4338
    ## 12553      3672
    ## 12554     53918
    ## 12555      3670
    ## 12556     79668
    ## 12557    133418
    ## 12558      <NA>
    ## 12559    348980
    ## 12560     10884
    ## 12561      <NA>
    ## 12562      2255
    ## 12563     23530
    ## 12564     10605
    ## 12565      <NA>
    ## 12566      <NA>
    ## 12567     64417
    ## 12568     56477
    ## 12569      3157
    ## 12570    167359
    ## 12571      <NA>
    ## 12572      <NA>
    ## 12573      <NA>
    ## 12574     10890
    ## 12575      3797
    ## 12576      <NA>
    ## 12577     55252
    ## 12578      1838
    ## 12579      1788
    ## 12580      5443
    ## 12581     22979
    ## 12582     51277
    ## 12583       109
    ## 12584     79172
    ## 12585    391356
    ## 12586      8648
    ## 12587      <NA>
    ## 12588      <NA>
    ## 12589     50618
    ## 12590      <NA>
    ## 12591    653140
    ## 12592    375190
    ## 12593    375189
    ## 12594     51639
    ## 12595      2281
    ## 12596     80304
    ## 12597    388931
    ## 12598    165324
    ## 12599     54454
    ## 12600    114818
    ## 12601      <NA>
    ## 12602      <NA>
    ## 12603     60526
    ## 12604      <NA>
    ## 12605     64342
    ## 12606       388
    ## 12607     23369
    ## 12608      6382
    ## 12609      9741
    ## 12610     57539
    ## 12611    130502
    ## 12612    130497
    ## 12613     57665
    ## 12614      3790
    ## 12615    348654
    ## 12616     79677
    ## 12617      7447
    ## 12618    729475
    ## 12619     81553
    ## 12620      4613
    ## 12621      1653
    ## 12622     51594
    ## 12623      <NA>
    ## 12624     28951
    ## 12625      <NA>
    ## 12626      <NA>
    ## 12627     23175
    ## 12628     23620
    ## 12629      9687
    ## 12630      1876
    ## 12631      <NA>
    ## 12632      9475
    ## 12633      <NA>
    ## 12634      <NA>
    ## 12635      3754
    ## 12636     10130
    ## 12637    245973
    ## 12638     79954
    ## 12639      4953
    ## 12640      3241
    ## 12641      <NA>
    ## 12642      <NA>
    ## 12643      <NA>
    ## 12644      8853
    ## 12645      <NA>
    ## 12646      9270
    ## 12647     51692
    ## 12648    285148
    ## 12649      6868
    ## 12650     10971
    ## 12651      <NA>
    ## 12652      <NA>
    ## 12653      <NA>
    ## 12654      9014
    ## 12655     29841
    ## 12656      8462
    ## 12657    192668
    ## 12658      6241
    ## 12659    129642
    ## 12660     57498
    ## 12661      3398
    ## 12662      <NA>
    ## 12663      9781
    ## 12664     91543
    ## 12665    129607
    ## 12666      <NA>
    ## 12667      6664
    ## 12668    728597
    ## 12669     55821
    ## 12670     78989
    ## 12671      6201
    ## 12672    246243
    ## 12673     55256
    ## 12674     51112
    ## 12675      7260
    ## 12676     23040
    ## 12677      7837
    ## 12678      7173
    ## 12679     54221
    ## 12680    129787
    ## 12681    285016
    ## 12682        52
    ## 12683     26751
    ## 12684    642273
    ## 12685      3912
    ## 12686      1738
    ## 12687     79872
    ## 12688      5172
    ## 12689     55973
    ## 12690     11062
    ## 12691     10466
    ## 12692      2845
    ## 12693     26959
    ## 12694      5577
    ## 12695      5294
    ## 12696    168455
    ## 12697     10135
    ## 12698      <NA>
    ## 12699      <NA>
    ## 12700    222256
    ## 12701    222255
    ## 12702 100130771
    ## 12703    221830
    ## 12704      7291
    ## 12705      9734
    ## 12706     23161
    ## 12707       196
    ## 12708     27075
    ## 12709     28969
    ## 12710     57037
    ## 12711     25928
    ## 12712      <NA>
    ## 12713    392636
    ## 12714      1607
    ## 12715      2115
    ## 12716     10124
    ## 12717    286006
    ## 12718      3475
    ## 12719      <NA>
    ## 12720      <NA>
    ## 12721      9732
    ## 12722     83943
    ## 12723     54674
    ## 12724      4189
    ## 12725      <NA>
    ## 12726     50640
    ## 12727      4897
    ## 12728     29091
    ## 12729      4857
    ## 12730      <NA>
    ## 12731      2290
    ## 12732      <NA>
    ## 12733      <NA>
    ## 12734      5587
    ## 12735     55632
    ## 12736     23256
    ## 12737      1690
    ## 12738     29966
    ## 12739     11154
    ## 12740     25831
    ## 12741     25938
    ## 12742    112487
    ## 12743     80224
    ## 12744       394
    ## 12745      9472
    ## 12746     64067
    ## 12747    112399
    ## 12748    171546
    ## 12749     55837
    ## 12750     58533
    ## 12751      1073
    ## 12752     11177
    ## 12753      <NA>
    ## 12754      <NA>
    ## 12755      <NA>
    ## 12756      <NA>
    ## 12757     55012
    ## 12758      <NA>
    ## 12759      5687
    ## 12760      4792
    ## 12761    253959
    ## 12762     84312
    ## 12763     51562
    ## 12764      <NA>
    ## 12765     89874
    ## 12766      <NA>
    ## 12767    145282
    ## 12768      6751
    ## 12769    161198
    ## 12770     10484
    ## 12771      8487
    ## 12772    122553
    ## 12773      5411
    ## 12774      4253
    ## 12775    254170
    ## 12776      <NA>
    ## 12777    145581
    ## 12778      <NA>
    ## 12779     54813
    ## 12780     23116
    ## 12781     55015
    ## 12782      2287
    ## 12783     57697
    ## 12784     55320
    ## 12785    161357
    ## 12786      6235
    ## 12787      6166
    ## 12788      4247
    ## 12789     55172
    ## 12790      5427
    ## 12791    122773
    ## 12792     23588
    ## 12793      9147
    ## 12794       382
    ## 12795     79609
    ## 12796      6655
    ## 12797     79944
    ## 12798      <NA>
    ## 12799      8814
    ## 12800     11183
    ## 12801     51062
    ## 12802     60485
    ## 12803     51199
    ## 12804      5836
    ## 12805    114088
    ## 12806     81542
    ## 12807      <NA>
    ## 12808    122786
    ## 12809      <NA>
    ## 12810     55860
    ## 12811      5684
    ## 12812      5926
    ## 12813     26520
    ## 12814      <NA>
    ## 12815     51339
    ## 12816      <NA>
    ## 12817     23002
    ## 12818     64582
    ## 12819    112849
    ## 12820     51528
    ## 12821      6252
    ## 12822    341883
    ## 12823     64430
    ## 12824     51635
    ## 12825      5494
    ## 12826      <NA>
    ## 12827     51804
    ## 12828      4331
    ## 12829     57570
    ## 12830    145389
    ## 12831      <NA>
    ## 12832      5583
    ## 12833      <NA>
    ## 12834      3091
    ## 12835      6617
    ## 12836     83851
    ## 12837      <NA>
    ## 12838      <NA>
    ## 12839     27133
    ## 12840     57381
    ## 12841      5529
    ## 12842    112840
    ## 12843     81537
    ## 12844     23224
    ## 12845      4522
    ## 12846      9495
    ## 12847      7597
    ## 12848     22890
    ## 12849      3306
    ## 12850    145376
    ## 12851      <NA>
    ## 12852     26030
    ## 12853      6710
    ## 12854     91612
    ## 12855      2342
    ## 12856    376267
    ## 12857      4149
    ## 12858      2530
    ## 12859     10243
    ## 12860     64398
    ## 12861     51382
    ## 12862      1965
    ## 12863    161145
    ## 12864     57475
    ## 12865      5283
    ## 12866       384
    ## 12867      <NA>
    ## 12868     10490
    ## 12869     51109
    ## 12870    145226
    ## 12871     23503
    ## 12872      5890
    ## 12873       677
    ## 12874      <NA>
    ## 12875      <NA>
    ## 12876        87
    ## 12877      8816
    ## 12878      <NA>
    ## 12879     55218
    ## 12880     57452
    ## 12881      2079
    ## 12882     55334
    ## 12883    400224
    ## 12884     56936
    ## 12885      9766
    ## 12886      6430
    ## 12887     64093
    ## 12888      6547
    ## 12889      <NA>
    ## 12890      <NA>
    ## 12891     51241
    ## 12892      <NA>
    ## 12893     55333
    ## 12894      8747
    ## 12895     10001
    ## 12896     23508
    ## 12897      <NA>
    ## 12898      4293
    ## 12899      <NA>
    ## 12900      <NA>
    ## 12901     26037
    ## 12902      9628
    ## 12903      8110
    ## 12904     26094
    ## 12905     53349
    ## 12906     58517
    ## 12907      5663
    ## 12908     89932
    ## 12909      8650
    ## 12910      <NA>
    ## 12911     79697
    ## 12912     10965
    ## 12913    641371
    ## 12914    122970
    ## 12915      <NA>
    ## 12916      <NA>
    ## 12917      <NA>
    ## 12918    641372
    ## 12919     83544
    ## 12920      9240
    ## 12921     91748
    ## 12922      <NA>
    ## 12923    145482
    ## 12924      <NA>
    ## 12925    145483
    ## 12926     51004
    ## 12927       957
    ## 12928     80127
    ## 12929      <NA>
    ## 12930      4329
    ## 12931     91750
    ## 12932    338917
    ## 12933      5826
    ## 12934    646658
    ## 12935     10577
    ## 12936    122961
    ## 12937      4053
    ## 12938      <NA>
    ## 12939      9870
    ## 12940     51077
    ## 12941     56252
    ## 12942      <NA>
    ## 12943      1743
    ## 12944     83694
    ## 12945      5228
    ## 12946      8892
    ## 12947     27030
    ## 12948        97
    ## 12949     79696
    ## 12950     91754
    ## 12951     10972
    ## 12952      2353
    ## 12953    122953
    ## 12954     55640
    ## 12955     11161
    ## 12956     23093
    ## 12957      7043
    ## 12958    112752
    ## 12959     55668
    ## 12960      2103
    ## 12961     22846
    ## 12962     23357
    ## 12963     64207
    ## 12964     85457
    ## 12965    283576
    ## 12966     57156
    ## 12967     58157
    ## 12968     29954
    ## 12969      2954
    ## 12970    283578
    ## 12971    161394
    ## 12972    122945
    ## 12973     63894
    ## 12974     10598
    ## 12975      9517
    ## 12976      8846
    ## 12977     81892
    ## 12978     22938
    ## 12979     57143
    ## 12980      9369
    ## 12981      1734
    ## 12982    145508
    ## 12983      7253
    ## 12984      2957
    ## 12985     85439
    ## 12986      6400
    ## 12987     23768
    ## 12988      <NA>
    ## 12989      2581
    ## 12990      <NA>
    ## 12991     54207
    ## 12992     55812
    ## 12993     11099
    ## 12994     79882
    ## 12995    161436
    ## 12996    123016
    ## 12997      1112
    ## 12998      <NA>
    ## 12999     90141
    ## 13000     55775
    ## 13001     56659
    ## 13002      5700
    ## 13003     55051
    ## 13004       801
    ## 13005      <NA>
    ## 13006      <NA>
    ## 13007    145567
    ## 13008      9252
    ## 13009     80017
    ## 13010      8111
    ## 13011    440193
    ## 13012      <NA>
    ## 13013     55671
    ## 13014      <NA>
    ## 13015     10516
    ## 13016      9321
    ## 13017      4287
    ## 13018     53981
    ## 13019      <NA>
    ## 13020    123041
    ## 13021     79890
    ## 13022      5641
    ## 13023      9950
    ## 13024      1113
    ## 13025      3705
    ## 13026      <NA>
    ## 13027     26175
    ## 13028     84520
    ## 13029     55148
    ## 13030     55727
    ## 13031     57578
    ## 13032    145270
    ## 13033     90050
    ## 13034      <NA>
    ## 13035     51676
    ## 13036     78990
    ## 13037      <NA>
    ## 13038     57062
    ## 13039      3429
    ## 13040      <NA>
    ## 13041     57718
    ## 13042    327657
    ## 13043      <NA>
    ## 13044     23405
    ## 13045     79789
    ## 13046      <NA>
    ## 13047    283596
    ## 13048     51218
    ## 13049 100507043
    ## 13050      <NA>
    ## 13051      <NA>
    ## 13052       624
    ## 13053     55102
    ## 13054     51527
    ## 13055     10914
    ## 13056      <NA>
    ## 13057      7443
    ## 13058      <NA>
    ## 13059      <NA>
    ## 13060     64919
    ## 13061     84193
    ## 13062      8812
    ## 13063    317762
    ## 13064     84439
    ## 13065     10858
    ## 13066      2009
    ## 13067     51466
    ## 13068    123099
    ## 13069      7528
    ## 13070    123096
    ## 13071    283600
    ## 13072      <NA>
    ## 13073     79446
    ## 13074     57596
    ## 13075      8788
    ## 13076     55384
    ## 13077      <NA>
    ## 13078      <NA>
    ## 13079      <NA>
    ## 13080      <NA>
    ## 13081      1735
    ## 13082      5527
    ## 13083      <NA>
    ## 13084      <NA>
    ## 13085      <NA>
    ## 13086      1778
    ## 13087      3320
    ## 13088     91833
    ## 13089      5891
    ## 13090      <NA>
    ## 13091     51550
    ## 13092      9895
    ## 13093    122416
    ## 13094     23186
    ## 13095      7187
    ## 13096     81693
    ## 13097      <NA>
    ## 13098      9578
    ## 13099 107984640
    ## 13100     91828
    ## 13101      7127
    ## 13102      <NA>
    ## 13103      <NA>
    ## 13104      1983
    ## 13105      4140
    ## 13106      1152
    ## 13107    115708
    ## 13108      9529
    ## 13109      <NA>
    ## 13110      3831
    ## 13111      7517
    ## 13112     79038
    ## 13113     23368
    ## 13114      <NA>
    ## 13115      <NA>
    ## 13116    122402
    ## 13117    647286
    ## 13118    374569
    ## 13119     26153
    ## 13120    388021
    ## 13121     64423
    ## 13122      <NA>
    ## 13123     10572
    ## 13124       207
    ## 13125 100128927
    ## 13126    283638
    ## 13127    122618
    ## 13128    113146
    ## 13129      <NA>
    ## 13130     55038
    ## 13131      3714
    ## 13132    256281
    ## 13133      2972
    ## 13134     90135
    ## 13135     23241
    ## 13136      9112
    ## 13137      1397
    ## 13138      1396
    ## 13139    283643
    ## 13140     80757
    ## 13141      <NA>
    ## 13142      <NA>
    ## 13143      <NA>
    ## 13144      3495
    ## 13145      3507
    ## 13146      <NA>
    ## 13147      7434
    ## 13148     55112
    ## 13149      <NA>
    ## 13150     57488
    ## 13151     54892
    ## 13152      5799
    ## 13153      9771
    ## 13154     55536
    ## 13155      8701
    ## 13156      <NA>
    ## 13157      6671
    ## 13158    221833
    ## 13159      3696
    ## 13160    346389
    ## 13161    256130
    ## 13162      6414
    ## 13163 100129792
    ## 13164      2690
    ## 13165     26272
    ## 13166      <NA>
    ## 13167      5019
    ## 13168    345557
    ## 13169     84674
    ## 13170      6167
    ## 13171      5562
    ## 13172     23548
    ## 13173      5734
    ## 13174      1601
    ## 13175      <NA>
    ## 13176    253260
    ## 13177      9180
    ## 13178      3977
    ## 13179    133584
    ## 13180      <NA>
    ## 13181     55100
    ## 13182      9631
    ## 13183      <NA>
    ## 13184     25836
    ## 13185      <NA>
    ## 13186      <NA>
    ## 13187      6507
    ## 13188    202151
    ## 13189    133686
    ## 13190      6502
    ## 13191     92255
    ## 13192    133690
    ## 13193      3575
    ## 13194      5618
    ## 13195    134218
    ## 13196     55299
    ## 13197      5810
    ## 13198     26064
    ## 13199     23600
    ## 13200     51289
    ## 13201      <NA>
    ## 13202      4883
    ## 13203     10923
    ## 13204     51663
    ## 13205     54545
    ## 13206     64083
    ## 13207     23037
    ## 13208      <NA>
    ## 13209     29102
    ## 13210      1004
    ## 13211      1007
    ## 13212      <NA>
    ## 13213      1008
    ## 13214      1010
    ## 13215      1016
    ## 13216      <NA>
    ## 13217     10409
    ## 13218      <NA>
    ## 13219      4651
    ## 13220     54463
    ## 13221      <NA>
    ## 13222      <NA>
    ## 13223     23194
    ## 13224      <NA>
    ## 13225      <NA>
    ## 13226     90268
    ## 13227      <NA>
    ## 13228      7204
    ## 13229      1767
    ## 13230      1501
    ## 13231      1611
    ## 13232    651746
    ## 13233     83853
    ## 13234      <NA>
    ## 13235    134147
    ## 13236     22948
    ## 13237      <NA>
    ## 13238 100505806
    ## 13239      9037
    ## 13240      <NA>
    ## 13241      6383
    ## 13242     10404
    ## 13243     85453
    ## 13244     92140
    ## 13245     55353
    ## 13246      4147
    ## 13247      6156
    ## 13248     10247
    ## 13249     10940
    ## 13250     79815
    ## 13251      3788
    ## 13252      6788
    ## 13253    157680
    ## 13254      1345
    ## 13255      5440
    ## 13256      6674
    ## 13257     25897
    ## 13258    157567
    ## 13259    169166
    ## 13260     26986
    ## 13261      7534
    ## 13262      <NA>
    ## 13263      <NA>
    ## 13264     83988
    ## 13265      <NA>
    ## 13266     50484
    ## 13267      <NA>
    ## 13268     51366
    ## 13269      7071
    ## 13270      <NA>
    ## 13271     51582
    ## 13272      <NA>
    ## 13273       528
    ## 13274     79870
    ## 13275      8323
    ## 13276    115908
    ## 13277     81034
    ## 13278     25879
    ## 13279      9699
    ## 13280      1807
    ## 13281      <NA>
    ## 13282     29967
    ## 13283      <NA>
    ## 13284     23414
    ## 13285     55074
    ## 13286    137735
    ## 13287       284
    ## 13288    340419
    ## 13289      3646
    ## 13290      9694
    ## 13291    157753
    ## 13292      7201
    ## 13293     84955
    ## 13294     56943
    ## 13295     93035
    ## 13296      9166
    ## 13297     55638
    ## 13298     27012
    ## 13299    114788
    ## 13300      7227
    ## 13301      8667
    ## 13302     84294
    ## 13303      5885
    ## 13304    441376
    ## 13305     90390
    ## 13306      2131
    ## 13307    401474
    ## 13308      4982
    ## 13309    114569
    ## 13310      <NA>
    ## 13311      5168
    ## 13312      6873
    ## 13313     64798
    ## 13314      7373
    ## 13315     28998
    ## 13316     27085
    ## 13317      6641
    ## 13318     22882
    ## 13319      <NA>
    ## 13320     79139
    ## 13321     93594
    ## 13322     84985
    ## 13323      <NA>
    ## 13324      <NA>
    ## 13325     11244
    ## 13326     29028
    ## 13327     55093
    ## 13328    114907
    ## 13329    157769
    ## 13330    654463
    ## 13331    157378
    ## 13332     55039
    ## 13333     11236
    ## 13334     83940
    ## 13335      4715
    ## 13336      9788
    ## 13337      6713
    ## 13338      9897
    ## 13339    286053
    ## 13340     10221
    ## 13341      <NA>
    ## 13342      4609
    ## 13343      5820
    ## 13344     51571
    ## 13345     50807
    ## 13346       114
    ## 13347     23167
    ## 13348      3786
    ## 13349     23639
    ## 13350    137835
    ## 13351     51105
    ## 13352      7038
    ## 13353      6503
    ## 13354      <NA>
    ## 13355     10397
    ## 13356      6482
    ## 13357     57623
    ## 13358     10656
    ## 13359     51059
    ## 13360    169044
    ## 13361     51305
    ## 13362     83696
    ## 13363     54108
    ## 13364     27161
    ## 13365      5747
    ## 13366     22898
    ## 13367     57210
    ## 13368      2843
    ## 13369     11156
    ## 13370    389690
    ## 13371      <NA>
    ## 13372       575
    ## 13373      <NA>
    ## 13374     23237
    ## 13375      8629
    ## 13376     51337
    ## 13377    137797
    ## 13378     66004
    ## 13379      8581
    ## 13380      <NA>
    ## 13381      4061
    ## 13382      <NA>
    ## 13383      <NA>
    ## 13384      <NA>
    ## 13385      4062
    ## 13386    286128
    ## 13387    116447
    ## 13388    114822
    ## 13389    389692
    ## 13390     23144
    ## 13391     79792
    ## 13392     93100
    ## 13393      1936
    ## 13394     84948
    ## 13395      <NA>
    ## 13396      7264
    ## 13397      <NA>
    ## 13398      <NA>
    ## 13399 100130274
    ## 13400    286077
    ## 13401      <NA>
    ## 13402     23513
    ## 13403     22827
    ## 13404    340371
    ## 13405      <NA>
    ## 13406      5339
    ## 13407     84875
    ## 13408      2907
    ## 13409    392275
    ## 13410     26873
    ## 13411     54512
    ## 13412      8733
    ## 13413      1537
    ## 13414     81858
    ## 13415     84232
    ## 13416      <NA>
    ## 13417     51236
    ## 13418    727957
    ## 13419     23246
    ## 13420    642658
    ## 13421      3297
    ## 13422      8694
    ## 13423     83482
    ## 13424     26233
    ## 13425     79581
    ## 13426    203054
    ## 13427     29894
    ## 13428     51160
    ## 13429      4796
    ## 13430     50626
    ## 13431     90990
    ## 13432     84988
    ## 13433      2875
    ## 13434    113655
    ## 13435      9401
    ## 13436      9684
    ## 13437    441381
    ## 13438      <NA>
    ## 13439     80728
    ## 13440      <NA>
    ## 13441      <NA>
    ## 13442     28991
    ## 13443      6132
    ## 13444      <NA>
    ## 13445      <NA>
    ## 13446     23543
    ## 13447      <NA>
    ## 13448      <NA>
    ## 13449      4627
    ## 13450     25828
    ## 13451     80020
    ## 13452      8664
    ## 13453     10369
    ## 13454      <NA>
    ## 13455     11020
    ## 13456      <NA>
    ## 13457      5816
    ## 13458      4689
    ## 13459      <NA>
    ## 13460      1439
    ## 13461      7263
    ## 13462      4357
    ## 13463     79734
    ## 13464    164656
    ## 13465      3560
    ## 13466    114904
    ## 13467      6753
    ## 13468      5880
    ## 13469      <NA>
    ## 13470     27128
    ## 13471    114794
    ## 13472      4242
    ## 13473     29775
    ## 13474     11135
    ## 13475      3957
    ## 13476     26088
    ## 13477     23616
    ## 13478      <NA>
    ## 13479     57026
    ## 13480      3956
    ## 13481     79159
    ## 13482     11078
    ## 13483      <NA>
    ## 13484      <NA>
    ## 13485      <NA>
    ## 13486     23464
    ## 13487    129138
    ## 13488     51386
    ## 13489     85377
    ## 13490      <NA>
    ## 13491      5435
    ## 13492      6663
    ## 13493      9463
    ## 13494      <NA>
    ## 13495     23539
    ## 13496     80115
    ## 13497      8398
    ## 13498     23764
    ## 13499     25829
    ## 13500      1454
    ## 13501      3761
    ## 13502     11015
    ## 13503     10521
    ## 13504     11144
    ## 13505    646851
    ## 13506     25776
    ## 13507     56993
    ## 13508      9929
    ## 13509      9567
    ## 13510     25777
    ## 13511     10126
    ## 13512     23467
    ## 13513     23466
    ## 13514      <NA>
    ## 13515     23492
    ## 13516      5155
    ## 13517      6122
    ## 13518      9145
    ## 13519     10454
    ## 13520      4248
    ## 13521     54471
    ## 13522       468
    ## 13523     91582
    ## 13524      8911
    ## 13525    113828
    ## 13526     23112
    ## 13527       158
    ## 13528     27352
    ## 13529      <NA>
    ## 13530      2847
    ## 13531     10478
    ## 13532      6767
    ## 13533     63929
    ## 13534      9978
    ## 13535      2033
    ## 13536     83746
    ## 13537    150356
    ## 13538      5905
    ## 13539     23264
    ## 13540      <NA>
    ## 13541      7008
    ## 13542     10766
    ## 13543      <NA>
    ## 13544     84844
    ## 13545        50
    ## 13546    171568
    ## 13547     27254
    ## 13548      5372
    ## 13549      2547
    ## 13550     27351
    ## 13551      4809
    ## 13552    150365
    ## 13553     79879
    ## 13554      6721
    ## 13555    440829
    ## 13556    115650
    ## 13557     79019
    ## 13558      <NA>
    ## 13559      <NA>
    ## 13560      4668
    ## 13561      <NA>
    ## 13562     91689
    ## 13563      4700
    ## 13564      <NA>
    ## 13565      6942
    ## 13566      <NA>
    ## 13567      <NA>
    ## 13568    150372
    ## 13569     94009
    ## 13570     27341
    ## 13571     84271
    ## 13572      1727
    ## 13573      <NA>
    ## 13574     53947
    ## 13575     26286
    ## 13576      <NA>
    ## 13577     11252
    ## 13578     25809
    ## 13579     27349
    ## 13580       706
    ## 13581     23170
    ## 13582     80274
    ## 13583       758
    ## 13584     64800
    ## 13585     25830
    ## 13586     25813
    ## 13587     29780
    ## 13588     64098
    ## 13589      <NA>
    ## 13590     84247
    ## 13591     55615
    ## 13592    112885
    ## 13593      <NA>
    ## 13594     10762
    ## 13595      <NA>
    ## 13596     55007
    ## 13597      2192
    ## 13598     25814
    ## 13599      7477
    ## 13600      <NA>
    ## 13601    150383
    ## 13602     10343
    ## 13603     55020
    ## 13604     51512
    ## 13605     55687
    ## 13606      9620
    ## 13607     23151
    ## 13608     64781
    ## 13609     25771
    ## 13610      <NA>
    ## 13611     23774
    ## 13612      9889
    ## 13613     79087
    ## 13614     79174
    ## 13615    415116
    ## 13616     23209
    ## 13617     56666
    ## 13618      <NA>
    ## 13619     80305
    ## 13620     83642
    ## 13621     85378
    ## 13622     83933
    ## 13623      6300
    ## 13624      5600
    ## 13625     23654
    ## 13626    414918
    ## 13627      9701
    ## 13628      6305
    ## 13629     91289
    ## 13630     29781
    ## 13631      9997
    ## 13632    440836
    ## 13633      1120
    ## 13634      <NA>
    ## 13635     23542
    ## 13636       410
    ## 13637      <NA>
    ## 13638     85358
    ## 13639      <NA>
    ## 13640        49
    ## 13641      <NA>
    ## 13642    341359
    ## 13643    144245
    ## 13644    144402
    ## 13645     55605
    ## 13646       225
    ## 13647    114134
    ## 13648    120892
    ## 13649      1272
    ## 13650     29951
    ## 13651    283464
    ## 13652     10138
    ## 13653     85437
    ## 13654     51535
    ## 13655    144165
    ## 13656      <NA>
    ## 13657     80070
    ## 13658     83448
    ## 13659     51135
    ## 13660      5756
    ## 13661     84216
    ## 13662      4753
    ## 13663    440097
    ## 13664      <NA>
    ## 13665    196527
    ## 13666      <NA>
    ## 13667    196528
    ## 13668      9169
    ## 13669     81539
    ## 13670     54407
    ## 13671     55089
    ## 13672    347902
    ## 13673     91523
    ## 13674     79657
    ## 13675      8909
    ## 13676     10411
    ## 13677      <NA>
    ## 13678     55652
    ## 13679     51564
    ## 13680      7421
    ## 13681     79022
    ## 13682      1280
    ## 13683     29843
    ## 13684      5213
    ## 13685    140461
    ## 13686    387856
    ## 13687      <NA>
    ## 13688      <NA>
    ## 13689     54934
    ## 13690       904
    ## 13691      <NA>
    ## 13692       112
    ## 13693       784
    ## 13694      9416
    ## 13695     27289
    ## 13696     85478
    ## 13697     51303
    ## 13698       377
    ## 13699      7480
    ## 13700      7471
    ## 13701     23109
    ## 13702      <NA>
    ## 13703      5571
    ## 13704      8085
    ## 13705      <NA>
    ## 13706    121268
    ## 13707     55716
    ## 13708     10376
    ## 13709      <NA>
    ## 13710      7846
    ## 13711     84790
    ## 13712      5630
    ## 13713     65244
    ## 13714     23416
    ## 13715     10445
    ## 13716     84070
    ## 13717     25766
    ## 13718     91010
    ## 13719      7009
    ## 13720     57701
    ## 13721    144233
    ## 13722     23017
    ## 13723       362
    ## 13724     29127
    ## 13725        41
    ## 13726      6602
    ## 13727      <NA>
    ## 13728      2819
    ## 13729     84987
    ## 13730     91012
    ## 13731     51474
    ## 13732    113251
    ## 13733      <NA>
    ## 13734     57609
    ## 13735      <NA>
    ## 13736       466
    ## 13737      <NA>
    ## 13738      4891
    ## 13739     25875
    ## 13740     81566
    ## 13741      7024
    ## 13742      5463
    ## 13743      <NA>
    ## 13744      9802
    ## 13745     57228
    ## 13746     51411
    ## 13747      1990
    ## 13748     11226
    ## 13749      9498
    ## 13750      <NA>
    ## 13751      6334
    ## 13752    401720
    ## 13753        94
    ## 13754        91
    ## 13755      <NA>
    ## 13756    160622
    ## 13757      3164
    ## 13758     60673
    ## 13759      <NA>
    ## 13760      <NA>
    ## 13761    144501
    ## 13762      3855
    ## 13763    319101
    ## 13764      3849
    ## 13765      3848
    ## 13766    374454
    ## 13767      1975
    ## 13768     23371
    ## 13769     84926
    ## 13770      3489
    ## 13771     51380
    ## 13772      <NA>
    ## 13773      3695
    ## 13774      5916
    ## 13775     84975
    ## 13776      9700
    ## 13777      5204
    ## 13778      <NA>
    ## 13779      8086
    ## 13780    121340
    ## 13781      6667
    ## 13782       269
    ## 13783     54458
    ## 13784      5094
    ## 13785      7786
    ## 13786      6895
    ## 13787      8620
    ## 13788     11016
    ## 13789      <NA>
    ## 13790      <NA>
    ## 13791     57658
    ## 13792      <NA>
    ## 13793     23583
    ## 13794     23468
    ## 13795      3178
    ## 13796     22818
    ## 13797     53831
    ## 13798      <NA>
    ## 13799      3678
    ## 13800    121355
    ## 13801      3071
    ## 13802      5153
    ## 13803      5502
    ## 13804      <NA>
    ## 13805      <NA>
    ## 13806      <NA>
    ## 13807     79903
    ## 13808      <NA>
    ## 13809     23059
    ## 13810    197358
    ## 13811     84464
    ## 13812      <NA>
    ## 13813     10131
    ## 13814      1387
    ## 13815       115
    ## 13816      6345
    ## 13817      7023
    ## 13818     84662
    ## 13819     51025
    ## 13820     79585
    ## 13821    114990
    ## 13822      9093
    ## 13823     57407
    ## 13824      3163
    ## 13825     29965
    ## 13826    124402
    ## 13827     23295
    ## 13828      <NA>
    ## 13829     84309
    ## 13830    124401
    ## 13831     79641
    ## 13832     84656
    ## 13833     29855
    ## 13834      9717
    ## 13835     51172
    ## 13836      <NA>
    ## 13837     56052
    ## 13838    196483
    ## 13839     54715
    ## 13840    283953
    ## 13841     79091
    ## 13842        18
    ## 13843     25880
    ## 13844      5373
    ## 13845     23589
    ## 13846      7874
    ## 13847      <NA>
    ## 13848      <NA>
    ## 13849      2903
    ## 13850      2013
    ## 13851    146279
    ## 13852      4682
    ## 13853    780776
    ## 13854     28955
    ## 13855     23274
    ## 13856      8651
    ## 13857    116028
    ## 13858      9516
    ## 13859      8303
    ## 13860     51061
    ## 13861     29066
    ## 13862      <NA>
    ## 13863     26156
    ## 13864      <NA>
    ## 13865      2935
    ## 13866       608
    ## 13867      <NA>
    ## 13868     55313
    ## 13869    729993
    ## 13870      2072
    ## 13871      <NA>
    ## 13872      5073
    ## 13873     51283
    ## 13874      <NA>
    ## 13875     54700
    ## 13876    123803
    ## 13877     23042
    ## 13878    255027
    ## 13879      <NA>
    ## 13880      <NA>
    ## 13881      9665
    ## 13882     54820
    ## 13883      4629
    ## 13884    123811
    ## 13885      4363
    ## 13886      6591
    ## 13887     79645
    ## 13888      7336
    ## 13889      4173
    ## 13890      5591
    ## 13891      <NA>
    ## 13892      1052
    ## 13893     23514
    ## 13894      5318
    ## 13895     51067
    ## 13896     10059
    ## 13897    121512
    ## 13898      <NA>
    ## 13899      8940
    ## 13900      9647
    ## 13901      5594
    ## 13902      <NA>
    ## 13903     29799
    ## 13904     23759
    ## 13905      <NA>
    ## 13906     23753
    ## 13907    164592
    ## 13908    150223
    ## 13909      7332
    ## 13910     85376
    ## 13911     23119
    ## 13912    645426
    ## 13913      5297
    ## 13914      3053
    ## 13915      9342
    ## 13916      1399
    ## 13917    150209
    ## 13918      8216
    ## 13919     80764
    ## 13920    400891
    ## 13921      9127
    ## 13922      6545
    ## 13923     55627
    ## 13924     90557
    ## 13925     51586
    ## 13926     84861
    ## 13927     91179
    ## 13928      <NA>
    ## 13929      <NA>
    ## 13930      9993
    ## 13931      <NA>
    ## 13932      <NA>
    ## 13933      6576
    ## 13934      8214
    ## 13935      5625
    ## 13936     65078
    ## 13937      <NA>
    ## 13938      <NA>
    ## 13939     29801
    ## 13940      5902
    ## 13941     27037
    ## 13942     54487
    ## 13943    128989
    ## 13944       421
    ## 13945      1312
    ## 13946     10587
    ## 13947     54584
    ## 13948      6899
    ## 13949      2812
    ## 13950      <NA>
    ## 13951      <NA>
    ## 13952      7122
    ## 13953      8318
    ## 13954      7353
    ## 13955      <NA>
    ## 13956     64976
    ## 13957      7290
    ## 13958      3538
    ## 13959     84002
    ## 13960      <NA>
    ## 13961     89857
    ## 13962     54800
    ## 13963     55689
    ## 13964     79929
    ## 13965     55486
    ## 13966     10057
    ## 13967      8893
    ## 13968      1857
    ## 13969      1173
    ## 13970     55324
    ## 13971     90113
    ## 13972     10195
    ## 13973 110599564
    ## 13974     94032
    ## 13975      9718
    ## 13976      5708
    ## 13977      1981
    ## 13978    131408
    ## 13979      1181
    ## 13980      5437
    ## 13981      7066
    ## 13982      8646
    ## 13983      <NA>
    ## 13984      2049
    ## 13985     23355
    ## 13986      <NA>
    ## 13987      1962
    ## 13988      <NA>
    ## 13989      9175
    ## 13990     90407
    ## 13991    200879
    ## 13992     59343
    ## 13993     10644
    ## 13994      6434
    ## 13995      2119
    ## 13996      1608
    ## 13997     55171
    ## 13998     51726
    ## 13999      1974
    ## 14000      5984
    ## 14001      6480
    ## 14002    132112
    ## 14003      5648
    ## 14004     64108
    ## 14005      6750
    ## 14006       604
    ## 14007      <NA>
    ## 14008      4026
    ## 14009      <NA>
    ## 14010     55214
    ## 14011      9076
    ## 14012      3556
    ## 14013    647309
    ## 14014    344901
    ## 14015    152137
    ## 14016      2257
    ## 14017    151963
    ## 14018      <NA>
    ## 14019    344905
    ## 14020     84239
    ## 14021      4976
    ## 14022      <NA>
    ## 14023      <NA>
    ## 14024      <NA>
    ## 14025      3280
    ## 14026      2814
    ## 14027     79572
    ## 14028     93109
    ## 14029     55341
    ## 14030    131583
    ## 14031    152002
    ## 14032     23527
    ## 14033      <NA>
    ## 14034      5504
    ## 14035       347
    ## 14036       622
    ## 14037      <NA>
    ## 14038      <NA>
    ## 14039      1739
    ## 14040      4241
    ## 14041     80235
    ## 14042      <NA>
    ## 14043     22916
    ## 14044    205564
    ## 14045      5062
    ## 14046     54965
    ## 14047     84984
    ## 14048    375387
    ## 14049    200933
    ## 14050    348793
    ## 14051    165918
    ## 14052     26043
    ## 14053    255758
    ## 14054      5130
    ## 14055      <NA>
    ## 14056      7037
    ## 14057     10188
    ## 14058      9711
    ## 14059     84248
    ## 14060     84859
    ## 14061     84223
    ## 14062      6165
    ## 14063     89782
    ## 14064    114885
    ## 14065      8723
    ## 14066      <NA>
    ## 14067     84561
    ## 14068     57493
    ## 14069      3693
    ## 14070      7372
    ## 14071      8997
    ## 14072     54763
    ## 14073     64770
    ## 14074      4638
    ## 14075      <NA>
    ## 14076    201562
    ## 14077       111
    ## 14078     26984
    ## 14079     10954
    ## 14080     54437
    ## 14081      <NA>
    ## 14082     79663
    ## 14083     54625
    ## 14084    151636
    ## 14085     83666
    ## 14086      <NA>
    ## 14087      3836
    ## 14088     54554
    ## 14089     26355
    ## 14090    131076
    ## 14091       942
    ## 14092      6565
    ## 14093     55840
    ## 14094      9657
    ## 14095      2804
    ## 14096      3059
    ## 14097     10721
    ## 14098      9515
    ## 14099      2960
    ## 14100    285282
    ## 14101      4710
    ## 14102     11167
    ## 14103    116064
    ## 14104    165829
    ## 14105      2932
    ## 14106     89876
    ## 14107     10063
    ## 14108     64091
    ## 14109     51365
    ## 14110       141
    ## 14111      <NA>
    ## 14112       941
    ## 14113     51300
    ## 14114     56983
    ## 14115     55254
    ## 14116     57514
    ## 14117      8702
    ## 14118      7348
    ## 14119    152404
    ## 14120      <NA>
    ## 14121      4045
    ## 14122      2596
    ## 14123     26137
    ## 14124      <NA>
    ## 14125     79691
    ## 14126     57577
    ## 14127    254887
    ## 14128     54762
    ## 14129       523
    ## 14130     80218
    ## 14131    205717
    ## 14132     54847
    ## 14133    152185
    ## 14134     55779
    ## 14135     91653
    ## 14136     25871
    ## 14137     29083
    ## 14138    131450
    ## 14139      <NA>
    ## 14140    151887
    ## 14141     55032
    ## 14142     64422
    ## 14143    151888
    ## 14144      4345
    ## 14145      <NA>
    ## 14146    344805
    ## 14147     29114
    ## 14148     55347
    ## 14149     90102
    ## 14150    257068
    ## 14151      <NA>
    ## 14152     25945
    ## 14153     50852
    ## 14154      <NA>
    ## 14155      9666
    ## 14156      <NA>
    ## 14157     55081
    ## 14158       961
    ## 14159     56987
    ## 14160    344595
    ## 14161       868
    ## 14162       214
    ## 14163     64332
    ## 14164     91775
    ## 14165     79598
    ## 14166      6152
    ## 14167      <NA>
    ## 14168     27107
    ## 14169     57092
    ## 14170     54931
    ## 14171     57337
    ## 14172     50939
    ## 14173     25890
    ## 14174     10342
    ## 14175     55076
    ## 14176      <NA>
    ## 14177     56954
    ## 14178     55773
    ## 14179      <NA>
    ## 14180     84319
    ## 14181     11259
    ## 14182      1295
    ## 14183      <NA>
    ## 14184    131566
    ## 14185     10402
    ## 14186      1371
    ## 14187     56650
    ## 14188     84864
    ## 14189    131544
    ## 14190     84100
    ## 14191    285220
    ## 14192     63899
    ## 14193    200894
    ## 14194      5627
    ## 14195      2042
    ## 14196      <NA>
    ## 14197      <NA>
    ## 14198      8545
    ## 14199      3355
    ## 14200     25978
    ## 14201    389136
    ## 14202    253559
    ## 14203      2632
    ## 14204      6091
    ## 14205      6092
    ## 14206     54033
    ## 14207      6782
    ## 14208     64092
    ## 14209      8204
    ## 14210     29761
    ## 14211      1525
    ## 14212     10950
    ## 14213      <NA>
    ## 14214    140578
    ## 14215      4685
    ## 14216     54148
    ## 14217     58494
    ## 14218      <NA>
    ## 14219      2551
    ## 14220       351
    ## 14221      <NA>
    ## 14222    116159
    ## 14223      9510
    ## 14224     11096
    ## 14225     29104
    ## 14226     26046
    ## 14227     10069
    ## 14228     10600
    ## 14229     10694
    ## 14230       571
    ## 14231      <NA>
    ## 14232      2897
    ## 14233      7074
    ## 14234      6647
    ## 14235     57466
    ## 14236     30811
    ## 14237      <NA>
    ## 14238     54069
    ## 14239      <NA>
    ## 14240     56246
    ## 14241      9875
    ## 14242     59271
    ## 14243      <NA>
    ## 14244      8867
    ## 14245     94104
    ## 14246      <NA>
    ## 14247     10215
    ## 14248    116448
    ## 14249      3455
    ## 14250      3588
    ## 14251      3454
    ## 14252      3460
    ## 14253       757
    ## 14254     54943
    ## 14255      2618
    ## 14256      6651
    ## 14257     29980
    ## 14258      9946
    ## 14259      6453
    ## 14260      <NA>
    ## 14261     64968
    ## 14262      6526
    ## 14263      9992
    ## 14264      <NA>
    ## 14265      1827
    ## 14266     54102
    ## 14267       861
    ## 14268     54093
    ## 14269       873
    ## 14270      <NA>
    ## 14271       874
    ## 14272      <NA>
    ## 14273     23515
    ## 14274      8208
    ## 14275     23562
    ## 14276      3141
    ## 14277     53820
    ## 14278     51227
    ## 14279      7267
    ## 14280      <NA>
    ## 14281      1859
    ## 14282      3763
    ## 14283      2114
    ## 14284      8624
    ## 14285     54014
    ## 14286      3150
    ## 14287      <NA>
    ## 14288    150082
    ## 14289      6450
    ## 14290     10317
    ## 14291    150084
    ## 14292      5121
    ## 14293      1826
    ## 14294     25825
    ## 14295      4599
    ## 14296     54101
    ## 14297     63977
    ## 14298     25966
    ## 14299     49854
    ## 14300      <NA>
    ## 14301      <NA>
    ## 14302      <NA>
    ## 14303     22828
    ## 14304      <NA>
    ## 14305     26230
    ## 14306     51106
    ## 14307     49861
    ## 14308      <NA>
    ## 14309     57492
    ## 14310    729515
    ## 14311     79683
    ## 14312     51429
    ## 14313      8871
    ## 14314     84947
    ## 14315    404672
    ## 14316     56995
    ## 14317      <NA>
    ## 14318      <NA>
    ## 14319      <NA>
    ## 14320      <NA>
    ## 14321      7430
    ## 14322      <NA>
    ## 14323      <NA>
    ## 14324      <NA>
    ## 14325      <NA>
    ## 14326      <NA>
    ## 14327      <NA>
    ## 14328      6196
    ## 14329     84624
    ## 14330      <NA>
    ## 14331    117289
    ## 14332      <NA>
    ## 14333      <NA>
    ## 14334     11116
    ## 14335     51660
    ## 14336    113402
    ## 14337    285800
    ## 14338      <NA>
    ## 14339     10846
    ## 14340      <NA>
    ## 14341      <NA>
    ## 14342      <NA>
    ## 14343    135138
    ## 14344      <NA>
    ## 14345      <NA>
    ## 14346      <NA>
    ## 14347     56895
    ## 14348      4216
    ## 14349      <NA>
    ## 14350      6581
    ## 14351      6582
    ## 14352      6580
    ## 14353      3482
    ## 14354 100271873
    ## 14355      4142
    ## 14356    154197
    ## 14357     29074
    ## 14358      6950
    ## 14359      <NA>
    ## 14360        39
    ## 14361      9589
    ## 14362      6648
    ## 14363      <NA>
    ## 14364      <NA>
    ## 14365      <NA>
    ## 14366      4301
    ## 14367    168002
    ## 14368     64094
    ## 14369      7058
    ## 14370    253769
    ## 14371      <NA>
    ## 14372     55274
    ## 14373      <NA>
    ## 14374      <NA>
    ## 14375      <NA>
    ## 14376      6991
    ## 14377     55780
    ## 14378     28514
    ## 14379      <NA>
    ## 14380     84498
    ## 14381      5689
    ## 14382      6908
    ## 14383      5134
    ## 14384     56979
    ## 14385      1105
    ## 14386    285704
    ## 14387      <NA>
    ## 14388      <NA>
    ## 14389      <NA>
    ## 14390     55781
    ## 14391      <NA>
    ## 14392    167410
    ## 14393      4012
    ## 14394    147650
    ## 14395      3036
    ## 14396      5518
    ## 14397      <NA>
    ## 14398      <NA>
    ## 14399      <NA>
    ## 14400      <NA>
    ## 14401      <NA>
    ## 14402      <NA>
    ## 14403      <NA>
    ## 14404      <NA>
    ## 14405      <NA>
    ## 14406      <NA>
    ## 14407      <NA>
    ## 14408      <NA>
    ## 14409      <NA>
    ## 14410      <NA>
    ## 14411      <NA>
    ## 14412      <NA>
    ## 14413      <NA>
    ## 14414      <NA>
    ## 14415      <NA>
    ## 14416      <NA>
    ## 14417      <NA>
    ## 14418      <NA>
    ## 14419      <NA>
    ## 14420      <NA>
    ## 14421     79228
    ## 14422     54985
    ## 14423     51330
    ## 14424      9080
    ## 14425      9088
    ## 14426    124222
    ## 14427     79412
    ## 14428      <NA>
    ## 14429     84256
    ## 14430    114984
    ## 14431     23524
    ## 14432      6923
    ## 14433      <NA>
    ## 14434     64063
    ## 14435     54442
    ## 14436      5170
    ## 14437     51005
    ## 14438      <NA>
    ## 14439       527
    ## 14440     57465
    ## 14441     80178
    ## 14442       899
    ## 14443      <NA>
    ## 14444        21
    ## 14445     10921
    ## 14446      1632
    ## 14447      1775
    ## 14448      1877
    ## 14449    283871
    ## 14450     64223
    ## 14451    283870
    ## 14452     57524
    ## 14453     84231
    ## 14454     25837
    ## 14455      <NA>
    ## 14456      5310
    ## 14457      7249
    ## 14458      4913
    ## 14459      9351
    ## 14460    283869
    ## 14461      <NA>
    ## 14462      9143
    ## 14463      2671
    ## 14464    124056
    ## 14465     10607
    ## 14466      6187
    ## 14467    735301
    ## 14468      4716
    ## 14469      6123
    ## 14470     51734
    ## 14471     64711
    ## 14472      3029
    ## 14473     81889
    ## 14474     10101
    ## 14475     90864
    ## 14476    197342
    ## 14477     23162
    ## 14478     65993
    ## 14479      4832
    ## 14480     90861
    ## 14481      <NA>
    ## 14482      9742
    ## 14483     79652
    ## 14484      9894
    ## 14485      1186
    ## 14486      <NA>
    ## 14487     64718
    ## 14488     84572
    ## 14489    115939
    ## 14490      8938
    ## 14491      7329
    ## 14492      8912
    ## 14493    150483
    ## 14494     30812
    ## 14495     64788
    ## 14496     51764
    ## 14497     63922
    ## 14498    113000
    ## 14499     10232
    ## 14500      <NA>
    ## 14501     84264
    ## 14502      <NA>
    ## 14503     79006
    ## 14504    146330
    ## 14505     84219
    ## 14506    339123
    ## 14507     10273
    ## 14508      9028
    ## 14509      <NA>
    ## 14510     89941
    ## 14511    197335
    ## 14512     84331
    ## 14513     84326
    ## 14514    117166
    ## 14515     57799
    ## 14516      9091
    ## 14517      <NA>
    ## 14518      6650
    ## 14519      9727
    ## 14520     26063
    ## 14521      4833
    ## 14522      <NA>
    ## 14523      <NA>
    ## 14524     10573
    ## 14525      8312
    ## 14526     64714
    ## 14527       398
    ## 14528      8786
    ## 14529     83986
    ## 14530     55692
    ## 14531     54492
    ## 14532      1843
    ## 14533     57222
    ## 14534      <NA>
    ## 14535    153222
    ## 14536       662
    ## 14537      <NA>
    ## 14538      <NA>
    ## 14539      5252
    ## 14540     51596
    ## 14541      8831
    ## 14542    449520
    ## 14543       578
    ## 14544      3710
    ## 14545     84300
    ## 14546    117283
    ## 14547    221496
    ## 14548      2914
    ## 14549      <NA>
    ## 14550      3159
    ## 14551      <NA>
    ## 14552     11165
    ## 14553      6204
    ## 14554      <NA>
    ## 14555     29993
    ## 14556     25803
    ## 14557      <NA>
    ## 14558      6631
    ## 14559     54887
    ## 14560      6882
    ## 14561      <NA>
    ## 14562    222663
    ## 14563      <NA>
    ## 14564     50619
    ## 14565      5467
    ## 14566      2178
    ## 14567      4736
    ## 14568      7005
    ## 14569      2289
    ## 14570      <NA>
    ## 14571    222662
    ## 14572      6732
    ## 14573    116369
    ## 14574      1432
    ## 14575      5603
    ## 14576     27154
    ## 14577    285848
    ## 14578      <NA>
    ## 14579    222658
    ## 14580     11329
    ## 14581      6428
    ## 14582      1026
    ## 14583     57699
    ## 14584     51645
    ## 14585      <NA>
    ## 14586    221476
    ## 14587     23787
    ## 14588    221472
    ## 14589      5292
    ## 14590     55633
    ## 14591      9025
    ## 14592     23070
    ## 14593    154467
    ## 14594    266727
    ## 14595     60685
    ## 14596    114781
    ## 14597      2739
    ## 14598      9619
    ## 14599     89765
    ## 14600     54020
    ## 14601      5152
    ## 14602     10785
    ## 14603      4731
    ## 14604      5316
    ## 14605       875
    ## 14606      7307
    ## 14607    150094
    ## 14608     23076
    ## 14609      <NA>
    ## 14610      4854
    ## 14611     79852
    ## 14612     23476
    ## 14613     10270
    ## 14614     26993
    ## 14615     58525
    ## 14616     64926
    ## 14617      <NA>
    ## 14618      <NA>
    ## 14619      <NA>
    ## 14620      <NA>
    ## 14621      <NA>
    ## 14622      <NA>
    ## 14623      <NA>
    ## 14624      <NA>
    ## 14625      <NA>
    ## 14626      <NA>
    ## 14627      <NA>
    ## 14628      <NA>
    ## 14629      <NA>
    ## 14630      <NA>
    ## 14631      <NA>
    ## 14632      <NA>
    ## 14633      <NA>
    ## 14634    284382
    ## 14635     81794
    ## 14636      4542
    ## 14637      <NA>
    ## 14638     84106
    ## 14639      4670
    ## 14640      <NA>
    ## 14641      9230
    ## 14642     51129
    ## 14643    256949
    ## 14644      6234
    ## 14645      4701
    ## 14646     51293
    ## 14647      3833
    ## 14648      <NA>
    ## 14649      <NA>
    ## 14650      1616
    ## 14651      6892
    ## 14652      9278
    ## 14653      <NA>
    ## 14654      5863
    ## 14655     10471
    ## 14656      9277
    ## 14657      8705
    ## 14658      6222
    ## 14659      6293
    ## 14660      <NA>
    ## 14661      6015
    ## 14662      <NA>
    ## 14663      7922
    ## 14664      6257
    ## 14665      1302
    ## 14666      <NA>
    ## 14667      6046
    ## 14668      <NA>
    ## 14669      <NA>
    ## 14670      <NA>
    ## 14671      5698
    ## 14672      6890
    ## 14673      5696
    ## 14674      6891
    ## 14675      <NA>
    ## 14676      <NA>
    ## 14677      <NA>
    ## 14678      <NA>
    ## 14679      <NA>
    ## 14680      4855
    ## 14681     63940
    ## 14682      5089
    ## 14683       177
    ## 14684      6048
    ## 14685     10554
    ## 14686     80864
    ## 14687      9374
    ## 14688     80863
    ## 14689     63943
    ## 14690      1388
    ## 14691      7148
    ## 14692       721
    ## 14693      8859
    ## 14694      1797
    ## 14695      6499
    ## 14696      7936
    ## 14697       717
    ## 14698    221527
    ## 14699     10919
    ## 14700     80736
    ## 14701      4758
    ## 14702      3304
    ## 14703      3303
    ## 14704      <NA>
    ## 14705     57819
    ## 14706      <NA>
    ## 14707     80737
    ## 14708    401251
    ## 14709      4439
    ## 14710      1192
    ## 14711     23564
    ## 14712     58530
    ## 14713     79136
    ## 14714    259215
    ## 14715      7920
    ## 14716      1460
    ## 14717      7918
    ## 14718      <NA>
    ## 14719      <NA>
    ## 14720      7917
    ## 14721      7916
    ## 14722      <NA>
    ## 14723       199
    ## 14724      7940
    ## 14725      4050
    ## 14726      4049
    ## 14727      4795
    ## 14728       534
    ## 14729      7919
    ## 14730      <NA>
    ## 14731      <NA>
    ## 14732      <NA>
    ## 14733      <NA>
    ## 14734      <NA>
    ## 14735      6941
    ## 14736     54535
    ## 14737    170680
    ## 14738      1041
    ## 14739     57176
    ## 14740      2968
    ## 14741       780
    ## 14742      8870
    ## 14743     10211
    ## 14744      <NA>
    ## 14745      9656
    ## 14746     11270
    ## 14747    170954
    ## 14748      8449
    ## 14749      <NA>
    ## 14750     79969
    ## 14751     28973
    ## 14752      5514
    ## 14753        23
    ## 14754     80742
    ## 14755      2794
    ## 14756      <NA>
    ## 14757      <NA>
    ## 14758      <NA>
    ## 14759     79897
    ## 14760     56658
    ## 14761      7726
    ## 14762     10107
    ## 14763      6992
    ## 14764     30834
    ## 14765      <NA>
    ## 14766      <NA>
    ## 14767      <NA>
    ## 14768    346171
    ## 14769      4340
    ## 14770      2550
    ## 14771      <NA>
    ## 14772      <NA>
    ## 14773      <NA>
    ## 14774      <NA>
    ## 14775     55166
    ## 14776      <NA>
    ## 14777    442213
    ## 14778    221393
    ## 14779     23607
    ## 14780     27242
    ## 14781    221395
    ## 14782 100287718
    ## 14783      7941
    ## 14784    221400
    ## 14785      9481
    ## 14786     51302
    ## 14787     10231
    ## 14788     59084
    ## 14789     22875
    ## 14790     53405
    ## 14791       860
    ## 14792      <NA>
    ## 14793      <NA>
    ## 14794       988
    ## 14795      <NA>
    ## 14796    221409
    ## 14797     57505
    ## 14798    202500
    ## 14799    441151
    ## 14800      4794
    ## 14801    347734
    ## 14802      3326
    ## 14803      2030
    ## 14804     11131
    ## 14805     55362
    ## 14806     64928
    ## 14807      <NA>
    ## 14808      <NA>
    ## 14809      7422
    ## 14810     55168
    ## 14811    221421
    ## 14812      9587
    ## 14813     54676
    ## 14814      5429
    ## 14815     57510
    ## 14816      9533
    ## 14817     25844
    ## 14818    221424
    ## 14819     93643
    ## 14820     65989
    ## 14821     89845
    ## 14822      <NA>
    ## 14823    401262
    ## 14824      <NA>
    ## 14825     84630
    ## 14826     10591
    ## 14827     23113
    ## 14828      6722
    ## 14829      5754
    ## 14830     89953
    ## 14831     51069
    ## 14832      9820
    ## 14833     88745
    ## 14834    116138
    ## 14835      4201
    ## 14836      5528
    ## 14837      5190
    ## 14838     27232
    ## 14839     10695
    ## 14840    171558
    ## 14841      <NA>
    ## 14842    285855
    ## 14843     23506
    ## 14844      6903
    ## 14845     23304
    ## 14846     55173
    ## 14847      2979
    ## 14848      2978
    ## 14849      <NA>
    ## 14850    129685
    ## 14851       896
    ## 14852       705
    ## 14853      9477
    ## 14854      <NA>
    ## 14855     25862
    ## 14856 100188893
    ## 14857      <NA>
    ## 14858     29964
    ## 14859     10817
    ## 14860      7942
    ## 14861      4188
    ## 14862    116113
    ## 14863     79865
    ## 14864     54209
    ## 14865    340205
    ## 14866      4800
    ## 14867    221443
    ## 14868    222642
    ## 14869     57497
    ## 14870      4337
    ## 14871     23500
    ## 14872    221458
    ## 14873     23180
    ## 14874      1618
    ## 14875     23228
    ## 14876      <NA>
    ## 14877      9779
    ## 14878      6304
    ## 14879    131096
    ## 14880      5868
    ## 14881      8850
    ## 14882    151648
    ## 14883     60482
    ## 14884     84620
    ## 14885      <NA>
    ## 14886      <NA>
    ## 14887      <NA>
    ## 14888      <NA>
    ## 14889      <NA>
    ## 14890     10148
    ## 14891      <NA>
    ## 14892     56961
    ## 14893     79187
    ## 14894     55620
    ## 14895     84954
    ## 14896      6455
    ## 14897     10036
    ## 14898     80700
    ## 14899     84717
    ## 14900    729359
    ## 14901    440503
    ## 14902    116844
    ## 14903     10501
    ## 14904    126282
    ## 14905     56005
    ## 14906     91039
    ## 14907     55527
    ## 14908    148022
    ## 14909     10226
    ## 14910     29128
    ## 14911     23030
    ## 14912      5802
    ## 14913      <NA>
    ## 14914      9667
    ## 14915      6294
    ## 14916      <NA>
    ## 14917     25873
    ## 14918      9361
    ## 14919    257062
    ## 14920      8498
    ## 14921    400673
    ## 14922    126328
    ## 14923      4902
    ## 14924     56931
    ## 14925    163154
    ## 14926      5990
    ## 14927      <NA>
    ## 14928      4298
    ## 14929      8192
    ## 14930     84266
    ## 14931      2962
    ## 14932      8570
    ## 14933     79085
    ## 14934     92359
    ## 14935     79958
    ## 14936     10382
    ## 14937      <NA>
    ## 14938      8744
    ## 14939       718
    ## 14940     56927
    ## 14941      9322
    ## 14942      7409
    ## 14943      2015
    ## 14944      <NA>
    ## 14945      <NA>
    ## 14946     83594
    ## 14947      1946
    ## 14948     64839
    ## 14949      2241
    ## 14950      9867
    ## 14951      4124
    ## 14952    642987
    ## 14953      9218
    ## 14954     11031
    ## 14955      9989
    ## 14956     10928
    ## 14957     57045
    ## 14958     23253
    ## 14959      4729
    ## 14960 100287171
    ## 14961      1663
    ## 14962     23255
    ## 14963    201475
    ## 14964      5797
    ## 14965    284217
    ## 14966    645369
    ## 14967     23136
    ## 14968      7541
    ## 14969      <NA>
    ## 14970    642597
    ## 14971      9229
    ## 14972      7050
    ## 14973      <NA>
    ## 14974    103910
    ## 14975     10627
    ## 14976      8736
    ## 14977      9663
    ## 14978     84034
    ## 14979     23347
    ## 14980     10403
    ## 14981    245711
    ## 14982     23160
    ## 14983    165186
    ## 14984     79745
    ## 14985     51646
    ## 14986     81606
    ## 14987    253558
    ## 14988     79623
    ## 14989     30845
    ## 14990      7498
    ## 14991      7795
    ## 14993     84661
    ## 14994      6683
    ## 14995     55676
    ## 14996     84272
    ## 14997     57448
    ## 14998     55622
    ## 14999      4052
    ## 15000     25780
    ## 15001     25940
    ## 15002     51232
    ## 15003      9637
    ## 15004      5212
    ## 15005      6801
    ## 15006     54497
    ## 15007    253635
    ## 15008      5610
    ## 15009 100505876
    ## 15010     10153
    ## 15011     55471
    ## 15012     23683
    ## 15013     25797
    ## 15014     10602
    ## 15015      <NA>
    ## 15016    151393
    ## 15017      1545
    ## 15018     64225
    ## 15019     92906
    ## 15020    130589
    ## 15021      6432
    ## 15022     79833
    ## 15023     90957
    ## 15024    729967
    ## 15025 100271715
    ## 15026      6654
    ## 15027    344387
    ## 15028      8491
    ## 15029      <NA>
    ## 15030     80745
    ## 15031      6546
    ## 15032     91461
    ## 15033     27436
    ## 15034      9167
    ## 15035     57504
    ## 15036     23498
    ## 15037       678
    ## 15038     63892
    ## 15039    130271
    ## 15040     51626
    ## 15041     10128
    ## 15042      <NA>
    ## 15043      5495
    ## 15044      6519
    ## 15045      9581
    ## 15046     79823
    ## 15047      6496
    ## 15048     55133
    ## 15049      5581
    ## 15050      2034
    ## 15051     23433
    ## 15052      5281
    ## 15053      9419
    ## 15054      9655
    ## 15055     90411
    ## 15056      <NA>
    ## 15057      <NA>
    ## 15058       805
    ## 15059      4072
    ## 15060      4436
    ## 15061     56660
    ## 15062      2956
    ## 15063     80204
    ## 15064      3344
    ## 15065    129285
    ## 15066      9378
    ## 15067       116
    ## 15068     64863
    ## 15069      <NA>
    ## 15070      <NA>
    ## 15071      <NA>
    ## 15072      <NA>
    ## 15073      8284
    ## 15074      <NA>
    ## 15075      <NA>
    ## 15076      7404
    ## 15077      8653
    ## 15078      <NA>
    ## 15079      1390
    ## 15080      <NA>
    ## 15081      8453
    ## 15082     25805
    ## 15083      1326
    ## 15084     55149
    ## 15085     57608
    ## 15086      6840
    ## 15087      <NA>
    ## 15088      6935
    ## 15089     94134
    ## 15090      3799
    ## 15091      <NA>
    ## 15092     80314
    ## 15093      <NA>
    ## 15094     22931
    ## 15095    283078
    ## 15096     55130
    ## 15097    143098
    ## 15098     51322
    ## 15099      8325
    ## 15100    219771
    ## 15101     81035
    ## 15102      9984
    ## 15103      9097
    ## 15104      6093
    ## 15105     80000
    ## 15106    114799
    ## 15107      6632
    ## 15108    171586
    ## 15109     57534
    ## 15110      5932
    ## 15111      <NA>
    ## 15112     91768
    ## 15113     85019
    ## 15114      8780
    ## 15115      <NA>
    ## 15116      4864
    ## 15117      <NA>
    ## 15118    147463
    ## 15119      3909
    ## 15120    125488
    ## 15121     26256
    ## 15122    114876
    ## 15123     55364
    ## 15124      <NA>
    ## 15125      6760
    ## 15126      6875
    ## 15127    284252
    ## 15128       361
    ## 15129      <NA>
    ## 15130      1000
    ## 15131      <NA>
    ## 15132      1824
    ## 15133      1829
    ## 15134      7276
    ## 15135      <NA>
    ## 15136      9331
    ## 15137     22878
    ## 15138     51444
    ## 15139     64762
    ## 15140     57565
    ## 15141      <NA>
    ## 15142     80816
    ## 15143      8715
    ## 15144      1837
    ## 15145      <NA>
    ## 15146     10982
    ## 15147      <NA>
    ## 15148      <NA>
    ## 15149      <NA>
    ## 15150    125476
    ## 15151      2589
    ## 15152      <NA>
    ## 15153     55197
    ## 15154     25800
    ## 15155     55250
    ## 15156     55034
    ## 15157     80206
    ## 15158     25941
    ## 15159      <NA>
    ## 15160     56853
    ## 15161      5289
    ## 15162      6014
    ## 15163      6860
    ## 15164     91137
    ## 15165      <NA>
    ## 15166     79595
    ## 15167     83607
    ## 15168      5433
    ## 15169     55339
    ## 15170     84826
    ## 15171     55679
    ## 15172      2840
    ## 15173     55677
    ## 15174     10746
    ## 15175      2071
    ## 15176      <NA>
    ## 15177       274
    ## 15178      2995
    ## 15179     85480
    ## 15180    134430
    ## 15181       814
    ## 15182    134429
    ## 15183      9315
    ## 15184      <NA>
    ## 15185     64097
    ## 15186       324
    ## 15187      6728
    ## 15188      7905
    ## 15189     51306
    ## 15190      8382
    ## 15191      <NA>
    ## 15192     10902
    ## 15193     10112
    ## 15194      8697
    ## 15195       995
    ## 15196     51307
    ## 15197     51780
    ## 15198     51308
    ## 15199      1958
    ## 15200      2107
    ## 15201      3313
    ## 15202      1495
    ## 15203     26045
    ## 15204     64374
    ## 15205      <NA>
    ## 15206      9782
    ## 15207     51247
    ## 15208     51237
    ## 15209    389333
    ## 15210    202051
    ## 15211    202052
    ## 15212    641700
    ## 15213      <NA>
    ## 15214      <NA>
    ## 15215     51523
    ## 15216     84249
    ## 15217      9542
    ## 15218      5813
    ## 15219    492311
    ## 15220     84418
    ## 15221      5201
    ## 15222      1839
    ## 15223     54882
    ## 15224     10011
    ## 15225    113829
    ## 15226       929
    ## 15227     55374
    ## 15228      4695
    ## 15229      3550
    ## 15230     54853
    ## 15231    373863
    ## 15232      <NA>
    ## 15233     23438
    ## 15234    153527
    ## 15235     56134
    ## 15236     56132
    ## 15237     56131
    ## 15238     26167
    ## 15239     56130
    ## 15240     56129
    ## 15241     56128
    ## 15242     56127
    ## 15243     56126
    ## 15244     56125
    ## 15245     56124
    ## 15246     56123
    ## 15247     56122
    ## 15248     56121
    ## 15249     57717
    ## 15250      <NA>
    ## 15251      <NA>
    ## 15252      <NA>
    ## 15253      <NA>
    ## 15254      <NA>
    ## 15255      <NA>
    ## 15256      6879
    ## 15257      <NA>
    ## 15258     56105
    ## 15259     26025
    ## 15260      <NA>
    ## 15261      1729
    ## 15262      8841
    ## 15263    285613
    ## 15264     89848
    ## 15265     64411
    ## 15266      5097
    ## 15267      <NA>
    ## 15268      <NA>
    ## 15269      9604
    ## 15270     10007
    ## 15271     80762
    ## 15272     81848
    ## 15273     23092
    ## 15274      2246
    ## 15275      <NA>
    ## 15276      2908
    ## 15277     81555
    ## 15278      <NA>
    ## 15279     57528
    ## 15280    153768
    ## 15281    153769
    ## 15282      <NA>
    ## 15283     54439
    ## 15284     10915
    ## 15285      5521
    ## 15286    202374
    ## 15287      1809
    ## 15288      <NA>
    ## 15289      9832
    ## 15290    167227
    ## 15291      4163
    ## 15292      <NA>
    ## 15293     64848
    ## 15294      3781
    ## 15295      <NA>
    ## 15296     55521
    ## 15297      5229
    ## 15298    153733
    ## 15299     56929
    ## 15300    353376
    ## 15301     51014
    ## 15302      <NA>
    ## 15303      1036
    ## 15304      9140
    ## 15305      1176
    ## 15306     51397
    ## 15307     57556
    ## 15308      <NA>
    ## 15309    285605
    ## 15310      1657
    ## 15311     25816
    ## 15312      3295
    ## 15313     51334
    ## 15314      <NA>
    ## 15315    153443
    ## 15316      4015
    ## 15317      9627
    ## 15318      6643
    ## 15319     28966
    ## 15320      5480
    ## 15321     93166
    ## 15322    153241
    ## 15323      1456
    ## 15324      <NA>
    ## 15325      <NA>
    ## 15326       501
    ## 15327     51808
    ## 15328      4001
    ## 15329      <NA>
    ## 15330     84466
    ## 15331    133619
    ## 15332    613212
    ## 15333      <NA>
    ## 15334      6558
    ## 15335      2201
    ## 15336     51015
    ## 15337    171019
    ## 15338      <NA>
    ## 15339      <NA>
    ## 15340     85027
    ## 15341     51164
    ## 15342     55696
    ## 15343     11346
    ## 15344      3340
    ## 15345      6208
    ## 15346       972
    ## 15347      6949
    ## 15348    340075
    ## 15349       815
    ## 15350      6534
    ## 15351      5159
    ## 15352      1436
    ## 15353     22993
    ## 15354      1836
    ## 15355    133522
    ## 15356      1452
    ## 15357      <NA>
    ## 15358     27190
    ## 15359     78991
    ## 15360    134266
    ## 15361      <NA>
    ## 15362    134265
    ## 15363     22885
    ## 15364     79628
    ## 15365       154
    ## 15366     81545
    ## 15367      <NA>
    ## 15368      <NA>
    ## 15369      <NA>
    ## 15370      8774
    ## 15371     63895
    ## 15372      9352
    ## 15373     23335
    ## 15374     51046
    ## 15375      9480
    ## 15376      2235
    ## 15377      <NA>
    ## 15378      5205
    ## 15379     23327
    ## 15380     10892
    ## 15381      <NA>
    ## 15382      <NA>
    ## 15383     90701
    ## 15384      2922
    ## 15385      3998
    ## 15386    147372
    ## 15387      5366
    ## 15388      4160
    ## 15389      2774
    ## 15390     57132
    ## 15391     65258
    ## 15392      3613
    ## 15393      1149
    ## 15394     84617
    ## 15395     10939
    ## 15396     10650
    ## 15397     56907
    ## 15398     79959
    ## 15399     56984
    ## 15400      5771
    ## 15401      <NA>
    ## 15402     81929
    ## 15403     55125
    ## 15404       753
    ## 15405    125228
    ## 15406      8731
    ## 15407      6925
    ## 15408      5874
    ## 15409      <NA>
    ## 15410    147323
    ## 15411     11201
    ## 15412      8932
    ## 15413     51320
    ## 15414      4089
    ## 15415     55520
    ## 15416      4200
    ## 15417     83876
    ## 15418      5596
    ## 15419      <NA>
    ## 15420     30827
    ## 15421      4152
    ## 15422      4645
    ## 15423     10449
    ## 15424      9388
    ## 15425      6139
    ## 15426      <NA>
    ## 15427     54808
    ## 15428      4092
    ## 15429      9811
    ## 15430    201501
    ## 15431      4087
    ## 15432     51124
    ## 15433     84064
    ## 15434     83473
    ## 15435      9063
    ## 15436      <NA>
    ## 15437     29906
    ## 15438    494470
    ## 15439      <NA>
    ## 15440    115106
    ## 15441      <NA>
    ## 15442     57724
    ## 15443      6563
    ## 15444      <NA>
    ## 15445     26040
    ## 15446     84552
    ## 15447     22850
    ## 15448     79863
    ## 15449      <NA>
    ## 15450     10907
    ## 15451    440498
    ## 15452      <NA>
    ## 15453     26251
    ## 15454      <NA>
    ## 15455      9150
    ## 15456      4772
    ## 15457    374868
    ## 15458     27164
    ## 15459      <NA>
    ## 15460      4155
    ## 15461      <NA>
    ## 15462      <NA>
    ## 15463     10194
    ## 15464    284273
    ## 15465      <NA>
    ## 15466     84735
    ## 15467     55748
    ## 15468      <NA>
    ## 15469      1528
    ## 15470     29090
    ## 15471     81832
    ## 15472    147381
    ## 15473      9306
    ## 15474     25914
    ## 15475    220164
    ## 15476     54495
    ## 15477      <NA>
    ## 15478      3508
    ## 15479    219927
    ## 15480      1374
    ## 15481     51083
    ## 15482     55291
    ## 15483      4041
    ## 15484      <NA>
    ## 15485     51111
    ## 15486      1119
    ## 15487     10312
    ## 15488      4728
    ## 15489       221
    ## 15490     81622
    ## 15491       222
    ## 15492     91703
    ## 15493      <NA>
    ## 15494      <NA>
    ## 15495      4723
    ## 15496      2950
    ## 15497     10263
    ## 15498      9600
    ## 15499      9049
    ## 15500     57175
    ## 15501      5790
    ## 15502      6199
    ## 15503     57571
    ## 15504    374403
    ## 15505      5499
    ## 15506      5883
    ## 15507     23529
    ## 15508     57804
    ## 15509     54961
    ## 15510    338692
    ## 15511       156
    ## 15512     22992
    ## 15513     29984
    ## 15514     91683
    ## 15515      <NA>
    ## 15516     78999
    ## 15517      9986
    ## 15518      <NA>
    ## 15519      6712
    ## 15520      <NA>
    ## 15521     83759
    ## 15522      5936
    ## 15523     10432
    ## 15524      9973
    ## 15525     55231
    ## 15526      8722
    ## 15527        89
    ## 15528    254359
    ## 15529       582
    ## 15530     10072
    ## 15531    246330
    ## 15532     65003
    ## 15533    266743
    ## 15534      3177
    ## 15535     11041
    ## 15536     25855
    ## 15537      9610
    ## 15538     57124
    ## 15539    256472
    ## 15540     10897
    ## 15541    254263
    ## 15542     81876
    ## 15543     64837
    ## 15544     55690
    ## 15545     10992
    ## 15546     89792
    ## 15547      1474
    ## 15548      8815
    ## 15549     84285
    ## 15550      9092
    ## 15551      <NA>
    ## 15552     10589
    ## 15553      <NA>
    ## 15554     11007
    ## 15555      9158
    ## 15556     30008
    ## 15557     80198
    ## 15558      1072
    ## 15559    254122
    ## 15560     91056
    ## 15561     84153
    ## 15562     10524
    ## 15563      5970
    ## 15564      6494
    ## 15565    399909
    ## 15566      4296
    ## 15567    254102
    ## 15568     23625
    ## 15569      <NA>
    ## 15570      4054
    ## 15571     57410
    ## 15572    378938
    ## 15573    283131
    ## 15574      <NA>
    ## 15575     83786
    ## 15576    283130
    ## 15577    220359
    ## 15578      5977
    ## 15579      <NA>
    ## 15580     10435
    ## 15581     23649
    ## 15582       823
    ## 15583      <NA>
    ## 15584     84447
    ## 15585       740
    ## 15586      2197
    ## 15587       741
    ## 15588      7108
    ## 15589       738
    ## 15590      7542
    ## 15591     29901
    ## 15592     29907
    ## 15593       402
    ## 15594      5526
    ## 15595     23130
    ## 15596     10938
    ## 15597     55561
    ## 15598      4221
    ## 15599      5871
    ## 15600      <NA>
    ## 15601      7536
    ## 15602      5837
    ## 15603     10235
    ## 15604      9379
    ## 15605      <NA>
    ## 15606      8986
    ## 15607    283234
    ## 15608     25824
    ## 15609     51504
    ## 15610      2101
    ## 15611     25858
    ## 15612     50801
    ## 15613     56834
    ## 15614       572
    ## 15615      5331
    ## 15616     26472
    ## 15617      2286
    ## 15618      7423
    ## 15619      3338
    ## 15620     84304
    ## 15621     83707
    ## 15622     83706
    ## 15623     10963
    ## 15624     28992
    ## 15625     23769
    ## 15626     55611
    ## 15627      1351
    ## 15628     79829
    ## 15629    283248
    ## 15630      2011
    ## 15631    144097
    ## 15632      <NA>
    ## 15633     10313
    ## 15634     25923
    ## 15635      <NA>
    ## 15636     85329
    ## 15637      9376
    ## 15638      9356
    ## 15639      1128
    ## 15640      6520
    ## 15641     54663
    ## 15642      <NA>
    ## 15643     10482
    ## 15644     10629
    ## 15645      5436
    ## 15646     79842
    ## 15647    283237
    ## 15648    221092
    ## 15649      2785
    ## 15650     26580
    ## 15651    221091
    ## 15652     51035
    ## 15653    790955
    ## 15654     79081
    ## 15655      <NA>
    ## 15656      <NA>
    ## 15657     23193
    ## 15658     26229
    ## 15659      6094
    ## 15660    256364
    ## 15661      <NA>
    ## 15662      9219
    ## 15663     64852
    ## 15664      1937
    ## 15665     79026
    ## 15666     80150
    ## 15667      3619
    ## 15668      2495
    ## 15669      7439
    ## 15670      5866
    ## 15671      3995
    ## 15672      9415
    ## 15673      3992
    ## 15674      2237
    ## 15675       746
    ## 15676       745
    ## 15677       747
    ## 15678      9066
    ## 15679    390205
    ## 15680    220004
    ## 15681     54949
    ## 15682     79869
    ## 15683     51259
    ## 15684     51524
    ## 15685     26007
    ## 15686      1642
    ## 15687     55048
    ## 15688       923
    ## 15689     51296
    ## 15690     54972
    ## 15691     79073
    ## 15692     27339
    ## 15693     79080
    ## 15694    219995
    ## 15695      <NA>
    ## 15696     58475
    ## 15697      <NA>
    ## 15698      <NA>
    ## 15699      <NA>
    ## 15700      <NA>
    ## 15701     54948
    ## 15702      6809
    ## 15703    219988
    ## 15704      5007
    ## 15705    219972
    ## 15706     23220
    ## 15707     63901
    ## 15708      <NA>
    ## 15709     80829
    ## 15710      9404
    ## 15711      7091
    ## 15712     29968
    ## 15713     84131
    ## 15714      2776
    ## 15715     23230
    ## 15716    158471
    ## 15717      2650
    ## 15718     55312
    ## 15719      5125
    ## 15720     26578
    ## 15721     54981
    ## 15722    138199
    ## 15723      <NA>
    ## 15724    140803
    ## 15725      6096
    ## 15726       301
    ## 15727       216
    ## 15728      <NA>
    ## 15729      <NA>
    ## 15730      7763
    ## 15731      9615
    ## 15732      <NA>
    ## 15733     51104
    ## 15734      <NA>
    ## 15735     80036
    ## 15736      <NA>
    ## 15737      <NA>
    ## 15738       687
    ## 15739     23137
    ## 15740    256691
    ## 15741      <NA>
    ## 15742      <NA>
    ## 15743      <NA>
    ## 15744    375743
    ## 15745       320
    ## 15746      9413
    ## 15747      9414
    ## 15748      2395
    ## 15749      8395
    ## 15750    116224
    ## 15751    169693
    ## 15752      5239
    ## 15753      <NA>
    ## 15754     55871
    ## 15755     81704
    ## 15756     23189
    ## 15757     10655
    ## 15758      6595
    ## 15759      7436
    ## 15760      9933
    ## 15761      5991
    ## 15762    169792
    ## 15763      6505
    ## 15764      <NA>
    ## 15765    403313
    ## 15766     55664
    ## 15767     50808
    ## 15768      <NA>
    ## 15769     10171
    ## 15770      3717
    ## 15771     11172
    ## 15772      6013
    ## 15773     55848
    ## 15774     29126
    ## 15775      <NA>
    ## 15776     57589
    ## 15777     79956
    ## 15778      <NA>
    ## 15779     26953
    ## 15780     90865
    ## 15781    115426
    ## 15782      2731
    ## 15783      5592
    ## 15784     23283
    ## 15785     56624
    ## 15786    259230
    ## 15787      <NA>
    ## 15788      <NA>
    ## 15789      9562
    ## 15790      9060
    ## 15791     84896
    ## 15792      5728
    ## 15793     55328
    ## 15794      <NA>
    ## 15795     57559
    ## 15796        59
    ## 15797       355
    ## 15798      9023
    ## 15799      <NA>
    ## 15800      3988
    ## 15801      3433
    ## 15802      3437
    ## 15803      <NA>
    ## 15804      <NA>
    ## 15805      <NA>
    ## 15806      3434
    ## 15807    387700
    ## 15808     53354
    ## 15809      9585
    ## 15810      <NA>
    ## 15811      3363
    ## 15812     10556
    ## 15813     84333
    ## 15814    143279
    ## 15815      5507
    ## 15816     80351
    ## 15817    143282
    ## 15818      9044
    ## 15819     22849
    ## 15820      <NA>
    ## 15821      3416
    ## 15822      3832
    ## 15823      3087
    ## 15824     54536
    ## 15825     26509
    ## 15826     55165
    ## 15827      5950
    ## 15828      5146
    ## 15829    118924
    ## 15830      9211
    ## 15831    159371
    ## 15832     51196
    ## 15833     64318
    ## 15834     23232
    ## 15835      3070
    ## 15836      9124
    ## 15837     10580
    ## 15838      5832
    ## 15839     26123
    ## 15840       953
    ## 15841     54619
    ## 15842      <NA>
    ## 15843      <NA>
    ## 15844     29760
    ## 15845     93377
    ## 15846     56889
    ## 15847     84458
    ## 15848      <NA>
    ## 15849      6585
    ## 15850     84986
    ## 15851     10023
    ## 15852     23401
    ## 15853     23223
    ## 15854      5223
    ## 15855     51013
    ## 15856     84287
    ## 15857     64210
    ## 15858     80019
    ## 15859     26287
    ## 15860    112817
    ## 15861    118812
    ## 15862     55361
    ## 15863     60370
    ## 15864     83742
    ## 15865    118813
    ## 15866      6425
    ## 15867    401647
    ## 15868     55118
    ## 15869     27291
    ## 15870     84795
    ## 15871      3257
    ## 15872      <NA>
    ## 15873     26507
    ## 15874      2805
    ## 15875     81894
    ## 15876     57089
    ## 15877      1355
    ## 15878     51076
    ## 15879     23268
    ## 15880      <NA>
    ## 15881     10613
    ## 15882      1147
    ## 15883     55280
    ## 15884    282991
    ## 15885      9033
    ## 15886      <NA>
    ## 15887      <NA>
    ## 15888      <NA>
    ## 15889      4714
    ## 15890     55662
    ## 15891     55719
    ## 15892     57715
    ## 15893     84545
    ## 15894     56652
    ## 15895     84445
    ## 15896     79955
    ## 15897     81855
    ## 15898     81621
    ## 15899      <NA>
    ## 15900      <NA>
    ## 15901      <NA>
    ## 15902      8945
    ## 15903     27343
    ## 15904      <NA>
    ## 15905      6468
    ## 15906     10360
    ## 15907      <NA>
    ## 15908     30819
    ## 15909      <NA>
    ## 15910     79803
    ## 15911      8861
    ## 15912     23082
    ## 15913      9221
    ## 15914      8729
    ## 15915      4791
    ## 15916      5662
    ## 15917     79176
    ## 15918     79004
    ## 15919     79847
    ## 15920      <NA>
    ## 15921     10121
    ## 15922     51684
    ## 15923     81603
    ## 15924       403
    ## 15925    118980
    ## 15926     54838
    ## 15927    119032
    ## 15928     57412
    ## 15929     54805
    ## 15930     22978
    ## 15931      9118
    ## 15932     84108
    ## 15933      6877
    ## 15934      <NA>
    ## 15935     22984
    ## 15936     51063
    ## 15937      <NA>
    ## 15938      9644
    ## 15939     79991
    ## 15940      9748
    ## 15941      1308
    ## 15942    119392
    ## 15943     80217
    ## 15944      9446
    ## 15945    119391
    ## 15946    159686
    ## 15947     22986
    ## 15948    114815
    ## 15949      7511
    ## 15950       120
    ## 15951      4601
    ## 15952     10285
    ## 15953      <NA>
    ## 15954      1847
    ## 15955      9126
    ## 15956    282996
    ## 15957     27250
    ## 15958     92482
    ## 15959      8036
    ## 15960       150
    ## 15961     57678
    ## 15962      <NA>
    ## 15963     51703
    ## 15964     64429
    ## 15965    143187
    ## 15966      6934
    ## 15967      4892
    ## 15968       840
    ## 15969     79949
    ## 15970      9937
    ## 15971    374354
    ## 15972       153
    ## 15973     55088
    ## 15974     56165
    ## 15975     84632
    ## 15976      3983
    ## 15977      <NA>
    ## 15978     57700
    ## 15979    142940
    ## 15980     26033
    ## 15981      2674
    ## 15982      5406
    ## 15983    259217
    ## 15984    387712
    ## 15985     57698
    ## 15986     11023
    ## 15987    118987
    ## 15988      <NA>
    ## 15989    196047
    ## 15990      2018
    ## 15991     22841
    ## 15992     63877
    ## 15993    143384
    ## 15994    340719
    ## 15995      8661
    ## 15996      <NA>
    ## 15997    119559
    ## 15998     10935
    ## 15999      2869
    ## 16000      <NA>
    ## 16001      1438
    ## 16002      <NA>
    ## 16003      <NA>
    ## 16004      <NA>
    ## 16005      <NA>
    ## 16006      <NA>
    ## 16007      <NA>
    ## 16008      <NA>
    ## 16009      <NA>
    ## 16010      <NA>
    ## 16011      <NA>
    ## 16012      <NA>
    ## 16013      <NA>
    ## 16014      <NA>
    ## 16015      6845
    ## 16016     10251
    ## 16017     55217
    ## 16018      <NA>
    ## 16019      <NA>
    ## 16020      <NA>

    test <- gene.list.1.0
    test[1:5]

    ##         Xkr4        Sox17       Mrpl15       Lypla1        Tcea1 
    ## -0.018046246 -0.010294148  0.006406841  0.206258596  0.014368925

    names(test)[1:5]

    ## [1] "Xkr4"   "Sox17"  "Mrpl15" "Lypla1" "Tcea1"

    typeof(unlist(annot['ENTREZID']))

    ## [1] "character"

    names(test) <- unlist(annot['ENTREZID'])
    test <- sort(test, decreasing=TRUE)
    test[1:5]

    ##      3776    149111     25907     57512      2675 
    ## 0.4946314 0.4868727 0.4770195 0.4432751 0.4220566

Older Spring 2019 code below.
=============================

    gene <- names(test)[abs(test) > 2]
    ggo <- groupGO(gene     = gene,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   level    = 2,
                   readable = TRUE)
    #order by count
    ggo@result = ggo@result[order(ggo@result[,"Count"],decreasing=TRUE),]
    ggo@result

    ##                    ID                      Description Count GeneRatio geneID
    ## GO:0003824 GO:0003824               catalytic activity     0       0/0       
    ## GO:0005198 GO:0005198     structural molecule activity     0       0/0       
    ## GO:0005215 GO:0005215             transporter activity     0       0/0       
    ## GO:0005488 GO:0005488                          binding     0       0/0       
    ## GO:0016209 GO:0016209             antioxidant activity     0       0/0       
    ## GO:0031386 GO:0031386                      protein tag     0       0/0       
    ## GO:0038024 GO:0038024          cargo receptor activity     0       0/0       
    ## GO:0044183 GO:0044183        protein folding chaperone     0       0/0       
    ## GO:0045182 GO:0045182   translation regulator activity     0       0/0       
    ## GO:0045735 GO:0045735      nutrient reservoir activity     0       0/0       
    ## GO:0060089 GO:0060089    molecular transducer activity     0       0/0       
    ## GO:0090729 GO:0090729                   toxin activity     0       0/0       
    ## GO:0098772 GO:0098772     molecular function regulator     0       0/0       
    ## GO:0140104 GO:0140104       molecular carrier activity     0       0/0       
    ## GO:0140110 GO:0140110 transcription regulator activity     0       0/0       
    ## GO:0140299 GO:0140299   small molecule sensor activity     0       0/0       
    ## GO:0140313 GO:0140313  molecular sequestering activity     0       0/0

    barplot(ggo,showCategory=8)

![](aug-5---enrichment_files/figure-markdown_strict/unnamed-chunk-18-1.png)

Further information about enirchment scores:
<a href="https://yulab-smu.github.io/clusterProfiler-book/chapter2.html#gene-set-enrichment-analysis" class="uri">https://yulab-smu.github.io/clusterProfiler-book/chapter2.html#gene-set-enrichment-analysis</a>

    #KEGG Gene Set Enrichment Analysis
    kk2 <- gseKEGG(geneList     = test,
                   organism     = 'hsa',
                   nPerm        = 1000,
                   minGSSize    = 100,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)

    ## Reading KEGG annotation online:
    ## 
    ## Reading KEGG annotation online:

    ## Warning in .GSEA(geneList = geneList, exponent = exponent, minGSSize =
    ## minGSSize, : We do not recommend using nPerm parameter incurrent and future
    ## releases

    ## Warning in fgsea(pathways = geneSets, stats = geneList, nperm = nPerm, minSize
    ## = minGSSize, : You are trying to run fgseaSimple. It is recommended to use
    ## fgseaMultilevel. To run fgseaMultilevel, you need to remove the nperm argument
    ## in the fgsea function call.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (11.71% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize,
    ## gseaParam, : There are duplicate gene names, fgsea may produce unexpected
    ## results.

    ## no term enriched under specific pvalueCutoff...

    kk2@result = kk2@result[order(kk2@result[,"setSize"],decreasing=TRUE),]
    kk2@result

    ## [1] ID              Description     setSize         enrichmentScore
    ## [5] NES             pvalue          p.adjust        qvalues        
    ## <0 rows> (or 0-length row.names)

    #KEGG over representation Analysis
    kk3 <- enrichKEGG(gene    = names(test),
                   organism     = 'hsa',
                   pvalueCutoff = 0.05
      )
    kk3@result[,c("Description","GeneRatio","p.adjust","Count")]

    ##                                                                      Description
    ## hsa04360                                                           Axon guidance
    ## hsa04010                                                  MAPK signaling pathway
    ## hsa04510                                                          Focal adhesion
    ## hsa04015                                                  Rap1 signaling pathway
    ## hsa05017                                                  Spinocerebellar ataxia
    ## hsa05205                                                 Proteoglycans in cancer
    ## hsa05132                                                    Salmonella infection
    ## hsa04140                                                      Autophagy - animal
    ## hsa05225                                                Hepatocellular carcinoma
    ## hsa04120                                          Ubiquitin mediated proteolysis
    ## hsa04144                                                             Endocytosis
    ## hsa04919                                       Thyroid hormone signaling pathway
    ## hsa05010                                                       Alzheimer disease
    ## hsa05165                                          Human papillomavirus infection
    ## hsa04150                                                  mTOR signaling pathway
    ## hsa04921                                              Oxytocin signaling pathway
    ## hsa04071                                          Sphingolipid signaling pathway
    ## hsa04722                                          Neurotrophin signaling pathway
    ## hsa04110                                                              Cell cycle
    ## hsa05210                                                       Colorectal cancer
    ## hsa05014                                           Amyotrophic lateral sclerosis
    ## hsa04810                                        Regulation of actin cytoskeleton
    ## hsa05212                                                       Pancreatic cancer
    ## hsa05220                                                Chronic myeloid leukemia
    ## hsa04725                                                     Cholinergic synapse
    ## hsa04261                                  Adrenergic signaling in cardiomyocytes
    ## hsa04390                                                 Hippo signaling pathway
    ## hsa01521                               EGFR tyrosine kinase inhibitor resistance
    ## hsa04012                                                  ErbB signaling pathway
    ## hsa04520                                                       Adherens junction
    ## hsa04141                             Protein processing in endoplasmic reticulum
    ## hsa04330                                                 Notch signaling pathway
    ## hsa05224                                                           Breast cancer
    ## hsa04728                                                    Dopaminergic synapse
    ## hsa04910                                               Insulin signaling pathway
    ## hsa04068                                                  FoxO signaling pathway
    ## hsa04020                                               Calcium signaling pathway
    ## hsa04934                                                        Cushing syndrome
    ## hsa04724                                                   Glutamatergic synapse
    ## hsa04926                                               Relaxin signaling pathway
    ## hsa00562                                           Inositol phosphate metabolism
    ## hsa04723                                    Retrograde endocannabinoid signaling
    ## hsa04142                                                                Lysosome
    ## hsa05213                                                      Endometrial cancer
    ## hsa04218                                                     Cellular senescence
    ## hsa05412                         Arrhythmogenic right ventricular cardiomyopathy
    ## hsa04022                                              cGMP-PKG signaling pathway
    ## hsa04014                                                   Ras signaling pathway
    ## hsa05211                                                    Renal cell carcinoma
    ## hsa05231                                            Choline metabolism in cancer
    ## hsa05222                                                  Small cell lung cancer
    ## hsa04070                                   Phosphatidylinositol signaling system
    ## hsa05223                                              Non-small cell lung cancer
    ## hsa04310                                                   Wnt signaling pathway
    ## hsa00310                                                      Lysine degradation
    ## hsa04024                                                  cAMP signaling pathway
    ## hsa05226                                                          Gastric cancer
    ## hsa05161                                                             Hepatitis B
    ## hsa05418                                  Fluid shear stress and atherosclerosis
    ## hsa04211                                            Longevity regulating pathway
    ## hsa04935                          Growth hormone synthesis, secretion and action
    ## hsa05016                                                      Huntington disease
    ## hsa04932                                       Non-alcoholic fatty liver disease
    ## hsa05166                                 Human T-cell leukemia virus 1 infection
    ## hsa04210                                                               Apoptosis
    ## hsa04713                                                   Circadian entrainment
    ## hsa05215                                                         Prostate cancer
    ## hsa05214                                                                  Glioma
    ## hsa05414                                                  Dilated cardiomyopathy
    ## hsa05135                                                      Yersinia infection
    ## hsa04072                                       Phospholipase D signaling pathway
    ## hsa04933                    AGE-RAGE signaling pathway in diabetic complications
    ## hsa04727                                                       GABAergic synapse
    ## hsa04720                                                  Long-term potentiation
    ## hsa01522                                                    Endocrine resistance
    ## hsa05130                                   Pathogenic Escherichia coli infection
    ## hsa04151                                              PI3K-Akt signaling pathway
    ## hsa04152                                                  AMPK signaling pathway
    ## hsa05032                                                      Morphine addiction
    ## hsa04371                                                Apelin signaling pathway
    ## hsa05410                                             Hypertrophic cardiomyopathy
    ## hsa04540                                                            Gap junction
    ## hsa05131                                                             Shigellosis
    ## hsa04066                                                 HIF-1 signaling pathway
    ## hsa03410                                                    Base excision repair
    ## hsa04530                                                          Tight junction
    ## hsa04929                                                          GnRH secretion
    ## hsa05163                                         Human cytomegalovirus infection
    ## hsa05235                  PD-L1 expression and PD-1 checkpoint pathway in cancer
    ## hsa03040                                                             Spliceosome
    ## hsa04137                                                      Mitophagy - animal
    ## hsa04146                                                              Peroxisome
    ## hsa04611                                                     Platelet activation
    ## hsa05221                                                  Acute myeloid leukemia
    ## hsa03015                                               mRNA surveillance pathway
    ## hsa00410                                                 beta-Alanine metabolism
    ## hsa00510                                                   N-Glycan biosynthesis
    ## hsa04550                Signaling pathways regulating pluripotency of stem cells
    ## hsa04625                                C-type lectin receptor signaling pathway
    ## hsa04931                                                      Insulin resistance
    ## hsa04730                                                    Long-term depression
    ## hsa05012                                                       Parkinson disease
    ## hsa04925                                     Aldosterone synthesis and secretion
    ## hsa03030                                                         DNA replication
    ## hsa04666                                        Fc gamma R-mediated phagocytosis
    ## hsa03018                                                         RNA degradation
    ## hsa04115                                                   p53 signaling pathway
    ## hsa05217                                                    Basal cell carcinoma
    ## hsa04670                                    Leukocyte transendothelial migration
    ## hsa04213                         Longevity regulating pathway - multiple species
    ## hsa03050                                                              Proteasome
    ## hsa04350                                              TGF-beta signaling pathway
    ## hsa04714                                                           Thermogenesis
    ## hsa05230                                     Central carbon metabolism in cancer
    ## hsa04916                                                           Melanogenesis
    ## hsa04136                                                       Autophagy - other
    ## hsa05031                                                   Amphetamine addiction
    ## hsa04370                                                  VEGF signaling pathway
    ## hsa05100                                  Bacterial invasion of epithelial cells
    ## hsa00280                              Valine, leucine and isoleucine degradation
    ## hsa05216                                                          Thyroid cancer
    ## hsa05218                                                                Melanoma
    ## hsa00600                                                 Sphingolipid metabolism
    ## hsa03420                                              Nucleotide excision repair
    ## hsa04912                                                  GnRH signaling pathway
    ## hsa04512                                                ECM-receptor interaction
    ## hsa01200                                                       Carbon metabolism
    ## hsa04260                                              Cardiac muscle contraction
    ## hsa04114                                                          Oocyte meiosis
    ## hsa01524                                                Platinum drug resistance
    ## hsa05170                                Human immunodeficiency virus 1 infection
    ## hsa03460                                                  Fanconi anemia pathway
    ## hsa05030                                                       Cocaine addiction
    ## hsa04928                     Parathyroid hormone synthesis, secretion and action
    ## hsa00230                                                       Purine metabolism
    ## hsa00520                             Amino sugar and nucleotide sugar metabolism
    ## hsa04668                                                   TNF signaling pathway
    ## hsa05167                         Kaposi sarcoma-associated herpesvirus infection
    ## hsa04917                                             Prolactin signaling pathway
    ## hsa05120              Epithelial cell signaling in Helicobacter pylori infection
    ## hsa03010                                                                Ribosome
    ## hsa04914                                 Progesterone-mediated oocyte maturation
    ## hsa04911                                                       Insulin secretion
    ## hsa04920                                         Adipocytokine signaling pathway
    ## hsa04924                                                         Renin secretion
    ## hsa03020                                                          RNA polymerase
    ## hsa00564                                          Glycerophospholipid metabolism
    ## hsa05142                                                          Chagas disease
    ## hsa00563                  Glycosylphosphatidylinositol (GPI)-anchor biosynthesis
    ## hsa03440                                                Homologous recombination
    ## hsa04340                                              Hedgehog signaling pathway
    ## hsa04380                                              Osteoclast differentiation
    ## hsa00534              Glycosaminoglycan biosynthesis - heparan sulfate / heparin
    ## hsa04064                                            NF-kappa B signaling pathway
    ## hsa04660                                       T cell receptor signaling pathway
    ## hsa01212                                                   Fatty acid metabolism
    ## hsa03430                                                         Mismatch repair
    ## hsa03013                                                           RNA transport
    ## hsa00514                                    Other types of O-glycan biosynthesis
    ## hsa00900                                         Terpenoid backbone biosynthesis
    ## hsa04971                                                  Gastric acid secretion
    ## hsa00062                                                   Fatty acid elongation
    ## hsa04930                                               Type II diabetes mellitus
    ## hsa04215                                            Apoptosis - multiple species
    ## hsa05110                                               Vibrio cholerae infection
    ## hsa00604                         Glycosphingolipid biosynthesis - ganglio series
    ## hsa03022                                             Basal transcription factors
    ## hsa04710                                                        Circadian rhythm
    ## hsa04721                                                  Synaptic vesicle cycle
    ## hsa05033                                                      Nicotine addiction
    ## hsa04961               Endocrine and other factor-regulated calcium reabsorption
    ## hsa04962                                Vasopressin-regulated water reabsorption
    ## hsa04927                                        Cortisol synthesis and secretion
    ## hsa00513                                  Various types of N-glycan biosynthesis
    ## hsa04923                                   Regulation of lipolysis in adipocytes
    ## hsa04922                                              Glucagon signaling pathway
    ## hsa04270                                      Vascular smooth muscle contraction
    ## hsa00051                                         Fructose and mannose metabolism
    ## hsa04130                               SNARE interactions in vesicular transport
    ## hsa05160                                                             Hepatitis C
    ## hsa00511                                                Other glycan degradation
    ## hsa03060                                                          Protein export
    ## hsa05219                                                          Bladder cancer
    ## hsa04750                        Inflammatory mediator regulation of TRP channels
    ## hsa00561                                                 Glycerolipid metabolism
    ## hsa00480                                                  Glutathione metabolism
    ## hsa00250                             Alanine, aspartate and glutamate metabolism
    ## hsa05162                                                                 Measles
    ## hsa00450                                               Selenocompound metabolism
    ## hsa05169                                            Epstein-Barr virus infection
    ## hsa00760                                  Nicotinate and nicotinamide metabolism
    ## hsa04662                                       B cell receptor signaling pathway
    ## hsa00020                                               Citrate cycle (TCA cycle)
    ## hsa00640                                                   Propanoate metabolism
    ## hsa00920                                                       Sulfur metabolism
    ## hsa04915                                              Estrogen signaling pathway
    ## hsa04392                              Hippo signaling pathway - multiple species
    ## hsa00532 Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate
    ## hsa00270                                      Cysteine and methionine metabolism
    ## hsa04216                                                             Ferroptosis
    ## hsa04664                                         Fc epsilon RI signaling pathway
    ## hsa05145                                                           Toxoplasmosis
    ## hsa04062                                             Chemokine signaling pathway
    ## hsa00260                                Glycine, serine and threonine metabolism
    ## hsa00330                                         Arginine and proline metabolism
    ## hsa00240                                                   Pyrimidine metabolism
    ## hsa00630                                 Glyoxylate and dicarboxylate metabolism
    ## hsa04918                                               Thyroid hormone synthesis
    ## hsa01230                                             Biosynthesis of amino acids
    ## hsa04726                                                    Serotonergic synapse
    ## hsa05133                                                               Pertussis
    ## hsa00670                                               One carbon pool by folate
    ## hsa05134                                                           Legionellosis
    ## hsa00380                                                   Tryptophan metabolism
    ## hsa01040                                 Biosynthesis of unsaturated fatty acids
    ## hsa00515                                      Mannose type O-glycan biosynthesis
    ## hsa00531                                           Glycosaminoglycan degradation
    ## hsa05203                                                    Viral carcinogenesis
    ## hsa04960                               Aldosterone-regulated sodium reabsorption
    ## hsa00340                                                    Histidine metabolism
    ## hsa00072                              Synthesis and degradation of ketone bodies
    ## hsa01210                                         2-Oxocarboxylic acid metabolism
    ## hsa04979                                                  Cholesterol metabolism
    ## hsa00770                                       Pantothenate and CoA biosynthesis
    ## hsa00620                                                     Pyruvate metabolism
    ## hsa04620                                    Toll-like receptor signaling pathway
    ## hsa03450                                              Non-homologous end-joining
    ## hsa00052                                                    Galactose metabolism
    ## hsa01523                                                   Antifolate resistance
    ## hsa02010                                                        ABC transporters
    ## hsa00071                                                  Fatty acid degradation
    ## hsa05020                                                          Prion diseases
    ## hsa00512                                        Mucin type O-glycan biosynthesis
    ## hsa00130                     Ubiquinone and other terpenoid-quinone biosynthesis
    ## hsa00430                                      Taurine and hypotaurine metabolism
    ## hsa00061                                                 Fatty acid biosynthesis
    ## hsa00190                                               Oxidative phosphorylation
    ## hsa04659                                               Th17 cell differentiation
    ## hsa00533                        Glycosaminoglycan biosynthesis - keratan sulfate
    ## hsa00010                                            Glycolysis / Gluconeogenesis
    ## hsa05202                                 Transcriptional misregulation in cancer
    ## hsa04966                                          Collecting duct acid secretion
    ## hsa00030                                               Pentose phosphate pathway
    ## hsa04964                                 Proximal tubule bicarbonate reclamation
    ## hsa05152                                                            Tuberculosis
    ## hsa00790                                                     Folate biosynthesis
    ## hsa00360                                                Phenylalanine metabolism
    ## hsa04621                                     NOD-like receptor signaling pathway
    ## hsa04658                                        Th1 and Th2 cell differentiation
    ## hsa00650                                                    Butanoate metabolism
    ## hsa00603              Glycosphingolipid biosynthesis - globo and isoglobo series
    ## hsa00220                                                   Arginine biosynthesis
    ## hsa04970                                                      Salivary secretion
    ## hsa00601              Glycosphingolipid biosynthesis - lacto and neolacto series
    ## hsa00100                                                    Steroid biosynthesis
    ## hsa05146                                                              Amoebiasis
    ## hsa00565                                                  Ether lipid metabolism
    ## hsa04977                                        Vitamin digestion and absorption
    ## hsa03320                                                  PPAR signaling pathway
    ## hsa04973                                   Carbohydrate digestion and absorption
    ## hsa00350                                                     Tyrosine metabolism
    ## hsa04614                                                Renin-angiotensin system
    ## hsa04622                                   RIG-I-like receptor signaling pathway
    ## hsa04623                                           Cytosolic DNA-sensing pathway
    ## hsa04974                                        Protein digestion and absorption
    ## hsa05140                                                           Leishmaniasis
    ## hsa00730                                                     Thiamine metabolism
    ## hsa04657                                                 IL-17 signaling pathway
    ## hsa03008                                       Ribosome biogenesis in eukaryotes
    ## hsa05340                                                Primary immunodeficiency
    ## hsa04145                                                               Phagosome
    ## hsa00120                                          Primary bile acid biosynthesis
    ## hsa05164                                                             Influenza A
    ## hsa04972                                                    Pancreatic secretion
    ## hsa04913                                                 Ovarian steroidogenesis
    ## hsa04978                                                      Mineral absorption
    ## hsa05144                                                                 Malaria
    ## hsa05416                                                       Viral myocarditis
    ## hsa04630                                              JAK-STAT signaling pathway
    ## hsa04217                                                             Necroptosis
    ## hsa05143                                                 African trypanosomiasis
    ## hsa00983                                         Drug metabolism - other enzymes
    ## hsa00500                                           Starch and sucrose metabolism
    ## hsa04514                                                 Cell adhesion molecules
    ## hsa05323                                                    Rheumatoid arthritis
    ## hsa00860                                    Porphyrin and chlorophyll metabolism
    ## hsa04975                                            Fat digestion and absorption
    ## hsa04080                                 Neuroactive ligand-receptor interaction
    ## hsa04744                                                       Phototransduction
    ## hsa00590                                             Arachidonic acid metabolism
    ## hsa00592                                         alpha-Linolenic acid metabolism
    ## hsa04650                               Natural killer cell mediated cytotoxicity
    ## hsa05321                                              Inflammatory bowel disease
    ## hsa04672                            Intestinal immune network for IgA production
    ## hsa04742                                                      Taste transduction
    ## hsa04610                                     Complement and coagulation cascades
    ## hsa04976                                                          Bile secretion
    ## hsa00040                                Pentose and glucuronate interconversions
    ## hsa04940                                                Type I diabetes mellitus
    ## hsa04950                                    Maturity onset diabetes of the young
    ## hsa00910                                                     Nitrogen metabolism
    ## hsa05034                                                              Alcoholism
    ## hsa00591                                                Linoleic acid metabolism
    ## hsa00053                                       Ascorbate and aldarate metabolism
    ## hsa04640                                              Hematopoietic cell lineage
    ## hsa00970                                             Aminoacyl-tRNA biosynthesis
    ## hsa04612                                     Antigen processing and presentation
    ## hsa04061           Viral protein interaction with cytokine and cytokine receptor
    ## hsa05204                                                 Chemical carcinogenesis
    ## hsa00982                                       Drug metabolism - cytochrome P450
    ## hsa00980                            Metabolism of xenobiotics by cytochrome P450
    ## hsa05330                                                     Allograft rejection
    ## hsa05310                                                                  Asthma
    ## hsa05332                                               Graft-versus-host disease
    ## hsa04060                                  Cytokine-cytokine receptor interaction
    ## hsa05320                                              Autoimmune thyroid disease
    ## hsa00140                                            Steroid hormone biosynthesis
    ## hsa05206                                                     MicroRNAs in cancer
    ## hsa00830                                                      Retinol metabolism
    ## hsa04740                                                  Olfactory transduction
    ## hsa05150                                         Staphylococcus aureus infection
    ## hsa05168                                        Herpes simplex virus 1 infection
    ## hsa05322                                            Systemic lupus erythematosus
    ##          GeneRatio     p.adjust Count
    ## hsa04360  174/5470 4.265201e-19   174
    ## hsa04010  267/5470 4.816900e-19   267
    ## hsa04510  184/5470 7.586440e-14   184
    ## hsa04015  191/5470 9.137233e-14   191
    ## hsa05017  127/5470 1.816115e-13   127
    ## hsa05205  186/5470 2.738751e-13   186
    ## hsa05132  192/5470 4.542971e-13   192
    ## hsa04140  129/5470 2.059265e-12   129
    ## hsa05225  154/5470 4.682862e-12   154
    ## hsa04120  130/5470 2.790368e-11   130
    ## hsa04144  216/5470 3.840718e-11   216
    ## hsa04919  114/5470 3.840718e-11   114
    ## hsa05010  309/5470 3.915561e-11   309
    ## hsa05165  279/5470 9.727228e-11   279
    ## hsa04150  141/5470 1.468890e-10   141
    ## hsa04921  140/5470 1.858634e-10   140
    ## hsa04071  111/5470 3.425516e-10   111
    ## hsa04722  111/5470 3.425516e-10   111
    ## hsa04110  115/5470 4.031046e-10   115
    ## hsa05210   83/5470 5.950156e-10    83
    ## hsa05014  301/5470 1.028114e-09   301
    ## hsa04810  186/5470 1.120444e-09   186
    ## hsa05212   74/5470 1.406378e-09    74
    ## hsa05220   74/5470 1.406378e-09    74
    ## hsa04725  105/5470 1.702443e-09   105
    ## hsa04261  134/5470 2.239704e-09   134
    ## hsa04390  140/5470 3.710159e-09   140
    ## hsa01521   76/5470 5.087139e-09    76
    ## hsa04012   81/5470 5.902198e-09    81
    ## hsa04520   69/5470 6.926471e-09    69
    ## hsa04141  147/5470 1.077941e-08   147
    ## hsa04330   53/5470 1.204156e-08    53
    ## hsa05224  131/5470 1.244610e-08   131
    ## hsa04728  119/5470 1.244610e-08   119
    ## hsa04910  123/5470 1.244610e-08   123
    ## hsa04068  118/5470 1.607205e-08   118
    ## hsa04020  173/5470 1.736030e-08   173
    ## hsa04934  137/5470 1.736030e-08   137
    ## hsa04724  104/5470 2.362597e-08   104
    ## hsa04926  116/5470 2.578014e-08   116
    ## hsa00562   70/5470 2.858044e-08    70
    ## hsa04723  131/5470 2.878309e-08   131
    ## hsa04142  115/5470 3.197385e-08   115
    ## hsa05213   57/5470 3.577177e-08    57
    ## hsa04218  137/5470 3.983663e-08   137
    ## hsa05412   73/5470 5.706578e-08    73
    ## hsa04022  145/5470 7.967030e-08   145
    ## hsa04014  195/5470 8.665390e-08   195
    ## hsa05211   66/5470 9.494061e-08    66
    ## hsa05231   90/5470 9.494061e-08    90
    ## hsa05222   85/5470 1.110732e-07    85
    ## hsa04070   89/5470 1.241726e-07    89
    ## hsa05223   65/5470 1.259253e-07    65
    ## hsa04310  139/5470 1.293594e-07   139
    ## hsa00310   59/5470 1.383562e-07    59
    ## hsa04024  182/5470 1.559527e-07   182
    ## hsa05226  130/5470 1.969030e-07   130
    ## hsa05161  140/5470 2.274699e-07   140
    ## hsa05418  122/5470 2.301252e-07   122
    ## hsa04211   82/5470 2.410793e-07    82
    ## hsa04935  106/5470 2.878179e-07   106
    ## hsa05016  249/5470 3.150037e-07   249
    ## hsa04932  130/5470 4.425572e-07   130
    ## hsa05166  183/5470 4.512829e-07   183
    ## hsa04210  119/5470 4.596777e-07   119
    ## hsa04713   88/5470 4.692212e-07    88
    ## hsa05215   88/5470 4.692212e-07    88
    ## hsa05214   70/5470 5.315059e-07    70
    ## hsa05414   87/5470 6.088716e-07    87
    ## hsa05135  114/5470 6.088716e-07   114
    ## hsa04072  128/5470 6.489053e-07   128
    ## hsa04933   90/5470 8.010026e-07    90
    ## hsa04727   81/5470 1.001138e-06    81
    ## hsa04720   63/5470 1.007886e-06    63
    ## hsa01522   88/5470 1.366216e-06    88
    ## hsa05130  161/5470 1.435851e-06   161
    ## hsa04151  282/5470 1.608055e-06   282
    ## hsa04152  105/5470 2.227197e-06   105
    ## hsa05032   82/5470 2.308007e-06    82
    ## hsa04371  118/5470 2.904609e-06   118
    ## hsa05410   81/5470 3.005387e-06    81
    ## hsa04540   79/5470 5.272929e-06    79
    ## hsa05131  197/5470 5.516264e-06   197
    ## hsa04066   95/5470 1.076934e-05    95
    ## hsa03410   33/5470 1.076934e-05    33
    ## hsa04530  135/5470 1.147211e-05   135
    ## hsa04929   59/5470 1.369216e-05    59
    ## hsa05163  183/5470 1.381735e-05   183
    ## hsa05235   79/5470 1.442248e-05    79
    ## hsa03040  126/5470 1.932532e-05   126
    ## hsa04137   62/5470 1.932532e-05    62
    ## hsa04146   74/5470 1.932532e-05    74
    ## hsa04611  106/5470 2.003110e-05   106
    ## hsa05221   61/5470 2.538139e-05    61
    ## hsa03015   84/5470 2.577555e-05    84
    ## hsa00410   30/5470 3.053622e-05    30
    ## hsa00510   47/5470 3.071742e-05    47
    ## hsa04550  120/5470 3.284686e-05   120
    ## hsa04625   90/5470 3.284686e-05    90
    ## hsa04931   93/5470 3.652882e-05    93
    ## hsa04730   55/5470 4.058152e-05    55
    ## hsa05012  199/5470 4.377298e-05   199
    ## hsa04925   85/5470 4.599307e-05    85
    ## hsa03030   35/5470 4.965547e-05    35
    ## hsa04666   81/5470 5.070773e-05    81
    ## hsa03018   70/5470 5.070773e-05    70
    ## hsa04115   65/5470 7.010672e-05    65
    ## hsa05217   57/5470 7.240444e-05    57
    ## hsa04670   96/5470 8.218542e-05    96
    ## hsa04213   56/5470 9.529113e-05    56
    ## hsa03050   43/5470 9.953490e-05    43
    ## hsa04350   81/5470 1.149245e-04    81
    ## hsa04714  184/5470 1.208752e-04   184
    ## hsa05230   62/5470 1.505989e-04    62
    ## hsa04916   86/5470 1.712561e-04    86
    ## hsa04136   31/5470 1.880664e-04    31
    ## hsa05031   61/5470 1.927321e-04    61
    ## hsa04370   53/5470 2.118830e-04    53
    ## hsa05100   64/5470 2.266413e-04    64
    ## hsa00280   44/5470 2.696187e-04    44
    ## hsa05216   35/5470 2.696187e-04    35
    ## hsa05218   63/5470 2.875413e-04    63
    ## hsa00600   43/5470 3.553366e-04    43
    ## hsa03420   43/5470 3.553366e-04    43
    ## hsa04912   79/5470 3.738266e-04    79
    ## hsa04512   75/5470 4.358416e-04    75
    ## hsa01200   97/5470 4.723947e-04    97
    ## hsa04260   74/5470 5.430607e-04    74
    ## hsa04114  105/5470 5.896365e-04   105
    ## hsa01524   63/5470 6.628563e-04    63
    ## hsa05170  167/5470 6.854990e-04   167
    ## hsa03460   48/5470 7.782045e-04    48
    ## hsa05030   44/5470 8.072640e-04    44
    ## hsa04928   88/5470 8.072640e-04    88
    ## hsa00230  106/5470 8.332614e-04   106
    ## hsa00520   43/5470 1.058063e-03    43
    ## hsa04668   92/5470 1.222815e-03    92
    ## hsa05167  149/5470 1.282229e-03   149
    ## hsa04917   60/5470 1.287594e-03    60
    ## hsa05120   60/5470 1.287594e-03    60
    ## hsa03010  126/5470 1.404199e-03   126
    ## hsa04914   82/5470 1.404199e-03    82
    ## hsa04911   72/5470 1.583327e-03    72
    ## hsa04920   59/5470 1.583327e-03    59
    ## hsa04924   59/5470 1.583327e-03    59
    ## hsa03020   29/5470 1.635555e-03    29
    ## hsa00564   81/5470 1.670801e-03    81
    ## hsa05142   84/5470 1.676912e-03    84
    ## hsa00563   24/5470 1.754793e-03    24
    ## hsa03440   37/5470 1.767923e-03    37
    ## hsa04340   44/5470 2.030800e-03    44
    ## hsa04380  103/5470 2.246031e-03   103
    ## hsa00534   23/5470 2.394510e-03    23
    ## hsa04064   85/5470 2.394510e-03    85
    ## hsa04660   85/5470 2.394510e-03    85
    ## hsa01212   49/5470 3.274764e-03    49
    ## hsa03430   22/5470 3.347247e-03    22
    ## hsa03013  144/5470 3.686792e-03   144
    ## hsa00514   41/5470 4.273131e-03    41
    ## hsa00900   21/5470 4.645256e-03    21
    ## hsa04971   63/5470 5.059165e-03    63
    ## hsa00062   25/5470 5.354706e-03    25
    ## hsa04930   40/5470 5.387549e-03    40
    ## hsa04215   29/5470 5.387549e-03    29
    ## hsa05110   43/5470 5.900793e-03    43
    ## hsa00604   15/5470 5.900793e-03    15
    ## hsa03022   39/5470 6.823994e-03    39
    ## hsa04710   28/5470 7.075259e-03    28
    ## hsa04721   64/5470 7.343719e-03    64
    ## hsa05033   35/5470 7.693221e-03    35
    ## hsa04961   45/5470 7.851707e-03    45
    ## hsa04962   38/5470 8.544287e-03    38
    ## hsa04927   54/5470 8.637914e-03    54
    ## hsa00513   34/5470 9.784347e-03    34
    ## hsa04923   47/5470 1.011688e-02    47
    ## hsa04922   84/5470 1.222768e-02    84
    ## hsa04270  105/5470 1.368496e-02   105
    ## hsa00051   29/5470 1.419284e-02    29
    ## hsa04130   29/5470 1.419284e-02    29
    ## hsa05160  120/5470 1.518834e-02   120
    ## hsa00511   17/5470 1.611336e-02    17
    ## hsa03060   21/5470 1.664469e-02    21
    ## hsa05219   35/5470 1.690563e-02    35
    ## hsa04750   79/5470 1.700128e-02    79
    ## hsa00561   50/5470 1.859492e-02    50
    ## hsa00480   47/5470 1.859492e-02    47
    ## hsa00250   31/5470 1.969674e-02    31
    ## hsa05162  107/5470 2.072835e-02   107
    ## hsa00450   16/5470 2.159338e-02    16
    ## hsa05169  152/5470 2.319516e-02   152
    ## hsa00760   30/5470 2.476738e-02    30
    ## hsa04662   65/5470 2.684166e-02    65
    ## hsa00020   26/5470 2.916846e-02    26
    ## hsa00640   29/5470 3.119883e-02    29
    ## hsa00920   10/5470 3.475489e-02    10
    ## hsa04915  105/5470 3.690150e-02   105
    ## hsa04392   25/5470 3.696771e-02    25
    ## hsa00532   18/5470 3.779238e-02    18
    ## hsa00270   40/5470 3.989194e-02    40
    ## hsa04216   34/5470 4.021439e-02    34
    ## hsa04664   54/5470 4.086349e-02    54
    ## hsa05145   86/5470 4.086349e-02    86
    ## hsa04062  141/5470 4.295596e-02   141
    ## hsa00260   33/5470 4.893932e-02    33
    ## hsa00330   40/5470 6.734647e-02    40
    ## hsa00240   45/5470 7.336348e-02    45
    ## hsa00630   25/5470 7.511240e-02    25
    ## hsa04918   58/5470 7.656955e-02    58
    ## hsa01230   57/5470 8.847168e-02    57
    ## hsa04726   86/5470 1.041140e-01    86
    ## hsa05133   58/5470 1.098573e-01    58
    ## hsa00670   17/5470 1.155739e-01    17
    ## hsa05134   44/5470 1.286890e-01    44
    ## hsa00380   33/5470 1.387611e-01    33
    ## hsa01040   22/5470 1.393705e-01    22
    ## hsa00515   19/5470 1.429647e-01    19
    ## hsa00531   16/5470 1.443847e-01    16
    ## hsa05203  147/5470 1.712561e-01   147
    ## hsa04960   29/5470 1.712561e-01    29
    ## hsa00340   18/5470 1.750531e-01    18
    ## hsa00072    9/5470 1.754329e-01     9
    ## hsa01210   15/5470 1.792648e-01    15
    ## hsa04979   38/5470 2.053163e-01    38
    ## hsa00770   17/5470 2.129537e-01    17
    ## hsa00620   30/5470 2.170111e-01    30
    ## hsa04620   76/5470 2.202421e-01    76
    ## hsa03450   11/5470 2.295187e-01    11
    ## hsa00052   24/5470 2.469326e-01    24
    ## hsa01523   24/5470 2.469326e-01    24
    ## hsa02010   34/5470 2.469326e-01    34
    ## hsa00071   33/5470 2.824863e-01    33
    ## hsa05020   30/5470 3.041976e-01    30
    ## hsa00512   24/5470 3.572316e-01    24
    ## hsa00130    9/5470 3.613955e-01     9
    ## hsa00430    9/5470 3.613955e-01     9
    ## hsa00061   14/5470 3.666021e-01    14
    ## hsa00190   94/5470 3.856204e-01    94
    ## hsa04659   76/5470 3.857250e-01    76
    ## hsa00533   11/5470 3.996662e-01    11
    ## hsa00010   48/5470 4.134285e-01    48
    ## hsa05202  134/5470 4.311563e-01   134
    ## hsa04966   20/5470 4.338099e-01    20
    ## hsa00030   22/5470 4.509175e-01    22
    ## hsa04964   17/5470 4.740819e-01    17
    ## hsa05152  125/5470 4.834962e-01   125
    ## hsa00790   19/5470 4.891952e-01    19
    ## hsa00360   12/5470 4.977540e-01    12
    ## hsa04621  125/5470 5.329923e-01   125
    ## hsa04658   64/5470 5.431435e-01    64
    ## hsa00650   20/5470 5.608483e-01    20
    ## hsa00603   11/5470 5.746135e-01    11
    ## hsa00220   15/5470 6.013099e-01    15
    ## hsa04970   64/5470 6.119835e-01    64
    ## hsa00601   19/5470 6.174098e-01    19
    ## hsa00100   14/5470 6.706144e-01    14
    ## hsa05146   69/5470 7.251794e-01    69
    ## hsa00565   33/5470 7.590887e-01    33
    ## hsa04977   16/5470 8.102756e-01    16
    ## hsa03320   51/5470 8.445196e-01    51
    ## hsa04973   31/5470 8.445196e-01    31
    ## hsa00350   23/5470 8.486908e-01    23
    ## hsa04614   15/5470 8.634762e-01    15
    ## hsa04622   46/5470 8.683621e-01    46
    ## hsa04623   41/5470 9.032477e-01    41
    ## hsa04974   62/5470 9.184136e-01    62
    ## hsa05140   50/5470 9.215486e-01    50
    ## hsa00730   10/5470 9.353526e-01    10
    ## hsa04657   61/5470 9.353526e-01    61
    ## hsa03008   72/5470 9.490323e-01    72
    ## hsa05340   24/5470 9.490323e-01    24
    ## hsa04145   98/5470 1.000000e+00    98
    ## hsa00120   10/5470 1.000000e+00    10
    ## hsa05164  110/5470 1.000000e+00   110
    ## hsa04972   64/5470 1.000000e+00    64
    ## hsa04913   31/5470 1.000000e+00    31
    ## hsa04978   36/5470 1.000000e+00    36
    ## hsa05144   30/5470 1.000000e+00    30
    ## hsa05416   36/5470 1.000000e+00    36
    ## hsa04630  102/5470 1.000000e+00   102
    ## hsa04217  100/5470 1.000000e+00   100
    ## hsa05143   21/5470 1.000000e+00    21
    ## hsa00983   46/5470 1.000000e+00    46
    ## hsa00500   19/5470 1.000000e+00    19
    ## hsa04514   89/5470 1.000000e+00    89
    ## hsa05323   54/5470 1.000000e+00    54
    ## hsa00860   22/5470 1.000000e+00    22
    ## hsa04975   22/5470 1.000000e+00    22
    ## hsa04080  210/5470 1.000000e+00   210
    ## hsa04744   12/5470 1.000000e+00    12
    ## hsa00590   32/5470 1.000000e+00    32
    ## hsa00592   10/5470 1.000000e+00    10
    ## hsa04650   70/5470 1.000000e+00    70
    ## hsa05321   30/5470 1.000000e+00    30
    ## hsa04672   21/5470 1.000000e+00    21
    ## hsa04742   42/5470 1.000000e+00    42
    ## hsa04610   41/5470 1.000000e+00    41
    ## hsa04976   43/5470 1.000000e+00    43
    ## hsa00040   12/5470 1.000000e+00    12
    ## hsa04940   16/5470 1.000000e+00    16
    ## hsa04950    7/5470 1.000000e+00     7
    ## hsa00910    3/5470 1.000000e+00     3
    ## hsa05034   96/5470 1.000000e+00    96
    ## hsa00591    7/5470 1.000000e+00     7
    ## hsa00053    6/5470 1.000000e+00     6
    ## hsa04640   43/5470 1.000000e+00    43
    ## hsa00970   25/5470 1.000000e+00    25
    ## hsa04612   31/5470 1.000000e+00    31
    ## hsa04061   43/5470 1.000000e+00    43
    ## hsa05204   32/5470 1.000000e+00    32
    ## hsa00982   26/5470 1.000000e+00    26
    ## hsa00980   27/5470 1.000000e+00    27
    ## hsa05330    7/5470 1.000000e+00     7
    ## hsa05310    3/5470 1.000000e+00     3
    ## hsa05332    7/5470 1.000000e+00     7
    ## hsa04060  141/5470 1.000000e+00   141
    ## hsa05320   10/5470 1.000000e+00    10
    ## hsa00140   13/5470 1.000000e+00    13
    ## hsa05206  147/5470 1.000000e+00   147
    ## hsa00830   16/5470 1.000000e+00    16
    ## hsa04740   35/5470 1.000000e+00    35
    ## hsa05150   26/5470 1.000000e+00    26
    ## hsa05168  135/5470 1.000000e+00   135
    ## hsa05322   20/5470 1.000000e+00    20

    #plot(kk2@result[,"setSize"])

    #library(GSEABase)
    #library(tidyr)

    #gene <- names(test)[abs(test) > 2]
    #wpgmtfile <- system.file("extdata", "wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
    #wp2gene <- read.gmt(wpgmtfile) #i had to install GSEAbase first ??
    #wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
    #wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
    #wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
    #ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
    #ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
    #ewp@result = ewp@result[order(ewp@result[,"Count"],decreasing=TRUE),]
    #ewp@result

    #barplot(ewp,showCategory=12)

    #ewp2 <- GSEA(test, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
    #ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
    #ewp2@result = ewp2@result[order(ewp2@result[,"setSize"],decreasing=TRUE),]
    #ewp2@result

    #gset@featureData@data[which(gset@featureData@data[,"Gene.ID"]==100128545),c("ID","Gene.ID")]
