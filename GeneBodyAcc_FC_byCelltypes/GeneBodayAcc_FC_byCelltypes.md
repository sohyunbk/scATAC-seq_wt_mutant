Getting FC (Bif3/WT) of gene body chromatin acc
================
Sohyun Bang
27 May, 2025

    ## Loading required package: limma

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: grid

    ## Loading required package: futile.logger

### 1) Get Gene body acc for all the genes & normalize the value & calculate logFC

``` r
### 1) Get Gene body acc for all the genes
## load gene*cell table.
WT_GeneXCT <- read.table("./A619_AnnV4.GeneBodyACC.byGeneXCT.txt",header=TRUE)
Bif3_GeneXCT <- read.table("./Bif3_AnnV4.GeneBodyACC.byGeneXCT.txt",header=TRUE)
Celltypes <- colnames(WT_GeneXCT)[-1]
colnames(Bif3_GeneXCT)[-1] <- paste("Bif3&", colnames(Bif3_GeneXCT)[-1], sep="")
colnames(WT_GeneXCT)[-1] <- paste("A619&", colnames(WT_GeneXCT)[-1], sep="")

GeneXCT <- merge(WT_GeneXCT, Bif3_GeneXCT, by = "gene")
head(GeneXCT)
```

    ##              gene A619&FloralMeristem_SuppressedBract A619&G2_M A619&IM.OC
    ## 1 Zm00001eb000010                                 178       188        232
    ## 2 Zm00001eb000020                                 463       715        547
    ## 3 Zm00001eb000050                                  18        30         31
    ## 4 Zm00001eb000060                                 210       215        256
    ## 5 Zm00001eb000070                                  53        68         85
    ## 6 Zm00001eb000080                                 535       484        507
    ##   A619&L1 A619&L1atFloralMeristem A619&PhloemPrecursor
    ## 1     168                     414                  448
    ## 2     402                     930                 1649
    ## 3      22                      58                   48
    ## 4     185                     364                  463
    ## 5      47                     137                  155
    ## 6     479                     927                 1175
    ##   A619&ProcambialMeristem_ProtoXylem_MetaXylem
    ## 1                                          207
    ## 2                                          771
    ## 3                                           34
    ## 4                                          254
    ## 5                                           61
    ## 6                                          620
    ##   A619&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
    ## 1                                                         71
    ## 2                                                        238
    ## 3                                                         29
    ## 4                                                         93
    ## 5                                                         40
    ## 6                                                        241
    ##   A619&SPM.base_SM.base A619&Unknown_lowFRiP A619&Unknown_Sclerenchyma
    ## 1                   179                  347                       287
    ## 2                   402                  669                       577
    ## 3                    36                   66                        26
    ## 4                   180                  251                       267
    ## 5                    45                  167                       100
    ## 6                   462                  593                       720
    ##   A619&Unknown1 A619&Unknown2 A619&XylemParenchyma_PithParenchyma
    ## 1           272           248                                 263
    ## 2           589           435                                 906
    ## 3            27            29                                  46
    ## 4           228           181                                 355
    ## 5            87            93                                  90
    ## 6           740           565                                 913
    ##   Bif3&FloralMeristem_SuppressedBract Bif3&G2_M Bif3&IM.OC Bif3&L1
    ## 1                                 160        89         40      95
    ## 2                                 322       307        158     297
    ## 3                                  21        11          6      12
    ## 4                                 136        86         52     115
    ## 5                                  37        25         39      56
    ## 6                                 394       271        146     268
    ##   Bif3&L1atFloralMeristem Bif3&PhloemPrecursor
    ## 1                     148                  260
    ## 2                     455                  946
    ## 3                      34                   30
    ## 4                     147                  237
    ## 5                      85                   87
    ## 6                     411                  621
    ##   Bif3&ProcambialMeristem_ProtoXylem_MetaXylem
    ## 1                                          133
    ## 2                                          379
    ## 3                                           17
    ## 4                                          131
    ## 5                                           44
    ## 6                                          281
    ##   Bif3&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
    ## 1                                                         50
    ## 2                                                        171
    ## 3                                                          8
    ## 4                                                         55
    ## 5                                                         22
    ## 6                                                        173
    ##   Bif3&SPM.base_SM.base Bif3&Unknown_lowFRiP Bif3&Unknown_Sclerenchyma
    ## 1                    74                  137                       159
    ## 2                   210                  334                       321
    ## 3                     4                   23                        11
    ## 4                    70                  106                       156
    ## 5                    24                   62                        66
    ## 6                   255                  272                       425
    ##   Bif3&Unknown1 Bif3&Unknown2 Bif3&XylemParenchyma_PithParenchyma
    ## 1           148           140                                 129
    ## 2           349           294                                 401
    ## 3            15             3                                  27
    ## 4           135           102                                 146
    ## 5            39            41                                  48
    ## 6           377           314                                 356

``` r
## normalization
GeneNames <- GeneXCT[, 1]
gene_counts_df <- GeneXCT[, -1]

## Set up the cut off to 50!
#head(GeneXCT)
#dim(GeneXCT)

GeneXCT_MoreThan50Tn5 <- GeneXCT %>%
  filter_all(all_vars(. > 50))
dim(GeneXCT_MoreThan50Tn5)
```

    ## [1] 18113    29

``` r
GeneNames <- GeneXCT[, 1]
gene_counts_df <- GeneXCT[, -1]

cpm_data <- cpm(DGEList(counts = gene_counts_df), log = FALSE)
QNorm <- normalize.quantiles(as.matrix(gene_counts_df))
rownames(QNorm) <- GeneNames
colnames(QNorm) <- colnames(GeneXCT[,-1])
#head(QNorm)
#dim(QNorm)

# Build fold change table
FCTable <- as.data.frame(sapply(Celltypes, function(Celltypes) {
  bif3_col <- paste0("Bif3&", Celltypes)
  A619_col <- paste0("A619&", Celltypes)
  log2(QNorm[, bif3_col] / QNorm[, A619_col])
}))
```

### 2) Matching with gene names

``` r
#### 2) Find ARF genes
GeneInfo <- read.delim("../maizev5_data/Zm00001eb.1.fulldata_Curated2.txt", stringsAsFactors = FALSE)

FCTable_with_id <- FCTable %>%
  tibble::rownames_to_column(var = "gene_model")

FCTable_annotated <- FCTable_with_id %>%
  left_join(GeneInfo[, c("gene_model", "locus_symbol")], by = "gene_model")
write.table(FCTable_annotated, file = "GeneBodyACC_FC_Bif3WT_byCelltypes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#FCTable_annotated[FCTable_annotated$gene_model == "Zm00001eb067310",]
#FCTable_annotated[FCTable_annotated$gene_model == "Zm00001eb168120",]
```

### 3) Check Homemodomain TFs

``` r
TFlist <- read.delim("./Zma_TF_list.txt", stringsAsFactors = FALSE) ## This data is from https://planttfdb.gao-lab.org/index.php?sp=Zma

## Filter TFs with HD
#HD_TFs <- TFlist %>%
#  filter(Family %in% c("ZF-HD", "HD-ZIP", "WOX", "HB0PHD", "HB-other"))

V3_V5 <- read.delim("B73v3_to_B73v5.tsv", stringsAsFactors = FALSE)

HD_TFs_merged <- merge(TFlist, V3_V5, by.x = "Gene_ID", by.y = "V3", all.x = TRUE)
HD_TFs_expanded <- HD_TFs_merged %>% ## This is because to make V5 gene id as the key
  separate_rows(V5, sep = ",") %>%
  rename(gene_model=V5) 

HD_TFs_V5Key <- HD_TFs_expanded %>%
  group_by(gene_model) %>%
  summarise(Family = paste(sort(unique(Family)), collapse = ","), .groups = "drop")

FCTable_annotated_TFInfo <- FCTable_annotated %>%
  left_join(HD_TFs_V5Key[, c("gene_model","Family")], by = "gene_model") %>%
  rename(TF_Family = Family)

FCTable_with_HD <- FCTable_annotated_TFInfo %>%
  filter(TF_Family %in% c("ZF-HD", "HD-ZIP", "WOX", "HB0PHD", "HB-other") )

FCTable_with_HD[FCTable_with_HD$gene_model=="Zm00001eb433010",] #ZmWUS2
```

    ##          gene_model FloralMeristem_SuppressedBract     G2_M      IM.OC
    ## 147 Zm00001eb433010                     -0.2269246 0.288721 -0.6049317
    ##            L1 L1atFloralMeristem PhloemPrecursor
    ## 147 0.2720308          0.1076469     -0.06974107
    ##     ProcambialMeristem_ProtoXylem_MetaXylem
    ## 147                             -0.08958879
    ##     ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM.base_SM.base
    ## 147                                            0.03539159        0.1504664
    ##     Unknown_lowFRiP Unknown_Sclerenchyma   Unknown1  Unknown2
    ## 147      -0.6020971            0.4151011 -0.2844753 0.1940037
    ##     XylemParenchyma_PithParenchyma locus_symbol TF_Family
    ## 147                     0.06770149        hb122       WOX

``` r
FCTable_with_HD[FCTable_with_HD$gene_model=="Zm00001eb067310",] #ZmWUS1
```

    ##         gene_model FloralMeristem_SuppressedBract     G2_M    IM.OC       L1
    ## 26 Zm00001eb067310                      0.8477587 1.204014 0.549141 1.105983
    ##    L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem
    ## 26           1.140301        0.720564                               0.5548946
    ##    ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM.base_SM.base
    ## 26                                             0.9063194        0.8340232
    ##    Unknown_lowFRiP Unknown_Sclerenchyma Unknown1 Unknown2
    ## 26        1.234771             1.266715 1.194005  1.14416
    ##    XylemParenchyma_PithParenchyma locus_symbol TF_Family
    ## 26                       1.079582         hb67       WOX

``` r
Up_inCentralZone <- FCTable_with_HD %>% filter(IM.OC > 0.5)
print(table(Up_inCentralZone$TF_Family))
```

    ## 
    ## HB-other   HD-ZIP      WOX    ZF-HD 
    ##        6       17        4        4

``` r
Up_inCentralZone[Up_inCentralZone$TF_Family=="WOX",] #ZmWUS2
```

    ##         gene_model FloralMeristem_SuppressedBract       G2_M     IM.OC
    ## 3  Zm00001eb015500                     0.05723111  0.8246329 1.2689584
    ## 9  Zm00001eb067310                     0.84775870  1.2040139 0.5491410
    ## 26 Zm00001eb395430                     0.23895627 -0.4362021 0.6982909
    ## 31 Zm00001eb414580                    -1.25096157 -0.9382282 0.5095710
    ##            L1 L1atFloralMeristem PhloemPrecursor
    ## 3  0.04484057         -0.2751072       0.1874470
    ## 9  1.10598340          1.1403009       0.7205640
    ## 26 0.14874899          0.3405043       0.2655402
    ## 31 0.22309566         -0.7628538       0.1715078
    ##    ProcambialMeristem_ProtoXylem_MetaXylem
    ## 3                                0.5819988
    ## 9                                0.5548946
    ## 26                              -0.2956852
    ## 31                               1.4364298
    ##    ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM.base_SM.base
    ## 3                                              0.3687242       0.07569172
    ## 9                                              0.9063194       0.83402323
    ## 26                                             0.4593562       0.47633726
    ## 31                                            -0.3251093      -0.29805212
    ##    Unknown_lowFRiP Unknown_Sclerenchyma   Unknown1   Unknown2
    ## 3      -0.03982848            0.9662018 0.07241811 -0.5539623
    ## 9       1.23477142            1.2667153 1.19400475  1.1441598
    ## 26      0.11905946            0.6087666 0.45239175  0.6134752
    ## 31      1.06504843            0.6551089 0.20256307 -0.2332745
    ##    XylemParenchyma_PithParenchyma locus_symbol TF_Family
    ## 3                      0.19961669        hb130       WOX
    ## 9                      1.07958199         hb67       WOX
    ## 26                     0.09038908        wox11       WOX
    ## 31                    -0.37417008         hb73       WOX

### 4) Check De novo marker gene list

``` r
DeNovo <- read.delim("A619_v4.IM.OC_deseq_2_results.tsv", stringsAsFactors = FALSE)
Sig_DeNovo <- DeNovo[!is.na(DeNovo$padj) & 
                     !is.na(DeNovo$log2FoldChange) & 
                     DeNovo$padj < 0.01 & 
                     DeNovo$log2FoldChange > 0, ]


Sig_DeNovo$gene_name
```

    ##   [1] "Zm00001eb197040" "Zm00001eb400120" "Zm00001eb005740" "Zm00001eb357220"
    ##   [5] "Zm00001eb257620" "Zm00001eb160120" "Zm00001eb316060" "Zm00001eb093680"
    ##   [9] "Zm00001eb103420" "Zm00001eb072020" "Zm00001eb250240" "Zm00001eb295920"
    ##  [13] "Zm00001eb001630" "Zm00001eb154380" "Zm00001eb430150" "Zm00001eb319560"
    ##  [17] "Zm00001eb062570" "Zm00001eb207800" "Zm00001eb210980" "Zm00001eb344540"
    ##  [21] "Zm00001eb101880" "Zm00001eb361790" "Zm00001eb290630" "Zm00001eb284010"
    ##  [25] "Zm00001eb275610" "Zm00001eb314970" "Zm00001eb109730" "Zm00001eb427970"
    ##  [29] "Zm00001eb243560" "Zm00001eb387910" "Zm00001eb197550" "Zm00001eb318500"
    ##  [33] "Zm00001eb375830" "Zm00001eb196470" "Zm00001eb157360" "Zm00001eb400100"
    ##  [37] "Zm00001eb196460" "Zm00001eb230580" "Zm00001eb425590" "Zm00001eb067310"
    ##  [41] "Zm00001eb104430" "Zm00001eb332090" "Zm00001eb315610" "Zm00001eb234860"
    ##  [45] "Zm00001eb152810" "Zm00001eb250710" "Zm00001eb122620" "Zm00001eb397140"
    ##  [49] "Zm00001eb276980" "Zm00001eb199780" "Zm00001eb398820" "Zm00001eb368350"
    ##  [53] "Zm00001eb117820" "Zm00001eb210000" "Zm00001eb159490" "Zm00001eb251540"
    ##  [57] "Zm00001eb008690" "Zm00001eb076810" "Zm00001eb411430" "Zm00001eb376530"
    ##  [61] "Zm00001eb081230" "Zm00001eb362780" "Zm00001eb354880" "Zm00001eb006480"
    ##  [65] "Zm00001eb161910" "Zm00001eb106770" "Zm00001eb331060" "Zm00001eb331070"
    ##  [69] "Zm00001eb359750" "Zm00001eb368050" "Zm00001eb161780" "Zm00001eb005600"
    ##  [73] "Zm00001eb297040" "Zm00001eb362740" "Zm00001eb202660" "Zm00001eb430440"
    ##  [77] "Zm00001eb317770" "Zm00001eb254530" "Zm00001eb366880" "Zm00001eb364700"
    ##  [81] "Zm00001eb073830" "Zm00001eb296480" "Zm00001eb035180" "Zm00001eb003910"
    ##  [85] "Zm00001eb175400" "Zm00001eb313660" "Zm00001eb257410" "Zm00001eb433010"
    ##  [89] "Zm00001eb005980" "Zm00001eb155610" "Zm00001eb258830" "Zm00001eb213820"
    ##  [93] "Zm00001eb256500" "Zm00001eb274800" "Zm00001eb337060" "Zm00001eb401070"
    ##  [97] "Zm00001eb254360" "Zm00001eb071800" "Zm00001eb156040" "Zm00001eb368180"
    ## [101] "Zm00001eb063300" "Zm00001eb395970" "Zm00001eb207690" "Zm00001eb368080"
    ## [105] "Zm00001eb013470" "Zm00001eb101640" "Zm00001eb337070" "Zm00001eb107130"
    ## [109] "Zm00001eb005110" "Zm00001eb246060" "Zm00001eb169360" "Zm00001eb426920"
    ## [113] "Zm00001eb365880" "Zm00001eb272500" "Zm00001eb192560" "Zm00001eb406950"
    ## [117] "Zm00001eb051010" "Zm00001eb069500" "Zm00001eb148390" "Zm00001eb401510"
    ## [121] "Zm00001eb016240" "Zm00001eb119540" "Zm00001eb359110" "Zm00001eb374420"
    ## [125] "Zm00001eb159400" "Zm00001eb374410" "Zm00001eb408280" "Zm00001eb103980"
    ## [129] "Zm00001eb152630" "Zm00001eb283680" "Zm00001eb396640" "Zm00001eb000450"
    ## [133] "Zm00001eb405340" "Zm00001eb434450" "Zm00001eb033960" "Zm00001eb051770"
    ## [137] "Zm00001eb266360" "Zm00001eb288270" "Zm00001eb294390" "Zm00001eb002310"
    ## [141] "Zm00001eb074370" "Zm00001eb209120" "Zm00001eb362790" "Zm00001eb000680"
    ## [145] "Zm00001eb000550" "Zm00001eb234600" "Zm00001eb359600" "Zm00001eb208730"
    ## [149] "Zm00001eb398190" "Zm00001eb338800" "Zm00001eb286760" "Zm00001eb238570"
    ## [153] "Zm00001eb419690" "Zm00001eb390520" "Zm00001eb367160" "Zm00001eb283520"
    ## [157] "Zm00001eb118150" "Zm00001eb384620" "Zm00001eb404830" "Zm00001eb434390"
    ## [161] "Zm00001eb319530" "Zm00001eb002950" "Zm00001eb318630" "Zm00001eb297050"
    ## [165] "Zm00001eb078070" "Zm00001eb293690" "Zm00001eb254070" "Zm00001eb360940"
    ## [169] "Zm00001eb405460" "Zm00001eb103460" "Zm00001eb376160" "Zm00001eb005890"
    ## [173] "Zm00001eb058020" "Zm00001eb149000" "Zm00001eb225480" "Zm00001eb065750"
    ## [177] "Zm00001eb331370" "Zm00001eb361450" "Zm00001eb073920" "Zm00001eb108270"
    ## [181] "Zm00001eb215460" "Zm00001eb363360" "Zm00001eb419240" "Zm00001eb157130"
    ## [185] "Zm00001eb195050" "Zm00001eb316720" "Zm00001eb088370" "Zm00001eb168120"
    ## [189] "Zm00001eb432460" "Zm00001eb323660" "Zm00001eb000430" "Zm00001eb295120"
    ## [193] "Zm00001eb387850" "Zm00001eb400420" "Zm00001eb294970" "Zm00001eb359580"
    ## [197] "Zm00001eb336880" "Zm00001eb127880" "Zm00001eb364310" "Zm00001eb167340"
    ## [201] "Zm00001eb227800" "Zm00001eb367750" "Zm00001eb402200" "Zm00001eb333620"
    ## [205] "Zm00001eb188100" "Zm00001eb207810" "Zm00001eb106300" "Zm00001eb108640"
    ## [209] "Zm00001eb367740" "Zm00001eb227630" "Zm00001eb068410" "Zm00001eb121150"
    ## [213] "Zm00001eb373110" "Zm00001eb224940" "Zm00001eb369590" "Zm00001eb401270"
    ## [217] "Zm00001eb359500" "Zm00001eb401150" "Zm00001eb220660" "Zm00001eb259340"
    ## [221] "Zm00001eb254920" "Zm00001eb336930" "Zm00001eb257770" "Zm00001eb002480"
    ## [225] "Zm00001eb147990" "Zm00001eb406640" "Zm00001eb053740" "Zm00001eb403460"
    ## [229] "Zm00001eb296770" "Zm00001eb167840" "Zm00001eb150420" "Zm00001eb008680"
    ## [233] "Zm00001eb294780" "Zm00001eb102450" "Zm00001eb156220" "Zm00001eb193850"
    ## [237] "Zm00001eb293310" "Zm00001eb072720" "Zm00001eb250800" "Zm00001eb257510"
    ## [241] "Zm00001eb068860" "Zm00001eb432330" "Zm00001eb369280" "Zm00001eb405320"
    ## [245] "Zm00001eb003630" "Zm00001eb291720" "Zm00001eb160780" "Zm00001eb266600"
    ## [249] "Zm00001eb294730" "Zm00001eb430130" "Zm00001eb375010" "Zm00001eb000440"
    ## [253] "Zm00001eb076800" "Zm00001eb107040" "Zm00001eb106390" "Zm00001eb157890"
    ## [257] "Zm00001eb076790" "Zm00001eb315600" "Zm00001eb276150" "Zm00001eb159310"
    ## [261] "Zm00001eb332960" "Zm00001eb101510" "Zm00001eb103470" "Zm00001eb424870"
    ## [265] "Zm00001eb325180" "Zm00001eb286430" "Zm00001eb327220" "Zm00001eb216020"
    ## [269] "Zm00001eb189550" "Zm00001eb002040" "Zm00001eb211640" "Zm00001eb333220"
    ## [273] "Zm00001eb258370" "Zm00001eb082140" "Zm00001eb330340" "Zm00001eb031250"
    ## [277] "Zm00001eb356760" "Zm00001eb174100" "Zm00001eb403250" "Zm00001eb209220"
    ## [281] "Zm00001eb403240" "Zm00001eb105840" "Zm00001eb158740" "Zm00001eb382650"
    ## [285] "Zm00001eb328750" "Zm00001eb428880" "Zm00001eb219190" "Zm00001eb188240"
    ## [289] "Zm00001eb401720" "Zm00001eb253260" "Zm00001eb257570" "Zm00001eb048560"

``` r
ven_list_all <- list(HD_TF_Up = Up_inCentralZone$gene_model, 
                     Central_Zone_Denovo_markergene = Sig_DeNovo$gene_name)
#print(ggVennDiagram(ven_list_all, label = "count",nintersects = 20))
intersect(Up_inCentralZone$gene_model, Sig_DeNovo$gene_name)
```

    ## [1] "Zm00001eb067310"

``` r
venn.diagram(ven_list_all,
        category.names = c("Set 1" , "Set 2"),
        filename = "VennDiagram_DeNovoCentralZone.pdf",
        imagetype="tiff" ,
        height = 150 , 
        width = 200 , 
        resolution = 300,
        compression = "lzw",
        lwd = 2,
        lty = 'blank',
        fill = c("red","blue"))
```

    ## Warning in tiff(filename = filename, height = height, width = width, units =
    ## units, : compression is not supported for type = "quartz"

    ## [1] 1
