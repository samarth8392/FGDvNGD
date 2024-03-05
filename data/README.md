# Datafiles and their descriptions

```bash
.
├── gene.500kb.updown.bed
├── geneRegions.bed
├── lower10_genes.bed.csv
├── sampleInfo.csv
├── sites
│   ├── geneDes.sites.csv
│   ├── lof.lower10.sites.csv
│   └── upper10.sites.csv
└── upper10_genes.bed.csv
```

1. **geneRegions.bed** : Genic regions i.e. regions that correspond to annotated genes in the *S. catenatus* genome.
2. **500kb.updown.genes.bed** : Genic regions +/- 500kb (i.e., geneic regions + 500kb upstream and downstream).
3. **lower10_geneInfo.csv** : Gene co-ordinates for genes in the lower 10% of the Direction of Selection (DoS) distribution.
4. **upper10_geneInfo.csv** : Gene co-ordinates for genes in the upper 10% of the Direction of Selection (DoS) distribution.
5. **sampleInfo.csv** : Location and mapping info for each sample analyzed in the study (DOC = Depth of coverage; BOC_x = Breadth of coverage at x depth of coverage).
6. **geneDes.sites.csv** : SNPs identified within gene desert regions of the *S. catenatus* genome.
7. **lof.lower10.sites.csv** : Non-synonymous (missense + non-sense i.e. Loss of Function LOF) SNPs identified within genes in the lower 10% of DOS distribution.
8. **upper10.sites.csv** : Non-synonymous SNPs identified within genes in the upper 10% of DOS distribution.