library(tidyverse)
library(readxl)
library(ape)
library(clusterProfiler)
library(geiger)

#phylogenies
tree <- ladderize(read.tree('BackboneTree.nh')) #backbone dataset
tree <- ladderize(read.tree('SaccharomycetalesTree.nh')) #Saccharomycetales dataset
tree <- ladderize(read.tree('DipodascalesTree.nh')) #Dipodascales dataset

#get branch lengths
b.length <- data.frame(Node=c(tree$tip.label, tree$node.label)[tree$edge[,2]], Length=tree$edge.length)

#duplication events from orthofinder output
duplication_raw <- read.delim('BackboneDuplications.tsv', sep='\t', header = TRUE) #backbone dataset
duplication_raw <- read.delim('SaccharomycetalesDuplications.tsv', sep='\t', header = TRUE)  #Saccharomycetales dataset
duplication_raw <- read.delim('DipodascalesDuplications.tsv', sep='\t', header = TRUE)  #Dipodascales dataset

#get rate of duplication for each branch
duplications <- duplication_raw %>% filter(Support>=0.1)
duplications <- aggregate(Orthogroup ~ Species.Tree.Node, duplications, length) %>% arrange(-Orthogroup) %>% rename(Duplications = Orthogroup, Node = Species.Tree.Node)
duplications <- merge(b.length, duplications)
duplications$Rate <- duplications$Duplications/duplications$Length

#function to pull and format wgd syn output
process_multiplicons <- function(species) {
    
    #read in scaffold sizes
    species_sequences <- read.delim(paste0("wgd_syn/", species, "/scaffold_length.tsv")) %>% rename(seq_id = scaffold)
    
    #read in and format gff files
    species_gff <- read.gff(paste0("gff/", species, ".gff3"))
    species_gff$ID <- gsub(".*ID=([^;]+).*", "\\1", species_gff$attributes) 
    species_gff <- species_gff %>% rename(seq_id = seqid)
    species_gff <- species_gff %>% filter(type=='mRNA') %>% select(ID, seq_id, start, end, strand)
    
    #read in and format multiplicon pairs
    multiplicons <- read.delim(paste0("wgd_syn/", species, "/iadhore-out/multiplicon_pairs.txt")) 
    species_multiplicons <- data.frame(multiplicon = c(multiplicons$multiplicon, multiplicons$multiplicon), ID = c(multiplicons$X, multiplicons$gene_x))
    species_multiplicons <- merge(species_gff, species_multiplicons)
    
    #fill in multiplicon segments
    species_multiplicons$paired <- TRUE
    multiplicon_spans <- species_multiplicons %>%
      group_by(multiplicon, seq_id) %>%
      summarise(
        multiplicon_start = min(start),
        multiplicon_end = max(end),
        .groups = "drop"
      )
    species_extra_genes <- species_gff %>%
      inner_join(multiplicon_spans, by = "seq_id", relationship = "many-to-many") %>%
      filter(
        start >= multiplicon_start,
        end <= multiplicon_end
      ) %>%
      mutate(paired = FALSE) %>%
      select(ID, seq_id, start, end, strand, multiplicon, paired)
    species_segments <-  bind_rows(
          species_multiplicons,
          species_extra_genes
        ) %>%
          distinct(seq_id, start, end, .keep_all = TRUE) %>%
          arrange(multiplicon, seq_id, start)
    species_segments$species <- species
    species_segments$ID <- gsub("-T1$", "", species_segments$ID)

    return(list(seq=species_sequences, multiplicons=species_multiplicons, segments=species_segments))
    
    }

#pull multiplicon segments for each target genome
Gcan_segments <- process_multiplicons("Geotrichum_candidum_CLIB_918")$segments
Dfer_segments <- process_multiplicons("Dipodascus_fermentans")$segments
Ssua_segments <- process_multiplicons("Saprochaete_suaveolens")$segments
Mmag_segments <- process_multiplicons("Magnusiomyces_magnusii")$segments
Mtet_segments <- process_multiplicons("Magnusiomyces_tetraspermus")$segments

#add OG info to multiplicons, keep names consistent
OGs <- read.delim("DipodascalesOrthogroups.txt")

all_segments <- rbind(Gcan_segments, Dfer_segments, Ssua_segments, Mmag_segments, Mtet_segments) %>% 
  left_join(., OGs, by = c("ID" = "gene", "species"))

#compare OG membership across  multiplicons (Data S2)
#precompute the OG sets for each (species, multiplicon)
multiplicon_OGs <- all_segments %>%
  group_by(species, multiplicon) %>%
  summarise(OGs = list(unique(OG)), .groups = "drop")

# all pairs of multiplicons from different species
pairs <- expand_grid(
  idx1 = 1:nrow(multiplicon_OGs),
  idx2 = 1:nrow(multiplicon_OGs)
) %>%
  filter(idx1 < idx2) %>%
  rowwise() %>%
  mutate(
    species1 = multiplicon_OGs$species[idx1],
    multiplicon1 = multiplicon_OGs$multiplicon[idx1],
    OGs1 = list(multiplicon_OGs$OGs[[idx1]]),
    species2 = multiplicon_OGs$species[idx2],
    multiplicon2 = multiplicon_OGs$multiplicon[idx2],
    OGs2 = list(multiplicon_OGs$OGs[[idx2]])
  ) %>%
  filter(species1 != species2) %>%
  mutate(
    intersection = length(intersect(unlist(OGs1), unlist(OGs2))),
    union = length(union(unlist(OGs1), unlist(OGs2))),
    jaccard = ifelse(union == 0, NA, intersection / union)
  ) %>%
  ungroup() %>%
  select(species1, multiplicon1, species2, multiplicon2, intersection, union, jaccard) %>%
  arrange(desc(jaccard))

#count multiplicon members from every contiguous genome
mult_stats <- function(x) {
mult <- read.delim(paste0("wgd_syn/", x, "/iadhore-out/multiplicon_pairs.txt"))
return(data.frame(species=x, n_multiplicon_pair=length(unique(c(mult$gene_x, mult$X)))))
    }
multiplicon_count <- lapply(list.files("wgd_syn/"), mult_stats)
multiplicon_count <- do.call(rbind, multiplicon_count)

out %>% arrange(-n_multiplicon_pair)

#pull all multiplicons from each post-WGD genome, by event (get ohnologs)
get_multiplicon_pairs <- function(sp) {
    multiplicon <- read.delim(paste0("wgd_syn/", sp, "/iadhore-out/multiplicon_pairs.txt"))
    gene <- paste0(c(multiplicon[[3]], multiplicon[[4]]), "|", sp)
    return(data.frame(gene))
    }

WGD1_syngenes <- lapply(c('Kazachstania_naganishii', 'Saccharomyces_arboricola', 'Saccharomyces_cerevisiae', 'Saccharomyces_eubayanus', 'Saccharomyces_jurei', 'Saccharomyces_kudriavzevii', 'Saccharomyces_mikatae', 'Saccharomyces_paradoxus', 'Saccharomyces_uvarum'), get_multiplicon_pairs)
WGD1_genes <- do.call(rbind, WGD1_syngenes)

WGD2_syngenes <- lapply(c("Geotrichum_candidum_CLIB_918", "Dipodascus_fermentans"), get_multiplicon_pairs)
WGD2_genes <- do.call(rbind, WGD2_syngenes)

WGD4_syngenes <- lapply(c('Magnusiomyces_magnusii'), get_multiplicon_pairs)
WGD4_genes <- do.call(rbind, WGD4_syngenes)

#for WGD4, only ohnologs that are also specific to M. magnusii
WGD4_orthogenes <- duplication_raw[duplication_raw$Species.Tree.Node == "Magnusiomyces_magnusii",]
WGD4_orthogenes <- gsub("Magnusiomyces_magnusii_", "", unlist(strsplit(c(WGD4_orthogenes$`Genes.1`, WGD4_orthogenes$`Genes.2`), ",")))
WGD4_orthogenes <- sub("(.+)\\|(FUN_\\d+)", "\\2-T1|\\1", WGD4_orthogenes)
WGD4_genes <- WGD4_syngenes %>% filter(WGD4_syngenes$gene %in% WGD4_orthogenes)

WGD3_syngenes <- lapply(c('Saprochaete_suaveolens', 'Magnusiomyces_magnusii', "Magnusiomyces_tetraspermus"), get_multiplicon_pairs)
WGD3_genes <- do.call(rbind, WGD3_genes)

#for WGD3, only ohnologs that are also NOT specific to M. magnusii
WGD3_genes <- WGD3_genes %>% filter(!(WGD3_genes$gene %in% WGD3_genes))

#KEGG annotations
KO <- read.delim('KOannotations.txt')

#KEGG enrichment analysis
kk <- enrichKEGG(gene = KO[KO$gene %in% WGD1_genes$gene,]$KO, #make sure you're using the right WGD gene set

                 #all post WGD1 genes
                 universe = KO[KO$species %in% c('Kazachstania_naganishii', 'Saccharomyces_arboricola', 'Saccharomyces_cerevisiae', 'Saccharomyces_eubayanus', 'Saccharomyces_jurei', 'Saccharomyces_kudriavzevii', 'Saccharomyces_mikatae', 'Saccharomyces_paradoxus', 'Saccharomyces_uvarum'),]$KO,

                 #all post WGD2 genes
                 #universe = KO[KO$species %in% c("Dipodascus_fermentans", "Geotrichum_candidum_CLIB_918"),]$KO,
                 
                 #all post WGD3 genes
                 #universe = KO[KO$species %in% c("Saprochaete_suaveolens", "Magnusiomyces_magnusii", "Magnusiomyces_tetraspermus"),]$KO,
                 
                 #all post WGD4 genes
                 #universe = KO[KO$species %in% c("Magnusiomyces_magnusii"),]$KO,
                 
                 organism = 'ko',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'fdr'
                    )

#InterPro annotations
IPR <- read.delim("IPRannotations.txt")

#InterPro enrichment analysis
ipr <- enricher(
                  gene = IPR[IPR$gene %in% WGD1_genes$gene,]$gene,
                    
                  #all post WGD1 genes
                  universe = IPR[IPR$species %in% c('Kazachstania_naganishii', 'Saccharomyces_arboricola', 'Saccharomyces_cerevisiae', 'Saccharomyces_eubayanus', 'Saccharomyces_jurei', 'Saccharomyces_kudriavzevii', 'Saccharomyces_mikatae', 'Saccharomyces_paradoxus', 'Saccharomyces_uvarum'),]$gene,
                
                  #all post WGD2 genes
                  #universe = IPR[IPR$species %in% c("Dipodascus_fermentans", "Geotrichum_candidum_CLIB_918"),]$gene,
                                 
                  #all post WGD3 genes
                  #universe = IPR[IPR$species %in% c("Saprochaete_suaveolens", "Magnusiomyces_magnusii", "Magnusiomyces_tetraspermus"),]$gene,
                                 
                  #all post WGD4 genes
                  #universe = IPR[IPR$species %in% c("Magnusiomyces_magnusii"),]$gene,
                    
                  TERM2GENE = IPR[, c("ID", "gene")],
                  TERM2NAME = IPR[, c("ID", "Name")],
                  pvalueCutoff = 0.05,
                  pAdjustMethod = 'fdr'
                    )

# Full species treea
# from https://figshare.com/ndownloader/files/48089638
full_tree <- read.tree("SpeciesTree.txt")

# Qualitative trait matrix
# from https://figshare.com/ndownloader/files/52745621 
niche_breadth <- read.table("TraitMatrix.txt", row.names=1)
niche_breadth <- data.frame(breadth=rowMeans(niche_breadth, na.rm = TRUE))

# Quantitative trait matrix
# from https://figshare.com/ndownloader/files/40705427
# Mannose here, pick your trait
growth <- read_excel("GrowthRates_Yeasts.xlsx") %>% select(Species, Mannose) %>% rename(species=Species) %>% drop_na() %>%
  group_by(species) %>%
  filter(n() == 1) %>%
  ungroup() %>% as.data.frame()

rownames(growth) <- gsub(" ", "_", growth$species)
growth <- growth %>% select(-species) %>% as.matrix()

#all post-WGD species
WGD_taxa <- c(extract.clade(full_tree, 1943)$tip.label, extract.clade(full_tree, 2106)$tip.label, extract.clade(full_tree, 2101)$tip.label)
group <- as.factor(ifelse(full_tree$tip.label %in% WGD_taxa, 1, 0))
names(group) <- full_tree$tip.label

#run ANOVAs
summary(aov.phylo(niche_breadth~group, full_tree))
summary(aov.phylo(growth~group, full_tree))
