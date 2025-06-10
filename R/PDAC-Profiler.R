# ##################################################
# Load library
# ##################################################
library(optparse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(openxlsx)
#library(pheatmap)
#library(GO.db)
#library(fgsea)
#library(viridis)
#library(msigdbr)
library(ggplot2)
library(limma)
library(tidyr)
library(clusterProfiler)
#library(colorspace)
#library(enrichplot)
library(gridExtra)
library(parallel)
#library(grid)
library(circlize)
library(ComplexHeatmap)
#library(RColorBrewer)
#library(ggrepel)
library(homologene)

NBCPU <- 8

# ##################################################
# Source functions
# ##################################################

source("./R/enrichment_plot_helper_fct.R")

# ##################################################
# FUNCTIONS
# ##################################################
unifyDF <- function(df, nameIn = "ENSEMBL", nameOut = "ENTREZID", org = org.Hs.eg.db, scoreFct = "MAD", count = FALSE) {
    
    # nameIn, nameOut can be ENSEMBL, ENTREZID, SYMBOL, ALIAS
    
    # create gene annotation matrix
    ann.gene <- AnnotationDbi::select(org, rownames(df), c(nameIn, nameOut), nameIn)
    
    # remove missing names
    keep <- !is.na(ann.gene[, nameOut])
    ann.gene <- ann.gene[keep, ]
    df <- df[rownames(df) %in% ann.gene[, nameIn], ]
    
    # calculate variability score (MAD, IQR, SD...)
    if(count){
        
        # norm lib size and get CPM
        y <- edgeR::DGEList(counts=df)
        y <- edgeR::calcNormFactors(y)
        df.cpm <- edgeR::cpm(y)
        
        if(scoreFct == "MAD") score <- apply(df.cpm, 1, mad)
        if(scoreFct == "IQR") score <- apply(df.cpm, 1, iqr)
        if(scoreFct == "SD") score <- apply(df.cpm, 1, sd)
    }else{
        if(scoreFct == "MAD") score <- apply(df, 1, mad)
        if(scoreFct == "IQR") score <- apply(df, 1, iqr)
        if(scoreFct == "SD") score <- apply(df, 1, sd)
    }
    
    ann.gene$SCORE <- score[match(ann.gene[, nameIn], names(score))]
    
    # re-order and select the first (best)
    newOrder <- order(-ann.gene$SCORE)
    ann.gene <- ann.gene[newOrder, ]
    
    ann.gene <- ann.gene[!duplicated(ann.gene[, nameIn]), ]
    ann.gene <- ann.gene[!duplicated(ann.gene[, nameOut]), ]
    
    # return unified data
    df <- df[match(ann.gene[, nameIn], rownames(df)), ]
    rownames(df) <- ann.gene[, nameOut]
    return(df)
}
convert_rownames_to_column <- function(matrix, colname) {
    df <- as.data.frame(matrix)
    df[[colname]] <- rownames(df)
    rownames(df) <- NULL
    shuffled_df <- df[, c(ncol(df), 1:(ncol(df)-1))]
    return(shuffled_df)
}
convert_column_to_rownames <- function(matrix, colname) {
    countMat <- as.data.frame(matrix)
    countMat <- na.omit(countMat)
    rownames(countMat) <- countMat[[colname]] 
    countMat <- countMat[, -which(colnames(countMat) == colname)]
    return(countMat)
}
human2mouseSymbol <- function(hg.symbol) {
    homog <- homologene(hg.symbol, inTax = 9606, outTax = 10090)
    mm.symbol <- homog[match(hg.symbol, as.character(homog[, 1])), 2]
    mm.symbol[is.na(mm.symbol)] <- "unknown"
    return(mm.symbol)
}

human2mouseEntrez <- function(hg.symbol) {
    homog <- homologene(hg.symbol, inTax = 9606, outTax = 10090)
    mm.symbol <- homog[match(hg.symbol, as.character(homog[, 3])), 4]
    mm.symbol[is.na(mm.symbol)] <- "unknown"
    return(mm.symbol)
}

my.read_xlsx <- function(inFile) {
    mysheets <- getSheetNames(inFile)
    mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
    names(mList) <- mysheets
    return(mList)
}

symbol2entrez <- function(symbol) {
    entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
    entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
}

entrez2symbol <- function(entrez) {
    symbol <- mget(as.character(entrez), org.Mm.egSYMBOL, ifnotfound=NA)
    symbol <- unique(unlist(lapply(symbol, function(i) return(i[1]))))
    return(symbol)
}

ensembl2entrez <- function(ensembl) {
    entrez <- mget(as.character(ensembl), org.Mm.egENSEMBL2EG, ifnotfound=NA)
    entrez <- lapply(entrez, function(i) return(i[1]))
    return(unlist(entrez))
}

refseq2entrez <- function(refseq) {
    entrez <- mget(as.character(refseq), org.Mm.egREFSEQ2EG, ifnotfound=NA)
    entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
    return(entrez)
}

mygsInfo <- function(object, geneSetID) {
    geneList <- object@geneList
    
    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]
    
    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")

getAvgNESmat <- function(nes, group){
    
    group.unique <- unique(group)
    
    nes.avg <- lapply(group.unique, function(i){
        idx <- which(group == i)
        return(rowMeans(nes[, idx]))
    })
    nes.avg <- do.call(cbind, nes.avg)
    colnames(nes.avg) <- group.unique
    return(nes.avg)
}

# ##################################################
# MAIN
# ##################################################

# ##################################################
# LOAD ARGUMENTS

option_list <- list(
    make_option(c("--input_file"), type = "character", default = NULL, 
            help = "path to the mRNA intensity file (.tsv)", metavar = "character"),
    make_option(c("--signature_file"), type = "character", default = NULL, 
                help = "path to the PDAc signature file", metavar = "character"),   
    make_option(c("--output_dir"), type = "character", default = NULL, 
                help = "path to the output directory", metavar = "character"),
    make_option(c("--species_name"), type = "character", default = "human", 
                help = "human or mouse", metavar = "character")
    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input_file
signature_file <- opt$signature_file
output_dir <- opt$output_dir
species_name <- opt$species_name

if(FALSE){
    # toy
    input_file <- "log2cpm.tsv"
    signature_file <- "PDAC_signatures.tsv"
    output_dir <- "output"
    species_name <- "human"
}

if (!is.null(species_name)) {
    if (!species_name %in% c("human", "mouse")) {
        stop("species_name must be either human or mouse")
    }
}else {
    species_name <- "human"
    print("species_name is missing; running human as default.")
}

# ##################################################
# BUILD DB FROM PDAC SUBTYPE SIGNATURES

sig_df <- readr::read_tsv(file.path(signature_file))

sig_list <- lapply(seq_len(ncol(sig_df)), function(i) {
    genes <- as.character(sig_df[, i])
    genes.symbol <- alias2Symbol(genes, species = "Hs")
    genes.entrez <- symbol2entrez(genes.symbol)
    return(sort(unique(genes.entrez)))
})
names(sig_list) <- colnames(sig_df)

sig_list$Basal_Union <- Union(sig_list[grepl("basal", names(sig_list), ignore.case = TRUE)])
sig_list$Classical_Union <- Union(sig_list[grepl("classi", names(sig_list), ignore.case = TRUE)])

dbList <- list("PDAC_subtype" = sig_list)

if (species_name == "mouse") {
    # Human 2 mouse
    dbList <- lapply(dbList, function(db){
        db <- lapply(db, human2mouseEntrez)
        db <- lapply(db, function(x) x[x != "unknown"])
        db
})
}

# SET GENE-SETS LIMITS
dbList <- lapply(dbList, function(db) {
    db[sapply(db, length) >= 5]
})

if (FALSE) {
    # ENTREZ 2 SYMBOL 
    dbList <- mclapply(dbList, function(db){
        db <- lapply(db, entrez2symbol)
        db <- lapply(db, unique)
        db
    }, mc.cores = NBCPU)
}

# LIST TO TIBBLE
dbTibbles <- lapply(dbList, function(db){
    df <- lapply(seq_along(db), function(i){
        return(data.frame(TERM = names(db)[i],
                        geneID = as.character(db[[i]])))
    })
    df <- do.call(rbind, df)
    tibble(df)
})

# ##################################################
# based on CPM (single sample)

expr <- readr::read_delim(input_file)

# ##################################################
# Ensembl 2 entrez
if (species_name == "mouse") {
    expr <- unifyDF(df = expr, nameIn = "SYMBOL", nameOut = "ENTREZID", org = org.Mm.eg.db, scoreFct = "MAD", count = FALSE)
} else {
    expr <- unifyDF(df = expr, nameIn = "SYMBOL", nameOut = "ENTREZID", org = org.Hs.eg.db, scoreFct = "MAD", count = FALSE)
}
# scale
expr <- t(scale(t(expr)))

ranked_genes_list <- lapply(seq_len(ncol(expr)), function(i){
    values <- as.numeric(expr[, i])
    names(values) <- rownames(expr)
    sort(values, decreasing = TRUE)
})
names(ranked_genes_list) <- colnames(expr)

# Define color
mycolor <- c("#845EC2", "#B0A8B9", "#C34A36", "#008D83")
if (length(mycolor) < length(ranked_genes_list)) {
    nb_cols <- length(ranked_genes_list)
    mycolor <- colorRampPalette(brewer.pal(8, "Set2"))(nb_cols)
}else {
    mycolor <- mycolor[seq_along(ranked_genes_list)]
}

#########################
# Perform fgsea analysis

fgsea_dir <- file.path(output_dir, "fgsea")
dir.create(fgsea_dir, recursive = TRUE, showWarnings = FALSE)

setwd(fgsea_dir)
lapply(seq_along(dbTibbles), function(i){
    mydb <- dbTibbles[[i]]
    mydb_name <- names(dbList)[i]
    
    dir.create(file.path(fgsea_dir, mydb_name), showWarnings = FALSE)
    setwd(file.path(fgsea_dir, mydb_name))   
    
    gsea_list <- mclapply(ranked_genes_list, GSEA, TERM2GENE = mydb, pvalueCutoff = 1, by = "fgsea",
                        nPermSimple = 10000, maxGSSize = 1000, mc.cores = NBCPU)
    
    #p <- cnetplot(gsea_list[[1]], foldChange = ranked_genes_list[[1]], colorEdge = TRUE)
    
    # Add symbols
    gsea_list_df <- lapply(gsea_list, function(i) as.data.frame(i@result))
    
    gsea_list_df <- lapply(gsea_list_df, function(k){
        k.entrez <- strsplit(k$core_enrichment, split = "\\/")
        k$NB <- sapply(k.entrez, length)
        k.symbol <- lapply(k.entrez, entrez2symbol)
        k.symbol <- lapply(k.symbol, function(j) unique(j[!is.na(j)]))
        k.symbol <- lapply(k.symbol, toString)
        k$leadingEdge.symbol <- unlist(k.symbol)
        k[, colnames(k) != "Description"]
    })
    
    gsea_list_df <- lapply(gsea_list_df, function(j){
        j$NES[is.na(j$NES)] <- 0
        j
    })
    
    # Divide UP and DOWN
    gsea_list_df_up <- lapply(gsea_list_df, function(j) j[j$NES > 0, ])
    names(gsea_list_df_up) <- paste0(names(gsea_list_df), ".UP")
    gsea_list_df_down <- lapply(gsea_list_df, function(j) j[j$NES < 0, ])
    names(gsea_list_df_down) <- paste0(names(gsea_list_df), ".DOWN")
    gsea_list_df <- c(gsea_list_df_up, gsea_list_df_down)
    
    # Save
    lapply(seq_along(ranked_genes_list), function(j) {
        write.xlsx(list(UP =  gsea_list_df_up[[j]], DOWN = gsea_list_df_down[[j]]),
                paste(names(ranked_genes_list)[j], "_", mydb_name, "_fgsea.xlsx", sep = ""),
                rowNames = FALSE, firstRow = TRUE, headerStyle = createStyle(textDecoration = 'bold'), overwrite = TRUE)  
    })	
    
    if (FALSE) {
        # PLOT
        dir.create(file.path(fgsea_dir, mydb_name, "plot"))
        
        nb <- 20
        gene_set_id <- lapply(gsea_list, function(x){
        gs <- x@result$ID
        if (length(gs) > nb) gs <- gs[1:nb]
        gs
        })
        gene_set_id <- unique(unlist(gene_set_id))
        if (mydb_name == "H") gene_set_id <- unique(mydb$TERM)
        
        lapply(gene_set_id, function(j){
        if (length(gsea_list) != 1 & length(gsea_list) <= 6) {
            p <- mygseaplot_multi(gsea_list, geneSetID = j, pvalue_table = TRUE,
                                color = mycolor, ES_geom = "line")
            ggsave(plot = p,
                filename = file.path(fgsea_dir, mydb_name, "plot", paste0(j, "_multiplot.pdf")),
                width = 8, height = 6)
        }else {
            dir.create(file.path(fgsea_dir, mydb_name, "plot", j))
            lapply(seq_along(gsea_list), function(x){
            
            p <- mygseaplot2(gsea_list[[x]], geneSetID = j, pvalue_table = TRUE,
                            color = "green4", ES_geom = "line")
            ggsave(plot = p,
                    filename = file.path(fgsea_dir, mydb_name, "plot", j, paste0(j, "_", names(gsea_list)[x], "_plot.pdf")),
                    width = 8, height = 6)
            })
        }
        })
        
    }
    
})

#######################################
# HEATMAPS
comp <- names(ranked_genes_list)

setwd(fgsea_dir)
lapply(seq_along(dbList), function(i){
    mydb_name <- names(dbList)[i]
    dir.create(file.path(mydb_name, "heatmaps"), showWarnings = FALSE)
    
    # load fisher results
    fh_files <- file.path(mydb_name, paste0(comp, "_", mydb_name, "_fgsea.xlsx"))
    
    fh_list <- lapply(fh_files, my.read_xlsx)
    names(fh_list) <- gsub(paste0("_", mydb_name, "_fgsea.xlsx"), "", fh_files)
    names(fh_list) <- gsub(paste0(mydb_name, "\\/"), "", names(fh_list))
    fh_list <- do.call(c, fh_list)
    
    heatmapPipe(fh_list, file.path(mydb_name, "heatmaps", paste0(mydb_name, "_fgsea_heatmap_PV.pdf")),
                nb = 5, TermIdx = 1, pvIdx = 5, qvIdx = NULL, colClust = TRUE)
    heatmapPipe(fh_list, file.path(mydb_name, "heatmaps", paste0(mydb_name, "_fgsea_heatmap_QV.pdf")),
                nb = 5, TermIdx = 1, pvIdx = 5, qvIdx = 6, colClust = TRUE)

})


################################################################################
# NES HEATMAP

setwd(file.path(fgsea_dir, "PDAC_subtype/"))
fh_files <- list.files(pattern = "*_PDAC_subtype_fgsea.xlsx")
fh_list <- lapply(fh_files, function(i) {
    rbind(read.xlsx(i, sheet = "UP"), read.xlsx(i, sheet = "DOWN"))
}) 
names(fh_list) <- gsub("_PDAC_subtype_fgsea.xlsx", "", fh_files)

nes_df <- getNESmat(fh_list)
colnames(nes_df) <- gsub("\\.", "-", colnames(nes_df))

nes_df <- nes_df[, match(ann.sample$FASTQ, colnames(nes_df))]

# save
toxlsx <- data.frame(SAMPLE = colnames(nes_df.avg),
                    BasalScore = nes_df.avg["Basal_Union", ] - nes_df.avg["Classical_Union", ],
                    ClassicalScore = nes_df.avg["Classical_Union", ] - nes_df.avg["Basal_Union", ]
                    )
toxlsx$BasalScore <- scales::rescale(toxlsx$BasalScore)
toxlsx$ClassicalScore <- scales::rescale(toxlsx$ClassicalScore)

toxlsx$Subtype <- NA
toxlsx$Subtype[toxlsx$BasalScore > toxlsx$ClassicalScore] <- "Basal"
toxlsx$Subtype[toxlsx$ClassicalScore > toxlsx$BasalScore] <- "Classical"

toxlsx <- toxlsx[order(toxlsx$SAMPLE), ]
write.xlsx(toxlsx, file.path(fgsea_dir, "PDAC_subtype.xlsx"))