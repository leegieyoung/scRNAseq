lapply(c("dplyr","Seurat","HGNChelper","openxlsx","scater","patchwork","metap","limma","multtest","ggplot2","cowplot","future.apply","ggraph","igraph","tidyverse","data.tree"), library, character.only = T)
pre_dir <- '/data/keeyoung/scRNA/iCD/'
work_dir <- paste0(pre_dir,"output/sctype/")

Change="4"

Plasma <- paste0(work_dir,"Plasma/",Change,"/","Plasma_",Change,".rds")
Tcell <- paste0("Tcell_",Change,".rds")
MNP <- paste0("MNP_",Change,".rds")
Stromal <- paste0("Stromal_",Change,".rds")

Cell_list <- c("Plasma","Tcell","MNP","Stromal")
#Cell_list="Plasma"
Change="4"

merge=""
for (i in 1:length(Cell_list)){
	assign(Cell_list[i], readRDS(paste0(work_dir,Cell_list[i],"/",Change,"/",Cell_list[i],"_",Change,".rds")))
	Cell <- eval(parse(text=Cell_list[i]))
	Idents(Cell) <- Cell$customclassif
	subtype = Cell$customclassif
	#2000 Genes on integrated.
	#Gene <- rownames(GetAssayData(object=Cell, slot="data", assay="integrated"))
	#30000 over Genes on RNA.
	Gene <- rownames(GetAssayData(object=Cell, slot="data", assay="RNA"))
	Cell <- GetAssayData(object=Cell, slot="counts",assay="RNA")[Gene,]
	Cell <- as.matrix(x=Cell)
	print("Is it equal that barcode and Cell name ?")
	print(all.equal(colnames(Cell), names(subtype)))
	colnames(Cell) = subtype
	print(dim(Cell))
	merge <- cbind(merge, Cell)
}

GIMATS <- c("IgG","INF.macs","Activated DC","ACKR1+EC","Activated Fibroblasts", "Highly Activated T cells")
GIMATS_counts <- merge[,(which(colnames(merge) %in% GIMATS))]
GIMATS_counts <- GIMATS_counts[,order(colnames(GIMATS_counts))]

write.table(GIMATS_counts,paste0(work_dir,"GIMATS_All_genes_counts.txt"), sep="\t", col.names=T)
