library(tidyverse)


#### Read in data ####

PBMC_scRNAseq_study_Berlin <- read.csv("path/to/PBMC_scRNAseq_study_Berlin.csv",sep=";")
PBMC_scRNAseq_study_Bonn <- read.csv("path/to/PBMC_scRNAseq_study_Bonn",sep=";")
respiratory_specimen_scRNAseq_study <- read.csv("path/to/respiratory_specimen_scRNAseq_study",sep=";")

cell_types <- c("B.cell", "Basal", "Ciliated.diff", "Ciliated", "CTL", "FOXN4", "Ionocyte", "IRC", "MC", "moDC", "moMa", "Neu","NK", "NKT.p",
                "NKT", "nrMa", "outliers_epithelial", "pDC", "rMa", "Secretory.diff", "Secretory", "Squamous", "Treg", "unknown_epithelial" )


#### Figure 4 ####

plot_expression_overtime <- function(cell_type="B.cell", data){
  data %>% 
    ggplot(aes_string(y=cell_type, x="days_post_symptoms", col="disease_status"))+
    geom_point(size=3)+
    geom_line(aes(group=Patient_ID))+
    geom_smooth()+
    theme_bw(base_size=20)+
    scale_color_manual(values=cols)+
    labs(x="days post symptom onset")+
    theme(legend.title=element_blank(), legend.position="top")
}

map(cell_types, plot_expression_overtime, data = PBMC_scRNAseq_study_Berlin)


#### Suppl. Figure 3 ####

map(cell_types, plot_expression_overtime, data = PBMC_scRNAseq_study_Bonn)


#### Suppl. Figure 4 ####

map(cell_types, plot_expression_overtime, data = respiratory_specimen_scRNAseq_study)
