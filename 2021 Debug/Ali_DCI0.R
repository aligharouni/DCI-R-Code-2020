
#  reading the fowchar of the functions,
library(tidyverse)
setwd("~/projects/DCI-R-Code-2020/2021 Debug/")
source("graph_fx.r")
source("graph_and_data_setup_for_DCI.r")
source("sum_fx.r")
source("get_adj_matrix_from_gis.r")
source("dci_calc_fx.r")
# adj_matrix <- read_csv("adjacency_matrix.csv", col_names = TRUE,col_types = cols())
adj_matrix<-get_adj_matrix_from_gis()

passability<-read.csv("segments_and_barriers.csv")
lengths<-read_csv("length.csv")

# output: graph of system
dirgraph<-graph_fx(edge_size=75,
         node_size=5,
         #adj_matrix=read.csv("adjacency_matrix.csv"),
         adj_matrix=adj_matrix,
         plot_it=F)
# debugonce(sum_fx)
sum_table<-sum_fx(adj_matrix,passability,lengths)

DCI_test <- dci_calc_fx(sum_table,lengths,all_sections=F)

# Note: The option all_sections=T was tried, it returns the same as "F"
# what is being calculated? Ans: see comment on the line 49 of dci_calc_f.r
# "#for diadromous fish we are only interested in the movement 
#     # from the segment which is closest to the ocean"


# Step 4? adding natural barrier;
summary_table_all <- graph_and_data_setup_for_DCI(adj_matrix,passability,lengths)
