# Ali: This is a practice to calculate the DCI by FIPEX machinary
#  reading the fowchar of the functions,
library(tidyverse)
## used in sum_fx.r for "graphAM" object
library(Rgraphviz)
library(RBGL)

setwd("~/projects/DCI-R-Code-2020/2021 Debug/")
# Load the required functions
source("get_adj_matrix_from_gis.r")
source("graph_fx.r")
source("graph_and_data_setup_for_DCI.r")
source("sum_fx.r")

source("dci_calc_fx.r")
# adj_matrix <- read_csv("adjacency_matrix.csv", col_names = TRUE,col_types = cols())

# the adj_matrix() reads segment_matrix.csv from the directory! and output the adj matrix
#The input is the segment matrix, which is really two vectors. For each 
#segment in the first column (Seg_ID), the second column (Seg) gives the id 
#of other sections it touches, including itself.
adj_matrix<-get_adj_matrix_from_gis()

# I will make my toy model's segment matrix and will work with that;
segment_matrix_ali <- read.csv(sep=",", strip.white=TRUE,
text="
         Seg_ID, Seg
         sink, sink
         sink, 2_s
         sink, 3_s
         2_s, 2_s
         2_s, sink
         3_s, 3_s
         3_s, sink
         ")
write.table(segment_matrix_ali,"segment_matrix_ali.csv",
            row.names=F, sep=",")

adj_matrix<-get_adj_matrix_from_gis(inputname="segment_matrix_ali.csv",
                                    outputname="adj_matrix_ali.csv")

## FIPEX_output_to_R_input.r
## Inputs: 
#   FIPEX_connectivity.csv
#   FIPEX_BarrierHabitatLine.csv
# Outputs: 
#   segment_matrix.csv
#   barrier.csv 

## Let's construct the barrier_ali.csv file
#  Q: I think, in the FIPEX_output_to_R_input.r (lines 100-110) when constructing the barrier.csv,
# in the column Seg_ID the upstream is first and the downstream seg is the 2nd, but the Pass is the same. Can we change this Pass?
# segments_and_barriers.csv is only set from upstream to downstream
#  segments_and_barriers.csv created in convert_gis_output_to_r_format.r

(0.5*0.7)

segments_and_barriers_ali <- read.csv(sep="|", strip.white=TRUE,
                                      text="
Bar_ID|Seg_1|Seg_2|Pass|nat_barrier|section1_2
2|2_s|sink|0.35|FALSE|2_s,sink
2|sink|2_s|0.35|FALSE|sink,2_s
3|3_s|sink|0.35|FALSE|3_s,sink
3|sink|3_s|0.35|FALSE|sink,3_s
") ## see metapop_params.R 



# passability<-read.csv("segments_and_barriers.csv")
passability<-segments_and_barriers_ali

# lengths<-read_csv("length.csv")

lengths <- read.csv(sep=",", strip.white=TRUE, text="
Seg_ID, Shape_Length
sink, 2
2_s, 2
3_s, 2
")


# # output: graph of system
# dirgraph<-graph_fx(edge_size=75,
#          node_size=5,
#          #adj_matrix=read.csv("adjacency_matrix.csv"),
#          adj_matrix=adj_matrix,
#          plot_it=F)

# # debugonce(sum_fx)
sum_table <- sum_fx(adj_matrix,passability,lengths)

# Ali: adding a column to amount for d_ij distance along a path

library(tidyverse)
# total_length <- as.numeric(lengths %>% filter(Seg_ID %in% path) %>% summarise(sum(Shape_Length)))

# This fun returns the cummulative usm of the lengths of a given path by looking up the segments from lengths df
cumlength <- function(lengths, path){
  # length is a df with columns: Seg_ID, Shape_Length
  # path is a string of "seg1,seg2,...,segn"
  o <-as.numeric(lengths 
                 %>% filter(Seg_ID %in% unlist(strsplit(path,","))  ) 
                 %>% summarise(sum(Shape_Length)))
  return(o)
}
# test it:
cumlength(lengths,sum_table$path2[6])

 
cum_length <- apply(sum_table %>% select(path2), 1, FUN = function(par) cumlength(lengths,par))

sum_table<- (sum_table  %>% mutate(cum_length=cum_length))

write.table(sum_table,file="sum_table_toy1.csv")
# write.csv(sum_table, file="sum_table_toy1.csv", row.names=F)

# Note, (i) not reading one of the directions, (is the input the product of passability? or in the process they) in the code the are just picking the first passability.
# 

DCI_test <- dci_calc_fx(sum_table,lengths,all_sections=F)
# result: DCIp=75.33333 DCId=80
# my indices_cor.R which uses dci_p_mod0(params) gives 67.77778
# TODO: find out why DCI_p's from FIPEX and mine are different?

# Note: The option all_sections=T was tried, it returns the same as "F"
# what is being calculated? Ans: see comment on the line 49 of dci_calc_f.r
# "#for diadromous fish we are only interested in the movement 
#     # from the segment which is closest to the ocean"


# Step 4? adding natural barrier;
summary_table_all <- graph_and_data_setup_for_DCI(adj_matrix,passability,lengths)
