# Ali: This is a practice to calculate the DCI by FIPEX machinary
#  reading the fowchar of the functions,
library(tidyverse)
library(graph)
library(grid)
## used in sum_fx.r for "graphAM" object
if (!require("BiocManager", quietly = TRUE))  install.packages("BiocManager")

# BiocManager::install("Rgraphviz")
# BiocManager::install("RBGL")
library(Rgraphviz)
library(RBGL)

setwd("~/projects/DCI-R-Code-2020/2021 Debug/")
# Load the required functions
source("get_adj_matrix_from_gis.r")
source("graph_fx.r")
source("graph_and_data_setup_for_DCI.r")
source("sum_fx.r")

source("dci_calc_fx.r")
source("dci_calc_fx_AG.r") ## Ali's version of DCI calc with directed passabilities
# adj_matrix <- read_csv("adjacency_matrix.csv", col_names = TRUE,col_types = cols())

# the adj_matrix() reads segment_matrix.csv from the directory! and output the adj matrix
#The input is the segment matrix, which is really two vectors. For each 
#segment in the first column (Seg_ID), the second column (Seg) gives the id 
#of other sections it touches, including itself.
# adj_matrix<-get_adj_matrix_from_gis()

################################################################################
# 3-node example:
###############################################################################
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

segments_and_barriers_ali <- read.csv(sep="|", strip.white=TRUE,
                                      text="
Bar_ID|Seg_1|Seg_2|Pass|nat_barrier|section1_2
2|2_s|sink|0.1|FALSE|2_s,sink
2|sink|2_s|0.1|FALSE|sink,2_s
3|3_s|sink|0.5|FALSE|3_s,sink
3|sink|3_s|0.5|FALSE|sink,3_s
") ## see metapop_params.R 

# lengths<-read_csv("length.csv")
lengths <- read.csv(sep=",", strip.white=TRUE, text="
Seg_ID, Shape_Length
sink, 1
2_s, 2
3_s, 3
")

passability<-segments_and_barriers_ali
sum_table <- sum_fx(adj_matrix=adj_matrix,passability=passability,lengths=lengths)
## Calculate the DCIs for the directed network:
# dci_calc_fx(sum_table = sum_table,lengths = lengths,all_sections=F)
dci_calc_fx(sum_table = sum_table,lengths = lengths,all_sections=T)


###############################################################################
## Example 2. Asymmetric network with 3 nodes
## AG:: the upstream, downstream should be appear here in the Bar_ID

segment_matrix_ali3 <- segment_matrix_ali

# write.table(segment_matrix_ali3,"segment_matrix_ali3.csv",
#             row.names=F, sep=",")

adj_matrix3<-get_adj_matrix_from_gis(inputname="segment_matrix_ali3.csv",
                                    outputname="adj_matrix_ali3.csv")

segments_and_barriers_ali_Asym3 <- read.csv(sep="|", strip.white=TRUE,
                                      text="
Bar_ID|Seg_1|Seg_2|Pass|nat_barrier|section1_2
2d|2_s|sink|1|FALSE|2_s,sink
2u|sink|2_s|0.1|FALSE|sink,2_s
3d|3_s|sink|1|FALSE|3_s,sink
3u|sink|3_s|0.5|FALSE|sink,3_s
") ## see metapop_params.R 

lengths3 <- lengths
passability_Asym3<-segments_and_barriers_ali_Asym3

sum_table3 <- sum_fx(adj_matrix=adj_matrix3,passability=passability_Asym3,lengths=lengths3)
## Calculate the DCIs for the directed network:
dci_calc_fx_AG(sum_table = sum_table3,lengths = lengths3,all_sections=F)
dci_calc_fx_AG(sum_table = sum_table3,lengths = lengths3,all_sections=T)
###############################################################################
## 4 node example:
################################################################################
segment_matrix_ali4 <- read.csv(sep=",", strip.white=TRUE,
                               text="
         Seg_ID, Seg
         sink, sink
         sink, 2_s
         sink, 3_s
         2_s, 2_s
         2_s, sink
         2_s, 4_s
         3_s, 3_s
         3_s, sink
         4_s, 4_s
         4_s, 2_s
         ")
write.table(segment_matrix_ali4,"segment_matrix_ali4.csv",
            row.names=F, sep=",")

adj_matrix4<-get_adj_matrix_from_gis(inputname="segment_matrix_ali4.csv",
                                    outputname="adj_matrix_ali4.csv")



segments_and_barriers_ali_Asym4 <- read.csv(sep="|", strip.white=TRUE,
                                           text="
Bar_ID|Seg_1|Seg_2|Pass|nat_barrier|section1_2
2d|2_s|sink|1|FALSE|2_s,sink
2u|sink|2_s|0.1|FALSE|sink,2_s
3d|3_s|sink|1|FALSE|3_s,sink
3u|sink|3_s|0.2|FALSE|sink,3_s
4d|4_s|2_s|1|FALSE|4_s,2_s
4u|2_s|4_s|0.3|FALSE|2_s,4_s
") ## see metapop_params.R 

lengths4 <- read.csv(sep=",", strip.white=TRUE, text="
Seg_ID, Shape_Length
sink, 1
2_s, 2
3_s, 3
4_s, 4
")

passability4<-segments_and_barriers_ali_Asym4
sum_table4 <- sum_fx(adj_matrix4,passability4,lengths4)
## Calculate the DCIs for the directed network:
dci_calc_fx(sum_table4,lengths4,all_sections=F)

###############################################################################

# passability<-read.csv("segments_and_barriers.csv")
# passability<-segments_and_barriers_ali








# # output: graph of system
# dirgraph<-graph_fx(edge_size=75,
#          node_size=5,
#          #adj_matrix=read.csv("adjacency_matrix.csv"),
#          adj_matrix=adj_matrix,
#          plot_it=F)

# # debugonce(sum_fx)
sum_table <- sum_fx(adj_matrix,passability,lengths)

# Ali: adding a column to amount for d_ij distance along a path

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

cumlength2 <- function(lengths, path){
  # length is a df with columns: Seg_ID, Shape_Length
  # path is a string of "seg1,seg2,...,segn"
  o <-as.numeric(lengths 
                 %>% filter(Seg_ID %in% unlist(strsplit(path,","))  ) 
                 ## the central lengths are added 
                 %>% summarise(sum(Shape_Length)-0.5*(Shape_Length[1]+Shape_Length[length(Shape_Length)]))
                 )
  return(o)
}

# test it:
cumlength(lengths,sum_table$path2[6])
# eg
path <- sum_table$path2[6]
lengths 
as.numeric(lengths 
  %>% filter(Seg_ID %in% unlist(strsplit(path,","))  )
  %>% summarise(sum(Shape_Length)-0.5*(Shape_Length[1]+Shape_Length[length(Shape_Length)]))
)


 
cum_length <- apply(sum_table %>% select(path2), 1, FUN = function(par) cumlength(lengths,par))
sum_table<- (sum_table  %>% mutate(cum_length=cum_length))

cum_length2 <- apply(sum_table %>% select(path2), 1, FUN = function(par) cumlength2(lengths,par))
sum_table<- (sum_table  %>% mutate(cum_length2=cum_length2))



write.table(sum_table,file="sum_table_toy1.csv")
# write.csv(sum_table, file="sum_table_toy1.csv", row.names=F)

# Note, (i) not reading one of the directions, (is the input the product of passability? or in the process they) in the code the are just picking the first passability.
# 

DCI_test <- dci_calc_fx(sum_table,lengths,all_sections=F)
DCI_test
# result: DCIp=75.33333 DCId=80
# my indices_cor.R which uses dci_p_mod0(params) gives 67.77778
# TODO: find out why DCI_p's from FIPEX and mine are different?

# Note: The option all_sections=T was tried, it returns the same as "F"
# what is being calculated? Ans: see comment on the line 49 of dci_calc_f.r
# "#for diadromous fish we are only interested in the movement 
#     # from the segment which is closest to the ocean"


# Step 4? adding natural barrier;
summary_table_all <- graph_and_data_setup_for_DCI(adj_matrix,passability,lengths)




# 7 node binary tree  -----------------------------------------------------

# Ali: I took the binary tree approach in Jul 2022, using igraph package
#  See "/home/ag/projects/connectivity/scripts/graphData.R"
adj_matrix2<-get_adj_matrix_from_gis(inputname="segment_matrix_igraph_ali.csv",
                                     outputname="adj_matrix_igraph_ali.csv")
lengths2 <- read.csv("length2_ali.csv", stringsAsFactors = FALSE)
passability2 <- read.csv("segments_and_barriers2_ali.csv", stringsAsFactors = FALSE)
sum_table2 <- sum_fx(adj_matrix2,passability2,lengths2)
dci <- dci_calc_fx(sum_table2,lengths2,all_sections=F)

graph_and_data_setup_for_DCI(adj_matrix2,passability2,lengths2)

#--------------------------------------------------------------------------
# Randimed sample from connectivity/scripts/netwSim_DataDir.R

## Directionless sample: linear_3_1
# lengths<-read_csv("length.csv")
lengths.in <- read.csv(sep=",", strip.white=TRUE, text="
Seg_ID, Shape_Length
sink, 0.2735782
2_s, 0.3994519
3_s, 0.3269699
")

sum_table.in<- read.csv(sep=" ", strip.white=TRUE, text= "
start end path2 barrier_id pathway_pass start_section_length finish_section_length cum_length
2_s 2_s 2_s NA 1.000000000 0.3994519 0.3994519 0.0000000
2_s 3_s 2_s,3_s 3 0.554699885 0.3994519 0.3269699 0.3632109
2_s sink 2_s,sink 2 0.009998980 0.3994519 0.2735782 0.3365151
3_s 2_s 3_s,2_s 3 0.554699885 0.3269699 0.3994519 0.3632109
3_s 3_s 3_s NA 1.000000000 0.3269699 0.3269699 0.0000000
3_s sink 3_s,2_s,sink 3,2 0.005546433 0.3269699 0.2735782 0.6997260
sink 2_s sink,2_s 2 0.009998980 0.2735782 0.3994519 0.3365151
sink 3_s sink,2_s,3_s 2,3 0.005546433 0.2735782 0.3269699 0.6997260
sink sink sink NA 1.000000000 0.2735782 0.2735782 0.0000000
")

dci_calc_fx(sum_table.in,lengths.in,all_sections=TRUE)

## Directed sample:
lengths.di <- lengths.in

sum_table.di<- read.csv(sep=" ", strip.white=TRUE, text= "
start end path2 barrier_id pathway_pass start_section_length finish_section_length cum_length
2_s 2_s 2_s NA 1.000000000 0.3994519 0.3994519 0.0000000
2_s 3_s 2_s,3_s 3u 0.554699885 0.3994519 0.3269699 0.3632109
2_s sink 2_s,sink 2d 1.000000000 0.3994519 0.2735782 0.3365151
3_s 2_s 3_s,2_s 3d 1.000000000 0.3269699 0.3994519 0.3632109
3_s 3_s 3_s NA 1.000000000 0.3269699 0.3269699 0.0000000
3_s sink 3_s,2_s,sink 3d,2d 1.000000000 0.3269699 0.2735782 0.6997260
sink 2_s sink,2_s 2u 0.009998980 0.2735782 0.3994519 0.3365151
sink 3_s sink,2_s,3_s 2u,3u 0.005546433 0.2735782 0.3269699 0.6997260
sink sink sink NA 1.000000000 0.2735782 0.2735782 0.0000000
")


dci_calc_fx_AG(sum_table.di,lengths.di,all_sections=TRUE)









