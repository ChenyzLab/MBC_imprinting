#This script is used to prepare files for circos plot
library(dplyr)
library(tidyverse)
order <- c("rNAV","aNAV","DN2","DN1","uMBC","sMBC","aMBC","PB")

#### 1. load dataframe ####
print(nrow(clone_df))
clone_df$celltype <- clone_df$subType

order <- intersect(order, unique(clone_df$celltype))

clone_df$celltype <- factor(clone_df$celltype, levels = order)
clone_df <- clone_df[order(clone_df$celltype), ]

#### 2. prepare table ####
# 2.1 Count the quantity according to clone_id and sample
step01 <- clone_df %>%
  group_by(clone_id, celltype) %>%
  summarise(n_seqs = n())

# 2.2 add repeated counts
step01 <- step01 %>%
  group_by(clone_id) %>%
  mutate(repeat_counts = n())

# 2.3 seperate df based on repeated counts
step02_2 <- step01[step01$repeat_counts == 2, ]
step02_n <- step01[step01$repeat_counts > 2, ]
step02_1 <- step01[step01$repeat_counts ==1, ]
rm(step01)

# 2.3.1 add shared info for repeat_counts == 2
dfs <- split(step02_2, step02_2$clone_id)
if (length(dfs) == 0) {
  stop("Terminating execution because no shared clone exists.")
} else {
  for (i in 1:length(dfs)) {
    dfs[[i]][1, 'shared'] <- dfs[[i]][2, 'celltype']
    dfs[[i]][2, 'shared'] <- dfs[[i]][1, 'celltype']
  }
}

rm(i)
step03_2 <- bind_rows(dfs)
rm(dfs, step02_2)
step03 <- rbind(step02_n, step03_2, step02_1)

# 2.3.2 add shared info for repeat_counts > 2

step03 <- step03[order(step03$celltype,-step03$repeat_counts,-step03$n_seqs,step03$shared), ]

predefined_orders <- step03 %>%
  group_by(celltype) %>%
  summarise(clone_order = list(clone_id)) %>%
  deframe()


# 2.3.3 add start & end
dfs <- split(step03, step03$celltype)
for (i in 1:length(dfs)) {
  start = 0
  for (j in 1:nrow(dfs[[i]])){
    dfs[[i]][j, 'start'] = start
    dfs[[i]][j, 'end'] = start + dfs[[i]][j, 'n_seqs'] 
    start = dfs[[i]][j, 'end']
  }
}
step04 <- bind_rows(dfs)
rm(dfs, start, i, j, step03,step02_1)

##### 2. create link #####
# 2.1 assign 2 link
step05_se <- step04[, c('clone_id', 'celltype', 'start', 'end')]
step05_2 <- step03_2

step05_2 <- merge(step05_2, step05_se, by = c('clone_id', 'celltype'), all.x = T)
step05_2$start1 <- step05_2$start
step05_2$end1 <- step05_2$end
step05_2$start <- NULL
step05_2$end <- NULL

step05_se$shared <- step05_se$celltype
step05_2 <- merge(step05_2, step05_se, by = c('clone_id', 'shared'), all.x = T)
step05_2$start2 <- step05_2$start
step05_2$end2 <- step05_2$end
step05_2$start <- NULL
step05_2$end <- NULL

step05_2 <- step05_2[, c('clone_id', 'celltype.x', 'start1', 'end1', 'shared', 'start2', 'end2')]
step05_2 <- step05_2[!duplicated(step05_2$clone_id), ]

step05_2 <- step05_2[, c('celltype.x', 'start1', 'end1', 'shared', 'start2', 'end2')]
colnames(step05_2) <- c("celltype1", "start1", "end1", "celltype2", "start2", "end2")

rm(step03_2)

# 2.2 assign multi-link
step05_n <- step02_n

step05_se$shared <- NULL

dfs <- split(step05_n, step05_n$clone_id)
if (length(dfs)>0) {
  for (i in 1:length(dfs)){
    df <- dfs[[i]][, c('clone_id', 'celltype')]
    df$celltype <- as.character(df$celltype)
    combinations <- t(combn(nrow(df), 2))
    combined_df <- data.frame(matrix(nrow = nrow(combinations), ncol = ncol(df) * 2))
    
    for(j in 1:nrow(combinations)) {
      combined_row <- c(df[combinations[j, 1], ], df[combinations[j, 2], ])
      combined_df[j, ] <- combined_row
    }
    dfs[[i]] <- combined_df
  }
  rm(combined_df, i, j, df, combinations, combined_row)
  step05_n <- bind_rows(dfs)
} 


colnames(step05_n) <- c('clone_id', 'celltype', 'clone', 'shared')

step05_n <- merge(step05_n, step05_se, by = c('clone_id', 'celltype'), all.x = T)
step05_n$start1 <- step05_n$start
step05_n$end1 <- step05_n$end
step05_n$start <- NULL
step05_n$end <- NULL

step05_se$shared <- step05_se$celltype
step05_n <- merge(step05_n, step05_se, by = c('clone_id', 'shared'), all.x = T)
step05_n$start2 <- step05_n$start
step05_n$end2 <- step05_n$end
step05_n$start <- NULL
step05_n$end <- NULL

step05_n <- step05_n[, c('clone_id', 'celltype.x', 'start1', 'end1', 'shared', 'start2', 'end2')]

step05_n <- step05_n[, c('celltype.x', 'start1', 'end1', 'shared', 'start2', 'end2')]
colnames(step05_n) <- c("celltype1", "start1", "end1", "celltype2", "start2", "end2")
rm(step02_n, step05_se)

#### 3. save ####
step06 <- rbind(step05_2, step05_n)
rm(step05_2, step05_n)

get_order <- function(celltype, clone_id) {
  order <- predefined_orders[[celltype]]
  if (is.null(order)) {
    return(Inf)
  }
  suppressWarnings(match(clone_id, order, nomatch = Inf))
}

clone_df_ordered <- clone_df %>%
  split(.$celltype) %>%
  lapply(function(subset) {
    cell_type <- as.character(subset$celltype[1])
    subset$order <- get_order(cell_type, subset$clone_id)
    subset %>% arrange(order) %>% dplyr::select(-order)
  }) %>%
  bind_rows() %>%
  arrange(celltype)

clone_df_ordered <- do.call(rbind, lapply(split(clone_df_ordered, clone_df_ordered$celltype), function(x) {
  n <- nrow(x)
  x$start <- 0:(n-1)
  x$end <- 1:n
  return(x)
}))

write.csv(step06, paste0(title,"_link.csv"), row.names = F)
write.table(step06, file = paste0(title,"_link.txt"),quote = F,row.names = F,col.names = F,sep = "\t")
write.csv(clone_df_ordered, paste0(title,"_cells.csv"), row.names = F)
rm(step04, step06, order, dfs,title, clone_df, clone_df_ordered, predefined_orders, get_order)
