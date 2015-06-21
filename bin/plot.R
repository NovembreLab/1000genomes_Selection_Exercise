library(dplyr)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
chr <- args[1]
start <- args[2]
end <- args[3]
clst1 <- args[4]
clst2 <- args[5]
clst3 <- args[6]
ihs_clst <- args[7]
region <- paste0("chr", chr, "_", start, "_", end)

# Plot Fst
fst_12_df <- read.table(paste0(clst1, "_", clst2, "_", region, "_filtered.fst"), header = TRUE) %>% 
                  rename(FST_12 = FST) %>% 
                  filter(is.nan(FST_12) == FALSE) %>% 
                  select(SNP, POS, FST_12)

p <- ggplot(fst_12_df, aes(x = POS, y = FST_12)) + geom_point() + xlab('Position on chr2') + 
     theme_bw() + ylab("Fst_12") 

# Plot PBS
fst_13_df <- read.table(paste0(clst1, "_", clst3, "_", region, "_filtered.fst"), header = TRUE) %>% 
             rename(FST_13 = FST) %>% 
             filter(is.nan(FST_13) == FALSE) %>% 
             select(SNP, POS, FST_13)

fst_23_df <- read.table(paste0(clst2, "_", clst3, "_", region, "_filtered.fst"), header = TRUE) %>% 
             rename(FST_23 = FST) %>% 
             filter(is.nan(FST_23) == FALSE) %>% 
             select(SNP, POS, FST_23)

pbs_df <- fst_12_df %>%
          inner_join(fst_13_df, by = c("SNP", "POS")) %>%
          inner_join(fst_23_df, by = c("SNP", "POS")) %>%
          mutate(T_12 = -log(1 - FST_12), T_13 = -log(1 - FST_13), T_23 = -log(1 - FST_23)) %>%
          mutate(PBS_1 = (T_12 + T_13 - T_23) / 2) %>%
          mutate(PBS_2 = (T_12 + T_13 - T_23) / 2) %>%
          select(POS, PBS_1, PBS_2)

melted_pbs_df <- melt(pbs_df, id = c("POS"))

p1 <- ggplot(melted_pbs_df, aes(x = POS, y = value)) + geom_point() + xlab('Position on chr2') + 
      theme_bw() + ylab("") + facet_grid(variable ~ .)

# Plot iHS
ihs_df <- read.table(paste0(ihs_clst, "_", region, "_filtered.", "ihs.out.100bins.norm"), header = FALSE)
colnames(ihs_df) <- c("SNP", "POS", "FREQ", "ihh1", "ihh0", "unstd_iHS", "std_iHS", "top")

# Take the absolute value of iHS
selection_df <- ihs_df %>% mutate(std_iHS = abs(std_iHS))

# plot Fst, PBS, and iHS
p2 <- ggplot(selection_df, aes(x = POS, y = std_iHS)) + geom_point() + xlab('Position on chr2') + 
      theme_bw() + ylab("|normalized iHS|") 

pdf("./test.pdf", width = 7, height = 4)
p
p1
p2
dev.off()