library(factoextra)
library(rFIA)
library(sf)
library(dplyr)
library(tidyr)
library(quantreg)


### LOAD DATA 

# load species data habitat amount, fragemntaiton traits, 
f <- readRDS("Data/post-nb-np_calcs-near-20260520.RDA")

# Puerto Rico outline 
pr <- st_read("Data/PR_Otline/PR_outline_Project.shp")

# Download FIA data for Puerto Rico (one time only)
# PR_Trees <- getFIA(states = 'PR', dir = 'Data/FIA/bob')
# Load data (excluding points off of mainland PR (e.g., Vieques, Culebra))

FIA <- readFIA("Data/FIA/")

# load the downloaded abundance data 
prfia <- clipFIA(FIA, mask=pr, mostRecent=F)

# Estimate total population size by species for trees with diameter of at least 2 inches
pr_tpa_diam2 <- tpa(prfia, bySpecies=TRUE, totals=TRUE, treeDomain=DIA>2)

# Select data from 2004 and 2014
pr_tpa_diam2 <- pr_tpa_diam2[pr_tpa_diam2$YEAR %in% c(2004, 2014),]

# Filter the data frame and use only the columns
pr_tpa_diam2_sel <- pr_tpa_diam2[,c("YEAR", "SPCD", "SCIENTIFIC_NAME", "TPA")]

A <- pr_tpa_diam2_sel %>%
  pivot_wider(
    names_from = YEAR,
    values_from = TPA,
    names_prefix = "TPA_"
  )

# Format species names 
f$X <- gsub("_", " ", f$X)

# Give spcode to FIA data
A$code <- f$sp[match(A$SCIENTIFIC_NAME, f$X)]

# Merge the abundance data to the habitat data 
H <- f[,c("X","sp","fcover_51_raw", "fcover_00_raw", "tot_change_raw", "enn_51", "enn_00", "enn_change",
            "ss.log.z", "thk.log.z", "la.log.z", "sla.log.z", "maxht.z", "wd.z")]

# Remove species with NA values for the traits 
df <- na.omit(H)

# Merge the abundance data 
df$TPA_2004 <- A$TPA_2004[match(df$sp, A$code)]
df$TPA_2014 <- A$TPA_2014[match(df$sp, A$code)]

sum(is.na(df$TPA_2004))
sum(is.na(df$TPA_2014))
                          
df <- df[!(is.na(df$TPA_2014) & is.na(df$TPA_2014)), ]

df$TPA_2004[is.na(df$TPA_2004)] <- 0
df$TPA_2014[is.na(df$TPA_2014)] <- 0

df$TPA_change <- df$TPA_2004 - df$TPA_2014

# Run the trait PCA
pcadata <- df[, c("ss.log.z","thk.log.z","la.log.z",      
                 "sla.log.z","maxht.z", "wd.z")]

pca <- prcomp(pcadata, center=TRUE, scale.=TRUE)
summary(pca)

# Plot the PCA
fviz_pca_biplot(pca, label = "var", col.var = "contrib", repel = FALSE, 
  pointshape = 16, pointsize = 2.4,col.ind = "gray40") +
  scale_color_gradient(low  = "seagreen", high = "sienna1", name = "Contribution (%)"
  ) +
  xlab("PC1 (39.9%)") +
  ylab("PC2 (22.4.9%)") +
  ggtitle(NULL) +   
   theme_classic() + 
    theme(
  axis.title = element_text(size = 16),
  axis.text  = element_text(size = 14),
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 14)
)

# Save PCA axis scores 
scores <- as.data.frame(pca$x)
scores$sp <- df$sp

# Add the loading on the axis on the data frame
df$PC1 <- scores$PC1[match(df$sp,scores$sp)]
df$PC1.z <- scale(df$PC1)

# scale the variables for the models 
df$fcover_51_raw.z <- scale(df$fcover_51_raw)
df$tot_change_raw.z <- scale(df$tot_change_raw) 
df$enn_51.z <- scale(df$enn_51)
df$enn_change.z <- scale(df$enn_change)


## How does the remaining habitat in 1951 relate to current species abundance when accounting for species’ life-history traits?

# Remove remove abundance values that are 0 
df_51 <- df[df$TPA_2014!= 0, ]

# run the quantile regression 
qr <- rq(TPA_2014 ~ fcover_51_raw.z + PC1.z + enn_51.z + fcover_51_raw.z * PC1.z, tau = .95, data = df_51)
summary(qr)

#interaction plot 
par(mar = c(5, 5, 1, 2))
plot(
  df_51$fcover_51_raw.z,
  df_51$TPA_2014,
  col = adjustcolor("grey30", alpha.f = 0.5),
  pch = 16,
  cex = 1.3,
  xlab = "Potential suitable Hhbitat ammount (1951)",
  ylab = "FIA abundnace (trees per acres, 2014)",
  cex.lab = 1.3,
)

line_cols <- c("indianred2", "steelblue3")
names(line_cols) <- c("PC1 = acquisitive", "PC1 = conservative")

fcover_seq <- seq(
  min(df_51$fcover_51_raw.z),
  max(df_51$fcover_51_raw.z),
  length.out = 100
)

#min and max pca values ---> -2.682475 , 2.298517
min(df_51$PC1.z)
max(df_51$PC1.z)

PC1_vals <- c(-2.6, 2.2)

for (i in seq_along(PC1_vals)) {
  
  A_hat <- 6.03570 +
    3.07685 * fcover_seq +
    (-1.51119) * PC1_vals[i] +
    0.50221 * 0 +                          # enn held at mean
    (-1.90568) * fcover_seq * PC1_vals[i]
  
  lines(
    fcover_seq,
    A_hat,
    col = line_cols[i],
    lwd = 2
  )
}

legend(
  "topleft",
  legend = c("Acquisitive species", "Conservative species"),
  col = line_cols,
  lwd = 2,
  bty = "n"
)

# How did past changes in species' potential suitable habitat (from 1951-2000) 
#relate to changes in contemporary species abundance when accounting for species 
#life history traits?  

lm<- lm(TPA_change ~ tot_change_raw.z * PC1.z + enn_change.z , data = df)
summary(lm)

par(mfrow = c(1, 3),        
    mar = c(5, 5, 3, 2),   
    oma = c(0, 0, 2, 0))

plot(df$tot_change_raw.z, df$TPA_change,
     pch = 16,
     col = rgb(0, 0, 0, 0.3),
     cex = 2,
     cex.axis=1.5,
     xlab = "sacled potencial suitable habitat change 1951 - 2000",
     ylab = "FIA abundnace change (2004 – 2014)",
     cex.lab = 1.8
)

abline(lm(TPA_change ~ df$tot_change_raw.z,
          data = df), lty = 2, lwd = 2)



plot(df$enn_change.z, df$TPA_change,
     pch = 16,
     col = rgb(0, 0, 0, 0.3),
     cex = 2,
     cex.axis=1.5,
     xlab = "scaled nearest neighbor distance change 1951 - 2000",
     ylab = "FIA abundnace change (2004 – 2014)",
     cex.lab = 1.8
)

abline(lm(TPA_change ~ df$enn_change.z,
          data = df), lty = 2, lwd = 2)

plot(df$PC1.z, df$TPA_change,
     pch = 16,
     col = rgb(0, 0, 0, 0.3),
     cex = 2,
     cex.axis=1.5,
     xlab = "PC1",
     ylab = "FIA abundnace change (2004 – 2014)",
     cex.lab = 1.8
)

abline(lm(TPA_change ~ df$PC1.z,
          data = df), lty = 2, lwd = 2)






