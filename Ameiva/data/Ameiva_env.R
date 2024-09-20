# Packages####
require(dartR)
require(geodata)

# Loading Variables####
setwd("/home/users/jnadaline/Ameiva/")
# gl.a<- gl.read.vcf("Lampro_ddRAD_filtered.recode.vcf")
gl.a<- gl.read.vcf("denovo_261_final.vcf")
gl.write.csv(gl.a,"gl_ameiva")
pop.a <- read.csv("pop261-dnv.csv")

pop(gl.a) <- 1:261
gl.a@other$latlongs <- pop.a[, 2:3]
head(gl.a@other$latlongs)
# making the pop column of the pop data frame into categorical data 
pop.a$pop <- as.factor(pop.a$POP)
# assigning pop info to sites 
gl.a@other$ind.metrics$sites <- pop.a$pop
# assigning site info to pop
gl.a@pop <- gl.a@other$ind.metrics$sites
nPop(gl.a)
gl.map.interactive(gl.a)
gl

# Download env layers ####
require(terra)
require(geodata)

setwd("/home/users/jnadaline/Coobirds/chelsa_current")
files=list.files()

climatic=rast(files)
occ=data.frame(pop.a.a$lon,pop.a.a$lat)

#Extracting raster values for each point and climatic layer
plot(occ)
clim.val=extract(climatic,occ)
row.names(x = clim.val)=pop.a.a$ind
#Changing the colnames to a friendly name using regex
colnames(clim.val)= gsub("CHELSA_", "", colnames(clim.val), fixed=T)
colnames(clim.val)= gsub("_1981.2010_V.2.1", "", colnames(clim.val), fixed=T)
colnames(clim.val)
write.csv(clim.val,"/home/users/jnadaline/Lampropholis/Clim-Lampro.csv")
clim.val=read.csv("/home/users/jnadaline/Lampropholis/Clim-Lampro.csv")

# Manipulating Net primary productivity variable
setwd("/home/users/jnadaline/Coobirds/NPP")
untar("NPP2018.tar.gz")
prod=rast(list.files(pattern="*tif"))
prodm=mean(prod) #mean of 2018 primary productivity time series, othe option is use mean in NPP of since the oldest NPP raster

#extracting values
prod.val=extract(prodm,occ)
row.names(x = prod.val)=pop.a.a$ind
colnames(prod.val)
write.csv(prod.val,"/home/users/jnadaline/Lampropholis/Prod-Lampro.csv")

prod.val=read.csv("/home/users/jnadaline/Lampropholis/Prod-Lampro.csv")
pred=data.frame(clim.val,prod.val$mean)

# RDA ####
require(psych)    # Used to investigate correlations among predictors
require(vegan)
require(caret)
require(zoo)

# Remove correlated environmental data
pairs.panels(pred, scale=T)
p32 <- cor(pred[,3:22])
re  <- findCorrelation(p32,cutoff = 0.9,verbose = FALSE,names = FALSE,exact = TRUE)
pred2 <- pred[,-re]
pred2
pred2[,-1]
colnames(pred2)= gsub("_1981.2010_V.2.1", "", colnames(pred2), fixed=T)
pairs.panels(pred2, scale=T)

# Format genetic data for RDA and interpolate missing data
#imputation of neighbour
gl.a<- gl.impute(gl,method="neighbour")
# Sequence Tag presence-absence data
tess_ES <- data.frame(gl.a)
sum(is.na(tess_ES))

# Replace missing genotype data by column (locus); missing genotypes replaced with most common genotype
#insert the 
find_mode <- function(x) {
  u   <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}
gen.imp  <- suppressWarnings(na.aggregate(tess_ES,FUN=find_mode))
sum(is.na(gen.imp))

## Other option is exclude NA?
#Running RDA

colnames(pred2) = c(
  "bio12" = "TempAnnRange",        # Temperature Annual Range
  "bio13" = "PrecipWetQtr",         # Precipitation of Wettest Quarter
  "bio14" = "PrecipDryQtr",         # Precipitation of Driest Quarter
  "bio17" = "PrecipWetMonth",       # Precipitation of Wettest Month
  "bio19" = "PrecipDryMonth",       # Precipitation of Driest Month
  "bio2"  = "MeanDiurnalRange",     # Mean Diurnal Temperature Range
  "bio3"  = "Isothermality",        # Isothermality
  "bio4"  = "TempSeasonality",      # Temperature Seasonality
  "bio6"  = "MinTempColdestMonth",  # Minimum Temperature of Coldest Month
  "bio7"  = "TempAnnualMean",       # Temperature Annual Mean
  "bio8"  = "TempWettestQuarter",   # Temperature of Wettest Quarter
  "bio9"  = "TempDriestQuarter",    # Temperature of Driest Quarter
  "prod.val.mean" = "ProdMean"      # Mean Production Value
)

lampro.rda <- rda(gl.a~ ., data=as.data.frame(pred2), scale=T)

RsquareAdj(lampro.rda)

summary(eigenvals(lampro.rda, model = "constrained"))
screeplot(lampro.rda)


# estimate significance
signif.full <- anova.cca(lampro.rda, parallel=getOption("mc.cores"))

plot(lampro.rda)

levels(env$ecotype) <- c("Western Forest","Boreal Forest","Arctic","High Arctic","British Columbia","Atlantic Forest")
eco <- env$ecotype
eco=gl.a$pop

library(RColorBrewer)

# Combine paletas para obter 23 cores
palette1 <- brewer.pal(n = 12, name = "Paired")
palette2 <- brewer.pal(n = 8, name = "Set2")
palette3 <- brewer.pal(n = 3, name = "Set1")

# Combine todas as paletas
combined_palette <- c(palette1, palette2, palette3)

# Ajuste o número de cores para 23
bg <- combined_palette[1:23]

# 6 nice colors for our ecotypes
plot(lampro.rda, type="n", scaling=3)
points(lampro.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(lampro.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
text(lampro.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

### #Super plot
require(cowplot)
require(ggplot2)
install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

# Criar o gráfico RDA (ajuste conforme seus dados reais)
x=pop.a.a$lon
y=pop.a.a$lat
rda_plot <- ggplot() +
  geom_point(data = lampro.rda$species, aes(x = x, y = y), color = "gray32", size = 1.2) + # SNPs
  geom_point(data = lampro.rda$sites, aes(x = x, y = y, fill = eco), shape = 21, size = 3) + # Amostras
  # geom_text(data = lampro.rda$biplot, aes(x = x, y = y, label = label), color = "#0868ac", size = 3) + # Preditores
  theme_minimal() +
  theme(legend.position = "bottom")

# Criar o mapa das amostras (ajuste conforme seus dados reais)
map_plot <- ggplot(pop[,3:4], aes(x = x, y = y, color = eco)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(legend.position = "none") # Remove a legenda do mapa

# Combine o gráfico RDA com o mapa
combined_plot <- grid.arrange(
  rda_plot,
  map_plot,
  widths = c(4, 1), # Ajuste as larguras conforme necessário
  layout_matrix = rbind(c(2, 1)) # Layout com o gráfico principal e o mapa
)

# Exibir o gráfico combinado
print(combined_plot)


####


plot(lampro.rda, type="n", scaling=3)
points(lampro.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # SNPs
points(lampro.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # Amostras
text(lampro.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # Preditores
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# Criar o mapa das amostras
map_plot <- ggplot(pop[,3:4], aes(x=x, y=y, color=eco)) +
  geom_point(size=3) +
  theme_minimal() +
  theme(legend.position = "none") # Remove a legenda do mapa

# Combine o gráfico RDA com o mapa
final_plot <- ggdraw() +
  draw_plot(plot(lampro.rda, type="n", scaling=3) + 
              points(lampro.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3) + 
              points(lampro.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) + 
              text(lampro.rda, scaling=3, display="bp", col="#0868ac", cex=1) + 
              legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg),
            x = 0, y = 0, width = 1, height = 1) +
  draw_plot(map_plot, x = 0.7, y = 0.7, width = 0.2, height = 0.2) # Ajuste a posição e tamanho do mapa conforme necessário

print(final_plot)