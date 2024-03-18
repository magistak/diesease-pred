library(tidyverse)
library(Hmisc)
library(purrr)
library(furrr)
library(patchwork)
library(tidymodels)
library(factoextra)
library(fastcluster)
library(irlba)
library(flashClust)
library(pvclust)
library(ggfortify)
library(dbscan)
library(diptest)

# read table for cox models
surv_all <- read_tsv("data/surv_all_log.tsv")
# create surv_all_c
met_list <- read_csv("data/biom_meta.csv") |> 
  filter(in_strict == 1) |> 
  select(met_name) |> 
  unlist(use.names = FALSE)
load("data/rdata/met_dis.RData")
# mets no need
nn_mets <- setdiff(mets_v, met_list)
csa <- colnames(surv_all)
surv_c_cols <- csa[!csa %in% nn_mets]
surv_c <- surv_all

# create table with no events
surv_c_events <- 
surv_c |> 
  rowwise() |> 
  mutate(t_events = sum(c_across(starts_with("event_"))))

met_table <- 
surv_c_events |> 
  filter(t_events == 0) |> 
  mutate(on_study=`time_all-cause mortality`) |> 
  select(-contains("event"), -contains("time_"), -on_study) |> 
  mutate(eid = as.character(eid), 
         sex = as.factor(sex),
         spectrometer = as.factor(spectrometer))

save(met_table, file = "data/met_prof_table_all.RData")
load("data/met_prof_table_all.RData")
load("data/filters/eids/no_disease_eids.RData")

# filter for healthy individuals
met_table <- 
surv_c |> 
  select(-contains("event"), -contains("time_")) |> 
  filter(eid %in% no_disease_eids) |> 
  mutate(eid = as.character(eid), 
         sex = as.factor(sex),
         spectrometer = as.factor(spectrometer))

data <- met_table
metabolite_names <- colnames(met_table)[5:ncol(met_table)]
# Create a function to perform regression for each metabolite
regress_metabolite <- function(metabolite_name, data=met_table) {
  # Filter the data for the current metabolite
  metabolite_data <- data %>% 
    select({{metabolite_name}}, sex, age, spectrometer) 
  colnames(metabolite_data)[1] <- "met"
  
  # Define the model specification
  metabolite_spec <- linear_reg() %>%
    set_engine("lm")
  
  # Create a recipe
  metabolite_recipe <- recipe(met ~ sex + age + spectrometer, data = metabolite_data)
  
  # Combine the recipe and model specification into a workflow
  metabolite_workflow <- workflow() %>%
    add_recipe(metabolite_recipe) %>%
    add_model(metabolite_spec)
  
  # Fit the model
  metabolite_fit <- metabolite_workflow %>% 
    fit(data = metabolite_data)
  
  # Predict the metabolite levels
  predicted_data <- predict(metabolite_fit, new_data = data) %>%
    rename({{metabolite_name}} := .pred)
  
  return(predicted_data)
}

# Apply the function to create a new dataset with predicted values
predicted_data <- map_dfc(metabolite_names, ~regress_metabolite(.x, met_table))
colnames(predicted_data) <- metabolite_names


# percent:
real_data <- met_table[ ,5:ncol(met_table)]
change_table <- (real_data-predicted_data)#/abs(predicted_data)
change_binary <- 
change_table |> 
  mutate_all(~ ifelse(. < 0, 0, 1)) |>  
  rowwise() |> 
  mutate(unique = paste(c_across(everything()), collapse = ""))
change_binary |> 
  group_by(unique) |> 
  summarise(n = n()) |> View()

pca_result <- prcomp(change_table, center = TRUE, scale. = TRUE)
pc_data <- as.data.frame(pca_result$x)

# Add original sample information if needed
pc_data <- bind_cols(metabolite_data, pc_data)

# Print or visualize the results
print(pc_data)
pc_data$SampleIdentifier <- rownames(pc_data)
ggplot(pc_data, aes(x = PC1, y = PC2, color = SampleIdentifier)) +
  geom_point() +
  labs(title = "PCA Scatter Plot", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()
# Calculate the proportion of variance explained by each principal component
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create a barplot
barplot(variance_explained, main = "Variance Explained by Principal Components", 
        xlab = "Principal Component", ylab = "Proportion of Variance Explained",
        col = "skyblue", names.arg = seq_along(variance_explained))
# Calculate cumulative variance explained
cumulative_variance <- cumsum(variance_explained)

# Create a line plot
plot(cumulative_variance, type = "b", pch = 16, col = "darkorange",
     main = "Cumulative Variance Explained by Principal Components",
     xlab = "Number of Principal Components", ylab = "Cumulative Proportion of Variance Explained")
biplot(pca_result, scale = 0)

colnames(pc_data)
# Create a scatter plot of PC1 vs PC2
ggplot(pc_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA Scatter Plot (PC1 vs PC2)", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()


# plot change table

pdf("plots/change_hist_logfiltered.pdf")
for(i in 1:length(colnames(change_table))){
plot_col <- colnames(change_table)[i]

p <- change_table |> 
  as_tibble() |> 
  ggplot(aes(!!sym(plot_col))) +
  geom_histogram(bins = 50)
  
print(p)
}
dev.off()

# examine ratios
change_stud <- 
change_table |> 
  as_tibble() |> 
  mutate(across(everything(), ~ (.-mean(.))/sd(.)))


combinations <- combn(colnames(change_stud), 2, simplify = TRUE)


# Create a new tibble with metabolite level ratios
ratio_tibble <- data.frame(matrix(ncol = 3, nrow = ncol(combinations)))
colnames(ratio_tibble) <- c("Metabolite1", "Metabolite2", "Ratio")

# Fill in the metabolite names
ratio_tibble$Metabolite1 <- combinations[1, ]
ratio_tibble$Metabolite2 <- combinations[2, ]

# Calculate ratios
ratio_tibble$Ratio <- apply(combinations, 2, function(pair) {
  change_stud[, pair[1]] / change_stud[, pair[2]]
})

# Print the resulting tibble
print(ratio_tibble)
original_table <- change_stud |> 
  mutate(eid = met_table$eid)

# Create all combinations of metabolites
metabolite_columns  <- colnames(change_stud)
ratios_table <- tibble(
  id = original_table$eid
)
# Generate all possible combinations of metabolite pairs
metabolite_combinations <- combn(metabolite_columns, 2, simplify = TRUE)

# Iterate through each combination and calculate the ratio
i = 1
for (i in 1:ncol(metabolite_combinations)) {
  metabolite_pair <- metabolite_combinations[, i]
  ratio_column_name <- paste(metabolite_pair, collapse = "_ratio")
  ratios_table[[ratio_column_name]] <- original_table[[metabolite_pair[1]]] / original_table[[metabolite_pair[2]]]
}
col <- colnames(ratios_table)[2]
ratios_table |> 
  ggplot(aes(x=!!sym(col))) +
  geom_histogram()
multimod <- tibble(pair=colnames(ratios_table)[-1])
multimod <- multimod |> 
  mutate(p=NA,
         stat=NA)
i=1
for(i in 1:nrow(multimod)){
dip=dip.test(ratios_table[[i+1]])
p=dip[2][[1]]
stat=dip[1][[1]]
multimod$p[i] <- p
multimod$stat[i <- stat]
}

write_tsv(multimod, "data/multimod.tsv")
write_tsv(ratios_table, file = "data/ratios_table.tsv")

multimod_signif <- 
multimod |> 
  filter(p < 0.05)

pdf("plots/multimod_pairs.pdf")
for(i in 1:nrow(multimod_signif)){
col <- multimod_signif$pair[i]
data <- select(ratios_table, !!sym(col))[[1]]
data <- data[data>-500 & data < 500]
hist(data, breaks = 1000, main = col)
}
dev.off()

ratios_table |> 
  ggplot(aes(x=!!sym(col))) +
  geom_histogram(bins = 200)
ratios_table |> 
  aes(x = )


pairs_t <- 
expand_grid(mets, mets) |> 
  filter(mets...1>mets...2) |> 
  mutate(pair=paste(mets...1,mets...2, by="_")) |> 
  rename(met1=mets...1,
         met2=mets...2)

metabolite_tibble <- 
  change_stud |> 
  mutate(eid = met_table$eid)

metabolite_long <- 
metabolite_tibble |> 
  pivot_longer(!eid)
  

all_pairs <- expand_grid(met_table$eid, pairs_t$pair)
colnames(all_pairs) <- c("eid", "pair")
all_pairs <- 
all_pairs |> 
  left_join(pairs_t)

all_pairs <- 
all_pairs |> 
  left_join(select(metabolite_long, eid, met1=name, value1=value)) |> 
  left_join(select(metabolite_long, eid, met2=name, value2=value)) |> 
  mutate(ratio=value1/value2)

all_pairs_w <- 
all_pairs |> 
  select(pair, ratio, eid) |> 
  pivot_wider(names_from = pair, values_from = ratio) |> 
  select(-eid)

i=2
pdf("plots/change_ratio_hist2.pdf")
for(i in 1:length(colnames(all_pairs_w))){
  plot_col <- colnames(all_pairs_w)[i]
  
  p <- all_pairs_w |> 
    ggplot(aes(!!sym(plot_col))) +
    geom_histogram(bins=50) +
    scale_x_continuous(trans = "log10")
  print(p)
}
dev.off()

# Use purrr to iterate over combinations and calculate ratios
ratio_tibble <- map_dfr(seq_along(combinations[1, ]), function(i) {
  pair <- combinations[, i]
  ratio <- metabolite_tibble[, pair[1]] / metabolite_tibble[, pair[2]]
  tibble(
    Individual = rownames(metabolite_tibble),
    Ratio = ratio
  )
})



ratio_tibble |> 
  head()

ratio_tibble |> 
  head()
# Print the resulting tibble
print(ratio_tibble)

# k means clustering

metabolite_data <- change_stud
metabolite_matrix <- as.matrix(change_stud)
k <- 4
# Step 3: Perform k-means clustering
kmeans_result <- kmeans(metabolite_matrix, centers = k, nstart = 30)

# Step 4: Assign cluster labels to your original tibble
metabolite_data$cluster <- as.factor(kmeans_result$cluster)

# Step 5: Explore the results
print(kmeans_result)
table(metabolite_data$cluster)
i=1
pdf("plots/hist_clust_6_logfiltered.pdf")
for(i in 1:(length(colnames(metabolite_data))-1)){
  plot_col <- colnames(metabolite_data)[i]
  
  p <- metabolite_data |> 
    as_tibble() |> 
    ggplot(aes(!!sym(plot_col),
               fill = cluster)) +
    geom_histogram(bins = 50, alpha = 0.5, position = "identity")
  
  print(p)
}
dev.off()

data <- 
metabolite_data |> 
  mutate(eid = met_table$eid) |> 
  pivot_longer(!c("eid", "cluster")) |> 
  mutate(cluster = as.character(cluster)) |> 
  group_by(cluster, name) |>
  summarise(mean_value= mean(value, na.rm = TRUE)) |> 
  ungroup() |> 
  pivot_wider(names_from = cluster, values_from = mean_value)
colnames(data) <- c("met", "cl1", "cl2", "cl3", "cl4")
hc <- hclust(dist(data[ , 2:5]), method = "complete")
cluster_assignments <- cutree(hc, k = 5) 
df_with_cluster <- data %>%
  mutate(cluster = as.factor(cluster_assignments))



df_with_cluster |> 
  pivot_longer(!c("met", "cluster")) |> 
  arrange(cluster) |> 
  mutate(cluster = as.character(cluster)) |> 
ggplot(aes(x=name, y =fct_reorder(met, cluster), fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab")
  

metabolite_data <- 
metabolite_data |> 
  mutate(eid = metabolite_tibble$eid,
         sex = met_table$sex,
         age = met_table$age)
metabolite_data |> 
  group_by(cluster, sex) |> 
  summarise(n = n())

disease_healthy <- 
surv_c |> 
  select(-contains("event"), -contains("time_"))

metabolite_matrix <- as.matrix(disease_healthy[ ,4:255])
k <- 2
# Step 3: Perform k-means clustering
kmeans_result <- kmeans(metabolite_matrix, centers = k, nstart = 30)

# Step 4: Assign cluster labels to your original tibble
disease_healthy$cluster <- as.factor(kmeans_result$cluster)

# Step 5: Explore the results
events_t <- 
surv_c |> 
  select(eid, contains("event"))

mean_columns <- colnames(events_t)[-1]
disease_healthy |> 
  left_join(events_t) |> 
  group_by(cluster,sex) |> 
  summarise(across(all_of(mean_columns), ~mean(.), .names = "av_{.col}")) |> 
  pivot_longer(!c("cluster", "sex")) |> 
  arrange(name) |> 
  pivot_wider(names_from = cluster, values_from = value) |> 
  rename(cluster1=`1`, cluster2=`2`) |> 
  mutate(ratio=cluster2/cluster1) |> 
  View()
  
disease_healthy |> 
  group_by(cluster, sex) |> 
  summarise(n=n())


# hclust
change_stud
metabolite_data <- change_stud
# Assuming your data frame is named 'metabolite_data' and has individuals as rows and metabolite levels as columns
# Example using PCA
pca_result <- prcomp(metabolite_data)
screeplot(pca_result, type = "lines")
biplot(pca_result, scale = 0)
# 3D Scatter plot with the first three principal components
scatterplot3d(pca_result$x[, 1:3], color = clusters, main = "3D Scatter Plot")
# Variable contributions for the first two principal components
fviz_contrib(pca_result, choice = "var", axes = 1:2)
# Observation contributions for the first two principal components
fviz_contrib(pca_result, choice = "ind", axes = 1:2)
# Pairs plot with the first few principal components
pairs(pca_result$x[, 1:3], col = clusters)
autoplot(pca_result)

reduced_data <- pca_result$x[, 1:10] |> 
  as_tibble()
# Adjust the number of principal components

distance_matrix <- dist(reduced_data, method = "euclidean")


# 1. Calculate the distance matrix
distance_matrix <- dist(metabolite_data, method = "euclidean")

# 2. Perform hierarchical clustering
hierarchical_cluster <- hclust(distance_matrix, method = "ward.D")

# 3. Plot the dendrogram
plot(hierarchical_cluster, main = "Hierarchical Clustering Dendrogram", xlab = "Individuals", sub = NULL)

# 4. Cut the dendrogram to form clusters (you can choose the number of clusters)
clusters <- cutree(hierarchical_cluster, k = 3)  # Change 'k' to the desired number of clusters

# 5. Add the cluster assignments to your original data frame
metabolite_data_with_clusters <- cbind(metabolite_data, Cluster = clusters)

# 6. View or analyze the clustered data
print(metabolite_data_with_clusters)

# Identify the most important variables (e.g., top 10) based on loadings
top_variables <- rownames(pca_result$rotation)[order(-abs(pca_result$rotation[, 1]))][1:10]

# Biplot with the first two principal components, highlighting top variables
biplot(pca_result, scale = 0, choices = c(1, 2), xlabs = top_variables)


# Extract the relevant columns from the original data
top_variables_data <- metabolite_data[, top_variables]

# Perform PCA on the subset of variables
pca_result_subset <- prcomp(top_variables_data)

# Load the factoextra library
svd_result <- irlba(t(metabolite_data), nv = 10)  # Adjust 'nv' as needed

hierarchical_cluster <- fastcluster::hclust(dist(svd_result$v), method = "ward.D")

hierarchical_cluster <- flashClust::flashClust(as.dist(dist(reduced_data)), method = "average")

clusters <- cutree(pca_result_subset, k = 3)
# Create a biplot with the first two principal components
fviz_pca_biplot(pca_result_subset, col.var = "blue", col.ind = clusters, repel = TRUE)

hc <- pvclust::pvclust(t(reduced_data), method.hclust = "average")
plot(hc)

clusters <- cutree(hc$hclust, k = 3)

svd_result_df <- as.data.frame(svd_result$v)

fviz_pca_biplot(svd_result_df, col.var = "blue", col.ind = clusters, repel = TRUE)
biplot <- ggbiplot(svd_result$v, labels = rownames(svd_result$v), groups = clusters, ellipse = TRUE, circle = TRUE)
par(mar = c(5, 4, 4, 2) + 0.1)
biplot(svd_result$v, svd_result$u, cex = 0.8, col = clusters)

scatterplot3d(svd_result$v[, 1:3], color = clusters, pch = 16, cex.symbols = 0.8)
clusters
dbscan_result <- dbscan(metabolite_data, eps = 0.5, MinPts = 5)

# pca on people:
change_stud
t_change <- t(change_stud) |> 
  as_tibble()
pca_result <- prcomp(t_change, center = TRUE, scale. = TRUE)
screeplot(pca_result, type = "lines")
pca_result
ve <- pca_result$sdev^2 / sum(pca_result$sdev^2)
options("scipen"=100, "digits"=4)
ve*100
sum(ve)                    
pca_result$x[,1:4] |> 
  as_tibble()

pca_res <- 
t(pca_result$x[,1:5]) |> 
  as_tibble()
colnames(pca_res) <- colnames(change_stud)
cor_matrix_t <- 
cor(pca_res) |> 
  as_tibble()
cor_matrix <- as.matrix(cor_matrix_t)
rownames(cor_matrix) <- colnames(change_stud)
heatmap(cor_matrix,
        symm = FALSE,       # Ensure the heatmap is symmetric
        margins = c(10,10), # Add margins to the plot
        col = colorRampPalette(c("blue", "white", "red"))(50))

rownames(cor_matrix_t) <- colnames(change_stud)





# k means clustering

metabolite_data <- change_stud
metabolite_matrix <- as.matrix(change_stud)
k <- 6
# Step 3: Perform k-means clustering
kmeans_result <- kmeans(metabolite_matrix, centers = k, nstart = 20)
?kmeans

wcss_values <- numeric()

for (i in 1:10) {
  kmeans_result <- kmeans(metabolite_matrix, centers = i, nstart = 5)
  wcss_values[i] <- kmeans_result$tot.withinss
}

# Plotting the Elbow Method
plot(1:10, wcss_values, type = "b", xlab = "Number of Clusters (k)", ylab = "WCSS")
