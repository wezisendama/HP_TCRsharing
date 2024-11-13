library(immunarch)
file_path = "."

# Load repertoires
kaminski_data <- repLoad(file_path)

# Include only TCRb chain data
kaminski_data <- repFilter(kaminski_data, .method = "by.clonotype", .query = list(V.name = include("TRBV")), .match = "startswith")

# Exclude incompletely sequenced CDR3 regions
kaminski_data <- repFilter(kaminski_data, .method = "by.clonotype", .query = list(CDR3.aa = exclude("?")), .match = "substring")
kaminski_data <- repFilter(kaminski_data, .method = "by.clonotype", .query = list(CDR3.aa = exclude("partial", "out_of_frame")))

# Subset BAL data
kaminski_data.bal <- repFilter(kaminski_data, .method = "by.meta", .query = list(tissue = include("BAL")))

# Subset PBMC data
kaminski_data.pbmc <- repFilter(kaminski_data, .method = "by.meta", .query = list(tissue = include("PBMC")))

# Repertoire size metrics
overall.volume <- repExplore(kaminski_data$data, .method = "volume")
overall.clones <- repExplore(kaminski_data$data, .method = "clones")
vis(overall.volume, .by = c("tissue"), .meta = kaminski_data$meta)
vis(overall.clones, .by = c("tissue"), .meta = kaminski_data$meta)

bal.volume <- repExplore(kaminski_data.bal$data, .method = "volume")
bal.clones <- repExplore(kaminski_data.bal$data, .method = "clones")
pbmc.volume <- repExplore(kaminski_data.pbmc$data, .method = "volume")
pbmc.clones <- repExplore(kaminski_data.pbmc$data, .method = "clones")

figure1A <- vis(bal.clones, .by = c("disease"), .meta = kaminski_data.bal$meta) # Figure 1A
figure1B <- vis(bal.volume, .by = c("disease"), .meta = kaminski_data.bal$meta) # Figure 1B
figure1C <- vis(pbmc.clones, .by = c("disease"), .meta = kaminski_data.pbmc$meta) # Figure 1C
figure1D <- vis(pbmc.volume, .by = c("disease"), .meta = kaminski_data.pbmc$meta) # Figure 1D

# used fixVis() to tidy up the charts. Charts render with Bonferroni-Holm adjusted p values
# which are reasonable to use in the context given three pairwise comparisons. Adjusted p
# values of 1 are therefore valid, but might unsettle some people. Straightforward p values
# therefore given with the below tables:

figure1A.stats <- figure1A$plot_env$p_df
figure1B.stats <- figure1B$plot_env$p_df
figure1C.stats <- figure1C$plot_env$p_df
figure1D.stats <- figure1D$plot_env$p_df

# Worth noting that the conclusions do not change even with unadjusted p values.
