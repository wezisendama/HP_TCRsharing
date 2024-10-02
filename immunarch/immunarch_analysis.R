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

vis(bal.volume, .by = c("disease"), .meta = kaminski_data.bal$meta)
vis(bal.clones, .by = c("disease"), .meta = kaminski_data.bal$meta)
vis(pbmc.volume, .by = c("disease"), .meta = kaminski_data.pbmc$meta)
vis(pbmc.clones, .by = c("disease"), .meta = kaminski_data.pbmc$meta)
