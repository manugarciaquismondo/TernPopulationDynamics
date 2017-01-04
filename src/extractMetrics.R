# Read the directory where the simulation files are and the file route where the output file will be stored

inputDirectory = commandArgs(T)[1]
outputDirectory = commandArgs(T)[2]
parametersRoute = commandArgs(T)[3]

medium.colony.lowerbound = 25
large.colony.lowerbound = 250

# Read the names of the colonies
setwd(parametersRoute)
colonyNames = read.csv("allpairs.csv", header = T, stringsAsFactors = F)[, 1]
colonyNames[length(colonyNames)] = "Young's I., Smithtown, NY"
largeColonies = read.csv("survivalestimates.csv", header = F, stringsAsFactors = F)[, 1]
# Read the birds table

setwd(inputDirectory)

total.birds.table = read.csv("TotalBirdsNextYear.csv", header = T)[, -1]
rownames(total.birds.table) = colonyNames

# Get maximum and minimum number of birds

max.min.num.pairs = range(colSums(total.birds.table))
max.min.num.pairs.matrix = matrix(nrow = 2, ncol = 1, data = max.min.num.pairs)
rownames(max.min.num.pairs.matrix) = c("MinPairs", "MaxPairs")
# Get proportion of population in large colonies

proportion.large.colonies = colSums(total.birds.table[largeColonies,]) / colSums(total.birds.table)

# Get colonization and decolonization statistics

colonization.events = decolonization.events = matrix(nrow = nrow(total.birds.table), ncol = ncol(total.birds.table), data = 0)
rownames(colonization.events) = rownames(decolonization.events) = rownames(total.birds.table)
for (colony.index in rownames(total.birds.table)) {
    for (year.index in 1:(ncol(total.birds.table) - 1)) {
        if (total.birds.table[colony.index, year.index] == 0 && total.birds.table[colony.index, year.index + 1] > 0) {
            colonization.events[colony.index, year.index + 1] = total.birds.table[colony.index, year.index + 1]
        }
        if (total.birds.table[colony.index, year.index] > 0 && total.birds.table[colony.index, year.index + 1] == 0) {
            decolonization.events[colony.index, year.index] = total.birds.table[colony.index, year.index]
        }
    }
}

# Get colonies that are initially empty

initially.empty.colonies = rownames(total.birds.table)[total.birds.table[, 1] == 0]

# Exclude large colonies
initially.empty.colonies = setdiff(initially.empty.colonies, largeColonies)

# Calculate initially empty colonies that become medium-sized and large
initially.empty.colonies.table = total.birds.table[initially.empty.colonies,]

last.year.populations = initially.empty.colonies.table[, ncol(total.birds.table)]
final.medium.colonies = initially.empty.colonies.table[last.year.populations > medium.colony.lowerbound,]
final.large.colonies = initially.empty.colonies.table[last.year.populations > large.colony.lowerbound,]

# Save statistics matrices

setwd(outputDirectory)

write.csv(max.min.num.pairs.matrix, "MinMaxPairs.csv")
write.csv(matrix(proportion.large.colonies, nrow = length(proportion.large.colonies), ncol = 1), "ProportionLargeColonies.csv")
write.csv(colonization.events, "ColonizationEvents.csv")
write.csv(decolonization.events, "DecolonizationEvents.csv")
write.csv(final.medium.colonies, "FinalMediumColonies.csv")
write.csv(final.large.colonies, "FinalLargeColonies.csv")