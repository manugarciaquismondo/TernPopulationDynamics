

get.values = function(file.name) {
    read.csv(paste(file.name, "csv",sep = "."), header = T)[, 2]
}

simulation.batch.id="juvenile"

# Read results
setwd(paste(simulation.batch.id, "metrics", sep="_"))
simulation.max.values = get.values("maxValues")
simulation.min.values = get.values("minValues")
simulation.max.values.large.colonies = get.values("maxValuesLargeColonies")
simulation.min.values.large.colonies = get.values("minValuesLargeColonies")
colonization.numbers = get.values("colonizations")
decolonization.numbers = get.values("decolonizations")
proportion.large.colonies = get.values("proportion_large_colonies")
proportion.large.colonies = proportion.large.colonies[seq(from = 50, to = length(proportion.large.colonies), by = 50)]
final.large.colonies=get.values("final_large_colonies")
final.medium.colonies = get.values("final_medium_colonies")

# Plot results
setwd("..")
setwd(paste(simulation.batch.id, "images", sep="_"))
jpeg("minValues.jpeg")
hist(simulation.min.values, ylab = "Minimum number of pairs in the simulation", xlab = "Number of simulations", main="Minimum number of pairs per simulation")
dev.off()
jpeg("maxValues.jpeg")
hist(simulation.max.values, ylab = "Maximum number of pairs in the simulation", xlab = "Number of simulations", main = "Maximum number of pairs per simulation")
dev.off()
jpeg("colonizations.jpeg")
hist(colonization.numbers, ylab = "Number of colonized sites", xlab = "Number of pairs", main = "Number of pairs after colonization", breaks = c(seq(from = 0, to = max(colonization.numbers), by = 2), max(colonization.numbers)))
dev.off()
jpeg("decolonizations.jpeg")
hist(decolonization.numbers, ylab = "Number of abandoned sites", xlab = "Number of pairs", main = "Number of pairs before abandonment")
dev.off()
jpeg("proportion_large_colonies.jpeg")
hist(proportion.large.colonies, ylab = "Number of simulations", xlab = "Proportion of pairs", main = "Proportion of pairs in large colonies in the last year")
dev.off()
jpeg("final_large_colonies.jpeg")
hist(final.large.colonies, ylab = "Number of simulations", xlab = "Number of colonies", main = "Small colonies that maintain large size", breaks = c(unique(final.large.colonies),2))
dev.off()
jpeg("final_medium_colonies.jpeg")
hist(final.medium.colonies, ylab = "Number of simulations", xlab = "Number of colonies", main = "Small colonies that maintain medium size", breaks = unique(final.medium.colonies))
dev.off()
jpeg("minValuesLargeColonies.jpeg")
hist(simulation.min.values.large.colonies, ylab = "Minimum number of pairs in the simulation", xlab = "Number of simulations", main="Minimum number of pairs per simulation in large colonies only")
dev.off()
jpeg("maxValuesLargeColonies.jpeg")
hist(simulation.max.values.large.colonies, ylab = "Maximum number of pairs in the simulation", xlab = "Number of simulations", main="Maximum number of pairs per simulation in large colonies only")
dev.off()


# Return to root directory
setwd("..")
