

get.values = function(file.name) {
    read.csv(paste(file.name, "csv",sep = "."), header = T)[, 2]
}

# Read results
setwd("downloaded_metrics")
simulation.max.values = get.values("maxValues")
simulation.min.values = get.values("minValues")
colonization.numbers = get.values("colonizations")
decolonization.numbers = get.values("decolonizations")
proportion.large.colonies = get.values("proportion_large_colonies")
proportion.large.colonies = proportion.large.colonies[seq(from = 50, to = length(proportion.large.colonies), by = 50)]
final.large.colonies=get.values("final_large_colonies")
final.medium.colonies = get.values("final_medium_colonies")

# Plot results
setwd("../downloaded_images")
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
hist(final.large.colonies, ylab = "Number of simulations", xlab = "Number of colonies", main = "Small colonies that maintain large size", breaks = unique(final.large.colonies))
dev.off()
jpeg("final_medium_colonies.jpeg")
hist(final.medium.colonies, ylab = "Number of simulations", xlab = "Number of colonies", main = "Small colonies that maintain medium size", breaks = unique(final.medium.colonies))
dev.off()

# Return to root directory
setwd("..")