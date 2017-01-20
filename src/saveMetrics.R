source.directory = commandArgs(T)[1]
destination.directory = commandArgs(T)[2]

# Set route to the stastistics from the simulations
setwd(source.directory)

# Set the number of years in which a small colony needs to remain medium or large to be considered medium or large

large.stable.years=10

# A function to check if a colony has stayed large for a long time
check.large.colony=function(input.table, input.large.stable.years=large.stable.years){
  
  simulated.years=ncol(input.table)
  large.colony.count=0
  #If there exist at least one colony that has maintained large size at the end of the simulation
  if(nrow(input.table)>0){
    for(large.colony in 1:nrow(input.table)){
      large.colony.vector=input.table[large.colony,]
      #If the colony has stayed large for the required number of years, add 1 to the count
      if(all(((simulated.years-input.large.stable.years):simulated.years) %in% which(large.colony.vector>0))){
        large.colony.count=large.colony.count+1
      }
      
    }
  }
  large.colony.count
}

# Initialize structures for simulation metrics
simulation.min.values = simulation.max.values = simulation.min.values.large.colonies = simulation.max.values.large.colonies = c()
proportion.large.colonies=c()
colonization.numbers=decolonization.numbers=c()
final.large.colonies=final.medium.colonies=c()

# For each simulation directory that is not the same or the parent, enter the directory
for(simulation in list.dirs(recursive = F)){
    if (!simulation %in% c(".", "..")) {
      setwd(file.path(source.directory, simulation))
    
    # Register maximum and minimum number of pairs in all colonies across all years
    minMaxTable=read.csv("MinMaxPairs.csv", header=T)
    simulation.min.values=c(simulation.min.values,minMaxTable[1,2])
    simulation.max.values=c(simulation.max.values,minMaxTable[2,2])
    
    # Register maximum and minimum number of pairs in large colonies across all years
    minMaxTableLargeColonies=read.csv("MinMaxPairsLargeColonies.csv", header=T)
    simulation.min.values.large.colonies=c(simulation.min.values.large.colonies,minMaxTableLargeColonies[1,2])
    simulation.max.values.large.colonies=c(simulation.max.values.large.colonies,minMaxTableLargeColonies[2,2])
    
    # Register number of pairs after colonization and before decolonization
    colonization.table=read.csv("ColonizationEvents.csv", header = T)[,-1]
    decolonization.table=read.csv("DecolonizationEvents.csv", header = T)[,-1]
    colonization.numbers=c(colonization.numbers,colonization.table[colonization.table>0])
    decolonization.numbers=c(decolonization.numbers,decolonization.table[decolonization.table>0])
    
    # Register the proportion of the population placed in large colonies
    proportion.large.colonies=c(proportion.large.colonies, read.csv("ProportionLargeColonies.csv", header = T)[,2])
    
    # Register the number of small colonies that become medium and large
    final.large.colony.table= read.csv("FinalLargeColonies.csv", header=T)[,-1]
    final.medium.colony.table = read.csv("FinalMediumColonies.csv", header=T)[,-1]
   
    final.large.colonies=c(final.large.colonies, check.large.colony(final.large.colony.table))
    final.medium.colonies=c(final.medium.colonies, check.large.colony(final.medium.colony.table))
  }
}

# Save results
setwd(destination.directory)
write.csv(simulation.min.values, "minValues.csv")
write.csv(simulation.max.values, "maxValues.csv")
write.csv(simulation.min.values.large.colonies, "minValuesLargeColonies.csv")
write.csv(simulation.max.values.large.colonies, "maxValuesLargeColonies.csv")
write.csv(colonization.numbers, "colonizations.csv")
write.csv(decolonization.numbers, "decolonizations.csv")
write.csv(proportion.large.colonies, "proportion_large_colonies.csv")
write.csv(final.large.colonies, "final_large_colonies.csv")
write.csv(final.medium.colonies, "final_medium_colonies.csv")