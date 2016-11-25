
# Define functions to extend year columns from year ranges

year.range = function(year.range.string) {
    min.year = as.numeric(sapply(strsplit(year.range.string, "\\."), function(x) {
        gsub("X", "", x)
    }))
    if (length(min.year) <= 1) {
        min.year[1]
    } else {
        min.year[1]:min.year[2]
    }

}

extend.year.range.with.value = function(year.range.string, value.vector) {
    year.range.values = year.range(year.range.string)
    year.range.matrix = matrix(nrow = length(value.vector), ncol = 0)
    for (year.in.range in year.range.values) {
        year.range.matrix = cbind(year.range.matrix, value.vector)
    }
    year.range.matrix
}

extend.year.matrix = function(input.matrix) {
    year.range.extended = matrix(nrow = nrow(input.matrix), ncol = 0)
    for (year.index in 3:ncol(input.matrix)) {

        year.range.extended = cbind(year.range.extended, extend.year.range.with.value(colnames(input.matrix)[year.index], input.matrix[, year.index]))
    }
    year.range.extended = cbind(input.matrix[, c("Origin", "Destination")], year.range.extended)
    colnames(year.range.extended)[3:ncol(year.range.extended)] = sapply(1989:2008, toString)

    year.range.extended
}

# Set the first column of each table as row name

get.named.table = function(table.file.route) {
    new.workspace = Sys.getenv("NEWWORKSPACE")
    setwd(file.path(new.workspace, "TernModel", "data"))
    input.table = read.csv(table.file.route, header = T)
    site.names = input.table[, 1]
    input.table = input.table[, 2:ncol(input.table)]
    rownames(input.table) = site.names
    return(input.table)
}


# Read survival probabilities per year and site

survival.table = get.named.table("survivalallestimates.csv")

# Average survival per site across all years

averaged.survival = rowSums(survival.table) / ncol(survival.table)

# Read productivity estimates

productivity.estimates = get.named.table("LogPredictedProductivityValues.csv")
colnames(productivity.estimates) = gsub("X", "", colnames(productivity.estimates))

# Generate random parameters for productivity

productivity.random.parameters = matrix(nrow = nrow(productivity.estimates), ncol = 2)
colnames(productivity.random.parameters) = c("Mean", "SD")
rownames(productivity.random.parameters) = rownames(productivity.estimates)
for (site.name in rownames(productivity.estimates)) {
    productivity.random.parameters[site.name,] = c(mean(as.matrix(productivity.estimates[site.name,])), sd(as.matrix(productivity.estimates[site.name,])))

}

# Read transition probabilities

new.workspace = Sys.getenv("NEWWORKSPACE")
setwd(file.path(new.workspace, "TernModel", "data"))
transition.probabilities.table = read.csv("transitionprobabilities.csv", header = T, stringsAsFactors = F)
transition.probabilities.table[, c("Origin")] = sapply(transition.probabilities.table[, c("Origin")], trimws)
transition.probabilities.table[, c("Destination")] = sapply(transition.probabilities.table[, c("Destination")], trimws)
colnames(transition.probabilities.table) = gsub("X", "", colnames(transition.probabilities.table))

# Extend transition probability table

transition.probabilities.extended = extend.year.matrix(transition.probabilities.table)


# Get phi_i,i for all sites

diagonal.sites = matrix(nrow = 0, ncol = ncol(transition.probabilities.table))
colnames(diagonal.sites) = colnames(transition.probabilities.table)
for (transition.index in 1:nrow(transition.probabilities.table)) {
    if (transition.probabilities.table[transition.index, "Origin"] == transition.probabilities.table[transition.index, "Destination"]) {
        diagonal.sites = rbind(diagonal.sites, matrix(nrow = 1, ncol = ncol(transition.probabilities.table), data = transition.probabilities.table[transition.index,]))
    }
}

# Extend phi_i,i values for all known years


# Obtain table with all phi_i,i values


year.range.extended = extend.year.matrix(diagonal.sites)
rownames(year.range.extended) = year.range.extended[, 1]
year.range.extended = year.range.extended[, 3:ncol(year.range.extended)]

# Obtain all attractiveness (A_i,i) values

all.pairs.table = read.csv("allpairs.csv", header = T)
rownames(all.pairs.table) = all.pairs.table[, 1]
all.pairs.table = all.pairs.table[, 2:ncol(all.pairs.table)]
colnames(all.pairs.table) = sapply(1988:2015, toString)
known.attractiveness.values = all.pairs.table / colSums(all.pairs.table)

# Multiply pool subtable by local number of pairs and normalize

pool.subtable = pool.subtable.multiplied = subset(transition.probabilities.extended, Origin != Destination)
for (pool.subtable.year in colnames(pool.subtable)[3:ncol(pool.subtable)]) {
    for (site.name in unique(pool.subtable$Origin)) {
        origin.name = pool.subtable$Origin == site.name && pool.subtable[, pool.subtable.year]>=0
        pool.subtable.multiplied[origin.name, pool.subtable.year] = pool.subtable[origin.name, pool.subtable.year] * all.pairs.table[site.name, pool.subtable.year]
    }
    pool.indexes.positive = pool.subtable.multiplied[, pool.subtable.year] >= 0
    pool.subtable.in.year = sum(pool.subtable.multiplied[pool.indexes.positive, pool.subtable.year])
    pool.subtable.multiplied[pool.indexes.positive, pool.subtable.year] = pool.subtable.multiplied[pool.indexes.positive, pool.subtable.year] / pool.subtable.in.year
    pool.subtable.multiplied[!pool.indexes.positive, pool.subtable.year] = 0
}

# Sum across all instances by destination variable
pool.subtable.multiplied.sum = pool.subtable.multiplied[, 2:ncol(pool.subtable.multiplied)]
pool.subtable.multiplied.sum = aggregate(pool.subtable.multiplied.sum[2:ncol(pool.subtable.multiplied.sum)], by = list(Destination = pool.subtable.multiplied.sum$Destination), FUN = function(x) {
    sum(x[x >= 0])
})
rownames(pool.subtable.multiplied.sum) = pool.subtable.multiplied.sum$Destination
pool.subtable.multiplied.sum = pool.subtable.multiplied.sum[, 2:ncol(pool.subtable.multiplied.sum)]
#Obtain pool indexes greater than 0

pool.valid.indexes = unlist(pool.subtable.multiplied.sum)>0
# Adjust model of pool proportion as a function of attractiveness and quality

# Obtain 3-year lagged A_i,i values

calculate.lagged.attractiveness = function(attractive.matrix, matrix.index) {
    (attractive.matrix[, matrix.index] + attractive.matrix[, matrix.index - 1] + attractive.matrix[, matrix.index - 2] * .5 + attractive.matrix[, matrix.index - 3] * .25) / 2.5
}

known.lagged.attractiveness.values = matrix(nrow = nrow(known.attractiveness.values), ncol = ncol(known.attractiveness.values) - 3)
rownames(known.lagged.attractiveness.values) = rownames(known.attractiveness.values)
colnames(known.lagged.attractiveness.values) = colnames(known.attractiveness.values)[4:ncol(known.attractiveness.values)]
for (lagged.year.index in 1:ncol(known.lagged.attractiveness.values)) {
    known.lagged.attractiveness.values[, lagged.year.index] = calculate.lagged.attractiveness(known.attractiveness.values, lagged.year.index+3)
}

# Obtain the intrinsic quality of each site

sum.known.attractiveness.values = rowSums(known.attractiveness.values)
intrinsic.quality = sum.known.attractiveness.values / sum(sum.known.attractiveness.values)

# Estimate the constant component of fidelity as a logistic function of quality

# First, build a matrix with quality_i, phi_i,i,t and TA_i,t

# First, build TA_i,t matrix
filtered.known.attractiveness.values = known.attractiveness.values[rownames(year.range.extended), colnames(year.range.extended)]

# Then, build quality_i matrix, repeating the quality_i vector per year
intrinsic.quality.calibration = intrinsic.quality[rownames(year.range.extended)]
intrinsic.quality.calibration.matrix = matrix(nrow = nrow(year.range.extended), ncol = ncol(year.range.extended), intrinsic.quality.calibration, byrow = F)
colnames(intrinsic.quality.calibration.matrix) = colnames(year.range.extended)
rownames(intrinsic.quality.calibration.matrix) = rownames(year.range.extended)

# Filter valid phi_i,i,t values (those that are greater or equal to 0)
selected.phi.rows = unlist(year.range.extended)
dated.phi.indexes = which(selected.phi.rows >= 0)


#Finally, build the matrix for model regression

create.adjustment.matrix = function(dependent.variable.indexes, dependent.variable.vector, dependent.variable.name) {

    # Create vectors for dependent and independent variables
    selected.dependent.vector = dependent.variable.vector[dependent.variable.indexes]
    intrinsic.quality.filtered.calibration = unlist(intrinsic.quality.calibration.matrix)[dependent.variable.indexes]
    temporal.attractiveness.vector = unlist(filtered.known.attractiveness.values)[dependent.variable.indexes]

    # Create matrix with these vectors
    adjustment.matrix = matrix(nrow = length(intrinsic.quality.filtered.calibration), ncol = 3)
    colnames(adjustment.matrix) = c("Quality", dependent.variable.name, "TemporalAttractiveness")
    adjustment.matrix[, "Quality"] = intrinsic.quality.filtered.calibration
    adjustment.matrix[, dependent.variable.name] = selected.dependent.vector
    adjustment.matrix[, "TemporalAttractiveness"] = temporal.attractiveness.vector
    adjustment.matrix
}

create.adjustment.model = function(dependent.variable.indexes, dependent.variable.vector, dependent.variable.name) {
    # Create adjustment matrix
    adjustment.matrix = create.adjustment.matrix(dependent.variable.indexes, dependent.variable.vector, dependent.variable.name)

    # Create logistic model from dependent and independent matrix
    glm(eval(parse(text=paste(dependent.variable.name,"~ Quality + TemporalAttractiveness"))), family = gaussian(link = 'logit'), data = as.data.frame(adjustment.matrix))

}

estimate.regression = function(input.data, input.model) {
    predict(input.model, newdata = as.data.frame(input.data), type = 'response')
}


# Next, build a logistic regression model
logistic.quality.model = create.adjustment.model(dated.phi.indexes, selected.phi.rows, "Phi")

pool.proportion.model = create.adjustment.model(pool.valid.indexes, unlist(pool.subtable.multiplied.sum), "PoolProportion")

# Calculate the random part of the proportion from the pool as the squared error in the model estimations
pool.adjustment.matrix = create.adjustment.matrix(pool.valid.indexes, unlist(pool.subtable.multiplied.sum), "PoolProportion")
pool.proportion.known.values = unlist(pool.subtable.multiplied.sum)[pool.valid.indexes]
pool.proportion.expectancy = mean(pool.proportion.known.values)
pool.proportion.sd = sd(pool.proportion.known.values)
pool.randomness = sum((estimate.regression(pool.adjustment.matrix, pool.proportion.model) - unlist(pool.subtable.multiplied.sum)[pool.valid.indexes]) ^ 2)
pool.randomness = pool.randomness / length(pool.valid.indexes)

# A function to calculate the proportion from the pool to each site
calculate.pool.proportion = function(model.parameters, input.number.of.sites = number.of.sites, input.pool.proportion.model = pool.proportion.model, input.pool.randomness = pool.randomness, input.pool.proportion.expectancy = pool.proportion.expectancy, input.pool.proportion.sd = pool.proportion.sd) {
    deterministic.component = estimate.regression(model.parameters, input.pool.proportion.model)
    stochastic.component = rnorm(n = input.number.of.sites, mean = input.pool.proportion.expectancy, sd = input.pool.proportion.sd)
    pool.proportion.vector = (1 - input.pool.randomness) * deterministic.component + input.pool.randomness * stochastic.component
    pool.proportion.vector / sum(pool.proportion.vector)
}


# Build a regression model for the proportion of birds drawn from the pool

intrinsic.quality.filtered.calibration = unlist(intrinsic.quality.calibration.matrix)[pool.valid.indexes]

known.sites = rownames(all.pairs.table)

# A random rounder of real numbers
randomRound = function(realNumbers) {
    fractionalPart = realNumbers %% 1
    selectionNumbers = runif(n = length(realNumbers))
    sapply(1:length(realNumbers), function(x) {
        if (selectionNumbers[x] > fractionalPart[x]) {
            floor(realNumbers[x])
        } else {
            ceiling(realNumbers[x])
        }

    })
}

# Read productivity and number of pairs
new.workspace = Sys.getenv("NEWWORKSPACE")
setwd(file.path(new.workspace, "TernModel", "data"))
productivity.matrix = read.csv("productivities6.csv", header = F)

# Separate in productivity and number of pairs
productivity.number.of.pairs = productivity.matrix[, seq(1, ncol(productivity.matrix), 2)]
productivity.known.values = productivity.matrix[, seq(2, ncol(productivity.matrix), 2)]

# Get valid productivity values
productivity.as.vector = unlist(productivity.known.values)
productivity.matrix.valid.indexes = productivity.as.vector >= 0
productivity.as.vector = productivity.as.vector[productivity.matrix.valid.indexes]
number.of.pairs.as.vector = unlist(productivity.number.of.pairs)[productivity.matrix.valid.indexes]

# Create regression model to predict productivity values
known.productivity.values.matrix = cbind(productivity.as.vector, number.of.pairs.as.vector)
colnames(known.productivity.values.matrix) = c("Productivity", "Pairs")
max.productivity = 2
known.productivity.values.matrix[, "Productivity"] = known.productivity.values.matrix[, "Productivity"] / max.productivity
known.productivity.values.matrix = known.productivity.values.matrix[known.productivity.values.matrix[, "Productivity"]>0,]
productivity.logistic.model = glm(Productivity ~ Pairs, family = gaussian(link = 'logit'), as.data.frame(known.productivity.values.matrix))

# Read juvenile survival distribution

juvenile.survival.data = read.csv("juvenile_survival.csv", header = F)
mean.juvenile.survival = mean(juvenile.survival.data[,1])
sd.juvenile.survival = sd(juvenile.survival.data[, 1])

# Read adult and juvenile fidelity

fidelity.for.juveniles = read.csv("fidelity_ratio.csv", header = T)
fidelity.ratio = fidelity.for.juveniles["Juvenile"] / fidelity.for.juveniles["Adult"]
fidelity.mean=mean(unlist(fidelity.ratio))

# Initialize data structures
simulation.years = 15
number.of.sites = length(known.sites)
yearly.fidelity = rep(0, times = length(known.sites))
TotalBirdsNextYear = simulation.attractiveness.values = matrix(nrow = number.of.sites, ncol = simulation.years)
BirdsStayingNextYear = BirdsLeavingNextYear = ImmigrantsNextYear  = JuvenileStaying = JuvenileLeaving = EstimatedProductivity = matrix(nrow = number.of.sites, ncol = simulation.years-1)
rownames(BirdsStayingNextYear) = rownames(BirdsLeavingNextYear) = rownames(TotalBirdsNextYear) = rownames(ImmigrantsNextYear) = rownames(simulation.attractiveness.values) = rownames(JuvenileStaying) = rownames(JuvenileLeaving) = rownames(EstimatedProductivity) = known.sites
names(yearly.fidelity) = known.sites
Pool=c()
# Simulate the system
averaged.overall.survival = mean(as.matrix(survival.table))
TotalBirdsNextYear[, 1] = all.pairs.table[, 1]
simulation.attractiveness.values[,1] =known.attractiveness.values[,1]
for (t in 1:(simulation.years-1)) {
    # Calculate the fidelity in the current year
    fidelity.predictor = cbind(intrinsic.quality, simulation.attractiveness.values[, t])
    colnames(fidelity.predictor)=c("Quality", "TemporalAttractiveness")
    yearly.fidelity = estimate.regression(fidelity.predictor, logistic.quality.model)
    temporal.survivors = TotalBirdsNextYear[, t] * averaged.overall.survival

    # Estimate productivity 
    total.birds.as.matrix = as.data.frame(TotalBirdsNextYear[, t])
    colnames(total.birds.as.matrix) = c("Pairs")
    EstimatedProductivity[, t] = pmax(0, estimate.regression(total.birds.as.matrix, productivity.logistic.model)*max.productivity)

    # Estimate juvenile survival
    yearly.juvenile.survival = pmax(0,pmin(1,rnorm(number.of.sites, mean = mean.juvenile.survival, sd=sd.juvenile.survival)))

    # Calculate juvenile product

    juvenile.generated = EstimatedProductivity[, t] * yearly.juvenile.survival * temporal.survivors
    juvenile.fidelity = yearly.fidelity*fidelity.mean
    JuvenileStaying[, t] = juvenile.generated * juvenile.fidelity
    JuvenileLeaving[, t] = juvenile.generated * (1 - juvenile.fidelity)
    # Calculate the number of birds that will remain and stay in the colony the next year
    BirdsStayingNextYear[, t] = temporal.survivors * yearly.fidelity + JuvenileStaying[, t]

    BirdsLeavingNextYear[, t] = temporal.survivors * (1 - yearly.fidelity) + JuvenileLeaving[, t]


    # Calculate the migrant pool and immigrants
    year.pool = sum(BirdsLeavingNextYear[, t])
    Pool = c(Pool, year.pool)
    ImmigrantsNextYear[, t] = year.pool * calculate.pool.proportion(fidelity.predictor)

    # Calculate total number of birds
    TotalBirdsNextYear[, t + 1] = BirdsStayingNextYear[, t] + ImmigrantsNextYear[, t]
    TotalBirdsNextYear[, t + 1] = randomRound(TotalBirdsNextYear[, t + 1])

    

    # Update attractiveness
    sumBirdsStaying = sum(BirdsStayingNextYear[, t])
    yearly.attractiveness = BirdsStayingNextYear[, t] / ifelse(sumBirdsStaying <= 0, 1, sumBirdsStaying)
    simulation.attractiveness.values[,t+1] = yearly.attractiveness
    #yearly.lagged.attractiveness = calculate.lagged.attractiveness(known.attractiveness.values, t)
        #known.lagged.attractiveness.values = cbind(known.lagged.attractiveness.values, yearly.lagged.attractiveness)
    #colnames(known.lagged.attractiveness.values)[ncol(known.lagged.attractiveness.values)] = toString(as.numeric(colnames(known.lagged.attractiveness.values)[ncol(known.lagged.attractiveness.values)]) + 1)
        

}
