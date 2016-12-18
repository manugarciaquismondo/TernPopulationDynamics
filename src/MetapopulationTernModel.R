# Get parameters from console

output.directory.route = commandArgs(T)[1]
data.directory.route = commandArgs(T)[2]
simulation.years = as.numeric(commandArgs(T)[3])

#output.directory.route = "C:/Users/manu_/localdata/workspaces/eclipse/newworkspace/TernModel/results"
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


# Set the data route according to the operating system
if (Sys.info()['sysname'] == "Windows") {
    move.to.data.route = function() {
        new.workspace = Sys.getenv("NEWWORKSPACE")
        setwd(file.path(new.workspace, "TernModel", "data"))
    }
} else {
    move.to.data.route = function() {
        setwd(data.directory.route)
    }
}
# Set the first column of each table as row name
get.named.table = function(table.file.route) {
    move.to.data.route()
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
site.distances = get.named.table("distances.csv")
colnames(site.distances) = rownames(site.distances)
colnames(productivity.estimates) = gsub("X", "", colnames(productivity.estimates))

# Generate random parameters for productivity

productivity.random.parameters = matrix(nrow = nrow(productivity.estimates), ncol = 2)
colnames(productivity.random.parameters) = c("Mean", "SD")
rownames(productivity.random.parameters) = rownames(productivity.estimates)
for (site.name in rownames(productivity.estimates)) {
    productivity.random.parameters[site.name,] = c(mean(as.matrix(productivity.estimates[site.name,])), sd(as.matrix(productivity.estimates[site.name,])))

}

# Read transition probabilities

move.to.data.route()
message("Sited in directory ", getwd())
transition.probabilities.table = read.csv("transitionprobabilities.csv", header = T, stringsAsFactors = F)
if (!exists("trimws")) {
    trimws = function(x, which = c("both", "left", "right")) {
        which <- match.arg(which)
        mysub <- function(re, x) sub(re, "", x, perl = TRUE)
        if (which == "left")
            return(mysub("^[ \t\r\n]+", x))
        if (which == "right")
            return(mysub("[ \t\r\n]+$", x))
        mysub("[ \t\r\n]+$", mysub("^[ \t\r\n]+", x))
    }
}
transition.probabilities.table[, c("Origin")] = sapply(transition.probabilities.table[, c("Origin")], trimws)
transition.probabilities.table[, c("Destination")] = sapply(transition.probabilities.table[, c("Destination")], trimws)
colnames(transition.probabilities.table) = gsub("X", "", colnames(transition.probabilities.table))
message("Transition probability table")
message(transition.probabilities.table)
# Extend transition probability table

transition.probabilities.extended = extend.year.matrix(transition.probabilities.table)


# Get phi_i,i for all sites
message("Diagonals processed")
message("Checking origins")
message(transition.probabilities.table[, "Origin"])
message("Checking destinations")
message(transition.probabilities.table[, "Destination"])
diagonal.sites = matrix(nrow = 0, ncol = ncol(transition.probabilities.table))
colnames(diagonal.sites) = colnames(transition.probabilities.table)
for (transition.index in 1:nrow(transition.probabilities.table)) {
    if (transition.probabilities.table[transition.index, "Origin"] == transition.probabilities.table[transition.index, "Destination"]) {
        message("Origin and destination coincide")
        diagonal.sites = rbind(diagonal.sites, matrix(nrow = 1, ncol = ncol(transition.probabilities.table), data = transition.probabilities.table[transition.index,]))
    }
}

# Extend phi_i,i values for all known years


# Obtain table with all phi_i,i values

message("Diagonal sites")
message(diagonal.sites)
year.range.extended = extend.year.matrix(diagonal.sites)
rownames(year.range.extended) = year.range.extended[, 1]
year.range.extended = year.range.extended[, 3:ncol(year.range.extended)]
message("Year range extended")
message(year.range.extended)
# Obtain all attractiveness (A_i,i) values
message("Calculating colonizations and decolonizations")
all.pairs.table = read.csv("allpairs.csv", header = T)
rownames(all.pairs.table) = all.pairs.table[, 1]
all.pairs.table = all.pairs.table[, 2:ncol(all.pairs.table)]
colnames(all.pairs.table) = sapply(1988:2015, toString)
known.attractiveness.values = all.pairs.table / colSums(all.pairs.table)

# Calculate site-specific carrying capacity
large.sites = read.csv("survivalestimates.csv", header = F, stringsAsFactors = F)[, 1]
carrying.capacity = apply(all.pairs.table, 1, max)
# Calculate colonization probability
colonization.opportunities = all.pairs.table
colonization.opportunities[all.pairs.table > 0] = -1
colonization.opportunities[all.pairs.table == 0] = 1
colonization.opportunities[colonization.opportunities == -1] = 0
decolonization.opportunities = 1 - colonization.opportunities
decolonization.opportunities[large.sites,] = 0

# A function to calculate colonization counts
calculate.colonizations = function(input.matrix, inverse.initial.position) {
    function(x) {
        locale.vector = input.matrix[x,]
        locale.colonizations = as.numeric(ifelse(inverse.initial.position, 1 - locale.vector[1], locale.vector[1]))
        #message("Locale col: ", locale.colonizations, ", Locale vector: ", locale.vector[1])
        for (y in 2:length(locale.vector)) {
            if (locale.vector[y - 1] == 1 && locale.vector[y] == 0) {
                locale.colonizations = locale.colonizations + 1
            }
        }
        if (rownames(input.matrix)[x] == "Breezy Point, Queens City, NY") {
            #message("Vector: ", locale.vector, ", Colonizations:", locale.colonizations)
        }
        locale.colonizations
    }
}
count.colonizations = sapply(1:nrow(colonization.opportunities), calculate.colonizations(colonization.opportunities, T))
count.decolonizations = sapply(1:nrow(decolonization.opportunities), calculate.colonizations(decolonization.opportunities, F))
names(count.colonizations) = names(count.decolonizations) = rownames(colonization.opportunities)
colonization.probability = count.colonizations / (rowSums(colonization.opportunities) + 1)
#colonization.probability = rep(sum(count.colonizations) / sum((rowSums(colonization.opportunities) + 1)), times = length(count.colonizations))
decolonization.probability = count.decolonizations / (rowSums(decolonization.opportunities) + 1)
names(colonization.probability) = names(decolonization.probability) = names(count.colonizations)
decolonization.probability = decolonization.probability * 2
#decolonization.probability = rep(sum(count.decolonizations) / sum((rowSums(decolonization.opportunities) + 1)), times = length(count.colonizations))*1.7
colonization.probability[is.infinite(colonization.probability)] = 1
colonization.probability[large.sites] = 1
decolonization.probability[colonization.probability == 1] = 0
# Create colonization matrix
colonization.probability.matrix = cbind(colonization.probability, decolonization.probability)
colnames(colonization.probability.matrix) = c("Colonization", "Decolonization")
# Multiply pool subtable by local number of pairs and normalize
message("Pool of suitable sites")
pool.subtable = pool.subtable.multiplied = subset(transition.probabilities.extended, Origin != Destination)
for (pool.subtable.year in colnames(pool.subtable)[3:ncol(pool.subtable)]) {
    for (site.name in unique(pool.subtable$Origin)) {
        origin.name = pool.subtable$Origin == site.name && pool.subtable[, pool.subtable.year] >= 0
        pool.subtable.multiplied[origin.name, pool.subtable.year] = pool.subtable[origin.name, pool.subtable.year] * all.pairs.table[site.name, pool.subtable.year]
    }
    pool.indexes.positive = pool.subtable.multiplied[, pool.subtable.year] >= 0
    pool.subtable.in.year = sum(pool.subtable.multiplied[pool.indexes.positive, pool.subtable.year])
    pool.subtable.multiplied[pool.indexes.positive, pool.subtable.year] = pool.subtable.multiplied[pool.indexes.positive, pool.subtable.year] / pool.subtable.in.year
    pool.subtable.multiplied[!pool.indexes.positive, pool.subtable.year] = 0
}
message("Summing of pool of suitable sites")
# Sum across all instances by destination variable
pool.subtable.multiplied.sum = pool.subtable.multiplied[, 2:ncol(pool.subtable.multiplied)]
pool.subtable.multiplied.sum = aggregate(pool.subtable.multiplied.sum[2:ncol(pool.subtable.multiplied.sum)], by = list(Destination = pool.subtable.multiplied.sum$Destination), FUN = function(x) {
    sum(x[x >= 0])
})
rownames(pool.subtable.multiplied.sum) = pool.subtable.multiplied.sum$Destination
pool.subtable.multiplied.sum = pool.subtable.multiplied.sum[, 2:ncol(pool.subtable.multiplied.sum)]
#Obtain pool indexes greater than 0

pool.valid.indexes = unlist(pool.subtable.multiplied.sum) > 0
# Adjust model of pool proportion as a function of attractiveness and quality

# Obtain 3-year lagged A_i,i values

calculate.lagged.attractiveness = function(attractive.matrix, matrix.index) {
    (attractive.matrix[, matrix.index] + attractive.matrix[, matrix.index - 1] + attractive.matrix[, matrix.index - 2] * .5 + attractive.matrix[, matrix.index - 3] * .25) / 2.5
}

known.lagged.attractiveness.values = matrix(nrow = nrow(known.attractiveness.values), ncol = ncol(known.attractiveness.values) - 3)
rownames(known.lagged.attractiveness.values) = rownames(known.attractiveness.values)
colnames(known.lagged.attractiveness.values) = colnames(known.attractiveness.values)[4:ncol(known.attractiveness.values)]
for (lagged.year.index in 1:ncol(known.lagged.attractiveness.values)) {
    known.lagged.attractiveness.values[, lagged.year.index] = calculate.lagged.attractiveness(known.attractiveness.values, lagged.year.index + 3)
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
message("Building matrix for model regression")
create.adjustment.matrix = function(dependent.variable.indexes, dependent.variable.vector, dependent.variable.name) {
    message("Dependent variable indexes")
    message(dependent.variable.indexes)
    # Create vectors for dependent and independent variables
    selected.dependent.vector = dependent.variable.vector[dependent.variable.indexes]
    intrinsic.quality.filtered.calibration = unlist(intrinsic.quality.calibration.matrix)[dependent.variable.indexes]
    temporal.attractiveness.vector = unlist(filtered.known.attractiveness.values)[dependent.variable.indexes]

    # Create matrix with these vectors
    adjustment.matrix = matrix(nrow = length(intrinsic.quality.filtered.calibration), ncol = 3)
    message("Internal adjustment matrix")
    colnames(adjustment.matrix) = c("Quality", dependent.variable.name, "TemporalAttractiveness")
    adjustment.matrix[, "Quality"] = intrinsic.quality.filtered.calibration
    adjustment.matrix[, dependent.variable.name] = selected.dependent.vector
    adjustment.matrix[, "TemporalAttractiveness"] = temporal.attractiveness.vector
    adjustment.matrix
}

create.adjustment.model = function(dependent.variable.indexes, dependent.variable.vector, dependent.variable.name) {
    # Create adjustment matrix
    adjustment.matrix = create.adjustment.matrix(dependent.variable.indexes, dependent.variable.vector, dependent.variable.name)
    message("Adjustment matrix")
    message(adjustment.matrix)
    # Create logistic model from dependent and independent matrix
    glm(eval(parse(text = paste(dependent.variable.name, "~ Quality + TemporalAttractiveness"))), family = gaussian(link = 'logit'), data = as.data.frame(adjustment.matrix))

}

estimate.regression = function(input.data, input.model) {
    predict(input.model, newdata = as.data.frame(input.data), type = 'response')
}

message("Quality model started")
message("Dated Phi indexes")
message(dated.phi.indexes)
# Next, build a logistic regression model
logistic.quality.model = create.adjustment.model(dated.phi.indexes, selected.phi.rows, "Phi")
message("Quality model finished")
pool.proportion.model = create.adjustment.model(pool.valid.indexes, unlist(pool.subtable.multiplied.sum), "PoolProportion")
message("Pool models finished")
message("Model regression matrices built")

# Calculate the random part of the proportion from the pool as the squared error in the model estimations
pool.adjustment.matrix = create.adjustment.matrix(pool.valid.indexes, unlist(pool.subtable.multiplied.sum), "PoolProportion")
pool.proportion.known.values = unlist(pool.subtable.multiplied.sum)[pool.valid.indexes]
pool.proportion.expectancy = mean(pool.proportion.known.values)
pool.proportion.sd = sd(pool.proportion.known.values)
pool.randomness = sum((estimate.regression(pool.adjustment.matrix, pool.proportion.model) - unlist(pool.subtable.multiplied.sum)[pool.valid.indexes]) ^ 2)
pool.randomness = pool.randomness / length(pool.valid.indexes)

# A function to calculate the proportion from the pool to each site
calculate.pool.proportion = function(model.parameters, number.of.pairs, colonized.sites, input.number.of.sites = number.of.sites, input.pool.proportion.model = pool.proportion.model, input.pool.randomness = pool.randomness, input.pool.proportion.expectancy = pool.proportion.expectancy, input.pool.proportion.sd = pool.proportion.sd) {
    deterministic.component = estimate.regression(model.parameters, input.pool.proportion.model)
    stochastic.component = rnorm(n = input.number.of.sites, mean = input.pool.proportion.expectancy, sd = input.pool.proportion.sd)

    # Penalize the existence of populated colonies according to their closeness to each site
    pool.proportion.vector = ((1 - input.pool.randomness) * deterministic.component + input.pool.randomness * stochastic.component) * colSums(number.of.pairs * site.distances)
    pool.proportion.vector = pool.proportion.vector * colonized.sites
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
move.to.data.route()
productivity.matrix = read.csv("productivities6.csv", header = F)

# Separate in productivity and number of pairs
productivity.number.of.pairs = productivity.matrix[, seq(1, ncol(productivity.matrix), 2)]
productivity.known.values = productivity.matrix[, seq(2, ncol(productivity.matrix), 2)]

# Get valid productivity values
productivity.as.vector = unlist(productivity.known.values)
productivity.matrix.valid.indexes = productivity.as.vector >= 0
productivity.as.vector = productivity.as.vector[productivity.matrix.valid.indexes]
number.of.pairs.as.vector = unlist(productivity.number.of.pairs)[productivity.matrix.valid.indexes]
intrinsic.quality.as.vector = unlist(matrix(data = intrinsic.quality, nrow = nrow(productivity.number.of.pairs), ncol = ncol(productivity.number.of.pairs)))[productivity.matrix.valid.indexes]

# Calculate productivity correlation per site
productivity.correlation.per.site = c()
for (site.index in 1:nrow(productivity.known.values)) {
    valid.indexes = productivity.known.values[site.index,] > 0
    if (length(valid.indexes[valid.indexes]) > 1) {
        latest.valid.indexes = valid.indexes
        locale.correlation = cor(unlist(productivity.known.values[site.index, latest.valid.indexes]), unlist(productivity.number.of.pairs[site.index, latest.valid.indexes]))
        productivity.correlation.per.site = c(productivity.correlation.per.site, locale.correlation)
    }
}
# Create regression model to predict productivity values
known.productivity.values.matrix = cbind(productivity.as.vector, number.of.pairs.as.vector, intrinsic.quality.as.vector)
colnames(known.productivity.values.matrix) = c("Productivity", "Pairs", "Quality")
max.productivity = 2
known.productivity.values.matrix[, "Productivity"] = known.productivity.values.matrix[, "Productivity"] / max.productivity
known.productivity.values.matrix = known.productivity.values.matrix[known.productivity.values.matrix[, "Productivity"] > 0,]
productivity.logistic.model = glm(Productivity ~ Pairs + Quality, family = gaussian(link = 'logit'), as.data.frame(known.productivity.values.matrix))

# Read adult and juvenile fidelity

fidelity.for.juveniles = read.csv("fidelity_ratio.csv", header = T)
juvenile.pairs.for.fidelity = get.named.table("juvenile_fidelity_pairs.csv")
juvenile.pairs.unlisted = unlist(juvenile.pairs.for.fidelity)
fidelity.ratio = fidelity.for.juveniles["Juvenile"] / fidelity.for.juveniles["Adult"]
fidelity.ratio.extended = rep(unlist(fidelity.ratio), each = nrow(juvenile.pairs.for.fidelity))
fidelity.data.frame = as.data.frame(cbind(juvenile.pairs.unlisted, fidelity.ratio.extended))
colnames(fidelity.data.frame) = c("Pairs", "Fidelity")
juvenile.fidelity.logistic.model = glm(Fidelity ~ Pairs, family = gaussian(link = 'logit'), fidelity.data.frame)

# Read juvenile survival distribution

juvenile.survival.data = read.csv("juvenile_survival.csv", header = F)
mean.juvenile.survival = mean(juvenile.survival.data[, 1])
sd.juvenile.survival = sd(juvenile.survival.data[, 1])

# Read juvenile survival by site

juvenile.ordered.survival.data = get.named.table("juvenile_ordered_survival.csv")
colnames(juvenile.ordered.survival.data) = sapply(colnames(juvenile.ordered.survival.data), function(x) gsub("X", "", x))
juvenile.ordered.survival.data = t(juvenile.ordered.survival.data)
colnames(juvenile.ordered.survival.data) = c("Falkner I., CT", "Bird I., Marion, MA", "Great Gull I., NY")
# Create logistic regression model of juvenile survival

filtered.survival.juvenile.pairs.for.fidelity = juvenile.pairs.for.fidelity[rownames(juvenile.ordered.survival.data),]
filtered.survival.juvenile.pairs.for.fidelity = filtered.survival.juvenile.pairs.for.fidelity[, c(1, 3, 2)]
colnames(filtered.survival.juvenile.pairs.for.fidelity) = colnames(juvenile.ordered.survival.data)
unlisted.survival.juvenile.pairs = unlist(filtered.survival.juvenile.pairs.for.fidelity)
intrinsic.quality.calibration.matrix.survival = intrinsic.quality.calibration.matrix[colnames(juvenile.ordered.survival.data), 1:nrow(juvenile.ordered.survival.data)]
unlisted.juvenile.survival = unlisted.quality.for.survival = c()
for (population.name in colnames(juvenile.ordered.survival.data)) {
    unlisted.quality.for.survival = c(unlisted.quality.for.survival, intrinsic.quality.calibration.matrix.survival[population.name,])
    unlisted.juvenile.survival = c(unlisted.juvenile.survival, juvenile.ordered.survival.data[, population.name])
}
juvenile.survival.matrix = cbind(unlisted.survival.juvenile.pairs, unlisted.juvenile.survival, unlisted.quality.for.survival)
colnames(juvenile.survival.matrix) = c("Pairs", "Survival", "Quality")
survival.logistic.model = glm(Survival ~ Pairs + Quality, family = gaussian(link = 'logit'), as.data.frame(juvenile.survival.matrix))

# Select matching adult survival and correlate it with juvenile survival

juvenile.corresponding.survival = t(survival.table[colnames(juvenile.ordered.survival.data), sapply(rownames(juvenile.ordered.survival.data), function(x) paste("X", x, sep = ""))])
unlisted.juvenile.corresponding.survival = c()
for (population.name in colnames(juvenile.corresponding.survival))
    unlisted.juvenile.corresponding.survival = c(unlisted.juvenile.corresponding.survival, juvenile.corresponding.survival[, population.name])
juvenile.and.adult.correlation = cor(unlisted.juvenile.corresponding.survival, unlisted.juvenile.survival)

# Initialize data structures

number.of.sites = length(known.sites)
yearly.fidelity = rep(0, times = length(known.sites))
TotalBirdsNextYear = simulation.attractiveness.values = matrix(nrow = number.of.sites, ncol = simulation.years)
NewDecoloners = NewColoners = DeNovoDecolonization = DeNovoColonization = JuvenileSurvival = TemporaryDecolonization = ColonizationVector = DecolonizationVector = BirdsStayingNextYear = BirdsLeavingNextYear = ImmigrantsNextYear = JuvenileStaying = JuvenileLeaving = EstimatedProductivity = matrix(nrow = number.of.sites, ncol = simulation.years - 1)
rownames(NewDecoloners) = rownames(NewColoners) = rownames(DeNovoDecolonization) = rownames(DeNovoColonization) = names(yearly.fidelity) = known.sites
Pool = c()

# Collect increases and decreases in colonization/decolonization

colonization.numbers = decolonization.numbers = c()
for (site.index in rownames(all.pairs.table)[!rownames(all.pairs.table) %in% large.sites]) {
    max.decolonization = 0
    for (year.index in 2:ncol(all.pairs.table)) {
        max.decolonization = max(max.decolonization, all.pairs.table[site.index, year.index - 1] - all.pairs.table[site.index, year.index])
        if (all.pairs.table[site.index, year.index - 1] == 0 && all.pairs.table[site.index, year.index] > 0)
            colonization.numbers = c(colonization.numbers, all.pairs.table[site.index, year.index])
        }
    if (max.decolonization > 0)
        decolonization.numbers = c(decolonization.numbers, max.decolonization)
    }

# Calculate parameters for colonization and decolonization distributions
max.possible.colonization = max(colonization.numbers)
colonization.probability = median(colonization.numbers) / max.possible.colonization

max.possible.decolonization = max(decolonization.numbers)
decolonization.probability = median(decolonization.numbers) / max.possible.decolonization
# Simulate the system

averaged.overall.survival = mean(as.matrix(survival.table))
TotalBirdsNextYear[, 1] = all.pairs.table[, 1]
simulation.attractiveness.values[, 1] = known.attractiveness.values[, 1]
TemporaryDecolonization[, 1] = F

# Test what happens when all sites are known sites
base.large.sites = large.sites
large.sites = known.sites

message("Simulation started")

for (t in 1:(simulation.years - 1)) {
    # Calculate the fidelity in the current year
    fidelity.predictor = cbind(intrinsic.quality, simulation.attractiveness.values[, t])
    colnames(fidelity.predictor) = c("Quality", "TemporalAttractiveness")
    yearly.fidelity = estimate.regression(fidelity.predictor, logistic.quality.model)
    temporal.survivors = TotalBirdsNextYear[, t] * averaged.overall.survival

    # Calculate colonized and decolonized sites
    colonized.sites = temporal.survivors >= 1
    colonization.probability.sampled.matrix = cbind(runif(number.of.sites), runif(number.of.sites))
    colnames(colonization.probability.sampled.matrix) = c("ColonizationSample", "DecolonizationSample")
    current.year.colonization = colonization.probability.sampled.matrix[, "ColonizationSample"] < colonization.probability.matrix[, "Colonization"]
    DeNovoColonization[, t] = ifelse(current.year.colonization & !colonized.sites, 1, 0)
    ColonizationVector[, t] = current.year.colonization | colonized.sites
    current.year.decolonization = colonization.probability.sampled.matrix[, "DecolonizationSample"] < colonization.probability.matrix[, "Decolonization"]
    DeNovoDecolonization[, t] = ifelse(current.year.decolonization & colonized.sites, 1, 0)
    DecolonizationVector[, t] = current.year.decolonization & colonized.sites
    colonization.status = ColonizationVector[, t] & !DecolonizationVector[, t]
    colonization.status = ifelse(colonization.status, 1, 0)
    # If the colony will be decolonized, set fidelity to 0
    #yearly.fidelity = ifelse(DecolonizationVector[, t] == 0, yearly.fidelity,0)
    # If decolonization occurs, flag temporary decolonization
    TemporaryDecolonization[, t] = DecolonizationVector[, t] == 1

    # Decolonization for large sites
    #carrying.capacity.threshold = pmin(1, exp( - (carrying.capacity[large.sites] - temporal.survivors[large.sites])))
    #TemporaryDecolonization[large.sites, t] = runif(length(large.sites)) < carrying.capacity.threshold
    if (t > 1) {

        # Decolonization finishes when the colony has no pairs

        TemporaryDecolonization[, t] = (TemporaryDecolonization[, t] | TemporaryDecolonization[, t - 1]) & TotalBirdsNextYear[, t - 1] > 0

        # Well-known large sites can rebound before total decolonization

        #for (large.site in base.large.sites) {
        #if (TemporaryDecolonization[large.site, t - 1]) {
        #site.carrying.capacity.threshold = 1-min(1, exp( - (carrying.capacity[large.site]*.6 - temporal.survivors[large.site])))
        #TemporaryDecolonization[, t] = !(runif(1) < site.carrying.capacity.threshold)
        #}
        #}

    }

    ## If all large sites are decolonized, stop decolonization to avoid the collapse of the metapopulation
    #if (all(TemporaryDecolonization[large.sites, t])) {
    #TemporaryDecolonization[large.sites, t][temporal.survivors[large.sites]<carrying.capacity[large.sites]] = F
    #}
    colonization.vector = colonization.status * ifelse(!TemporaryDecolonization[, t], 1, 0)
    # Estimate productivity 
    total.birds.as.matrix = cbind(as.data.frame(temporal.survivors), intrinsic.quality)
    colnames(total.birds.as.matrix) = c("Pairs", "Quality")
    EstimatedProductivity[, t] = pmax(0, estimate.regression(total.birds.as.matrix, productivity.logistic.model) * max.productivity)
    # If the colony will be decolonized, do not generate new juveniles
    EstimatedProductivity[, t] = ifelse(!TemporaryDecolonization[, t], EstimatedProductivity[, t], 0)
    # Estimate juvenile survival
    total.birds.plus.quality = cbind(total.birds.as.matrix, intrinsic.quality)
    colnames(total.birds.plus.quality) = c("Pairs", "Quality")
    JuvenileSurvival[, t] = estimate.regression(total.birds.plus.quality, survival.logistic.model)
    #yearly.juvenile.survival = pmax(0,pmin(1,rnorm(number.of.sites, mean = mean.juvenile.survival, sd=sd.juvenile.survival)))

    # Calculate juvenile product

    juvenile.generated = EstimatedProductivity[, t] * JuvenileSurvival[, t] * temporal.survivors
    yearly.fidelity.juvenile.ratio = estimate.regression(total.birds.as.matrix, juvenile.fidelity.logistic.model)
    juvenile.fidelity = yearly.fidelity * yearly.fidelity.juvenile.ratio
    JuvenileStaying[, t] = juvenile.generated * juvenile.fidelity
    JuvenileLeaving[, t] = juvenile.generated * (1 - juvenile.fidelity)
    # Calculate the number of birds that will remain and stay in the colony the next year
    BirdsStayingNextYear[, t] = temporal.survivors * yearly.fidelity + JuvenileStaying[, t]

    BirdsLeavingNextYear[, t] = temporal.survivors * (1 - yearly.fidelity) + JuvenileLeaving[, t]
    NewDecoloners[, t] = pmin(BirdsStayingNextYear[, t], rbinom(number.of.sites, max.possible.decolonization, decolonization.probability) * DeNovoDecolonization[, t])
    BirdsLeavingNextYear[, t] = BirdsLeavingNextYear[, t] + NewDecoloners[, t]
    BirdsStayingNextYear[, t] = BirdsStayingNextYear[, t] - NewDecoloners[, t]
    # Calculate the migrant pool and immigrants
    year.pool = sum(BirdsLeavingNextYear[, t])
    NewColoners[, t] = rbinom(number.of.sites, max.possible.colonization, colonization.probability) * DeNovoColonization[, t]
    year.pool = year.pool - sum(NewColoners[, t])
    Pool = c(Pool, year.pool)

    # If any colony receives immigrants, distribute them among the receiving colonies
    if (any(colonization.vector == 1)) {
        ImmigrantsNextYear[, t] = year.pool * calculate.pool.proportion(fidelity.predictor, temporal.survivors, colonization.vector * (1 - DeNovoColonization[, t]))
    } else {
        ImmigrantsNextYear[, t] = 0
    }
    ImmigrantsNextYear[, t] = ImmigrantsNextYear[, t] + NewColoners[, t]
    # Calculate total number of birds
    TotalBirdsNextYear[, t + 1] = BirdsStayingNextYear[, t] + ImmigrantsNextYear[, t]
    TotalBirdsNextYear[, t + 1] = randomRound(TotalBirdsNextYear[, t + 1])



    # Update attractiveness
    sumBirdsStaying = sum(BirdsStayingNextYear[, t])
    yearly.attractiveness = BirdsStayingNextYear[, t] / ifelse(sumBirdsStaying <= 0, 1, sumBirdsStaying)
    simulation.attractiveness.values[, t + 1] = yearly.attractiveness

}

# Create directory to store the results

dir.create(output.directory.route, showWarnings = FALSE)
setwd(output.directory.route)

# Store simulation results

write.csv(TotalBirdsNextYear, "TotalBirdsNextYear.csv")
write.csv(simulation.attractiveness.values, "SiteAttractiveness.csv")
write.csv(NewDecoloners, "NewDecoloners.csv")
write.csv(DeNovoColonization, "DeNovoColonization.csv")
write.csv(DeNovoDecolonization, "DeNovoDecolonization.csv")
write.csv(JuvenileSurvival, "JuvenileSurvival.csv")
write.csv(TemporaryDecolonization, "TemporaryDecolonization.csv")
write.csv(ColonizationVector, "ColonizationVector.csv")
write.csv(DecolonizationVector, "DecolonizationVector.csv")
write.csv(BirdsStayingNextYear, "BirdsStayingNextYear.csv")
write.csv(BirdsLeavingNextYear, "BirdsLeavingNextYear.csv")
write.csv(ImmigrantsNextYear, "Immigrants.csv")
write.csv(JuvenileStaying, "JuvenileStaying.csv")
write.csv(JuvenileLeaving, "JuvenileLeaving.csv")
write.csv(EstimatedProductivity, "EstimatedProductivity.csv")
