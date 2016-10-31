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

# Read transition probabilities


new.workspace = Sys.getenv("NEWWORKSPACE")
setwd(file.path(new.workspace, "TernModel", "data"))
transition.probabilities.table = read.csv("transitionprobabilities.csv", header = T, stringsAsFactors = F)
transition.probabilities.table[, c("Origin")] = sapply(transition.probabilities.table[, c("Origin")], trimws)
transition.probabilities.table[, c("Destination")] = sapply(transition.probabilities.table[, c("Destination")], trimws)

# Get phi_i,i for all sites

diagonal.sites = matrix(nrow = 0, ncol = ncol(transition.probabilities.table))
colnames(diagonal.sites) = colnames(transition.probabilities.table)
for (transition.index in 1:nrow(transition.probabilities.table)) {
    if (transition.probabilities.table[transition.index, "Origin"] == transition.probabilities.table[transition.index, "Destination"]) {
        diagonal.sites = rbind(diagonal.sites, matrix(nrow = 1, ncol = ncol(transition.probabilities.table), data = transition.probabilities.table[transition.index,]))
    }
}

# Extend phi_i,i values for all known years

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

# Obtain table with all phi_i,i values

extend.year.range.with.value = function(year.range.string, value.vector) {
    year.range.values = year.range(year.range.string)
    year.range.matrix = matrix(nrow = length(value.vector), ncol = 0)
    for (year.in.range in year.range.values) {
        year.range.matrix = cbind(year.range.matrix, value.vector)
    }
    year.range.matrix
}

year.range.extended = matrix(nrow = nrow(diagonal.sites), ncol=0)
for (year.index in 3:ncol(diagonal.sites)) {

    year.range.extended = cbind(year.range.extended, extend.year.range.with.value(colnames(diagonal.sites)[year.index], diagonal.sites[, year.index]))
}
year.range.extended = cbind(diagonal.sites[, c("Origin", "Destination")], year.range.extended)
colnames(year.range.extended)[3:ncol(year.range.extended)] = sapply(1989:2008, toString)
rownames(year.range.extended) = year.range.extended[, 1]
year.range.extended = year.range.extended[, 3:ncol(year.range.extended)]

# Obtain all attractiveness (A_i,i) values

all.pairs.table = read.csv("allpairs.csv", header = T)
rownames(all.pairs.table) = all.pairs.table[, 1]
all.pairs.table = all.pairs.table[, 2:ncol(all.pairs.table)]
colnames(all.pairs.table) = sapply(1988:2015, toString)
known.attractiveness.values = all.pairs.table / colSums(all.pairs.table)

# Obtain 3-year lagged A_i,i values

calculate.attractiveness = function(attractive.matrix, matrix.index) {
    (attractive.matrix[, matrix.index] + attractive.matrix[, matrix.index - 1] + attractive.matrix[, matrix.index - 2] * .5 + attractive.matrix[, matrix.index - 3] * .25) / 2.5
}

known.lagged.attractiveness.values = matrix(nrow = nrow(known.attractiveness.values), ncol = ncol(known.attractiveness.values) - 3)
rownames(known.lagged.attractiveness.values) = rownames(known.attractiveness.values)
colnames(known.lagged.attractiveness.values) = colnames(known.attractiveness.values)[4:ncol(known.attractiveness.values)]
for (lagged.year.index in 1:ncol(known.lagged.attractiveness.values)) {
    known.lagged.attractiveness.values[, lagged.year.index] = calculate.attractiveness(known.attractiveness.values, lagged.year.index+3)
}

# Obtain the intrinsic quality of each site

sum.known.attractiveness.values = rowSums(known.attractiveness.values)
intrinsic.quality = sum.known.attractiveness.values / sum(sum.known.attractiveness.values)

# Estimate the constant component of fidelity as a logistic function of quality

# First, build a matrix with quality and phi_i,i

intrinsic.quality.calibration = intrinsic.quality[rownames(year.range.extended)]
intrinsic.quality.calibration.matrix = matrix(nrow = nrow(year.range.extended), ncol = ncol(year.range.extended), intrinsic.quality.calibration, byrow = F)
selected.phi.rows = unlist(year.range.extended)
dated.phi.indexes = which(selected.phi.rows >= 0)
selected.phi.filtered.rows = selected.phi.rows[dated.phi.indexes]
intrinsic.quality.filtered.calibration = unlist(intrinsic.quality.calibration.matrix)[dated.phi.indexes]
phi.adjustment.matrix = matrix(nrow = length(intrinsic.quality.filtered.calibration), ncol = 2)
colnames(phi.adjustment.matrix) = c("Quality", "Phi")
phi.adjustment.matrix[, 1] = intrinsic.quality.filtered.calibration
phi.adjustment.matrix[, 2] = selected.phi.filtered.rows

# Next, build a logistic regression model

logistic.quality = glm(Phi ~ Quality, family = binomial(link = 'logit'), data = as.data.frame(phi.adjustment.matrix))

# Finally, estimate the fidelity

estimate.phi = function(input.quality) {
    predict(logistic.quality, newdata = input.quality, type = 'response')
}
phi.weight = .9
lagged.attractiveness.weight = 1-phi.weight

estimate.fidelity = function(input.quality, ta.value, age.weight) {
    max(0, min(1, age.weight*(lagged.attractiveness.weight * (rnorm(length(input.quality)) + 1) * ta.value + phi.weight * estimate.phi(input.quality))))
}

known.sites = rownames(all.pairs.table)

# Initialize data structures
simulation.years = 50

yearly.fidelity = rep(0, times = length(known.sites))
BirdsStayingNextYear = BirdsLeavingNextYear = matrix(ncol = length(known.sites), nrow = simulation.years)
colnames(BirdsStayingNextYear) = colnames(BirdsLeavingNextYear) = known.sites
names(yearly.fidelity) = known.sites

# Simulate the system

for (t in 1:(simulation.years-1)) {
    for (site.name in known.sites) {
        # Calculate the fidelity in the current year
        yearly.fidelity[site.name] = estimate.fidelity(intrinsic.quality[site.name], known.lagged.attractiveness.values[site.name, t], age.weight)
        # Calculate the number of birds that will remain and stay in the colony the next year
        BirdsStayingNextYear[site.name, t + 1] = BirdsStayingNextYear[site.name, t] * averaged.survival[site.name] * yearly.fidelity[site.name] # + Immigrants
        BirdsLeavingNextYear[site.name, t + 1] = BirdsNextYear[site.name, t] * averaged.survival[site.name] * (1 - yearly.fidelity[site.name])
        yearly.attractiveness = BirdsStayingNextYear[, t + 1] / sum(BirdsStayingNextYear[, t + 1])

        # Update attractiveness
        known.attractiveness.values = cbind(known.attractiveness.values, yearly.attractiveness)
        yearly.lagged.attractiveness = calculate.attractiveness(known.attractiveness.values, t)
        known.lagged.attractiveness.values = cbind(known.lagged.attractiveness.values, yearly.lagged.attractiveness)
        colnames(known.attractiveness.values)[ncol(known.attractiveness.values)] = colnames(known.lagged.attractiveness.values)[ncol(known.lagged.attractiveness.values)] = toString(as.numeric(colnames(known.lagged.attractiveness.values)[ncol(known.lagged.attractiveness.values)])+1)
        

    }
}
