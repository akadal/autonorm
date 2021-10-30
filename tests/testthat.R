#library("testthat")


# including the library
source("R/index.R")


# first example
dataset_digimon <- read.csv("tests/DigiDB_digimonlist.csv", header = FALSE, sep = ",")
autonormResult <- autonorm(dataset_digimon)
View(dataset_digimon)

# second example
dataset_bites <- read.csv("tests/Health_AnimalBites.csv", header = FALSE, sep = ",")
autonormResult <- autonorm(dataset_bites)
View(dataset_bites)

# third example
dataset_payments <- read.csv("tests/online_payments.csv", header = FALSE, sep = ",")
autonormResult <- autonorm(dataset_payments)
View(dataset_payments)

dataset_access <- read.csv("tests/SNAP_Online_State_EBT_Accounts.csv", header = FALSE, sep = ";")
autonormResult <- autonorm(dataset_access)
View(dataset_access)

dataset_transactions <- read.csv("tests/new_transaction.csv", header = T, sep = ",")
autonormResult <- autonorm(dataset_access, maxiter = 100, w = c(35,35,15,15))
View(dataset_access)
