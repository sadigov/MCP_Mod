library(MCPModPack)

#Simulation-based evaluation of dose-finding trials with a count endpoint
#Plot od an emax model


# Endpoint type
endpoint_type = "Binary"

# Select the candidate dose-response models and initial values 
# of the non-linear model parameters (linear, quadratic, exponential, 
# emax, logistic and sigemax)
models = list(exponential = 0.2, 
              emax = 0.25,
              emax = 0.5,
              sigemax = c(0.25, 3),
              sigemax = c(0.5, 4),
              logistic = c(0.5, 0.1))


doses <- c(0, 1, 2, 3)/3 # normalized equidistant doses

# mods5 <- Mods(emax = c(0.25, 0.5), 
#               sigEmax = rbind(c(0.25, 3), c(0.5, 4)), 
#               betaMod = c(1.5, 0.1),
#               placEff = logit(p.ocr), 
#               maxEff = logit(p.magli) - logit(p.ocr), 
#               doses = doses)



# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship
direction = "decreasing"

# Model selection criterion
model_selection = "aveAIC"

# The treatment effect for identifying the target dose 
# (this effect is defined relative to the placebo effect)

p.ocr <- 0.43 # Based on oratorio at W96

Delta = -p.ocr*0.5  # optimistic TPP = 50% relative risk reduction


# Select the assumed dose-response model and values of the non-linear model parameters
sim_models = list(emax = 0.5, 
                  placebo_effect = p.ocr, 
                  max_effect = seq(from = -0.25, to = -0.1, by = 0.005))

# Simulation parameters 
# (go threshold is defined relative to the placebo effect)
sim_parameters1 = list(n = c(85, 85, 85, 85),
                      doses = c(0, 1, 2, 3)/3,
                      dropout_rate = 0.1,
                      nsims = 1000,
                      go_threshold = -p.ocr*0.3) # Target TPP

sim_parameters2 = list(n = c(100, 100, 100, 100),
                      doses = c(0, 1, 2, 3)/3,
                      dropout_rate = 0.05,
                      nsims = 1000,
                      go_threshold = -p.ocr*0.3) # Target TPP



# Perform an MCP-Mod-based simulation
results = MCPModSimulation(endpoint_type = endpoint_type, 
                           models = models, 
                           alpha = alpha, 
                           direction = direction, 
                           model_selection = model_selection, 
                           Delta = Delta,
                           sim_models = sim_models,
                           sim_parameters = sim_parameters1)

# Simple summary of the MCPMod simulation results
results

# Detailed summary of the MCPMod simulation results (remove tempfile)

setwd("/Users/sadikhos/Library/CloudStorage/GoogleDrive-shamil.sadikhov@roche.com/My Drive/My Documents/CNS/MS/MAGLi/CDP/2023/Sample_Size/Alun/")

SimulationReport(results, 
                 "MCPMod simulation summary_cCDP12.docx", 
                 "MCPMod simulation summary_cCDP12.docx")





