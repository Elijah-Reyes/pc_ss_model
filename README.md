# Parental care and sexual selection model (pc_ss_model)
This individual-based simulation model of parental care and sexual selection is from the study "Interactions between parental care and sexual selection as a driver of reproductive isolation" (Reyes et al., In prep 2024)

Any feedback on this code is welcome!

If there are any questions or issues with the code, feel free to contact me at `elijahr@sfu.ca`
\- Elijah Reyes

## To run a single replicate\:

`run_test.R` sources the model and runs a single replicate for the parameters you have put into the arguements (default paramter values are set in src/prms.R).
Before running, you will need to set the working directory `setwd()` found in the first line to wherever you have stored the files (`run_test.R` and `src`).

`src` stores the model. Everything in there is sourced via `src/init.R`, which is sourced in `run_test.R`.
Any adjustments you would like to make to the model itself, not just to the parameter space, should be done in here.

Inside `src` you will find the following: 
  - `init.R`  
  - `lifecycle.R`
  - `meta_fxns.R`
  - `pop_fxns.R`
  - `prms.R`

Next, we'll go through what each of these does:

#### `init.R`
Initializes the model by sourcing the other files. Additionally, it makes functions for creating the initial population and other information that will need to be stored throughout the model.

#### `prms.R`
Sets the default parameters for the model. Anything listed here can be added to the arguments of `run.model()` to change the parameters that you are running for the replicate. However, if you would like to change what the default values are (i.e., so you don't have to include in the arguements of `run.model()`), you can change them here.

#### `pop_fxns.R`
Assigns functions for `lifecycle.R` to use. Includes events such as mating, mortality, mutation, and maturity.

#### `lifecycle.R`
Calls the functions made in `pop_fxns.R`. Importantly, sets the order for life events that individuals experience. e.g., if death or birth occurs first in the lifecycle. 

#### `meta_fxns.R`
Packages the creation of all the inital objects and then loops through generations of the lifecycle. Additionally, packages the data for output. (Optional: plot code for viewing a single replicate run as it progresses.) `run_test.R` calls the function made here to run the model.

## Looking at output
After running a replicate, the output will be in a list called `out`. In its current state, the following is stored there:
1. whether the popualtion went extinct
2. the population in the final generation
3. the fitness measures over time
4. the distributions of the trait values over time
5. measurements of popualtion growth
6. the initial population
7. the summary statistics of the population over time
8. the parameters used to run this output 

What is stored here can be adjusted in `meta_fxns.R`


## For running multiple parameter combinations:
While I have not included it here, I created a dataframe with each row being the parameter values of what to run. I then ran a loop through the number of rows for the parameter index I made. This method very easily can be run in parallel instead of a loop.  

