
#install.packages("gdata")
#install.packages("ggplot2")
#install.packages("ggplot")
#install.packages("scales")



library(datasets) # 
library(ggplot2)
source(file="Chemical_Kinetics_Functions_M1.r")


##################################

#. Part I: Model M1

        
# 2)

lbound = 2
ubound = 4260
n_k = 100

priors = produce_all_priors_M1(lbound,ubound,n_k)

test_produce_priors_M1(lbound,ubound,n_k)

plot(priors[1:100, 2])
plot(priors[1:100, 3])
plot(priors[1:100, 4])

integrate_3densities_M1(priors)

####

# 3) 

k=1357
R0=1200
n_t=1000
step = (1000-0)/(n_t-1)
t_model <- seq(0,1000,step)
R_model = Compute_R_profile_M1(t_model, R0, k)
plot_profile_M1(t_model, R_model)



####

# 4) 

epsilon = 0.065
result = Compute_likelihood_all_M1(R0,lbound,ubound,n_k,epsilon)
Examine_likelihood_M1(R0,lbound,ubound,n_k,epsilon)


####

#5) 

lb_plot = 1000
ub_plot = 1500
Compute_all_posteriors_M1(R0,lbound,ubound,n_k,epsilon=0.065,lb_plot,ub_plot)


####

# 6) 

integrate_3densities_M1(priors, 1000, 1500)

####

# 7) 

posterior_k <- read.table("posterior_k_M1.csv", header=FALSE)
posterior_u <- read.table("posterior_u_M1.csv", header=FALSE)
posterior_w <- read.table("posterior_w_M1.csv", header=FALSE)


integrate_density_M1(priors[1:100, 1], posterior_k[1:100, 2], 1000, 1500 )
integrate_density_M1(priors[1:100, 1], posterior_u[1:100, 2], 1000, 1500 )
integrate_density_M1(priors[1:100, 1], posterior_w[1:100, 2], 1000, 1500 )


####

# 8) 

lb_plot = 2
ub_plot = 4260
Compute_all_posteriors_M1(R0,lbound,ubound,n_k,epsilon=0.3,lb_plot,ub_plot)

posterior_k <- read.table("posterior_k_M1.csv", header=FALSE)
posterior_u <- read.table("posterior_u_M1.csv", header=FALSE)
posterior_w <- read.table("posterior_w_M1.csv", header=FALSE)

integrate_density_M1(priors[1:100, 1], posterior_k[1:100, 2], 2, 4260 )
integrate_density_M1(priors[1:100, 1], posterior_u[1:100, 2], 2, 4260 )
integrate_density_M1(priors[1:100, 1], posterior_w[1:100, 2], 2, 4260 )


####

# 9)

integrate_density_M1(priors[1:100, 1], posterior_k[1:100, 2], 1000, 1500 )
integrate_density_M1(priors[1:100, 1], posterior_u[1:100, 2], 1000, 1500 )
integrate_density_M1(priors[1:100, 1], posterior_w[1:100, 2], 1000, 1500 )

epsilon = 0.3
result = Compute_likelihood_all_M1(R0,lbound,ubound,n_k,epsilon)
Examine_likelihood_M1(R0,lbound,ubound,n_k,epsilon)

epsilon = 0.065
result = Compute_likelihood_all_M1(R0,lbound,ubound,n_k,epsilon)
Examine_likelihood_M1(R0,lbound,ubound,n_k,epsilon)

###################################################################"

# Part II: Model M2


# 2) 

lbound = 5*10^(4)
ubound = 5*10^(9)
n_k = 100

test_produce_priors_M2(lbound,ubound,n_k)

#integrate_3densities_M2(priors)

###


# 3) 

R0=1200
n_t=1000
step = (1000000-0)/(n_t-1)
t_model <- seq(0,1000,step)
lbound = 5*10^(4)
ubound = 5*10^(9)
epsilon = 0.065
result = Compute_likelihood_all_M2(R0,lbound,ubound,n_k,epsilon)
Examine_likelihood_M2(R0,lbound,ubound,n_k,epsilon)

###

# 4) 

k=50000
R0=1200
n_t=1000
step = (1000-0)/(n_t-1)
t_model <- seq(0,1000,step)
R_model = Compute_R_profile_M2(t_model, R0, k)
plot_profile_M2(t_model, R_model)

###

# 5) 

lb_plot = 5*10^(4)
ub_plot = 5*10^(6)
Compute_all_posteriors_M2(R0,lbound,ubound,n_k,epsilon=0.065,lb_plot,ub_plot)


###################################################################"

#  III. Bayes factor


# 1)


#  2)   




###

#  3) 








































