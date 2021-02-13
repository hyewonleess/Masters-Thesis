# SVMDLB R code

# ver.1 by Hyewon Lee
# SVM Deep Learning for Bivariate Tensor B-Splines

SVMDLB = function(response,
                  predictors,
                  number_nodes = 10,
                  number_lambdas = 100,
                  degree = 1,
                  lambda_max = 10.0,
                  lambda_max_min_ratio = 1e-10,
                  mu = 0.0,
                  sd = 0.1,
                  max_iterations = 1000,
                  epsilon_iterations = 1e-4,
                  verbose = T)

{
   cat("=================================================\n")
   cat("Single layer Deep Learning Univariate SVM \n")
   cat("by HW from SDMLAB (January 20, 2021)\n")
   cat("Department of Statistics, Korea University, Korea\n")
   cat("=================================================\n")
   svm = list()
   sample_size = length(response)
   dimension = 0
   x1 = predictors[,1]
   x2 = predictors[,2]
   knots1 = runif(number_nodes, min(x1), max(x1))
   knots2 = runif(number_nodes, min(x2), max(x2))
   beta1 = rnorm(number_nodes, mu, sd)
   beta2 = rnorm(number_nodes, mu, sd)
   gamma = rnorm(number_nodes, mu, sd)
   alpha1 = - knots1 * beta1
   alpha2 = - knots2 * beta2
   fitted_values = rep(0, sample_size)
   lambda_list = lambdas_all(number_lambdas, lambda_max, lambda_max_min_ratio)
   partial_fitted_values = rep(0, sample_size)
   univariate_lasso_x = rep(0, sample_size)
   univariate_lasso_y = rep(0, sample_size)
   active_nodes = 1:number_nodes
   candidate_nodes = 1:number_nodes
   iteration_list = rep(0, number_lambdas)
   number_nodes_list = rep(0, number_lambdas)
   active_planes = matrix(0, sample_size, number_nodes)
   active_planes1 = matrix(0, sample_size, number_nodes)
   active_planes2 = matrix(0, sample_size, number_nodes)
   derivative_active_planes1 = matrix(0, sample_size, number_nodes)
   derivative_active_planes2 = matrix(0, sample_size, number_nodes)
   planes1 = matrix(0, sample_size, number_nodes)
   planes2 = matrix(0, sample_size, number_nodes)
   
   lambda =  penalty = H = Hlambda = store_Hlambda = Inf
   a = m = iter = lambda_index = number_active_nodes = number_nodes
   for (m in 1:number_nodes)
   {
      planes1[,m] = alpha1[m] + beta1[m] * x1
      planes2[,m] = alpha2[m] + beta2[m] * x2
      active_planes1[,m] = bspline(planes1[,m], m)
      active_planes2[,m] = bspline(planes2[,m], m)
      active_planes[,m] = active_planes1[,m] * active_planes2[,m]
      derivative_active_planes1[,m] = bspline_derivative(planes1[,m], m)
      derivative_active_planes2[,m] = bspline_derivative(planes2[,m], m)
      fitted_values = fitted_values + gamma[m] * active_planes[,m]
   }
   
   derivative_f = rep(0, sample_size)
   R_m = rep(0, sample_size)
   z_m = rep(0, sample_size)
   old_active_planes_m = rep(0, sample_size)
   
   # Model update
   for (lambda_index in 1:number_lambdas)
   {  
      if (verbose)
         cat("\n", lambda_index, "th lambda runs \n")
      lambda = lambda_list[lambda_index]
      # iteration starts
      for (iter in 1:max_iterations)
      {
         
         # 1. gamma
         
         
         for (a in from_to(1, number_active_nodes))
         {
            m = active_nodes[a]
            partial_fitted_values = fitted_values - active_planes[,m] * gamma[m]
            R_m = 1.0 - response * partial_fitted_values 
            z_m = response * active_planes[,m]
            gamma[m] = svm_lasso(R_m, z_m, lambda)
            fitted_values = partial_fitted_values + active_planes[,m] * gamma[m]
         }
         if (lambda_index == 4)
            cat("\n", gamma)
         active_nodes = candidate_nodes[gamma != 0] # pruning
         number_active_nodes = length(active_nodes) 
         # 2. alpha1
         for (a in from_to(1, number_active_nodes))
         {
            m = active_nodes[a]
            old_active_planes_m = active_planes[,m]
            derivative_f = gamma[m] * derivative_active_planes1[,m] * active_planes2[,m]
            partial_fitted_values = fitted_values - derivative_f * alpha1[m]
            R_m = 1.0 - response * partial_fitted_values
            z_m = response * derivative_f
            alpha1[m] = svm_lasso(R_m, z_m, lambda)
            planes1[,m] = alpha1[m] + beta1[m] * x1
            active_planes1[,m] = bspline(planes1[,m], m)
            derivative_active_planes1[,m] = bspline_derivative(planes1[,m], m)
            active_planes[,m] = active_planes1[,m] * active_planes2[,m]
            fitted_values = fitted_values - gamma[m] * (old_active_planes_m - active_planes[,m])
         }
         # 3. beta 1
         for (a in from_to(1, number_active_nodes))
         {
            m = active_nodes[a]
            old_active_planes_m = active_planes[,m]
            derivative_f = gamma[m] * x1 * derivative_active_planes1[,m] * active_planes2[,m]
            partial_fitted_values = fitted_values - derivative_f * beta1[m]
            R_m = 1.0 - response * partial_fitted_values
            z_m = response * derivative_f
            beta1[m] = svm_lasso(R_m, z_m, lambda)
            planes1[,m] = alpha1[m] + beta1[m] * x1
            active_planes1[,m] = bspline(planes1[,m], m)
            derivative_active_planes1[,m] = bspline_derivative(planes1[,m], m)
            active_planes[,m] = active_planes1[,m] * active_planes2[,m]
            fitted_values = fitted_values - gamma[m] * (old_active_planes_m - active_planes[,m])
         }
         # 4. alpha2
         for (a in from_to(1, number_active_nodes))
         {
            m = active_nodes[a]
            old_active_planes_m = active_planes[,m]
            derivative_f = gamma[m] * active_planes1[,m] * derivative_active_planes2[,m]
            partial_fitted_values = fitted_values - derivative_f * alpha2[m]
            R_m = 1.0 - response * partial_fitted_values
            z_m = response * derivative_f
            alpha2[m] = svm_lasso(R_m, z_m, lambda)
            planes2[,m] = alpha2[m] + beta2[m] * x2
            active_planes2[,m] = bspline(planes2[,m], m)
            derivative_active_planes2[,m] = bspline_derivative(planes2[,m], m)
            active_planes[,m] = active_planes1[,m] * active_planes2[,m]
            fitted_values = fitted_values - gamma[m] * (old_active_planes_m - active_planes[,m])
         }
         # 5. beta2
         
         for (a in from_to(1, number_active_nodes))
         {
            m = active_nodes[a]
            old_active_planes_m = active_planes[,m]
            derivative_f = gamma[m] * x2 * active_planes1[,m] * derivative_active_planes2[,m]
            partial_fitted_values = fitted_values - derivative_f * beta2[m]
            R_m = 1.0 - response * partial_fitted_values
            z_m = response * derivative_f
            beta2[m] = svm_lasso(R_m, z_m, lambda)
            planes2[,m] = alpha2[m] + beta2[m] * x2;
            active_planes2[,m] = bspline(planes2[,m], m)
            derivative_active_planes2[,m] = bspline_derivative(planes2[,m], m)
            active_planes[,m] = active_planes1[,m] * active_planes2[,m]
            fitted_values = fitted_values - gamma[m] * (old_active_planes_m - active_planes[,m])
         }
         
         # if (lambda_index == 4 & iter == 2)
         # {
         #    cat("\n", "derivative_f ", derivative_f , "\n")
         #    
         #    
         # }
         
         # check model convergence (H: hinge loss function)
         H = mean(svm_plus(1.0 - response * fitted_values))
         Hlambda = H
         if (abs(Hlambda - store_Hlambda) < epsilon_iterations)
            break
         store_Hlambda = Hlambda
         if (verbose)
            cat("\n", "Number of iterations = ", iter, ", Hlambda = ", Hlambda)
      }
      # results
      pruned_alpha1 = alpha1[active_nodes]
      pruned_alpha2 = alpha2[active_nodes]
      pruned_beta1 = beta1[active_nodes]
      pruned_beta2 = beta2[active_nodes]
      pruned_gamma = gamma[active_nodes]
      dimension = number_active_nodes + sum(pruned_alpha1 != 0) + sum(pruned_alpha2 != 0) +
         sum(pruned_beta1 != 0) + sum(pruned_beta2 != 0)
      number_nodes_list[lambda_index] = number_active_nodes
      iteration_list[lambda_index] = iter
      svm[[lambda_index]] = list(lambda = lambda,
                                 number_nodes = number_active_nodes,
                                 alpha1 = pruned_alpha1,
                                 alpha2 = pruned_alpha2,
                                 beta1 = pruned_beta1,
                                 beta2 = pruned_beta2,
                                 gamma = pruned_gamma,
                                 fitted_values = fitted_values,
                                 active_nodes = active_nodes)
   }
   svm$lambda_list = lambda_list
   svm$iteration_list = iteration_list
   svm$number_nodes_list = number_nodes_list
   cat("\n\n", "Done!", "\n\n")
   return(svm)
}


from_to = function(from, to)
{
   if (from > to)
      return(NULL)
   else
      return(from : to)
}





