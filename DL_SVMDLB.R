# DL_SVMDLB

# 1-2: minimize hinge loss function
svm_lasso = function(r, z, lambda)
{
   z_nonzero = z[abs(z) > 0]
   if (length(z_nonzero) < 1)
      return(0)
   r_nonzero = r[abs(z) > 0]
   candidates = r_nonzero / z_nonzero
   candidates = append(candidates, 0)
   return(find_which_min(r_nonzero, z_nonzero, candidates, lambda))
}

# 2
find_which_min = function(r, z, b, lambda)
{
   n = length(b)
   hr = rep(0,n)
   for (i in 1:n)
      hr[i] = mean(svm_plus(r - b[i] * z)) + lambda * abs(b[i]);
   return (b[which.min(hr)])
}

# 3
lambdas_all = function(number_lambdas, lambda_max, epsilon_lambda)
{
   lambdas_all = rep(0, number_lambdas)
   lambda_min = epsilon_lambda * lambda_max
   ratio_max_min = 1.0 / epsilon_lambda
   div = exponent = 0
   div = number_lambdas - 1
   for (lambda_index in 1:number_lambdas-1)
   {
      exponent = lambda_index / div
      lambdas_all[lambda_index + 1] = lambda_min * (ratio_max_min^exponent)
   }
   lambdas_all[1] = 0
   return(lambdas_all)
}

# 4
svm_plus = function(h)
{
   h_plus = h
   h_plus[h < 0] = 0
   return(h_plus)
}

# 5
bspline2 = function(h, m) # µ¹¸±¶§ 
{
   if (m <= 1)
   {
      a = h
   }
   else
   {
      a = rep(0, length(h))
      left = -1 < h & h <= 0
      a[left] = h[left] + 1
      right = 0 < h & h < +1
      a[right] = 1 - h[right]
   }
   return(a)
}


bspline = function(h,m)
{
   sigma = rep(0, length(h))
   abs_h = 0
   for (i in 1:length(h))
   {
      abs_h = abs(h[i])
      if (abs_h < 1)
         sigma[i] = 1 - abs_h
   }
   return(sigma)
}



# 6
bspline_derivative2 = function(h, m)
{
   n = length(h)
   if (m <= 1)
   {
      da = rep(1, n)
   }
   else
   {
      da = rep(0, n)
      da[-1 < h & h <= 0] = +1
      da[0 < h & h < +1] = -1
   }
   return(da)
}


# 7
bspline_derivative = function(h, m, degree = 1)
{
   {
      i = n = length(h)
      dsigma = rep(0, n)
      abs_h = 0
      if (degree == 1)
      {
         for (i in 1:n)
            if (abs(h[i]) < 1)
               dsigma[i] = - sign(h[i])
      }
      else if (degree == 2)
      {
         for (i in 1:n)
         {
            abs_h = abs(h[i])
            if (abs_h < 0.5)
               dsigma[i] = - 2.0 * h[i]
            else if (abs_h < 1.5)
               dsigma[i] = sign(h[i]) * (abs_h - 1.5)
         }
      }
      else
         stop("unknown degree in bspline")
      return(dsigma)
   }
}
