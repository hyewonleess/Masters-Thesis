# fitted b-spline activation functions
k = lambda_index
res = fit[[k]]

act_1 = matrix(0, nrow = length(xgrid[,1]), ncol = res$number_nodes)
act_2 = matrix(0, nrow = length(xgrid[,2]), ncol = res$number_nodes)

act = matrix(0, nrow = length(xgrid[,1]), ncol = res$number_nodes)

for(i in 1:res$number_nodes)
  {
   act_1[, i] = bspline(res$alpha1[i] + res$beta1[i] * xgrid[ ,1], res$active_nodes[i])
   act_2[, i] = bspline(res$alpha1[i] + res$beta1[i] * xgrid[ ,2], res$active_nodes[i])
   act[, i] =  bspline(res$alpha1[i] + res$beta1[i] * xgrid[ ,1], res$active_nodes[i]) *
                bspline(res$alpha2[i] + res$beta2[i] * xgrid[ ,2], res$active_nodes[i])
   
}

par(mfrow = c(1, 2))
matplot(xgrid[,1][1:50], act_1[1:50, ], type = "l", lty = 1, main = "x1 activation function", xlab = 'x1', ylab = "")
matplot(xgrid[,2], act_2, type = "l", lty = 1, main = "x2 activation function", xlab = 'x2', ylab = "")     
