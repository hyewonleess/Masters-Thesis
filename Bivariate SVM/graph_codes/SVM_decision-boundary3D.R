# func: decision boundary function

## define grids on x1 and x2
n_grid = 50
px1 = seq(min(x[, 1]), max(x[, 1]), length = n_grid)
px2 = seq(min(x[, 2]), max(x[, 2]), length = n_grid)
xgrid = expand.grid(px1, px2)

# func: decision boundary function on SVM
func = rep(0, n_grid^2)
for(m in 1 : res$number_nodes)
{
   func = func + res$gamma[m] * 
         bspline(res$alpha1[m] + res$beta1[m] * xgrid[ ,1], res$active_nodes[m]) *
         bspline(res$alpha2[m] + res$beta2[m] * xgrid[ ,2], res$active_nodes[m])
}

func = func + res$bias

# 3d plot1: class(1, -1) seperation
func_1 = func; func_2 = func
func_1[which(func<0)]  = 0
func_2[which(func>0)]  = 0
func_1 = matrix(func_1, n_grid, n_grid)
func_2 = matrix(func_2, n_grid, n_grid)

x_1 = x[1:100, ]; x_2 = x[101:200, ]
y_1 = y[1:100]; y_2 = y[101:200]

persp3d(px1, px2, func_1, col = 'orange', cex = 0.2, back = 'lines', front = 'lines', theta = 30, pi = 30)
persp3d(px1, px2, func_2, col = 'blue', cex = 0.2, back = 'lines', front = 'lines', add = TRUE, theta = 30, pi = 30, contour = T)
points3d(x_1[, 1], x_1[, 2], y_1, col = "blue", size = 10, add=T)
points3d(x_2[, 1], x_2[, 2], y_2, col = "orange", size = 10, add=T)

# 3d plot2: overall decision boundary plot
func = matrix(func, n_grid, n_grid)
persp3d(px1, px2, func, col = 'navy', cex = 0.2, back = 'lines', front = 'lines', expand = 0.5, contour = T)
