## Process simulation outputs
setwd("C:/Users/smajumdar/Documents/Depth-scatter1/Codes")
Required.Packages <- c("data.table","tidyverse","patchwork","reshape2","latex2exp")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

# load files
normt4 = readRDS("elliptic_results.rds")
normcont4 = readRDS("nonelliptic_results.rds")

# FSE of elliptical dists
n = seq(50,500,by=50)
method_names = c(
    "SCM","Tyler",
    sprintf(r'($\tilde{\Sigma}-H$)'), sprintf(r'($\Sigma_*-H$)'),
    sprintf(r'($\tilde{\Sigma}-P$)'), sprintf(r'($\Sigma_*-P$)'),
    sprintf(r'($\tilde{\Sigma}-M$)'), sprintf(r'($\Sigma_*-M$)')
)
err_list = lapply(normt4, function(x) x[[1]])
eff_list = lapply(err_list, function(x){
    xmat = matrix(x[,1], nrow=10, ncol=8, byrow=F)/x[,-1]
    colnames(xmat) = method_names
    cbind(n, xmat)
})

# FSE in presence of contamination: divide by uncontaminated error
eff_list2 = lapply(normcont4[1:2], function(x) cbind(n, matrix(err_list[[1]][,1], nrow=10, ncol=8, byrow=F)/x[[1]][,-1]))
cop.table = normcont4[[3]][[1]]
eff_list2[[3]] = cbind(n, matrix(cop.table[,1], nrow=10, ncol=8, byrow=F)/cop.table[,-1])
for(i in 1:3){
    colnames(eff_list2[[i]])[-1] = method_names
}

# plot
plot_table = function(plotmat, line_width=.5, point_size=1){
    # plot(plotmat[,1], type='b', ylim=c(0,1))
    # for(i in 2:8){
    #     lines(plotmat[,i], type='b', col=i)
    # }
    # legend("bottomright", col=1:6, lty=1, legend=TeX(method_names))
    # ggplot
    eff_melt = melt(data.table(plotmat), id.vars="n", measure.vars=colnames(plotmat)[-1])
    names(eff_melt) = c("n","Method","FSE")
    p1 = ggplot(eff_melt, aes(x=n, y=FSE, colour=Method)) +
        geom_line(size=line_width) +
        geom_point(size=point_size) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90)) +
        # theme(legend.position = c(.82,.4)) +
        xlab(expression(n)) +
        scale_x_continuous(breaks=n) +
        scale_color_manual(labels = TeX(method_names), values=1:8)
    p1
}
plot_list = lapply(eff_list, plot_table)
plot_list2 = lapply(eff_list2, plot_table)
plot_list2[[2]] = NULL
dist_names = c("Normal", sprintf(r'($t_3$)'),sprintf(r'($t_{10}$)'),sprintf(r'($t_{20}$)'),
               "Contaminated normal","Copula")
plot_list = lapply(1:4, function(i) plot_list[[i]]+ggtitle(TeX(dist_names[i])))
plot_list2 = lapply(1:2, function(i) plot_list2[[i]]+ggtitle(TeX(dist_names[i+4])))
combined <- plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[1]] +
    plot_list2[[1]] + plot_list2[[2]] & theme(legend.position = "bottom")

# save plot
pdf("simplot.pdf", width=8, height=6)
combined + plot_layout(ncol=3, nrow=2, guides = "collect")
dev.off()