# Fig 4 (main text) from the v6 PRIMARY all-data A-vs-CTD model. Posterior MEDIAN
# throughout (phi per family A vs CTD at the Summer/offshore/avg-site baseline + OR).
suppressMessages({library(ggplot2); library(patchwork)})
load("outputs/v6/MOSAIC_mammal_families_twolevel_stage1_AvsD_v6_results.RData")
n_families <- length(focal_families); S <- samples_combined; d_method <- S[,"d_method"]
q3 <- function(x) c(mid=median(x), lo=quantile(x,.025), hi=quantile(x,.975))
pred <- do.call(rbind, lapply(seq_len(n_families), function(f){
  a <- S[, sprintf("alpha[%d]",f)]
  rbind(data.frame(family=focal_families[f], method="CTD",      t(q3(plogis(a)))),
        data.frame(family=focal_families[f], method="MOSAIC-A", t(q3(plogis(a+d_method)))))}))
names(pred) <- c("family","method","mid","lo","hi")
ord <- pred$family[pred$method=="MOSAIC-A"][order(pred$mid[pred$method=="MOSAIC-A"])]
pred$family <- factor(pred$family, levels=ord)
wong <- c("MOSAIC-A"="#0072B2","CTD"="#D55E00")
pa <- ggplot(pred, aes(mid, family, colour=method)) +
  geom_linerange(aes(xmin=lo,xmax=hi), position=position_dodge(.5), linewidth=.6) +
  geom_point(position=position_dodge(.5), size=2.4) +
  scale_colour_manual(values=wong, name=NULL) + xlim(0,1) +
  labs(x=expression(phi~"(capture probability)"), y=NULL) + theme_bw() + theme(legend.position="top")
orv <- exp(d_method)
or <- data.frame(mid=median(orv), lo=quantile(orv,.025), hi=quantile(orv,.975))
pb <- ggplot(or, aes(mid, "MOSAIC-A vs CTD")) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  geom_linerange(aes(xmin=lo,xmax=hi), linewidth=.7) +
  geom_point(size=3, colour="#0072B2") + scale_x_log10() +
  labs(x="odds ratio (log scale)", y=NULL) + theme_bw()
fig <- pa/pb + plot_layout(heights=c(3,1)) + plot_annotation(tag_levels="A")
ggsave("outputs/v6/figures/fig4_v6primary.pdf", fig, width=7.5, height=7, device=grDevices::cairo_pdf)
cat(sprintf("v6 primary Fig4: OR median=%.2f (95%%CrI %.2f-%.2f) P(A>CTD)=%.3f\n",
    median(orv), quantile(orv,.025), quantile(orv,.975), mean(d_method>0)))
