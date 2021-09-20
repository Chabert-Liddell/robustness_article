################################################################################
# Classif avec degree as edge covariates
################################################################################
library(cropanalysis)
library(tidyverse)
library(blockmodels)
library(robber)

param_degcov <-
  readRDS(file = "param_degcov.rds")


if (! file.exists("rob_cov_dc.rds")) {

  rob_cov_dc <- purrr::map_dbl(
    .x = seq_along(param_degcov),
    .f = function(i) {
      if(is.null(param_degcov[[i]])) return(0)
      nr <- web_of_life[[i]]$nr
      nc <- web_of_life[[i]]$nc
      param <- param_degcov[[i]]
      rob.fun <- matrix(1, nr+1, max(300, nr+nc))
      rob.fun[nr+1,] <- 0

      for (b in seq(ncol(rob.fun))) {
        for (q in seq_along(param$rho)) {
          #        browser()
          bq <- 0
          for (k in seq_along(param$pi)) {
            bq <- bq + (param$pi[k]/(1+exp(-param$m[k,q] -
                                             param$thetanu[,,drop=FALSE])))
          }
          rob.fun[1, b] <- rob.fun[1,b] -
            mean(param$rho[q] *
                   (exp(colSums(log(
                     1 - bq)))))
        }
        s <- sample(nr, nr)
        for (m in seq(nr-1)) {
          for (q in seq_along(param$rho)) {
            bq <- 0
            for (k in seq_along(param$pi)) {
              bq <- bq + (param$pi[k]/(1+exp(-param$m[k,q] -
                                               param$thetanu[s[(m+1):nr],,drop=FALSE])))
            }
            rob.fun[m+1, b] <- rob.fun[m+1,b] -
              mean(param$rho[q] *
                     (exp(colSums(log(
                       1 - bq)))))
          }
        }
      }
      rob.auc <- sum(rowMeans(rob.fun))/nr
      print(paste0(web_of_life[[i]]$id, ": ", round(rob.auc, digits = 3),
                   " vs ", round(tb_emp[i,]$uniform.emp.lin, 3), " vs ",
                   round(tb_robustness_wol[i, ]$bm.unif, 3)))
      return(rob.auc)
    }
  )
  saveRDS(object = rob_cov_dc,
          file = "rob_cov_dc.rds")
} else {
  rob_cov_dc <- readRDS(file = "rob_cov_dc.rds")
}




tb_bm <- readRDS(file = "tb_robustness_wol.rds")
tb_emp <- readRDS(file = "empirical_robustness.rds")

type <- stringr::str_sub(tb_bm$id, 1, 4)

dplyr::inner_join(tb_bm, tb_emp) %>%
  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::mutate("# Std Deviation" = abs(uniform.emp.lin - bm.unif)/sqrt(Variance),
                "Empty Column" = (1-dens)**nr) %>%
  dplyr::filter(nr>=10 & nc>=10) %>%
  ggplot2::ggplot(ggplot2::aes(x = uniform.emp.lin, y = bm.unif,#/(1-`Empty Column`),
                               size = `# Std Deviation`,
                               col = `Empty Column`)) +
  ggplot2::geom_abline(slope = 1, intercept = 0) +
  ggplot2::geom_point(alpha = .5) +
  ggplot2::scale_color_gradient( low = "grey50", high = "red") +
  ggplot2::xlim(c(0.5,1)) +
  ggplot2::ylim(c(0.5,1)) +
  ggplot2::xlab(label = latex2exp::TeX("\\hat{\\bar{R}} Uniform")) +
  ggplot2::ylab(label = latex2exp::TeX("\\bar{R} Uniform")) +
  ggplot2::scale_size(name = latex2exp::TeX("Z_R(A)")) +
  #  ggplot2::coord_fixed() +
  ggplot2::theme_minimal(base_size = 15, base_rect_size = .25)+
  ggplot2::theme(legend.box = "vertical") +
  #  ggplot2::annotate("point", x = tb_emp$uniform.emp.lin, y = rob_dc) +
  #  ggplot2::annotate("point", x = tb_emp$uniform.emp.lin, y = rob_edd, col = "blue") +
  # ggplot2::annotate("point", x = tb_emp$uniform.emp.lin[tb_bm$nr >=10 & tb_bm$nc >= 10],
  #                   y = rob_cov_edd[tb_bm$nr >=10 & tb_bm$nc >= 10], shape = as.factor(type[tb_bm$nr >=10 & tb_bm$nc >= 10]),
  #                   col = "green") +
  ggplot2::annotate("point", x = tb_emp$uniform.emp.lin[tb_bm$nr >=10 & tb_bm$nc >= 10], #shape = as.factor(type[tb_bm$nr >=10 & tb_bm$nc >= 10]), #stringr::str_sub(tb_bm$id, 1, 4),
                    y = rob_cov_dc[tb_bm$nr >=10 & tb_bm$nc >= 10],#*purrr::map_lgl(param_degcov, function(p) length(p$pi) > 1)*,
                    col = "purple", shape = 17)



# Calcul de la robustesse sur données wol avec SBM cov avec Q = 1
#

library(cropanalysis)
library(blockmodels)
if (! file.exists("xhat_edd.rds")) {
  xhat_edd <- pbmcapply::pbmclapply(
    X = web_of_life,
    FUN = function(wol) {
      if (wol$nr + wol$nc > 1000) return(NULL)
      covdc <- covariatesDegreeCorrectedLBM(wol$nr, wol$nc)
      lbm_cov <- blockmodels::BM_bernoulli_covariates_fast("LBM",
                                                           wol$net,
                                                           covariates=covdc,
                                                           verbosity=0,
                                                           plotting="",
                                                           exploration_direction = c(1,1))

      lbm_cov$estimate()
      return (list(edd = list(id = wol$id,
                              pred = lbm_cov$prediction(Q = 2))))
    }, mc.cores = 6L)
  saveRDS(xhat_edd, file = "article/xhat_edd.rds")
} else {
  readRDS(xhat_edd, file = "xhat_edd.rds")
}




rob_cov_edd <- purrr::map_dbl(
  .x = seq_along(xhat_edd),
  .f = function(i) {
    if(is.null(xhat_edd[[i]])) return(0)
    nr <- web_of_life[[i]]$nr
    nc <- web_of_life[[i]]$nc
    rob.fun <- matrix(1, nr+1, max(300, nr+nc))
    rob.fun[nr+1,] <- 0
    rob.fun[1, ] <- rob.fun[1,] -
      mean((exp(colSums(log(
        1 - xhat_edd[[i]]$edd$pred)))))
    for (b in seq(ncol(rob.fun))) {
      s <- sample(nr, nr)
      for (m in seq(nr-1)) {
        rob.fun[m+1, b] <- rob.fun[m+1,b] -
          mean((exp(colSums(log(
            1 - xhat_edd[[i]]$edd$pred[s[(m+1):nr],,drop=FALSE])))))
      }
    }
    rob.auc <- sum(rowMeans(rob.fun))/nr
    print(paste0(web_of_life[[i]]$id, ": ", round(rob.auc, digits = 3)))
    return(rob.auc)
  }
)



tb_bm <- readRDS(file = "tb_robustness_wol.rds")
tb_emp <- readRDS(file = "empirical_robustness.rds")

dplyr::inner_join(tb_bm, tb_emp) %>%
  #  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::mutate("# Std Deviation" = abs(uniform.emp.lin - bm.unif)/sqrt(Variance),
                "Empty Column" = (1-dens)**nr) %>%
  #  dplyr::filter(nr>=10 & nc>=10) %>%
  ggplot2::ggplot(ggplot2::aes(x = uniform.emp.lin, y = bm.unif,#/(1-`Empty Column`),
                               size = `# Std Deviation`,
                               col = `Empty Column`)) +
  ggplot2::annotate("point", x = tb_emp$uniform.emp.lin[],#tb_bm$nr >=10 & tb_bm$nc >= 10],
                    y = rob_cov_edd[],#tb_bm$nr >=10 & tb_bm$nc >= 10],
                    #shape = as.factor(type[tb_bm$nr >=10 & tb_bm$nc >= 10]),
                    col = "green", shape = "triangle", alpha = .5, size = 3)  +
  ggplot2::geom_abline(slope = 1, intercept = 0) +
  ggplot2::geom_point(alpha = .5, ggplot2::aes())+#shape = type)) +
  ggplot2::scale_color_gradient( low = "grey50", high = "red") +
  ggplot2::xlim(c(0.5,1)) +
  ggplot2::ylim(c(0.5,1)) +
  ggplot2::xlab(label = latex2exp::TeX("\\hat{R} Uniform")) +
  ggplot2::ylab(label = latex2exp::TeX("\\bar{R} Uniform")) +
  #  ggplot2::coord_fixed() +
  ggplot2::theme_minimal(base_size = 15, base_rect_size = .25)#+
