# Win ratio and win odds of prioritized OS, PFS and response for one trial
#
# Every experimental patient is compared with every control patient. The pair
# is decided on OS (Gehan scores) if possible, otherwise on PFS (Gehan scores),
# otherwise on response (a responder wins against a non-responder); pairs not
# decided on any endpoint are ties. With p_win and p_loss the proportions of
# pairs won and lost by the experimental patient and p_tie = 1 - p_win -
# p_loss,
#   win ratio WR = p_win / p_loss, net benefit NB = p_win - p_loss,
#   win odds WO = (p_win + p_tie / 2) / (p_loss + p_tie / 2) = (1 + NB) / (1 - NB).
# Standard errors are first-order two-sample U-statistic (Hajek projection)
# estimates: Var(p_win) = s1^2 / n1 + s0^2 / n0, where s1^2 is the variance
# (divisor n1) over experimental patients of their proportions of pairs won,
# and similarly for the control patients, the losses and the covariance. The
# delta method gives the standard errors of log(WR) and log(WO).
# The function also returns the contribution of each endpoint to the wins and
# losses, and the marginal proportions of pairs won and lost when each
# endpoint is compared alone.
#
# Arguments
#   trt       logical, TRUE for experimental patients
#   os_tte, os_event, pfs_tte, pfs_event observed OS and PFS
#   resp      0/1 response observed by the analysis cutoff
# Value
#   named numeric vector
win_statistics <- function(trt, os_tte, os_event, pfs_tte, pfs_event, resp) {
  trt <- as.logical(trt)
  n1 <- sum(trt)
  n0 <- sum(!trt)
  os <- gehan_scores(os_tte[trt], os_event[trt], os_tte[!trt], os_event[!trt])
  pf <- gehan_scores(pfs_tte[trt], pfs_event[trt], pfs_tte[!trt], pfs_event[!trt])
  r1 <- resp[trt] == 1
  r0 <- resp[!trt] == 1
  rw <- outer(r1, r0, ">")
  rl <- outer(r1, r0, "<")
  und1 <- !(os$win | os$loss)
  w2 <- und1 & pf$win
  l2 <- und1 & pf$loss
  und2 <- und1 & !(pf$win | pf$loss)
  w3 <- und2 & rw
  l3 <- und2 & rl
  w <- os$win | w2 | w3
  l <- os$loss | l2 | l3
  p_win <- mean(w)
  p_loss <- mean(l)
  p_tie <- 1 - p_win - p_loss
  # first-order U-statistic variance (Hajek projection)
  aw <- rowMeans(w)
  al <- rowMeans(l)
  bw <- colMeans(w)
  bl <- colMeans(l)
  mv <- function(x, y) mean((x - mean(x)) * (y - mean(y)))
  v_w <- mv(aw, aw) / n1 + mv(bw, bw) / n0
  v_l <- mv(al, al) / n1 + mv(bl, bl) / n0
  c_wl <- mv(aw, al) / n1 + mv(bw, bl) / n0
  nb <- p_win - p_loss
  v_nb <- v_w + v_l - 2 * c_wl
  v_log_wr <- v_w / p_win ^ 2 + v_l / p_loss ^ 2 - 2 * c_wl / (p_win * p_loss)
  v_log_wo <- v_nb * (2 / (1 - nb ^ 2)) ^ 2
  wr <- p_win / p_loss
  wo <- (1 + nb) / (1 - nb)
  c(n1 = n1, n0 = n0, p_win = p_win, p_loss = p_loss, p_tie = p_tie,
    wr = wr, wo = wo, nb = nb,
    se_log_wr = sqrt(v_log_wr), se_log_wo = sqrt(v_log_wo), se_nb = sqrt(v_nb),
    z_wr = log(wr) / sqrt(v_log_wr), z_wo = log(wo) / sqrt(v_log_wo),
    win_os = mean(os$win), loss_os = mean(os$loss),
    win_pfs = mean(w2), loss_pfs = mean(l2),
    win_resp = mean(w3), loss_resp = mean(l3),
    mwin_os = mean(os$win), mloss_os = mean(os$loss),
    mwin_pfs = mean(pf$win), mloss_pfs = mean(pf$loss),
    mwin_resp = mean(rw), mloss_resp = mean(rl))
}
