
MREM.plotParamTraj = function(MCEM_fit) {
  graphics::par(mfrow = c(3,1), mar = c(2, 4.1, 1, 2.1), xaxs = "i")
  graphics::matplot(t(MCEM_fit$paramTraj$pis), type = 'l', lty = 1, ylab = "pis", xaxt = 'n', bty = "L")
  graphics::matplot(t(MCEM_fit$paramTraj$mus), type = 'l', lty = 1, ylab = "mus", xaxt = 'n', bty = "L")
  graphics::matplot(t(MCEM_fit$paramTraj$sds), type = 'l', lty = 1, ylab = "sds", bty = "L")
  graphics::par(mfrow = c(1,1))
  ## TODO: add legend 
}
