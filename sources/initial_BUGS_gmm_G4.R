
x1 = c(0:(Time-1))
lm.cf.est = apply(y, 1, function(x) lm(x~x1)$coefficients)
lm.rs.est = apply(y, 1, function(x) summary(lm(x~x1))$sigma^2)
I.quant = quantile(lm.cf.est[1,], probs= c(0.80, 0.60, 0.40, 0.20))
S.quant = quantile(lm.cf.est[2,], probs= c(0.80, 0.60, 0.40, 0.20))
muLS.init = cbind(I.quant, S.quant)
