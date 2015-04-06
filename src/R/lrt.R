lrt <- function(pop_lhc, pop_hn, rest_lhc, rest_hn) {
  pop_lhc <- pop_lhc + 1
  rest_lhc <- rest_lhc + 1
  pop_lhf <- pop_lhc/pop_hn
  rest_lhf <- rest_lhc/rest_hn
  global_lhc <- pop_lhc + rest_lhc
  global_hn <- pop_hn + rest_hn
  global_lhf <- global_lhc/global_hn
  l_null <- global_lhc*log(global_lhf) + (global_hn-global_lhc)*log(1-global_lhf)
  l_mle <- pop_lhc*log(pop_lhf) + (pop_hn-pop_lhc)*log(1-pop_lhf) + rest_lhc*log(rest_lhf) + (rest_hn-rest_lhc)*log(1-rest_lhf)
  chisq <- 2*(l_mle-l_null)
  return(chisq)
}