if (length(commandArgs(trailingOnly = T)) == 0) {
  # Random seed for replicability
  s_val <- 50
  # 1: Sample randomly or 0: take highest N values
  r_val <- 0
  # Number of samples per iteration
  n_val <- 1
} else {
  args <- commandArgs(trailingOnly = T)
  arg.l <- do.call(rbind, strsplit(args, '\\s*=\\s*'))
  formatted.args <- arg.l[,2]
  names(formatted.args) <- arg.l[,1]
  s_val <- as.integer(formatted.args['s_val'])
  m_val <- formatted.args['m_val']
  n_val <- as.integer(formatted.args['n_val'])
}

## Tuning parameters
pi.xi <- 1
ucb.delta <- 0.1

theme_set(theme_bw())

## Helpers
logit <- function(x) {
  log(x/(1-x))
}

inv.logit <- function(x) {
  1/(1+exp(-x))
}

link <- logit
inv.link <- inv.logit

## "True" prevalence
## MODIFY FILE PATH!
xy.dt <- readRDS(file = 'PATH TO mwi_2017.Rds')
xy.dt[, i := .I]
n.x <- length(unique(xy.dt[, x]))
n.y <- length(unique(xy.dt[, y]))
xy.dt <- xy.dt[!is.na(prev)]

## Length of resulting gif
gif.time <- 10

## Include noise in GP?
noise <- 0.0
n.obs <- n_val

xy.dt[, true.val := link(prev)]
xy.dt[, true.prev := prev]
xy.dt[, obs.val := true.prev]
mu <- xy.dt[, mean(true.val)]

## Kernel for GP
k.obs <- rbfdot(25)
plot.dt <- data.table()
mse.dt <- data.table()
n.iter <- 50

## Iterate over methods
## SRS: simple random sampling
## EI: expected improvement
## UCB: upper confidence bound
## PI: probability of improvement
## VAR: highest variance
for (m_val in c('SRS', 'UCB')) {
  i <- 1
  curr.max <- Inf
  curr.sum <- Inf
  set.seed(s_val)
  xy.dt[, is.new := F]
  obs.i <- sample(1:nrow(xy.dt), n.obs)
  xy.dt[, in.sample := i %in% obs.i]


  while (i <= n.iter & n.obs > 0) {
    ## Covariance of observed points
    Sigma.obs <- kernelMatrix(k.obs, as.matrix(xy.dt[(in.sample), .(x, y)]))
    diag(Sigma.obs) <- diag(Sigma.obs) + pmax(noise, 1e-4)
    Sigma.obs.inv <- solve(Sigma.obs)

    ## Covariance of all points
    prop.x <- as.matrix(xy.dt[, .(x,y)])
    Sigma.x.x <- kernelMatrix(k.obs, prop.x)
    diag(Sigma.x.x) <- diag(Sigma.x.x)

    ## Covariance between observed and all points
    Sigma.x.obs <- kernelMatrix(k.obs, prop.x, as.matrix(xy.dt[(in.sample), .(x, y)]))

    ## Predicted mean
    pred.mean <- Sigma.x.obs %*% Sigma.obs.inv %*% xy.dt[(in.sample), link(obs.val) - mu] + mu

    ## Predicted SD
    pred.sd <- sqrt(diag(Sigma.x.x - Sigma.x.obs %*% (Sigma.obs.inv %*% t(Sigma.x.obs))))

    ## Fill everything in
    xy.dt[, curr.mean := pred.mean]
    xy.dt[, curr.sd := pred.sd]
    xy.dt[, curr.pred := inv.link(curr.mean)]
    xy.dt[, curr.lower := inv.link(curr.mean - qnorm(0.975) * curr.sd)]
    xy.dt[, curr.upper := inv.link(curr.mean + qnorm(0.975) * curr.sd)]

    ## Append to rolling data set of MSEs
    tmp.i <- i
    mse.dt <- rbind(mse.dt, data.table(iter = i,
                                       mse = xy.dt[!(in.sample), mean((true.val - curr.mean) ^ 2)],
                                       method = m_val,
                                       n.obs = n_val,
                                       n.samples = xy.dt[, sum(in.sample)]))

    ## Predict at unobserved points
    prop.x <- as.matrix(xy.dt[(!in.sample), .(x, y)])
    Sigma.x.x <- kernelMatrix(k.obs, prop.x)
    Sigma.x.obs <- kernelMatrix(k.obs, prop.x, as.matrix(xy.dt[(in.sample), .(x, y)]))
    test.mean <- Sigma.x.obs %*% Sigma.obs.inv %*% xy.dt[(in.sample), link(obs.val) - mu] + mu
    mean.diff <- test.mean - max(xy.dt[(in.sample), link(obs.val) - mu])
    pred.sd <- sqrt(diag(Sigma.x.x - Sigma.x.obs %*% (Sigma.obs.inv %*% t(Sigma.x.obs))))

    ## Calculate acquisition at proposal points
    if (m_val == 'EI') {
      prop.a <- t(pmax(0, mean.diff) + dnorm(mean.diff/pred.sd) - abs(mean.diff) * pnorm(mean.diff/pred.sd))
    } else if (m_val == 'PI') {
      prop.a <- pnorm((mean.diff-pi.xi) / pred.sd)
    } else if (m_val == 'UCB') {
      b <- 2 * log(nrow(xy.dt) * i^2 * pi ^ 2 / 6 * ucb.delta)
      prop.a <- test.mean + b * pred.sd
    } else if (m_val == 'SRS') {
      prop.a <- rep(1/length(test.mean), length(test.mean))
    } else if (m_val == 'VAR') {
      prop.a <- pred.sd^2
    }

    ## Select new points
    xy.dt[!(in.sample), a.val := as.vector(prop.a)]
    xy.dt[is.na(a.val), a.val := 0]
    plot.dt <- rbind(plot.dt, xy.dt[, .(method = m_val, n.obs = n_val, x, y, curr.pred, iter = tmp.i, is.new, curr.lower, curr.upper, obs.val, in.sample, curr.sd, a.val)])
    if (xy.dt[!(in.sample), sum(a.val > 0)] <= n.obs) {
      break()
    }
    if (r_val == 0 & m_val != 'SRS') {
      new.x <- xy.dt[!(in.sample)][order(a.val, decreasing = T), i][1:n.obs]
    } else if (r_val == 1 | m_val == 'SRS') {
      new.x <- xy.dt[!(in.sample)][, sample(i, n.obs, prob = a.val)]
    }
    best.a <- xy.dt[i %in% new.x, a.val]
    curr.max <- xy.dt[!(in.sample), max(a.val)]
    curr.sum <- xy.dt[!(in.sample), sum(a.val)]

    print(sprintf('Iter %i (%s): %0.3f (%i samples)', i,m_val, max(best.a), n.obs))

    ## Take new points out of circulation
    xy.dt[i %in% new.x, in.sample := T]
    xy.dt[, is.new := i %in% new.x]

    if (i < n.iter) {
      xy.dt[, a.val := NULL]
    }
    i <- i + 1
  }
}
ggplot(mse.dt, aes(x=iter, y=mse, color=method)) + geom_line()

ggplot(plot.dt[iter==max(iter)], aes(x=x, y=y, fill = curr.pred)) +
  geom_tile() +
  geom_point(data = plot.dt[(in.sample)], color='red', shape=3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_equal() +
  scale_fill_viridis_c(na.value = NA) +
  facet_wrap('method', nrow=1) +
  theme(legend.position = 'bottom')

if (make.gif) {
  anim <- ggplot(plot.dt, aes(x=x, y=y, fill = curr.pred)) +
    geom_tile() +
    geom_tile(data = plot.dt[(is.new)], aes(fill=NaN), color='red', size=1) +
    geom_point(data = plot.dt[(in.sample)], color='red', shape=3) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_equal() +
    scale_fill_viridis_c(na.value = NA) +
    transition_states(iter, transition_length = 0) +
    facet_wrap('method', nrow=1) +
    labs(title = 'Iteration: {closest_state}') +
    theme(legend.position = 'bottom')
  animate(anim, duration = gif.time, fps = ceiling(plot.dt[,max(iter)]/gif.time),
          width=600, height=500)
}
