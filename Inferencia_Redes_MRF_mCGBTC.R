###############################################################################
#  OBJETIVO 8  |  f.1.8  Inferencia de redes probabilísticas a partir de
#  registros neuronales y su relación con la organización estructural del
#  circuito motor corteza-ganglios basales-tálamo-corteza (mCGBTC)
#
#  Propuesta de implementación en R.  Uso:  source("objetivo8_redes_MRF_mCGBTC.R",
#  encoding = "UTF-8");  demo_objetivo8()   # corrida de prueba con datos simulados
#
#  FLUJO (corresponde párrafo por párrafo a f.1.8):
#   1. Binning de espigas -> matriz binaria (ensayos x bins) por sesión.
#   2. MRF binario de pares (modelo de Ising) por sesión, estimado con regresión
#      logística L1 nodo-a-nodo (Ravikumar et al. 2010), selección de lambda por
#      EBIC y estabilidad por submuestreo de ensayos.  El PSTH por condición entra
#      como offset: las aristas son acoplamiento funcional más allá de la
#      modulación común por el estímulo/condición.
#   3. MRF bipartito temporal (X[t-lag] -> X[t]) = Fig. 8 (capas t-1 y t).
#   4. Ensambles funcionales (módulos) y su asociación con razón, duración total
#      y tarea.
#   5. Topología de las redes por sesión y tendencia a lo largo de las 10 etapas
#      de aprendizaje (prueba de permutación).
#   6. Grafos estructurales (información mutua entre parcelas/vértices y/o
#      densidad de tractos) y comparación con los grafos funcionales:
#      similitud de grafos (Mantel), alineamiento espectral, correspondencia
#      multiescala, co-localización de cambios y rezago estructura->función.
#
#  SUPUESTO CLAVE: las neuronas de sesiones distintas NO son las mismas (registro
#  semicrónico), así que las redes neuronales no se comparan arista a arista entre
#  etapas.  Se comparan (a) métricas topológicas y (b) matrices acoplamiento
#  AGREGADAS POR PARCELA, que sí están alineadas entre sesiones y con la
#  tractografía.  Cada neurona debe asignarse a una parcela (coordenadas de f.1.5).
###############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(glmnet)
  library(igraph)
  library(parallel)
})

# ---- 0. Configuración --------------------------------------------------------
cfg_default <- function(...) {
  cfg <- list(
    bw          = 20,      # ancho de bin (ms)
    min_rate    = 0.5,     # tasa mínima (Hz) para incluir una unidad
    min_events  = 30,      # mínimo de bins con espiga (y sin espiga)
    cond_vars   = c("task", "ratio_class"),  # columnas de trl que definen el PSTH-offset
    psth_sd     = 2,       # suavizado gaussiano del PSTH (en bins)
    gamma       = 0.25,    # hiperparámetro EBIC (0 = BIC; 0.5 = muy conservador)
    nlambda     = 30,
    rule        = "AND",   # regla de simetrización nodo-a-nodo: "AND" | "OR"
    B           = 20,      # submuestreos para selección por estabilidad
    frac        = 0.5,     # fracción de ensayos por submuestra
    thr         = 0.7,     # frecuencia mínima de selección de una arista
    lag         = 1,       # rezago (bins) del MRF bipartito temporal
    min_ens     = 4,       # tamaño mínimo de ensamble
    sw_nrand    = 10,      # grafos aleatorios para el índice small-world
    nperm       = 999,     # permutaciones (tendencias y comparaciones)
    epochs      = list(ensayo = c(0, Inf)),  # ventanas (ms desde t_start) p/ ensambles
    macro_pair  = c("AMS", "Putamen"),       # par de áreas para el análisis de rezago
    ncores      = 1,       # >1 usa mclapply (no disponible en Windows)
    outdir      = "resultados_obj8",
    seed        = 123
  )
  modifyList(cfg, list(...))
}

###############################################################################
# FORMATO DE ENTRADA
#   data   : lista de sesiones; cada sesión es una lista con
#              $session (id), $animal, $stage (1..10),
#              $neu : data.frame(neuron_id, area, parcel)
#              $spk : data.frame(trial, neuron_id, t_ms)  # t relativo al inicio de sesión de ensayo
#              $trl : data.frame(trial, task, ratio, total_dur, iter, t_start, t_end,
#                                ratio_class)      # t_start/t_end en ms
#   struct : lista con
#              $parcels : data.frame(parcel, area)
#              $stages  : vector de etapas (p. ej. 1:10)
#              $feat    : lista (una por etapa) de matrices parcelas x rasgos
#                         (grosor, T1w/T2w, FA, NDI, densidad de fibras...) -> grafo MI
#              $tract   : (opcional) lista (una por etapa) de matrices parcelas x parcelas
#                         (densidad de fibras / streamlines, p. ej. de MRtrix3)
###############################################################################

# ---- 1. Datos simulados (para probar el flujo completo) ------------------------
simulate_objetivo8 <- function(n_stage = 10, n_sess = 1, n_per_parcel = 5,
                               n_trials = 60, bw = 20, T_ms = 1600, seed = 1) {
  set.seed(seed)
  parcels <- data.frame(
    parcel = c(paste0("AMS_", 1:3), paste0("PUT_", 1:3), paste0("TAL_", 1:3)),
    area   = rep(c("AMS", "Putamen", "Talamo"), each = 3), stringsAsFactors = FALSE)
  nb <- floor(T_ms / bw)
  sg <- function(s) plogis((s - n_stage / 2 - 0.5) / 1.2)          # crecimiento sigmoide con el aprendizaje
  g_fun <- function(s) (sg(s) - sg(1)) / (sg(n_stage) - sg(1))
  sessions <- list(); sid <- 0
  for (s in seq_len(n_stage)) for (k in seq_len(n_sess)) {
    sid <- sid + 1
    neu <- parcels[rep(seq_len(nrow(parcels)), each = n_per_parcel), ]
    neu$neuron_id <- sprintf("s%02d_n%03d", sid, seq_len(nrow(neu)))
    rownames(neu) <- NULL
    p <- nrow(neu)
    trl <- data.frame(trial = seq_len(n_trials), task = sample(1:3, n_trials, TRUE),
                      ratio = 2^runif(n_trials, -2, 2),
                      total_dur = sample(c(1300, 1600, 1900), n_trials, TRUE),
                      iter = NA_integer_, t_start = 0, t_end = T_ms)
    trl$ratio[trl$task == 1] <- 1
    trl$ratio_class <- ifelse(abs(log2(trl$ratio)) < 0.6, "simple", "compleja")
    ia <- which(neu$area == "AMS"); ip <- which(neu$area == "Putamen"); it <- which(neu$area == "Talamo")
    load <- matrix(0, p, 4)
    load[ia[1:5], 1] <- 1; load[ip[1:5], 2] <- 1; load[it[1:7], 3] <- 1
    load[c(ia[6:15], ip[6:15]), 4] <- 0.2 + 1.3 * g_fun(s)   # ensamble AMS-putamen (parcelas 2-3) crece con el aprendizaje
    base <- rnorm(p, -1.6, 0.3); amp <- runif(p, 0.2, 0.8); ph <- runif(p, 0, 2 * pi)
    tt <- (seq_len(nb) - 0.5) / nb
    rows <- list()
    for (tr in seq_len(n_trials)) {
      cond_gain <- if (trl$ratio_class[tr] == "compleja") 1.5 else 1
      L <- t(apply(matrix(rnorm(4 * nb), 4), 1, function(e) {
        for (i in 2:nb) e[i] <- 0.8 * e[i - 1] + sqrt(1 - 0.64) * e[i]; e }))
      eta <- base + sapply(tt, function(t) amp * cond_gain * sin(2 * pi * 2 * t + ph)) +
        0.8 * (load %*% L)
      hit <- which(matrix(rbinom(p * nb, 1, plogis(eta)), p, nb) == 1, arr.ind = TRUE)
      if (nrow(hit))
        rows[[tr]] <- data.frame(trial = tr, neuron_id = neu$neuron_id[hit[, 1]],
                                 t_ms = (hit[, 2] - 1 + runif(nrow(hit))) * bw)
    }
    sessions[[sid]] <- list(session = sid, animal = "M1", stage = s, neu = neu,
                            spk = do.call(rbind, rows), trl = trl)
  }
  # estructura: la conectividad AMS-putamen "se adelanta" una etapa a la función
  q <- 12; z <- rnorm(q); mu <- matrix(rnorm(3 * q), 3); rownames(mu) <- c("AMS", "Putamen", "Talamo")
  feat <- tract <- vector("list", n_stage)
  for (s in seq_len(n_stage)) {
    gs <- g_fun(min(s + 1, n_stage))
    w <- ifelse(parcels$parcel %in% c("AMS_2", "AMS_3", "PUT_2", "PUT_3"), 1.5 * gs, 0)
    f <- mu[parcels$area, ] + outer(w, z) + matrix(rnorm(9 * q, 0, 0.6), 9)
    dimnames(f) <- list(parcels$parcel, paste0("rasgo", seq_len(q)))
    feat[[s]] <- f
    tr <- matrix(0.1, 9, 9, dimnames = list(parcels$parcel, parcels$parcel))
    for (k in 2:3) tr[paste0("AMS_", k), paste0("PUT_", k)] <- 0.2 + 0.8 * gs
    tr <- pmax(tr + t(tr) + matrix(rnorm(81, 0, 0.02), 9), 0); diag(tr) <- 0
    tract[[s]] <- tr
  }
  list(data = sessions, struct = list(parcels = parcels, stages = seq_len(n_stage),
                                      feat = feat, tract = tract))
}

# ---- 2. Preparación: binning, PSTH-offset -------------------------------------
gauss_smooth <- function(x, sd) {
  if (sd <= 0 || length(x) < 3) return(x)
  h <- ceiling(3 * sd); k <- dnorm(-h:h, 0, sd); k <- k / sum(k)
  xp <- c(rep(x[1], h), x, rep(x[length(x)], h))
  as.numeric(stats::filter(xp, k, sides = 2))[(h + 1):(h + length(x))]
}

prepare_session <- function(sess, cfg) {
  bw <- cfg$bw; trl <- sess$trl; neu <- sess$neu; spk <- sess$spk
  ntr <- nrow(trl); nb <- floor(max(trl$t_end - trl$t_start) / bw)
  ti <- match(spk$trial, trl$trial); ni <- match(spk$neuron_id, neu$neuron_id)
  tb <- floor((spk$t_ms - trl$t_start[ti]) / bw) + 1
  ok <- !is.na(ti) & !is.na(ni) & !is.na(tb) & tb >= 1 & tb <= nb
  X <- sparseMatrix(i = (ti[ok] - 1) * nb + tb[ok], j = ni[ok], x = 1,
                    dims = c(ntr * nb, nrow(neu)))
  X@x <- pmin(X@x, 1)                                  # binariza (espiga/no espiga por bin)
  n_ev <- colSums(X); rate <- n_ev / (ntr * nb * bw / 1000)
  keep <- rate >= cfg$min_rate & n_ev >= cfg$min_events & (nrow(X) - n_ev) >= cfg$min_events
  X <- X[, keep, drop = FALSE]; neu <- neu[keep, , drop = FALSE]
  trial <- rep(seq_len(ntr), each = nb); tbin <- rep(seq_len(nb), times = ntr)
  # condiciones para el PSTH (offset)
  cv <- intersect(cfg$cond_vars, names(trl))
  key <- if (length(cv)) do.call(paste, c(trl[cv], sep = "|")) else rep("all", ntr)
  cid <- as.integer(factor(key)); ncond <- max(cid)
  gidx <- (cid[trial] - 1) * nb + tbin
  G <- sparseMatrix(i = seq_along(gidx), j = gidx, x = 1, dims = c(length(gidx), ncond * nb))
  den <- pmax(colSums(G), 1)
  psth <- as.matrix(crossprod(G, X)) / den           # (cond*bin) x neuronas
  for (cc in seq_len(ncond)) {
    r <- ((cc - 1) * nb + 1):(cc * nb)
    psth[r, ] <- apply(psth[r, , drop = FALSE], 2, gauss_smooth, sd = cfg$psth_sd)
  }
  offs <- qlogis(pmin(pmax(psth, 0.002), 0.998))
  list(X = X, neu = neu, trial = trial, tbin = tbin, gidx = gidx, offs = offs,
       nb = nb, ntr = ntr, trl = trl, bw = bw)
}

# ---- 3. Inferencia del MRF (Ising) nodo-a-nodo ------------------------------------
# Para cada neurona j:  logit P(x_j=1 | resto) = offset_j(cond,t) + b_j + sum_k theta_jk x_k
# (lasso, lambda por EBIC).  Los theta_jk son los parámetros de pares del MRF.
fit_node_glmnet <- function(y, Z, off, gamma, nlambda) {
  if (sum(y) == 0 || sum(y) == length(y)) return(rep(0, ncol(Z)))
  fit <- tryCatch(glmnet(Z, y, family = "binomial", offset = off, nlambda = nlambda,
                         lambda.min.ratio = 0.01, standardize = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(0, ncol(Z)))
  n <- length(y)
  dev <- fit$nulldev * (1 - fit$dev.ratio)               # = -2 log-verosimilitud (saturado = 0)
  ebic <- dev + fit$df * log(n) + 2 * gamma * fit$df * log(ncol(Z))
  as.numeric(fit$beta[, which.min(ebic)])
}

par_lapply <- function(X, FUN, cfg) {
  if (cfg$ncores > 1) mclapply(X, FUN, mc.cores = cfg$ncores) else lapply(X, FUN)
}

# MRF contemporáneo (no dirigido).  'rows' permite submuestrear (estabilidad).
infer_contemp <- function(prep, cfg, rows = NULL) {
  X <- prep$X; g <- prep$gidx
  if (!is.null(rows)) { X <- X[rows, , drop = FALSE]; g <- g[rows] }
  p <- ncol(X)
  co <- par_lapply(seq_len(p), function(j)
    fit_node_glmnet(as.numeric(X[, j]), X[, -j, drop = FALSE], prep$offs[g, j],
                    cfg$gamma, cfg$nlambda), cfg)
  Bm <- matrix(0, p, p)
  for (j in seq_len(p)) Bm[j, -j] <- co[[j]]
  sym_rule(Bm, cfg$rule)
}

sym_rule <- function(Bm, rule = "AND") {
  W <- (Bm + t(Bm)) / 2
  keep <- if (rule == "AND") (Bm != 0 & t(Bm) != 0) else (Bm != 0 | t(Bm) != 0)
  W * keep
}

# MRF bipartito temporal (Fig. 8): capa X[t-lag] (pre) -> capa X[t] (post).  W[i, j]: i -> j
infer_lagged <- function(prep, cfg, rows = NULL) {
  lag <- cfg$lag
  idx <- which(prep$tbin > lag)                       # el pasado existe dentro del mismo ensayo
  if (!is.null(rows)) idx <- intersect(idx, rows)
  Xr <- prep$X[idx, , drop = FALSE]; Xp <- prep$X[idx - lag, , drop = FALSE]; g <- prep$gidx[idx]
  p <- ncol(Xr)
  co <- par_lapply(seq_len(p), function(j)
    fit_node_glmnet(as.numeric(Xr[, j]), Xp, prep$offs[g, j], cfg$gamma, cfg$nlambda), cfg)
  W <- do.call(cbind, co)                              # columna j = coeficientes sobre j
  diag(W) <- 0                                         # se descarta autohistoria (refractariedad)
  W
}

# Selección por estabilidad: frecuencia de cada arista en submuestras de ensayos
stability_select <- function(prep, cfg, fun) {
  if (cfg$B <= 0) return(NULL)
  p <- ncol(prep$X); freq <- matrix(0, p, p)
  for (b in seq_len(cfg$B)) {
    tr <- sample.int(prep$ntr, max(2, floor(cfg$frac * prep$ntr)))
    Wb <- fun(prep, cfg, which(prep$trial %in% tr))
    freq <- freq + (Wb != 0)
  }
  freq / cfg$B
}

# ---- 4. Topología y ensambles --------------------------------------------------------
net_metrics <- function(W, area, n_rand = 10) {
  A <- abs(W); diag(A) <- 0; p <- nrow(A)
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  m <- ecount(g)
  empty <- data.frame(n_nodes = p, n_edges = m, density = m / (p * (p - 1) / 2), frac_pos = NA,
                      mean_degree = 0, clustering = NA, eff_weighted = 0, eff_binary = 0,
                      modularity = NA, n_modules = NA, frac_interarea = NA,
                      sigma_sw = NA, economy = NA)
  if (m < 3) return(empty)
  ed <- as_edgelist(g, names = FALSE)
  D <- distances(g, weights = 1 / E(g)$weight); iD <- 1 / D; diag(iD) <- 0
  Db <- distances(g); iDb <- 1 / Db; diag(iDb) <- 0
  eff_w <- mean(iD[upper.tri(iD)]); eff_b <- mean(iDb[upper.tri(iDb)])
  lc <- cluster_louvain(g)
  C <- transitivity(g, "global"); L <- mean_distance(g, directed = FALSE, unconnected = TRUE)
  rr <- replicate(n_rand, {
    gr <- rewire(g, keeping_degseq(niter = 10 * m))
    c(transitivity(gr, "global"), mean_distance(gr, directed = FALSE, unconnected = TRUE)) })
  sigma <- (C / mean(rr[1, ], na.rm = TRUE)) / (L / mean(rr[2, ], na.rm = TRUE))
  dens <- m / (p * (p - 1) / 2)
  data.frame(n_nodes = p, n_edges = m, density = dens,
             frac_pos = mean(W[upper.tri(W)][W[upper.tri(W)] != 0] > 0),
             mean_degree = 2 * m / p, clustering = C, eff_weighted = eff_w, eff_binary = eff_b,
             modularity = modularity(lc), n_modules = length(lc),
             frac_interarea = mean(area[ed[, 1]] != area[ed[, 2]]),
             sigma_sw = sigma, economy = eff_b / dens)
}

# Acoplamiento medio |w| entre grupos (áreas o parcelas); sirve para redes dirigidas
group_coupling <- function(W, grp, levels) {
  A <- abs(W); diag(A) <- 0
  out <- matrix(NA_real_, length(levels), length(levels), dimnames = list(levels, levels))
  for (a in seq_along(levels)) for (b in seq_along(levels)) {
    ia <- which(grp == levels[a]); ib <- which(grp == levels[b])
    if (!length(ia) || !length(ib)) next
    blk <- A[ia, ib, drop = FALSE]
    if (a == b) { if (length(ia) < 2) next; out[a, b] <- sum(blk) / (length(ia) * (length(ia) - 1)) }
    else out[a, b] <- mean(blk)
  }
  out
}

# Ensambles funcionales = módulos de la red de coactivación (pesos positivos)
get_ensembles <- function(W, min_size = 4) {
  A <- pmax(W, 0); diag(A) <- 0
  memb <- rep(0L, nrow(A))
  if (sum(A > 0) < 2) return(memb)
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  cl <- membership(cluster_louvain(g))
  tab <- table(cl); big <- as.integer(names(tab)[tab >= min_size])
  for (k in seq_along(big)) memb[cl == big[k]] <- k
  memb
}

ensemble_tests <- function(prep, memb, cfg) {
  nens <- max(memb); if (nens == 0) return(NULL)
  trl <- prep$trl; out <- list()
  for (e in seq_len(nens)) for (ep in names(cfg$epochs)) {
    w <- cfg$epochs[[ep]]
    tms <- (seq_len(prep$nb) - 0.5) * prep$bw
    sel <- which(tms >= w[1] & tms <= w[2])
    act <- Matrix::rowMeans(prep$X[, memb == e, drop = FALSE])
    a <- tapply(act[prep$tbin %in% sel], prep$trial[prep$tbin %in% sel], mean)
    d <- data.frame(act = as.numeric(a), task = factor(trl$task), lr = log2(trl$ratio),
                    tot = trl$total_dur)
    an <- tryCatch(anova(lm(act ~ task + lr + I(lr^2) + tot, data = d)), error = function(e) NULL)
    if (is.null(an)) next
    an <- an[rownames(an) != "Residuals", ]
    out[[length(out) + 1]] <- data.frame(ensemble = e, epoch = ep, term = rownames(an), F = an$`F value`,
                                         p = an$`Pr(>F)`)
  }
  res <- do.call(rbind, out)
  if (!is.null(res)) res$p_adj <- p.adjust(res$p, "BH")
  res
}

# ---- 5. Pipeline por sesión ------------------------------------------------------------
session_pipeline <- function(sess, cfg, parcels, areas) {
  prep <- prepare_session(sess, cfg); p <- ncol(prep$X)
  W_full <- infer_contemp(prep, cfg)
  fr <- stability_select(prep, cfg, infer_contemp)
  W <- if (is.null(fr)) W_full else W_full * (fr >= cfg$thr)
  Wl_full <- infer_lagged(prep, cfg)
  frl <- stability_select(prep, cfg, infer_lagged)
  Wl <- if (is.null(frl)) Wl_full else Wl_full * (frl >= cfg$thr)
  area <- prep$neu$area; par <- prep$neu$parcel
  met <- net_metrics(W, area, cfg$sw_nrand)
  memb <- get_ensembles(W, cfg$min_ens)
  ens_multi <- if (max(memb) > 0)
    mean(sapply(seq_len(max(memb)), function(e) length(unique(area[memb == e])) > 1)) else NA
  cpl <- group_coupling(W, area, areas); cpl_lag <- group_coupling(Wl, area, areas)
  flat <- function(M, pre, dir = FALSE) {
    idx <- which(if (dir) !is.na(M) else upper.tri(M, diag = TRUE), arr.ind = TRUE)
    v <- M[idx]; names(v) <- paste0(pre, rownames(M)[idx[, 1]], if (dir) "_to_" else "_", colnames(M)[idx[, 2]]); v }
  extra <- as.data.frame(t(c(flat(cpl, "cpl_"), flat(cpl_lag, "lag_", TRUE))))
  metrics <- cbind(data.frame(session = sess$session, animal = sess$animal, stage = sess$stage),
                   met, frac_ens_multiarea = ens_multi, n_ensembles = max(memb), extra)
  list(metrics = metrics, W = W, W_lag = Wl, freq = fr, neu = prep$neu, memb = memb,
       ens_tests = ensemble_tests(prep, memb, cfg),
       pc = group_coupling(W, par, parcels$parcel), pc_lag = group_coupling(Wl, par, parcels$parcel),
       session = sess$session, stage = sess$stage, animal = sess$animal)
}

# ---- 6. Tendencia a lo largo de las etapas de aprendizaje -------------------------------
trend_tests <- function(df, vars, nperm = 999) {
  res <- lapply(vars, function(v) {
    d <- df[is.finite(df[[v]]), ]                         # descarta NA/NaN/Inf
    if (nrow(d) < 5 || length(unique(d$stage)) < 3 || !isTRUE(sd(d[[v]]) > 0)) return(NULL)
    y <- d[[v]]; st <- as.numeric(d$stage); an <- factor(d$animal)
    slope <- function(st) {
      X <- if (nlevels(an) > 1) model.matrix(~ st + an) else model.matrix(~ st)
      qr.coef(qr(X), y)[2] }
    obs <- slope(st)
    null <- replicate(nperm, slope(ave(st, an, FUN = function(x) x[sample.int(length(x))])))
    data.frame(var = v, slope_per_stage = unname(obs), kendall_tau = cor(y, st, method = "kendall"),
               p_perm = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm))
  })
  res <- do.call(rbind, res)
  if (!is.null(res)) res$p_adj <- p.adjust(res$p_perm, "BH")
  res
}

# ---- 7. Grafos estructurales ---------------------------------------------------------------
# Información mutua gaussiana (cópula) entre perfiles de rasgos de dos parcelas/vértices:
#   MI = -0.5 log(1 - r^2), r = correlación de los perfiles tras normalizar cada rasgo.
# Supuesto: los rasgos (grosor, T1w/T2w, FA, NDI, profundidades cortico-subcorticales...) son
# muestras intercambiables del perfil de cada nodo; usar >= 8-10 rasgos para MI estable.
struct_mi_graph <- function(feat) {
  Z <- apply(feat, 2, function(x) qnorm((rank(x) - 0.5) / length(x)))   # parcelas x rasgos
  R <- cor(t(Z)); MI <- -0.5 * log(1 - pmin(R^2, 0.999))
  diag(MI) <- 0; dimnames(MI) <- list(rownames(feat), rownames(feat)); MI
}

# ---- 8. Comparación estructura <-> función ----------------------------------------------------
norm_laplacian <- function(A) {
  A <- abs(A); diag(A) <- 0; d <- rowSums(A); d[d == 0] <- 1
  Dm <- diag(1 / sqrt(d)); diag(nrow(A)) - Dm %*% A %*% Dm
}

mantel_sp <- function(A, B, nperm = 999) {          # unilateral: correspondencia positiva
  ut <- upper.tri(A); n <- nrow(A)
  r0 <- suppressWarnings(cor(A[ut], B[ut], method = "spearman"))
  if (is.na(r0)) return(list(rho = NA_real_, p = NA_real_))
  null <- replicate(nperm, { pm <- sample.int(n); suppressWarnings(cor(A[ut], B[pm, pm][ut], method = "spearman")) })
  list(rho = r0, p = (1 + sum(null >= r0, na.rm = TRUE)) / (1 + nperm))
}

# Alineamiento espectral: distancia entre los k primeros autovalores no triviales del
# Laplaciano normalizado y coseno medio de los ángulos principales entre subespacios.
spectral_alignment <- function(A, B, k = 3, nperm = 999) {
  sp <- function(M) { e <- eigen(norm_laplacian(M), symmetric = TRUE); o <- order(e$values)
                      list(val = e$values[o], vec = e$vectors[, o, drop = FALSE]) }
  sa <- sp(A); sb <- sp(B); idx <- 2:(k + 1)
  al <- function(Ua, Ub) mean(svd(crossprod(Ua, Ub))$d)
  a0 <- al(sa$vec[, idx], sb$vec[, idx])
  null <- replicate(nperm, al(sa$vec[, idx], sb$vec[sample.int(nrow(B)), idx]))
  list(spec_dist = sqrt(sum((sa$val[idx] - sb$val[idx])^2)), align = a0,
       p_align = (1 + sum(null >= a0)) / (1 + nperm))
}

sym_na <- function(M) { S <- (M + t(M)) / 2; S[is.na(S)] <- M[is.na(S)]; S[is.na(S)] <- t(M)[is.na(S)]; S }

coarse <- function(M, grp) {
  ug <- sort(unique(grp)); K <- length(ug); out <- matrix(0, K, K)
  for (a in seq_len(K)) for (b in seq_len(K)) if (a != b)
    out[a, b] <- mean(M[grp == ug[a], grp == ug[b]])
  out
}

# Correspondencia multiescala: se agrupan las parcelas con el dendrograma del grafo estructural
# y se compara función-estructura a cada escala de agregación (k grupos).
multiscale_corr <- function(Fm, Sm, kmin = 5, nperm = 499) {
  n <- nrow(Sm); hc <- hclust(as.dist(max(Sm) - Sm), method = "average")
  do.call(rbind, lapply(n:min(kmin, n), function(k) {
    grp <- if (k == n) seq_len(n) else cutree(hc, k)
    mt <- mantel_sp(coarse(Fm, grp), coarse(Sm, grp), nperm)
    data.frame(k = k, rho = mt$rho, p = mt$p)
  }))
}

edge_slopes <- function(mats, stages) {
  arr <- simplify2array(mats); n <- dim(arr)[1]
  sl <- matrix(NA_real_, n, n, dimnames = dimnames(mats[[1]]))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    y <- arr[i, j, ]; ok <- !is.na(y)
    if (sum(ok) >= 3) sl[i, j] <- if (sd(y[ok]) > 0) cov(stages[ok], y[ok]) / var(stages[ok]) else 0
  }
  sl
}

# x[t] vs y[t+lag]: lag > 0 significa que x (estructura) precede a y (función).
# Se usan PRIMERAS DIFERENCIAS (cambio etapa a etapa): dos series con tendencia monótona
# correlacionan a cualquier rezago y el nulo por desplazamiento circular no lo corrige.
# Con 10 etapas (9 diferencias) la potencia es baja: resultado exploratorio.
lagged_crosscor <- function(x, y, lags = -2:2, nperm = 2000, use_diff = TRUE) {
  if (use_diff) { x <- diff(x); y <- diff(y) }
  n <- length(x)
  cc <- function(x, y, l) suppressWarnings(
    if (l >= 0) cor(x[1:(n - l)], y[(1 + l):n]) else cor(x[(1 - l):n], y[1:(n + l)]))
  obs <- sapply(lags, function(l) cc(x, y, l))
  nm <- replicate(nperm, { sh <- sample.int(n - 1, 1); ys <- y[c((sh + 1):n, 1:sh)]
                           max(abs(sapply(lags, function(l) cc(x, ys, l))), na.rm = TRUE) })
  data.frame(lag = lags, r = obs, p_maxabs = sapply(abs(obs), function(a) (1 + sum(nm >= a)) / (1 + nperm)))
}

block_mean <- function(M, grp, a, b) {
  ia <- which(grp == a); ib <- which(grp == b); v <- M[ia, ib]
  if (a == b) v <- v[upper.tri(v)]; mean(v, na.rm = TRUE)
}

# Análisis completo para un tipo de grafo estructural (MI o tractos)
structure_function_analysis <- function(Fst, Sst, stages, parcels, cfg) {
  valid <- function(M) { d <- !is.na(M); diag(d) <- FALSE; colSums(d) > 0 }
  keep <- Reduce(intersect, lapply(Fst, function(M) names(which(valid(M)))))
  if (length(keep) < 5) stop("Menos de 5 parcelas muestreadas en todas las etapas: ampliar parcelación.")
  Fk <- lapply(Fst, function(M) { M <- sym_na(M)[keep, keep]; diag(M) <- 0; M })
  Sk <- lapply(Sst, function(M) { M <- M[keep, keep]; diag(M) <- 0; M })
  ke <- min(3, length(keep) - 2)
  per_stage <- do.call(rbind, lapply(seq_along(stages), function(i) {
    mt <- mantel_sp(Fk[[i]], Sk[[i]], cfg$nperm)
    sa <- spectral_alignment(Fk[[i]], Sk[[i]], ke, cfg$nperm)
    data.frame(stage = stages[i], mantel_rho = mt$rho, mantel_p = mt$p,
               spec_dist = sa$spec_dist, spec_align = sa$align, spec_p = sa$p_align)
  }))
  trends <- do.call(rbind, lapply(c("mantel_rho", "spec_align", "spec_dist"), function(v) {
    ct <- suppressWarnings(cor.test(per_stage$stage, per_stage[[v]], method = "kendall"))
    data.frame(var = v, kendall_tau = unname(ct$estimate), p = ct$p.value) }))
  multi <- do.call(rbind, lapply(seq_along(stages), function(i)
    cbind(stage = stages[i], multiscale_corr(Fk[[i]], Sk[[i]], 5, 499))))
  # co-localización del cambio: ¿las aristas que más se refuerzan estructuralmente son las que
  # más aumentan su acoplamiento funcional durante el aprendizaje?
  sF <- edge_slopes(Fk, stages); sS <- edge_slopes(Sk, stages)
  co <- mantel_sp(sF, sS, cfg$nperm)
  # rezago estructura -> función para el par macro (p. ej. AMS-putamen por la cápsula externa)
  area <- setNames(parcels$area, parcels$parcel)[keep]; pr <- cfg$macro_pair
  lag <- NULL
  if (all(pr %in% area)) {
    xs <- sapply(Sk, block_mean, grp = area, a = pr[1], b = pr[2])
    ys <- sapply(Fk, block_mean, grp = area, a = pr[1], b = pr[2])
    lag <- lagged_crosscor(xs, ys, nperm = cfg$nperm)
  }
  list(per_stage = per_stage, trends = trends, multiscale = multi,
       change_colocalization = data.frame(rho = co$rho, p = co$p),
       slopes_function = sF, slopes_structure = sS, lag_structure_to_function = lag, parcels = keep)
}

# ---- 9. Pipeline completo -----------------------------------------------------------------------
run_objetivo8 <- function(data, struct, cfg = cfg_default()) {
  set.seed(cfg$seed); dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)
  parcels <- struct$parcels; areas <- unique(parcels$area)
  nets <- lapply(seq_along(data), function(i) {
    message(sprintf("[%d/%d] sesion %s (etapa %s)", i, length(data), data[[i]]$session, data[[i]]$stage))
    session_pipeline(data[[i]], cfg, parcels, areas)
  })
  metrics <- do.call(rbind, lapply(nets, `[[`, "metrics"))
  vars <- setdiff(names(metrics)[sapply(metrics, is.numeric)], c("session", "stage", "n_nodes"))
  trends <- trend_tests(metrics, vars, cfg$nperm)
  ens <- do.call(rbind, lapply(nets, function(n) if (!is.null(n$ens_tests)) cbind(session = n$session, stage = n$stage, n$ens_tests)))
  # acoplamiento por parcela promediado por etapa (alineado entre sesiones y con la estructura)
  stages <- struct$stages
  avg <- function(key) lapply(stages, function(s) {
    ms <- lapply(nets[sapply(nets, `[[`, "stage") == s], `[[`, key)
    a <- simplify2array(ms); if (length(dim(a)) == 2) a <- array(a, c(dim(a), 1))
    out <- apply(a, c(1, 2), function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
    dimnames(out) <- dimnames(ms[[1]]); out })
  Fst <- avg("pc")
  Sst_mi <- lapply(struct$feat, struct_mi_graph)
  sf <- list(mi = structure_function_analysis(Fst, Sst_mi, stages, parcels, cfg))
  if (!is.null(struct$tract)) sf$tract <- structure_function_analysis(Fst, struct$tract, stages, parcels, cfg)
  res <- list(cfg = cfg, nets = nets, metrics = metrics, trends = trends, ensembles = ens,
              parcel_coupling_stage = Fst, structure_function = sf)
  saveRDS(res, file.path(cfg$outdir, "resultados_obj8.rds"))
  write.csv(metrics, file.path(cfg$outdir, "metricas_topologicas_por_sesion.csv"), row.names = FALSE)
  write.csv(trends, file.path(cfg$outdir, "tendencias_aprendizaje.csv"), row.names = FALSE)
  for (nm in names(sf)) {
    write.csv(sf[[nm]]$per_stage, file.path(cfg$outdir, paste0("estructura_funcion_", nm, ".csv")), row.names = FALSE)
    write.csv(sf[[nm]]$multiscale, file.path(cfg$outdir, paste0("multiescala_", nm, ".csv")), row.names = FALSE)
  }
  edges <- do.call(rbind, lapply(nets, function(n) {
    ix <- which(n$W != 0 & upper.tri(n$W), arr.ind = TRUE)
    if (!nrow(ix)) return(NULL)
    data.frame(session = n$session, stage = n$stage, from = n$neu$neuron_id[ix[, 1]], to = n$neu$neuron_id[ix[, 2]],
               area_from = n$neu$area[ix[, 1]], area_to = n$neu$area[ix[, 2]], weight = n$W[ix],
               sel_freq = if (is.null(n$freq)) NA else n$freq[ix]) }))
  write.csv(edges, file.path(cfg$outdir, "aristas_funcionales.csv"), row.names = FALSE)
  plot_results(res, cfg$outdir)
  res
}

# ---- 10. Figuras ------------------------------------------------------------------------------------
plot_results <- function(res, outdir) {
  m <- res$metrics
  vv <- intersect(c("density", "clustering", "eff_binary", "modularity", "frac_interarea",
                    "sigma_sw", "economy", "frac_ens_multiarea"), names(m))
  png(file.path(outdir, "fig_topologia_vs_etapa.png"), 1400, 900, res = 130)
  par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
  for (v in vv) { ok <- is.finite(m[[v]])
    plot(m$stage[ok], m[[v]][ok], pch = 19, xlab = "Etapa de aprendizaje", ylab = v, main = v)
    if (sum(ok) > 3) abline(lm(m[[v]][ok] ~ m$stage[ok]), col = "red", lwd = 2) }
  dev.off()
  png(file.path(outdir, "fig_estructura_funcion.png"), 1400, 500, res = 130)
  par(mfrow = c(1, length(res$structure_function)), mar = c(4, 4, 2, 1))
  for (nm in names(res$structure_function)) {
    ps <- res$structure_function[[nm]]$per_stage
    plot(ps$stage, ps$mantel_rho, type = "b", pch = 19, ylim = c(-0.2, 1), xlab = "Etapa",
         ylab = "Mantel (Spearman)", main = paste("Funcion vs estructura:", nm)); abline(h = 0, lty = 2)
  }
  dev.off()
  st <- sort(unique(m$stage)); pick <- c(min(st), max(st))
  png(file.path(outdir, "fig_redes_primera_ultima_etapa.png"), 1600, 800, res = 130)
  par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
  cols <- c(AMS = "#E69F00", Putamen = "#56B4E9", Talamo = "#009E73")
  for (s in pick) {
    n <- res$nets[[which(sapply(res$nets, `[[`, "stage") == s)[1]]]
    A <- abs(n$W); g <- graph_from_adjacency_matrix(A, "undirected", weighted = TRUE, diag = FALSE)
    ed <- as_edgelist(g, names = FALSE)
    E(g)$color <- ifelse(n$W[ed] > 0, "#3B6FB6", "#C0392B"); E(g)$width <- 0.5 + 3 * E(g)$weight / max(E(g)$weight, 1e-9)
    V(g)$color <- cols[n$neu$area]; V(g)$size <- 6
    plot(g, layout = layout_with_fr(g), vertex.label = NA, main = paste("Etapa", s))
  }
  dev.off()
}

# ---- 11. Carga de datos reales (esqueleto) -----------------------------------------------------------
# Convertir la salida de Kilosort2 (spike_times.npy, spike_clusters.npy, cluster_info.tsv) a las
# tablas spk/neu/trl (p. ej. con RcppCNPy o reticulate::import("numpy")) y guardarlas como CSV.
load_session_csv <- function(spk_csv, neu_csv, trl_csv, session, animal, stage) {
  list(session = session, animal = animal, stage = stage,
       spk = read.csv(spk_csv), neu = read.csv(neu_csv, stringsAsFactors = FALSE),
       trl = read.csv(trl_csv))
}

# ---- 12. Demo con datos simulados --------------------------------------------------------------------
demo_objetivo8 <- function(B = 10, ncores = 1, outdir = "resultados_obj8_demo") {
  sim <- simulate_objetivo8(n_stage = 10, n_sess = 1)
  res <- run_objetivo8(sim$data, sim$struct, cfg_default(B = B, ncores = ncores, outdir = outdir, nperm = 499))
  cat("\n== Tendencias topologicas durante el aprendizaje (permutacion) ==\n")
  print(head(res$trends[order(res$trends$p_perm), ], 8), digits = 3, row.names = FALSE)
  for (nm in names(res$structure_function)) {
    cat("\n== Estructura-funcion:", nm, "==\n")
    print(res$structure_function[[nm]]$trends, digits = 3, row.names = FALSE)
    cat("Co-localizacion de cambios (rho, p):",
        round(unlist(res$structure_function[[nm]]$change_colocalization), 3), "\n")
    cat("Rezago estructura->funcion (par macro):\n")
    print(res$structure_function[[nm]]$lag_structure_to_function, digits = 3, row.names = FALSE)
  }
  invisible(res)
}

if (!interactive() && sys.nframe() == 0) demo_objetivo8()
