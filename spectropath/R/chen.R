# ---------- Signature increments for one linear segment ----------
.sig_inc_d2 <- function(v) {
  list(L1 = v, L2 = 0.5 * .kron2(v, v))
}

.sig_inc_d3 <- function(v) {
  list(L1 = v, L2 = 0.5 * .kron2(v, v), L3 = (1/6) * .kron3(v, v, v))
}

.sig_inc_d4 <- function(v) {
  list(
    L1 = v,
    L2 = 0.5 * .kron2(v, v),
    L3 = (1/6) * .kron3(v, v, v),
    L4 = (1/24) * .kron4(v, v, v, v)
  )
}

# ---------- Chen products (concatenation) ----------
.chen_prod_d2 <- function(A, B) {
  list(L1 = A$L1 + B$L1,
       L2 = A$L2 + .kron2(A$L1, B$L1) + B$L2)
}

.chen_prod_d3 <- function(A, B) {
  e <- diag(2)
  C1 <- A$L1 + B$L1
  C2 <- A$L2 + .kron2(A$L1, B$L1) + B$L2
  A2B1 <- numeric(8)
  for (i in 1:2) for (j in 1:2) {
    A2B1 <- A2B1 + A$L2[i, j] * .kron3(e[, i], e[, j], B$L1)
  }
  A1B2 <- numeric(8)
  for (i in 1:2) {
    tmp <- numeric(8)
    for (j in 1:2) for (k in 1:2) {
      tmp <- tmp + B$L2[j, k] * .kron3(e[, i], e[, j], e[, k])
    }
    A1B2 <- A1B2 + A$L1[i] * tmp
  }
  C3 <- A$L3 + A2B1 + A1B2 + B$L3
  list(L1 = C1, L2 = C2, L3 = C3)
}

.chen_prod_d4 <- function(A, B) {
  e <- diag(2)
  C1 <- A$L1 + B$L1
  C2 <- A$L2 + .kron2(A$L1, B$L1) + B$L2

  # degree-3 block
  A2B1 <- numeric(8)
  for (i in 1:2) for (j in 1:2)
    A2B1 <- A2B1 + A$L2[i, j] * .kron3(e[, i], e[, j], B$L1)

  A1B2 <- numeric(8)
  for (i in 1:2) {
    tmp <- numeric(8)
    for (j in 1:2) for (k in 1:2)
      tmp <- tmp + B$L2[j, k] * .kron3(e[, i], e[, j], e[, k])
    A1B2 <- A1B2 + A$L1[i] * tmp
  }
  C3 <- A$L3 + A2B1 + A1B2 + B$L3
  idx3 <- rbind(
    c(1,1,1), c(1,1,2), c(1,2,1), c(1,2,2),
    c(2,1,1), c(2,1,2), c(2,2,1), c(2,2,2)
  )

  A3B1 <- numeric(16)
  for (p in 1:8) {
    i <- idx3[p,1]; j <- idx3[p,2]; k <- idx3[p,3]
    A3B1 <- A3B1 + A$L3[p] * .kron4(e[, i], e[, j], e[, k], B$L1)
  }

  A2B2 <- numeric(16)
  for (i in 1:2) for (j in 1:2) for (k in 1:2) for (l in 1:2)
    A2B2 <- A2B2 + A$L2[i, j] * B$L2[k, l] * .kron4(e[, i], e[, j], e[, k], e[, l])

  A1B3 <- numeric(16)
  for (i in 1:2) {
    tmp <- numeric(16)
    for (p in 1:8) {
      j <- idx3[p,1]; k <- idx3[p,2]; l <- idx3[p,3]
      tmp <- tmp + B$L3[p] * .kron4(e[, i], e[, j], e[, k], e[, l])
    }
    A1B3 <- A1B3 + A$L1[i] * tmp
  }

  C4 <- A$L4 + A3B1 + A2B2 + A1B3 + B$L4
  list(L1 = C1, L2 = C2, L3 = C3, L4 = C4)
}
