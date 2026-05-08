# ---------- Log-signature (tensor form) up to degree 4 ----------
#' @keywords internal
.logsig_tensor_d4 <- function(sig) {
  S1 <- sig$L1; S2 <- sig$L2; S3 <- sig$L3; S4 <- sig$L4

  # degree-2
  L1 <- S1
  L2 <- S2 - 0.5 * .kron2(S1, S1)

  # helpers for explicit (i,j,k) positions in 8-vector:
  .pos3 <- function(a,b,c) switch(paste0(a,b,c),
                                  "111"=1, "112"=2, "121"=3, "122"=4,
                                  "211"=5, "212"=6, "221"=7, "222"=8)

  # degree-3 pieces built explicitly to match (ttt,ttf,tft,tff,ftt,ftf,fft,fff)
  S1S2 <- numeric(8)
  for (i in 1:2) for (j in 1:2) for (k in 1:2)
    S1S2[.pos3(i,j,k)] <- S1S2[.pos3(i,j,k)] + S1[i] * S2[j,k]

  S2S1 <- numeric(8)
  for (j in 1:2) for (k in 1:2) for (i in 1:2)
    S2S1[.pos3(j,k,i)] <- S2S1[.pos3(j,k,i)] + S2[j,k] * S1[i]

  L3 <- S3 - 0.5 * (S1S2 + S2S1) + (1/6) * .kron3(S1, S1, S1) * 2  # (1/3) = 2*(1/6)

  # degree-4 (layouts already consistent with .kron4)
  S1S3 <- as.vector(kronecker(S1, S3))
  S3S1 <- as.vector(kronecker(S3, S1))
  S1S1S2 <- as.vector(kronecker(kronecker(S1, S1), as.vector(S2)))
  S1S2S1 <- as.vector(kronecker(kronecker(S1, as.vector(S2)), S1))
  S2S1S1 <- as.vector(kronecker(as.vector(S2), kronecker(S1, S1)))
  S2S2   <- as.vector(kronecker(as.vector(S2), as.vector(S2)))
  S1S1S1S1 <- as.vector(kronecker(kronecker(S1, S1), kronecker(S1, S1)))

  L4 <- S4
  L4 <- L4 - 0.5 * (S1S3 + S3S1)
  L4 <- L4 + (1/3) * (S1S1S2 + S1S2S1 + S2S1S1)
  L4 <- L4 - 0.25 * S1S1S1S1
  L4 <- L4 + (1/12) * S2S2

  list(L1 = L1, L2 = L2, L3 = L3, L4 = L4)
}
