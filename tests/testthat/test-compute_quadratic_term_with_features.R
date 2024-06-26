test_that("computing a quadratic term with multiple features works", {
  set.seed(334)
  # Number of nodes
  N <- 6
  # Number of clusters
  K <- 3

  # Create an adjacency matrix
  edgelist <-
    tibble::tibble(
      tail = 1:N,
      head = 1:N
    ) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail < head) %>%
    dplyr::mutate(connect = as.integer(unlist(rbinom(size = 1,prob = 0.5,n = nrow(.))))) %>%
    dplyr::filter(connect == 1)

  net <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE)
  adj <- network::as.matrix.network.adjacency(net)
  adj <- as(adj, "dgCMatrix")

  # Create feature matrices
  x <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  S <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  S <- as(S, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        s_ij <- ifelse(x[i] == x[j], 1, 0)
        S[i, j] <- s_ij
      }
    }
  }

  y <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  V <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  V <- as(V, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        v_ij <- ifelse(y[i] == y[j], 1, 0)
        V[i, j] <- v_ij
      }
    }
  }

  z <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  W <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  W <- as(W, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        w_ij <- ifelse(z[i] == z[j], 1, 0)
        W[i, j] <- w_ij
      }
    }
  }

  # Create a N x K matrix whose (i, k) element represents the probability that node i belongs to block k.
  tau <-
    matrix(c(
      0.2, 0.5, 0.3,
      0.4, 0.4, 0.2,
      0.1, 0.4, 0.5,
      0.4, 0.4, 0.2,
      0.1, 0.1, 0.8,
      0.05, 0.05, 0.9
    ),
    nrow = K, ncol = N
    )
  tau <- t(tau)

  ###########################################################
  # Compute the true quadratic term in a naive way
  ###########################################################
  # Compute pi for D_ij = 1
  minPi <- 1e-4
  list_pi <- list()
  for (w in 0:1) {
    for (v in 0:1) {
      for (s in 0:1) {
        print(glue::glue("Compute pi for pi_s{s}v{v}w{w}"))
        denom <- matrix(0, nrow = K, ncol = K)
        num <- matrix(0, nrow = K, ncol = K)
        index <- s + 2 * v + 4 * w + 1
        print(index)
        for (k in 1:K) {
          for (l in 1:K) {
            for (i in 1:N) {
              for (j in 1:N) {
                if (i != j & S[i, j] == s & V[i, j] == v & W[i, j] == w) {
                  denom[k, l] <- denom[k, l] + tau[i, k] * tau[j, l]
                }
                if (i != j & adj[i, j] == 1 & S[i, j] == s & V[i, j] == v & W[i, j] == w) {
                  num[k, l] <- num[k, l] + tau[i, k] * tau[j, l]
                }
              }
            }
          }
        }
        pi <- num / denom
        # Remove extremely small elements in pi
        for (k in 1:K) {
          for (l in 1:K) {
            if (pi[k, l] < minPi) {
              pi[k, l] <- minPi
            }
          }
        }
        list_pi[[index]] <- pi
      }
    }
  }
  # Compute the quadratic term
  A_true <- matrix(0, nrow = N, ncol = K)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        # For each ij, determine which pi must be used.
        index_ij <- S[i, j] + 2 * V[i, j] + 4 * W[i, j] + 1
        pi_ij <- list_pi[[index_ij]]
        for (k in 1:K) {
          for (l in 1:K) {
            a_ij <- tau[j, l] * (log(pi_ij[k, l]*adj[i,j]+(1-pi_ij[k, l])*(1-adj[i,j])))
            A_true[i, k] <- A_true[i, k] + a_ij
          }
        }
      }
    }
  }

  A_true <- 1 - A_true / 2

  # Divide A by alpha_{ik}
  A_true <- A_true / tau

  ###########################################################
  # Compute the quadratic term using the cpp function
  ###########################################################
  list_feature_adjmat <- list(S, V, W)
  list_multiplied_feature_adjmat <- get_elementwise_multiplied_matrices_R(adj, list_feature_adjmat)
  list_multiplied_feature_adjmat <- lapply(list_multiplied_feature_adjmat, function(x) x*1)
  A <- compute_quadratic_term_with_features(N, K, list_multiplied_feature_adjmat, tau, LB = 0, verbose = 0)

  expect_equal(A, A_true, check.attributes = FALSE, tolerance = 1e-10)
})


test_that("computing a quadratic term with multiple features works and a directed network", {

  set.seed(334)
  # Number of nodes
  N <- 12
  # Number of clusters
  K <- 3
  
  # Create an adjacency matrix
  edgelist <-
    tibble::tibble(
      tail = 1:N,
      head = 1:N
    ) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail != head) %>%
    dplyr::mutate(connect = as.integer(unlist(rbinom(size = 1,n = nrow(.), prob = 0.5)))) %>%
    dplyr::filter(connect == 1)
  
  net <- network::network(edgelist, matrix.type = "edgelist", directed = TRUE)
  adj <- network::as.matrix.network.adjacency(net)
  adj <- as(adj, "dMatrix")
  # adj[,] = 0
  
  # Create feature matrices
  x <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  S <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  S <- as(S, "dMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        s_ij <- ifelse(x[i] == x[j], 1, 0)
        S[i, j] <- s_ij
      }
    }
  }
  
  y <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  V <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  V <- as(V, "dMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        v_ij <- ifelse(y[i] == y[j], 1, 0)
        V[i, j] <- v_ij
      }
    }
  }
  
  z <- as.integer(unlist(rbinom(size = 1,prob = 0.5,n = N)))
  W <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  W <- as(W, "dMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        w_ij <- ifelse(z[i] == z[j], 1, 0)
        W[i, j] <- w_ij
      }
    }
  }
  
  # Create a N x K matrix whose (i, k) element represents the probability that node i belongs to block k.
  tau <-
    matrix(c(
      0.2, 0.5, 0.3,
      0.4, 0.4, 0.2,
      0.1, 0.4, 0.5,
      0.4, 0.4, 0.2,
      0.1, 0.1, 0.8,
      0.05, 0.05, 0.9,
      0.8, 0.1, 0.1,
      0.3, 0.4, 0.3,
      0.1, 0.8, 0.1,
      0.5, 0.4, 0.1,
      0.3, 0.3, 0.4,
      0.8, 0.1, 0.1
    ),
    nrow = K, ncol = N
    )
  tau <- t(tau)
  
  ###########################################################
  # Compute the true quadratic term in a naive way
  ###########################################################
  # Compute pi for D_ij = 1
  minPi <- 1e-4
  list_pi <- list()
  
  S[,] = 0
  V[,] = 0
  W[,] = 0
  for (w in 0:1) {
    for (v in 0:1) {
      for (s in 0:1) {
        print(glue::glue("Compute pi for pi_s{s}v{v}w{w}"))
        denom <- matrix(0, nrow = K, ncol = K)
        num <- matrix(0, nrow = K, ncol = K)
        index <- s + 2 * v + 4 * w + 1
        print(index)
        for (k in 1:K) {
          for (l in 1:K) {
            for (i in 1:N) {
              for (j in 1:N) {
                if (i != j & S[i, j] == s & V[i, j] == v & W[i, j] == w) {
                  denom[k, l] <- denom[k, l] + tau[i, k] * tau[j, l]
                }
                if (i != j & adj[i, j] == 1 & S[i, j] == s & V[i, j] == v & W[i, j] == w) {
                  num[k, l] <- num[k, l] + tau[i, k] * tau[j, l]
                }
              }
            }
          }
        }
        
        
        if(sum(denom) != 0){
          pi <- num / denom
        } else {
          pi <- num
          }
        # Remove extremely small elements in pi
        for (k in 1:K) {
          for (l in 1:K) {
            if (pi[k, l] < minPi) {
              pi[k, l] <- minPi
            }
          }
        }
        list_pi[[index]] <- pi
      }
    }
  }
  
  
  # Compute the quadratic term
  A_true <- matrix(0, nrow = N, ncol = K)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        # For each ij, determine which pi must be used.
        index_ij <- S[i, j] + 2 * V[i, j] + 4 * W[i, j] + 1
        pi_ij <- list_pi[[index_ij]]
        for (k in 1:K) {
          for (l in 1:K) {
            a_ij <- tau[j, l] *((log(pi_ij[k, l]*adj[i,j]+(1-pi_ij[k, l])*(1-adj[i,j]))) + 
                                   (log(pi_ij[l, k]*adj[j,i]+(1-pi_ij[l, k])*(1-adj[j,i]))))
            A_true[i, k] <- A_true[i, k] + a_ij
          }
        }
      }
    }
  }
  
  A_true <- 1 - A_true / 2

  # Divide A by alpha_{ik}
  A_true <- A_true / tau
  
  ###########################################################
  # Compute the quadratic term using the cpp function
  ###########################################################
  list_feature_adjmat <- list(S, V, W)
  
  
  g <- as(adj, "TsparseMatrix")
  list_multiplied_feature_adjmat <- get_elementwise_multiplied_matrices_R(g, list_feature_adjmat)
  list_multiplied_feature_adjmat <- lapply(list_multiplied_feature_adjmat, FUN = function(x) x*1)
  A <- compute_quadratic_term_with_features_directed(N, K,list_multiplied_feature_adjmat = list_multiplied_feature_adjmat,
                                                     tau = tau,LB =  0, verbose = 0)
  expect_equal(A, A_true, check.attributes = FALSE, tolerance = 1e-10)
  })

