test_that("getting within networks works", {
  # Which node belongs to which block
  df <-
    tibble::tibble(
      node_id = c("A", "B", "C", "D", "E", "F", "G", "H"), 
      block = c(1, 1, 1, 1, 2, 2, 2, 2)
    )

  # Edgelist
  edgelist <-
    data.frame(
      source_id = c("A", "A", "A", "B"),
      target_id = c("B", "C", "E", "F")
    )

  network <- as.network(edgelist, vertices = df, directed = FALSE)
  network %v% "vertex.names" <- 1:length(network %v% "vertex.names")
  # When not all nodes are isolated in a block
  subgraphs <- get_within_networks(network =network,block = network%v% "block")
  adj_true <- matrix(c(0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 4)
  expect_equal(adj_true, network::as.matrix.network.adjacency(subgraphs)[1:4,1:4], check.attributes = FALSE)
  adj_true <- matrix(0, nrow = 4, ncol = 4)
  expect_equal(adj_true, network::as.matrix.network.adjacency(subgraphs)[1:4+4,1:4+4], check.attributes = FALSE)
})
