# This source file is based on the file R/pcalg.R in the CRAN source package
# pcalg_2.6-10.tar.gz (the source code for the PCAlg R package, version 2.6-10).
#
# It is licensed as GPL (>= 2).
#
# The modifications are copyright (C) Joris M. Mooij, 2020.
#
# The hope is that these modifications will end up in the upstream PCAlg package.


###############################################################################
## FCI-JCI
###############################################################################

fcijci <- function(suffStat, indepTest, alpha, labels, p,
                   skel.method = c("stable", "original", "stable.fast"),
                   type = c("normal", "anytime", "adaptive"),
                   fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                   m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                   doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                   maj.rule = FALSE, numCores = 1, contextVars = NULL,
                   selectionBias = TRUE, jci = c("0","1","123"), verbose = FALSE)
{
  ## Purpose: Perform FCI-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximum size of conditioning set
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - rules: array of length 10 wich contains TRUE or FALSE corresponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - doPdsep: compute possible dsep
  ## - biCC: TRUE or FALSE variable containing if biconnected components are
  ##         used to compute pdsep
  ## - conservative: TRUE or FALSE defining if
  ##          the v-structures after the pdsep
  ##          have to be oriented conservatively or nor
  ## - maj.rule: TRUE or FALSE variable containing if the majority rule is
  ##             used instead of the normal conservative
  ## - labels: names of the variables or nodes
  ## - type: it specifies the version of the FCI that has to be used.
  ##         Per default it is normal, the normal FCI algorithm. It can also be
  ##         anytime for the Anytime FCI and in this cas m.max must be specified;
  ##         or it can be adaptive for Adaptive Anytime FCI and in this case
  ##         m.max must not be specified.
  ## - numCores: handed to skeleton(), used for parallelization
  ## - contextVars: subset of variables that are treated as context variables
  ## - selectionBias: allow for selection bias (default: TRUE)
  ## - jci: specifies the JCI background knowledge that is used; can be either:
  ##     "0"   no JCI background knowledge (default),
  ##     "1"   JCI assumption 1 only (i.e., no system variable causes any context variable),
  ##     "123" all JCI assumptions 1, 2 and 3
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Dec 2009; update: Diego Colombo, 2012; Martin Maechler, 2013; Joris Mooij, 2020

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")

  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")

  ## Check that jci background knowledge is valid
  jci <- match.arg(jci)
  ## Set fixed edges from JCI assumption 3 if asked for
  if( jci == "123" ) {
    if( any(is.null(fixedEdges)) )
      fixedEdges<-matrix(FALSE,p,p)
    if( length(contextVars) > 0 )
      fixedEdges[contextVars,contextVars]<-TRUE
  }

  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")

  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                     NAdelete=NAdelete, m.max=m.max, numCores=numCores, verbose=verbose)
  skel@call <- cl # so that makes it into result
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL

  if (doPdsep) {
    if (verbose) cat("\nCompute PDSEP\n=============\n")
    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
                            alpha = alpha, version.unf = c(1,1),
                            maj.rule = FALSE, verbose = verbose)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- pdsep2(skel@graph, suffStat, indepTest = indepTest, p = p,
                       sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax,
                       m.max = if (type == "adaptive") max.ordSKEL else m.max,
                       pdsep.max = pdsep.max, NAdelete = NAdelete,
                       unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                       biCC = biCC, fixedEdges = fixedEdges, verbose = verbose)

    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                            verbose = verbose, version.unf = c(1, 1),
                            maj.rule = maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                                verbose = verbose, version.unf = c(2, 1),
                                maj.rule = maj.rule)
      tripleList <- nopdsep$unfTripl
      ## update the sepsets
      sepset <- nopdsep$sk@sepset
    }
  }
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <- udag2pag2(pag = G, sepset, rules = rules, unfVect = tripleList,
                   contextVars = contextVars, selectionBias = selectionBias, jci = jci, verbose = verbose)
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
} ## {fci}

## only called in fci() [by default:  doPdsep=TRUE]
pdsep2 <- function (skel, suffStat, indepTest, p, sepset, alpha, pMax, m.max = Inf,
                   pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL,
                   biCC = FALSE, fixedEdges = NULL, verbose = FALSE) ## FIXME: verbose : 2 --> qreach(verbose)
{
  ## Purpose: Compute Possible-D-SEP for each node, perform the condittional
  ##          independent tests and adapt graph accordingly
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - skel: Graph object returned by function skeleton
  ## - suffStat, indepTest: infofor the independence tests
  ## - p: number of nodes in the graph
  ## - sepset: Sepset that was used for finding the skeleton
  ## - alpha: niveau for the tests
  ## - pMax: Maximal p-values during estimation of skeleton
  ## - m.max: maximal size of the conditioning sets
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - unfVect: vector containing the unfaithful triples, used for the
  ##   conservative orientation of the v-structures
  ## - biCC: if the biconnected components have to be used
  ## - fixedEdges: Edges marked here are not changed (logical)  JORIS
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G: Updated boolean adjacency matrix
  ## - sepset: Updated sepsets
  ## - pMax: Updated pMax
  ## - allPdsep: Possible d-sep for each node [list]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  9 Dec 2009
  ## Modification: Diego Colombo; Martin Maechler; Joris Mooij (2020)

  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  ord <- 0L
  allPdsep.tmp <- vector("list", p)
  if(biCC)
    conn.comp <- lapply(biConnComp(skel), as.numeric)
  if (any(is.null(fixedEdges))) { ## MM: could be sparse          JORIS
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  if (any(G)) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    if (verbose) cat("\nCompute collider:\n")
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      for (z in allZ) {
        if (amat[x, z] == 0 &&
            !(y %in% sepset[[x]][[z]] ||
              y %in% sepset[[z]][[x]])) {

          if (length(unfVect) == 0) { ## normal version -------------------
            amat[x, y] <- amat[z, y] <- 2
            if (verbose) cat("\n",x,"*->", y, "<-*", z, "\n")
          }
          else { ## conservative version : check if x-y-z is faithful
            if (!any(unfVect == triple2numb(p,x,y,z), na.rm = TRUE) &&
                !any(unfVect == triple2numb(p,z,y,x), na.rm = TRUE)) {
              amat[x, y] <- amat[z, y] <- 2
              if (verbose)
                cat("\n",x,"*->", y, "<-*", z, "\n")
            }
          }
        }
      } ## for( z )
    } ## for( i  )
    allPdsep <- lapply(1:p, qreach, amat = amat)# verbose = (verbose >= 2)
    allPdsep.tmp <- vector("list", p)
    for(x in 1:p) {
      if(verbose) cat("\nPossible D-Sep of", x, "is:", allPdsep[[x]], "\n")
      if (any(an0 <- amat[x, ] != 0)) {
        tf1 <- setdiff(allPdsep[[x]], x)
        adj.x <- which(an0)
        for (y in adj.x)
         if( !fixedEdges[x,y] ) {
          if(verbose) cat(sprintf("\ny = %3d\n.........\n", y))
          tf <- setdiff(tf1, y)
          diff.set <- setdiff(tf, adj.x)
          ## bi-connected components
          if (biCC) {
            for(cci in conn.comp) {
              if (x %in% cci && y %in% cci)
                break ## found it
            }
            bi.conn.comp <- setdiff(cci, c(x,y))
            tf <- intersect(tf, bi.conn.comp)
            if (verbose) {
              cat("There is an edge between",x,"and",y,"\n")
              cat("Possible D-Sep of", x,
                  "intersected with the biconnected component of",x,"and",y,
                  "is:", tf, "\n")
            }
          } ## if(biCC)
          allPdsep.tmp[[x]] <- c(tf,y) ## you must add y to the set
          ## for the large scale simulations, we need to stop the algorithm if
          ## it takes to much time, i.e. sepset>25
          if (length(tf) > pdsep.max) {
            if(verbose)
              cat("Size of Possible-D-SEP bigger than",pdsep.max,
                  ". Break the search for the edge between", x,"and",y,"\n")
          } else if (length(diff.set) > 0) {
            done <- FALSE
            ord <- 0L
            while (!done && ord < min(length(tf), m.max)) {
              ord <- ord + 1L
              if(verbose) cat("ord = ", ord, "\n")
              if (ord == 1) {
                for (S in diff.set) {
                  pval <- indepTest(x, y, S, suffStat)
                  n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                  if (is.na(pval))
                    pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                  if (pval > pMax[x, y])
                    pMax[x, y] <- pval
                  if (pval >= alpha) {
                    amat[x, y] <- amat[y, x] <- 0
                    sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                    done <- TRUE
                    if (verbose)
                      cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                    break
                  }
                }
              }
              else { ## ord > 1
                tmp.combn <- combn(tf, ord) ## has  choose( |tf|, ord ) columns
                if (ord <= length(adj.x)) {
                  for (k in seq_len(ncol(tmp.combn))) {
                    S <- tmp.combn[, k]
                    if (!all(S %in% adj.x)) {
                      n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                      pval <- indepTest(x, y, S, suffStat)
                      if (is.na(pval))
                        pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                      if(pMax[x, y] < pval)
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        amat[x, y] <- amat[y, x] <- 0
                        sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                        done <- TRUE
                        if (verbose)
                          cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                        break
                      }
                    }
                  } ## for(k ..)
                }
                else { ## ord > |adj.x| :
                  ## check all combinations; no combination has been tested before
                  for (k in seq_len(ncol(tmp.combn))) {
                    S <- tmp.combn[, k]
                    n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                    pval <- indepTest(x, y, S, suffStat)
                    if (is.na(pval))
                      pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                    if(pMax[x, y] < pval)
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                      done <- TRUE
                      if (verbose)
                        cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                      break
                    }
                  } ## for(k ..)
                } ## else: { ord > |adj.x| }
              } ## else

            } ## while(!done ..)
          }

        } ## for(y ..)

      } ## if(any( . ))

    } ## for(x ..)
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE

  } ## if(any(G))

  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep.tmp,
       max.ord = ord, n.edgetests = n.edgetests[1:(ord + 1)])
} ## {pdsep2}


udag2pag2 <- function(pag, sepset, rules = rep(TRUE,10), unfVect = NULL, contextVars = NULL, selectionBias = TRUE, jci = c("0","1","123"), verbose = FALSE, orientCollider = TRUE)
{
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PAG using
  ## the rules of Zhang. The output is an adjacency matrix.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - pag: adjacency matrix of size pxp
  ## - sepset: list of all separation sets
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - unfVect: Vector with ambiguous triples (coded as number using triple2numb)
  ## - contextVars: subset of variables to treat as context variables   JORIS
  ## - selectionBias: allow for selection bias (default: TRUE)
  ## - jci: specifies the JCI background knowledge that is used; can be either:
  ##     "0"   no JCI background knowledge (default),
  ##     "1"   JCI assumption 1 only (i.e., no system variable causes any context variable),
  ##     "123" all JCI assumptions 1, 2 and 3
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 6 Mar 2009; cleanup: Martin Maechler, 2010
  ## update: Diego Colombo, 2012; Joris Mooij (2020)

  ## Notation:
  ## ----------------------------------------------------------------------
  ## 0: no edge
  ## 1: -o
  ## 2: -> (arrowhead)
  ## 3: - (tail)
  ## a=alpha
  ## b=beta
  ## c=gamma
  ## d=theta

  stopifnot(is.logical(rules), length(rules) == 10)
  if(!is.numeric(pag))  storage.mode(pag) <- "numeric"

  jci <- match.arg(jci)

  if( !selectionBias ) {
    rules[5:7]<-FALSE
  }
  if( !is.null(contextVars) && selectionBias ) {
    stop('The combination selection bias and context is not implemented yet.')
  }

  if (any(pag != 0)) {
    p <- as.numeric(dim(pag)[1])

    ## special treatment for context variables
    if (!is.null(contextVars) && (jci == "123" || jci == "1")) {
      for (x in 1:p) {
        for (y in 1:p) {
          if (jci == "123" && x != y && x %in% contextVars && y %in% contextVars ) {
            if (verbose)
              cat("\nRule 0 (JCI Assumption 3)",
                  "\nOrient:", x, "*->", y, " because they are contextVars\n")
            pag[x,y]<-2
          } 
          if (pag[x,y] && (x %in% contextVars) && !(y %in% contextVars)) {
            if( jci == "123" ) {
              if (verbose)
                cat("\nRule 0 (JCI Assumptions 1,2)",
                    "\nOrient:", x, "-->", y, " because the first is in contextVars, the second isn't\n")
              pag[y,x]<-3
              pag[x,y]<-2
            } else if( jci == "1" ) {
              if (verbose)
                cat("\nRule 0 (JCI Assumption 1)",
                    "\nOrient:", x, "*->", y, " because the first is in contextVars, the second isn't\n")
              pag[x,y]<-2
            }
          }
        }
      }
    }

    ## orient collider
    if (orientCollider) {
      ind <- which(pag == 1, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        allZ <- setdiff(which(pag[y, ] != 0), x)
        if( !(jci == "123" && y %in% contextVars) ) {
        for (z in allZ) {
          if (pag[x, z] == 0 && !((y %in% sepset[[x]][[z]]) ||
                   (y %in% sepset[[z]][[x]]))) {
            if (length(unfVect) == 0) {
              if (verbose) {
                cat("\n", x, "*->", y, "<-*", z, "\n")
                cat("Sxz=", sepset[[z]][[x]], "and",
                    "Szx=", sepset[[x]][[z]], "\n")
              }
              pag[x, y] <- pag[z, y] <- 2
            }
            else {
              if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
                if (verbose) {
                  cat("\n", x, "*->", y, "<-*", z, "\n")
                  cat("Sxz=", sepset[[z]][[x]], "and",
                      "Szx=", sepset[[x]][[z]], "\n")
                }
                pag[x, y] <- pag[z, y] <- 2
              }
            }
          }
        }
        }
      }
    } ## end: Orient collider

    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) { ## R1--R4
      old_pag1 <- pag
      ##-- R1 ------------------------------------------------------------------
      if (rules[1]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2] # pag[a,b] == 2, pag[b,a] != 0
          indC <- which((pag[b, ] != 0 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          indC <- setdiff(indC, a) # pac[b,C] != 0, pag[C,b] ==1, pag[a,C] == 0, pag[C,a] == 0
          if (length(indC) > 0) {
            if (length(unfVect) == 0) {
              if( any(pag[b, indC] == 3) )
                stop('Contradiction in Rule 1')
              pag[b, indC] <- 2
              pag[indC, b] <- 3
              if (verbose)
                cat("\nRule 1",
                    "\nOrient:", a, "*->", b, "o-*", indC,
                    "as:", b, "->", indC, "\n")
            }
            else {
              for (c in indC) {
                if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                    !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                  if( pag[b, c] == 3 )
                    stop('Contradiction in Rule 1')
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                  if (verbose)
                    cat("\nRule 1",
                        "\nConservatively orient:", a, "*->", b, "o-*",
                        c, "as:", b, "->", c, "\n")
                }
              } ## for( c )
            }
          }
        } ## for( i )
      }
      ##-- R2 ------------------------------------------------------------------
      if (rules[2]) {
        ind <- which((pag == 1 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2] # pag[a,c] == 1, pag[c,a] != 0
          indB <- which((pag[a, ] == 2 & pag[, a] == 3 & pag[c, ] != 0 & pag[, c] == 2) | (pag[a, ] == 2 & pag[, a] != 0 & pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) > 0) {
            pag[a, c] <- 2
            if (verbose) {
              cat("\nRule 2","\n")
              cat("Orient:", a, "->", indB, "*->", c, "or", a, "*->", indB, "->", c, "with", a, "*-o", c, "as:", a, "*->", c, "\n")
            }
          }
        }
      }
      ##-- R3 ------------------------------------------------------------------
      if (rules[3]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          d <- ind[i, 2] # pag[b,d] != 0, pag[d,b] == 1
          indAC <- which((pag[b, ] != 0 & pag[, b] == 2) & (pag[, d] == 1 & pag[d, ] != 0))
          if (length(indAC) >= 2) {
            if (length(unfVect) == 0) {
              counter <- 0
              while ((counter < (length(indAC) - 1)) && (pag[d, b] != 2)) {
                counter <- counter + 1
                ii <- counter
                while (ii < length(indAC) && pag[d, b] != 2) {
                  ii <- ii + 1
                  if (pag[indAC[counter], indAC[ii]] == 0 && pag[indAC[ii], indAC[counter]] == 0) {
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Orient:", d, "*->", b, "\n")
                    }
                    pag[d, b] <- 2
                  }
                }
              }
            }
            else {
              comb.indAC <- combn(indAC, 2)
              for (j in 1:dim(comb.indAC)[2]) {
                a <- comb.indAC[1, j]
                c <- comb.indAC[2, j]
                if (pag[a, c] == 0 && pag[c, a] == 0 && c != a) {
                  if (!any(unfVect == triple2numb(p, a, d, c), na.rm = TRUE) &&
                      !any(unfVect == triple2numb(p, c, d, a), na.rm = TRUE)) {
                    pag[d, b] <- 2
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Conservatively orient:", d, "*->", b, "\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
      ##-- R4 ------------------------------------------------------------------
      if (rules[4]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        while (length(ind) > 0) {
          b <- ind[1, 1]
          c <- ind[1, 2] # pag[b,c] != 0, pag[c,b] == 1
          ind <- ind[-1,, drop = FALSE]
          ## find all a s.t. a -> c and a <-* b
          indA <- which((pag[b, ] == 2 & pag[, b] != 0) &
                        (pag[c, ] == 3 & pag[, c] == 2))
          # pag[b,a] == 2, pag[a,b] != 0
          # pag[c,a] == 3, pag[a,c] == 2
          ## chose one a s.t. the initial triangle structure exists and the edge hasn't been oriented yet
          while (length(indA) > 0 && pag[c,b] == 1) {
            a <- indA[1]
            indA <- indA[-1]
            ## path is the initial triangle
            ## abc <- c(a, b, c)
            ## Done is TRUE if either we found a minimal path or no path exists for this triangle
            Done <- FALSE
### MM: FIXME?? Isn't  Done  set to TRUE in *any* case inside the following
### while(.), the very first time already ??????????
            while (!Done && pag[a,b] != 0 && pag[a,c] != 0 && pag[b,c] != 0) {
              ## find a minimal discriminating path for a,b,c
              md.path <- minDiscrPath(pag, a,b,c, verbose = verbose)
              ## if the path doesn't exists, we are done with this triangle
              if ((N.md <- length(md.path)) == 1) {
                Done <- TRUE
              }
              else {
                ## a path exists
                ## if b is in sepset
                if ((b %in% sepset[[md.path[1]]][[md.path[N.md]]]) ||
                    (b %in% sepset[[md.path[N.md]]][[md.path[1]]])) {
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is in Sepset of",
                        c, "and", md.path[1], ". Orient:", b, "->", c, "\n")
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                }
                else {
                  ## if b is not in sepset
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is not in Sepset of",
                        c, "and", md.path[1], ". Orient:", a, "<->", b, "<->",
                        c, "\n")
                  pag[b,c] <- pag[c,b] <- 2
                  if( pag[a,b] ==  3 ) {
                    warning('Contradiction in Rule 4b, not overwriting...')
                  } else {
                    pag[a,b] <- 2
                  }
                }
                Done <- TRUE
              }
            }
          }
        }
      }
    } ## R1-R4
    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) { ## R5-R7
      old_pag1 <- pag
      ##-- R5 ------------------------------------------------------------------
      if (rules[5]) {
        ind <- which((pag == 1 & t(pag) == 1), arr.ind = TRUE) ## a o-o b
        while (length(ind) > 0) {
          a <- ind[1, 1]
          b <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all c s.t. a o-o c and c is not connected to b
          indC <- which((pag[a, ] == 1 & pag[, a] == 1) & (pag[b, ] == 0 & pag[, b] == 0))
          ## delete b since it is surely in indC
          indC <- setdiff(indC, b)
          ## find all d s.t. b o-o d and d is not connected to a
          indD <- which((pag[b, ] == 1 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          ## delete a since it is surely in indD
          indD <- setdiff(indD, a)
          if (length(indC) > 0 && length(indD) > 0) {
            counterC <- 0
            while ((counterC < length(indC)) && pag[a, b] == 1) {
              counterC <- counterC + 1
              c <- indC[counterC]
              counterD <- 0
              while ((counterD < length(indD)) && pag[a, b] == 1) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if (pag[c, d] == 1 && pag[d, c] == 1) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[a, b] <- pag[b, a] <- 3
                    pag[a, c] <- pag[c, a] <- 3
                    pag[c, d] <- pag[d, c] <- 3
                    pag[d, b] <- pag[b, d] <- 3
                    if (verbose)
                      cat("\nRule 5",
                          "\nThere exists an uncovered circle path between", a, "and", b,
                          ". Orient:", a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                  }
                  else { ## conservative: check that every triple on the circle is faithful
                    path2check <- c(a,c,d,b)
                    if (faith.check(path2check, unfVect, p)) {
                      pag[a, b] <- pag[b, a] <- 3
                      pag[a, c] <- pag[c, a] <- 3
                      pag[c, d] <- pag[d, c] <- 3
                      pag[d, b] <- pag[b, d] <- 3
                      if (verbose)
                        cat("\nRule 5",
                            "\nThere exists a faithful uncovered circle path between",
                            a, "and", b, ". Conservatively orient:",
                            a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                    }
                  }
                }
                ## search with a breitensuche a minimal uncovered circle path
                else {
                  ## Find a minimal uncovered circle path for these a,b,c, and d.
                  ## This path has already been checked to be uncovered and
                  ## to be faithful for the conservative case
		  ucp <- minUncovCircPath(p, pag = pag, path = c(a,c,d,b),
					  unfVect = unfVect, verbose = verbose)
                  ## there is a path ---> orient
                  if (length(ucp) > 1) {
                    ## orient every edge on the path as --
                    n <- length(ucp)
                    pag[ucp[1], ucp[n]] <- pag[ucp[n], ucp[1]] <- 3 ## a--b
                    for (j in 1:(length(ucp)-1)) ## each edge on the path --
                      pag[ucp[j], ucp[j + 1]] <- pag[ucp[j + 1], ucp[j]] <- 3
                  }
                }
              }
            }
          }
        }
      }
      ##-- R6 ------------------------------------------------------------------
      if (rules[6]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          if (any(pag[b, ] == 3 & pag[, b] == 3)) {
            pag[c, b] <- 3
            if (verbose)
              cat("\nRule 6",
                  "\nOrient:", b, "o-*", c, "as", b, "-*", c, "\n")
          }
        }
      }
      ##-- R7 ------------------------------------------------------------------
      if (rules[7]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          indA <- which((pag[b, ] == 3 & pag[, b] == 1) & (pag[c, ] == 0 & pag[, c] == 0))
          indA <- setdiff(indA, c)
          if (length(indA) > 0) {
            if (length(unfVect) == 0) {
              pag[c, b] <- 3
              if (verbose)
                cat("\nRule 7",
                    "\nOrient:", indA, "-o", b, "o-*",
                    c, "as", b, "-*", c, "\n")
            }
            else for (a in indA)
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pag[c, b] <- 3
                if (verbose)
                  cat("\nRule 7",
                      "\nConservatively orient:", a, "-o", b, "o-*",
                      c, "as", b, "-*", c, "\n")
              }
          }
        }
      }
    } ## R5-R7
    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) { ## R8-R10
      old_pag1 <- pag
      ##-- R8 ------------------------------------------------------------------
      if (rules[8]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          # pag[a,c] == 2, pag[c,a] == 1
          indB <- which(pag[, a] == 3 & (pag[a, ] == 2 | pag[a, ] == 1) &
                        pag[c, ] == 3 & pag[, c] == 2)
          if (length(indB) > 0) {
            pag[c, a] <- 3
            if (verbose)
              cat("\nRule 8",
                  "\nOrient:", a, "->", indB, "->", c,
                  "or", a, "-o", indB, "->", c, "with", a,
                  "o->", c, "as", a, "->", c, "\n")
          }
        }
      }
      ##-- R9 ------------------------------------------------------------------
      if (rules[9]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          # pag[a,c] == 2, pag[c,a] == 1
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. a (o-)--(o>) b and b and c are not connected
          indB <- which((pag[a, ] == 2 | pag[a, ] == 1) &
                        (pag[, a] == 1 | pag[, a] == 3) &
                        (pag[c, ] == 0 & pag[, c] == 0))
          ## delete c from indB since it is surely inside
          indB <- setdiff(indB, c)
          ## chose one b s.t. the initial structure exists and the edge hasn't been oriented yet
          while ((length(indB) > 0) && (pag[c,a] == 1)) {
            b <- indB[1]
            indB <- indB[-1]
            ## find a minimal uncovered pd path from initial (a,b,c) :
            upd <- minUncovPdPath(p, pag, a,b,c,
                unfVect = unfVect, verbose = verbose)
            ## there is a path ---> orient it
            if (length(upd) > 1) {
              pag[c, a] <- 3
              if (verbose)
                cat("\nRule 9",
                    "\nThere exists an uncovered potentially directed path between", a, "and", c,
                    ". Orient:", a, " ->",c, "\n")
            }
          }
        }
      }
      ##-- R10 ------------------------------------------------------------------
      if (rules[10]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          # pag[a,c] == 2, pac[c,a] == 1
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. b --> c
          indB <- which((pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) >= 2) {
            counterB <- 0
            while (counterB < length(indB) && (pag[c, a] == 1)) {
              counterB <- counterB + 1
              b <- indB[counterB]
              indD <- setdiff(indB, b)
              counterD <- 0
              while ((counterD < length(indD)) && (pag[c, a] == 1)) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if ((pag[a, b] == 1 || pag[a, b] == 2) &&
                    (pag[b, a] == 1 || pag[b, a] == 3) &&
                    (pag[a, d] == 1 || pag[a, d] == 2) &&
                    (pag[d, a] == 1 || pag[d, a] == 3) && pag[d, b] == 0 && pag[b, d] == 0) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[c, a] <- 3
                    if (verbose)
                      cat("\nRule 10 [easy]",
                          "\nOrient:", a, "->", c, "\n")
                  }
                  else ## conservative version: check faithfulness of b-a-d
                    if (!any(unfVect == triple2numb(p,b,a,d), na.rm = TRUE) &&
                        !any(unfVect == triple2numb(p,d,a,b), na.rm = TRUE)) {
                      pag[c, a] <- 3
                      if (verbose)
                        cat("\nRule 10 [easy]",
                            "\nConservatively orient:", a, "->", c, "\n")
                    }
                }
                ## search with a breitensuche two minimal uncovered circle paths
                else {
                  ## find all x s.t. a (o-)--(o>) x
                  indX <- which((pag[a, ] == 1 | pag[a, ] == 2) &
                                (pag[, a] == 1 | pag[, a] == 3), arr.ind = TRUE)
                  indX <- setdiff(indX, c)
                  if (length(indX >= 2)) {
                    counterX1 <- 0
                    while (counterX1 < length(indX) && pag[c, a] == 1) {
                      counterX1 <- counterX1 + 1
                      first.pos <- indA[counterX1]
                      indX2 <- setdiff(indX, first.pos)
                      counterX2 <- 0
                      while (counterX2 < length(indX2) && pag[c, a] == 1) {
                        counterX2 <- counterX2 + 1
                        sec.pos <- indX2[counterX2]
                        t1 <- minUncovPdPath(p, pag, a, first.pos, b,
                                             unfVect = unfVect, verbose = verbose)
                        if (length(t1) > 1) { # otherwise, can skip next minUnc..()
                          t2 <- minUncovPdPath(p, pag, a, sec.pos, d,
                                               unfVect = unfVect, verbose = verbose)
                          if (length(t2) > 1 &&
                              first.pos != sec.pos && pag[first.pos, sec.pos] == 0) {
                            ## we found 2 uncovered pd paths
                            if (length(unfVect) == 0) { ## normal version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10", "\nOrient:", a, "->", c, "\n")
                            }
                            else if(!any(unfVect == triple2numb(p,first.pos, a, sec.pos), na.rm = TRUE) &&
                                    !any(unfVect == triple2numb(p,sec.pos, a, first.pos), na.rm = TRUE)) {
                              ## conservative version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10",
                                    "\nConservatively orient:", a, "->", c, "\n")
                            }
                          }
                        }
                      } #  # while ( counterX2 .. )
                    }
                  }
                } # else
              } # while ( counterD .. )
            } # while ( counterB .. )
          } # if (length(indB) .)
        }
      } ## if (rules[10] ..)
    } ## R8-R10
  }
  pag
} ## udag2pag2()




######################################################################################################
# everything below this line is a literal copy of parts of pcalg.R

updateList <- function(path, set, old.list)
{
  ## Purpose: update the list of all paths in the iterative functions
  ## minDiscrPath, minUncovCircPath and minUncovPdPath
  ## ----------------------------------------------------------------------
  ## Arguments: - path: the path under investigation
  ##            - set: (integer) index set of variables to be added to path
  ##            - old.list: the list to update
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2011; Without for() by Martin Maechler
  c(old.list, lapply(set, function(s) c(path,s)))
}

## R9-R10
minUncovPdPath <- function(p, pag, a,b,c, unfVect, verbose = FALSE)
{
  ## Purpose: find a minimal uncovered pd path for a,b,c saved in path.
  ## Check also for the conservative case that it is unambiguous
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - p: number of nodes in the graph
  ##            - pag: adjacency matrix
  ##            - a,b,c : nodes under interest
  ##            - unfVect: vector containing the ambiguous triples
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 19 Oct 2011; small changes: Martin Maechler

  visited <- rep(FALSE, p)
  visited[c(a,b,c)] <- TRUE
  min.upd.path <- NA
  ## find all neighbours of b not visited yet
  indD <- which((pag[b,] == 1 | pag[b,] == 2) &
                (pag[,b] == 1 | pag[,b] == 3) &
                (pag[,a] == 0) & !visited)
  if (length(indD) > 0) {
    path.list <- updateList(b, indD, NULL)
    done <- FALSE
    while ((length(path.list) > 0) && (!done)) {
      ## next element in the queue
      mpath <- path.list[[1]]
      m <- length(mpath)
      d <- mpath[m]
      path.list[[1]] <- NULL
      visited[d] <- TRUE
      if (any(pag[d,c] == 1:2) && any(pag[c,d] == c(1,3))) {
        ## pd path found
        mpath <- c(a, mpath, c)
        n <- length(mpath)
        ## check the path to be uncovered
        uncov <- TRUE
	for (l in seq_len(n - 2)) {
	  if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
		pag[mpath[l + 2], mpath[l]] == 0)) {

	    uncov <- FALSE
	    break ## speed up!
	  }
	}
        ## if it is uncovered
        if (uncov)
          if (length(unfVect) == 0 || ## <<- normal version: just save
              ## conservative version, check the path to be faithful:
              faith.check(mpath, unfVect, p)) {
            ## save the path to be returned
            min.upd.path <- mpath
            done <- TRUE
          }
      }
      else {
        ## d and c are either not connected or connected with a "wrong" edge -----> search iteratively
        ## find all neighbours of d not visited yet
        indR <- which((pag[d,] == 1 | pag[d,] == 2) &
                      (pag[,d] == 1 | pag[,d] == 3) & !visited)
        if (length(indR) > 0) {
          ## update the queues
          path.list <- updateList(mpath, indR, path.list)
        }
      }
    } ## {while}
  }
  min.upd.path
} ## {minUncovPdPath}

## R5
minUncovCircPath <- function(p, pag, path, unfVect, verbose = FALSE)
{
  ## Purpose: find a minimal uncovered circle path for a,c,d,b saved in path.
  ## Check also for the conservative case that it is unambiguous
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - p: number of nodes in the graph
  ##            - pag: adjacency matrix
  ##            - path: a,c,d,b under interest
  ##            - unfVect: vector containing the unfaithful triples
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 19 Oct 2011, 13:11

  visited <- rep(FALSE, p)
  visited[path] <- TRUE # (a,b,c,d) all 'visited'
  a <- path[1]
  c <- path[2]
  d <- path[3]
  b <- path[4]
  min.ucp.path <- NA
  ## find all neighbours of c not visited yet
  indX <- which(pag[c,] == 1 & pag[,c] == 1 & !visited) ## c o-o x
  if (length(indX) > 0) {
    path.list <- updateList(c, indX, NULL)
    done <- FALSE
    while (!done && length(path.list) > 0) {
      ## next element in the queue
      mpath <- path.list[[1]]
      x <- mpath[length(mpath)]
      path.list[[1]] <- NULL
      visited[x] <- TRUE
      if (pag[x,d] == 1 && pag[d,x] == 1) {
        ## circle path found
        mpath <- c(a, mpath, d, b)
        n <- length(mpath)
        ## check the path to be uncovered
        uncov <- TRUE
	for (l in seq_len(n - 2)) {
	  if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
		pag[mpath[l + 2], mpath[l]] == 0)) {

	    uncov <- FALSE
	    break ## speed up!
	  }
	}
        ## if it is uncovered
        if (uncov)
          if (length(unfVect) == 0 || ## <<- normal version: just save
              ## conservative version, check the path to be faithful:
              faith.check(mpath, unfVect, p)) {
            ## save the path to be returned
            min.ucp.path <- mpath
            done <- TRUE
          }
      }
      else {
        ## x and d are either not connected or connected with an edge which is not o--o -----> search iteratively
        ## find all neighbours of x not visited yet
        indR <- which(pag[x,] == 1 & pag[,x] == 1 & !visited) ## x o--o r
        if (length(indR) > 0) {
          ## update the queues
          path.list <- updateList(mpath, indR, path.list)
        }
      }
    }## {while}
  }
  return(min.ucp.path)
}

## R4
minDiscrPath <- function(pag, a,b,c, verbose = FALSE)
{
  ## Purpose: find a minimal discriminating path for a,b,c.
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - pag: adjacency matrix
  ##            - a,b,c: node positions under interest
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 Jan 2011; speedup: Martin Maechler

  p <- as.numeric(dim(pag)[1])
  visited <- rep(FALSE, p)
  visited[c(a,b,c)] <- TRUE # {a,b,c} "visited"
  ## find all neighbours of a  not visited yet
  indD <- which(pag[a,] != 0 & pag[,a] == 2 & !visited) ## d *-> a
  if (length(indD) > 0) {
    path.list <- updateList(a, indD, NULL)
    while (length(path.list) > 0) {
      ## next element in the queue
      mpath <- path.list[[1]]
      m <- length(mpath)
      d <- mpath[m]
      if (pag[c,d] == 0 & pag[d,c] == 0)
	## minimal discriminating path found :
	return( c(rev(mpath), b,c) )

      ## else :
      pred <- mpath[m-1]
      path.list[[1]] <- NULL


      ## d is connected to c -----> search iteratively
      if (pag[d,c] == 2 && pag[c,d] == 3 && pag[pred,d] == 2) {
        visited[d] <- TRUE
        ## find all neighbours of d not visited yet
        indR <- which(pag[d,] != 0 & pag[,d] == 2 & !visited) ## r *-> d
        if (length(indR) > 0)
          ## update the queues
          path.list <- updateList(mpath[-1], indR, path.list)
      }
    } ## {while}
  }
  ## nothing found:  return
  NA
} ## {minDiscrPath}
