# (c) Douglas G. Scofield, McGill University
#
# determine effect of SNP, given GFF annotation

SNP.ref.must.match.genome <- TRUE


# rows of GFF features that overlap with this window
posFeaturesOverlap <-
function(x, gff)
{
  if (length(x) == 1)
    x <- c(x, x)  # how to give a one-bp window
  x <- range(sort(x))
  sel <- x[1] >= start(gff) & x[2] <= end(gff)
  ###
  gff[sel, ]
}

SNP_EFFECT <- c("CODING",    # within protein-coding region of gene
                "TRNA",      # within tRNA gene
                "RRNA",       # within rRNA gene
                "INTRAGENIC_NON_CODING_RNA",
                "5PRIME_UTR",       # within the 5' UTR
                "3PRIME_UTR",       # within the 3' UTR
                "EXTRAGENIC",  # outside CDS and RNA genes
                "SYNONYMOUS_CODING",   # no change in AA
                "NON_SYNONYMOUS_CODING",   # change in AA
                "STOP_GAINED",         # new Stop codon
                "STOP_LOST",         # conversion of Stop codon
                "FRAME_SHIFT",       # via indel
                "INTRON") 
names(SNP_EFFECT) <- SNP_EFFECT



snpEffect <-
function(snp, gff, genome, nHaplotypes=20, report.interval=10, ...)
{

  # snp is a (currently) VarScan-formatted data frame; gff is a BioC
  # import.gff-formatted file; genome is a BioC read.DNAStringSet-formatted
  # XStringSet.

  snp.ans <- NULL

  #          data.frame(chromosome=I(""),
  #                     position=0,
  #                     ref=I(""),
  #                     change=I(""),
  #                     changetype=I(""),
  #                     Homozygous=I(""),
  #                     changefreq=0,
  #                     gene_id=I(""),
  #                     gene_name=I(""),
  #                     strand=I(""),
  #                     effect=I(""),
  #                     oldAA_newAA=I(""),
  #                     oldcodon_newcodon=I(""),
  #                     codon_num=0,
  #                     CDS_size=0,
  #                     window1010bp=I(""),
  #                     note=I(""))[0, ]
  cat(nrow(snp), "SNPs\n")
  for (i in 1:nrow(snp)) {
    pos <- snp[i, "Position"]
    feat <- posFeaturesOverlap(pos, gff)

    # first, check to see if we have overlapping genes (+ and -, + and +, etc.)
    locus_tags <- unique.no.NA(feat$locus_tag[feat$codon_start == 1])
    if (length(locus_tags) > 0) {
      for (tag in locus_tags) {
        #thisfeat <- subset(feat, locus_tag == tag)
        rw <- feat$locus_tag == tag
        rw[is.na(rw)] <- FALSE
        thisfeat <- feat[rw, ]
        thisfeat <- rbind(thisfeat, subset(feat, type == "source"))
        thissnpEffect <- singleSNPEffect(snp[i, ], thisfeat, genome, nHaplotypes, ...)
        snp.ans <- rbind(snp.ans, thissnpEffect)
      }
    } else {
      thissnpEffect <- singleSNPEffect(snp[i, ], feat, genome, nHaplotypes, ...)
      snp.ans <- rbind(snp.ans, thissnpEffect)
    }
    if (! (i %% report.interval)) 
      cat(i, "SNPs processed, last at position", pos, "\n")
  }
  rownames(snp.ans) <- 1:nrow(snp.ans)
  ####
  snp.ans
}


singleSNPEffect <-
function(thissnp, feat, genome, nHaplotypes=20, debug=FALSE, verbose=FALSE)  
# thissnp is position etc. of snp, features is overlap of gff
{
  pos <- thissnp$Position
  ref <- as.character(thissnp$ReferenceAllele)
  refref <- as.character(subseq(genome, start=pos, width=1))
  if (ref != refref && SNP.ref.must.match.genome) {
    stop("SNP @", pos, ": SNP reference allele and reference genome do not match")
  }
  change <- as.character(thissnp$SNPAlleleQ)
  changefreq <- thissnp$NumStrainsQ / nHaplotypes
  effect.feat <- feat
  if ("CDS" %in% feat$type) {
    # at least one hit is coding
    effect <- SNP_EFFECT["CODING"]
    effect.feat <- subset(feat, type == "CDS")
  } else if ("tRNA" %in% feat$type) {
    effect <- SNP_EFFECT["TRNA"]
    effect.feat <- subset(feat, type == "tRNA")
  } else if ("rRNA" %in% feat$type) {
    effect <- SNP_EFFECT["RRNA"]
    effect.feat <- subset(feat, type == "rRNA")
  } else if ("gene" %in% feat$type) {
    effect <- SNP_EFFECT["INTRAGENIC_NON_CODING_RNA"]
    effect.feat <- subset(feat, type == "gene")
  } else if (nrow(feat) == 1 && feat$type == "source") {
    effect <- SNP_EFFECT["EXTRAGENIC"]
    effect.feat <- subset(feat, type == "source")
  }
  strand <- unique.no.NA(effect.feat$strand)
  # find the gene; this may be a group entry so we have to find the
  # first of the group
  gene_name <- ""
  if (effect != SNP_EFFECT["EXTRAGENIC"]) {
    if (! all(is.na(effect.feat$gene))) {
      # the matching feature has a gene name
      gene_name <- unique.no.NA(effect.feat$gene)
    } else {
      # we may be part of a group, find the gene via the group
      group.name <- unique.no.NA(effect.feat$group)
      if (group.name != "") {
        gfeat <- subset(feat, group == group.name)
        gene_name <- unique.no.NA(gfeat$gene)
      }
    }
  }
  oldAA_newAA <- ""
  oldcodon_newcodon <- ""
  codon_num <- 0
  CDS_size <- 0
  types <- unique.no.NA(feat$type)
  types <- types[! types %in% c("gene", "source")]
  if (verbose) 
    cat("SNP", "@", pos, paste(sep="/", ref, change), paste(collapse=":", types), " ")

  if (effect == "CODING") {
    while (TRUE) {
      CDS <- subset(feat, type == "CDS")
      if (nrow(CDS) > 1) {
        cat("SNP @", pos, ": can't handle SNPs in genes with >1 CDS\n")
        break
      }
      #if (CDS$strand != "+") {
      #  cat("  : can only handle + strand\n")
      #  break
      #}
      if (debug) 
        cat(" ", CDS$strand, "start", start(CDS), "end", end(CDS), "width", width(CDS), "\n")
      translation <- unique.no.NA(CDS$translation)
      if (translation == "") {
        cat("SNP @", pos, ": no annotated translation available\n")
      } else {  # use the existing translation

        if (CDS$strand == "+") {

          seq.CDS <- subseq(genome, start=start(CDS), width=width(CDS))
          CDS.pos <- pos - start(CDS) # 0-based
          CDS.codon <- CDS.pos %/% 3  # 0-based
          codon_num <- CDS.codon + 1 # 1-based
          CDS.phase <- CDS.pos %% 3  # 0, 1, 2
          if (debug) 
            cat("  + CDS.pos", CDS.pos, "CDS.codon", CDS.codon, "CDS.phase", CDS.phase, "\n")
          seq.codon.pos <- start(CDS) + (CDS.codon * 3)
          ref.codon <- as.character(subseq(genome, start=seq.codon.pos, width=3))
          snp.codon <- ref.codon
          substr(snp.codon, start=(CDS.phase + 1), stop=(CDS.phase + 1)) <- change
          if (debug)
            cat("  + seq.codon.pos", seq.codon.pos, "ref.codon", ref.codon, "snp.codon", snp.codon, "\n")
          ref.AA <- GENETIC_CODE[ref.codon]
          snp.AA <- GENETIC_CODE[snp.codon]
          trans.AA <- "#"
          if (nchar(translation) < codon_num) { # translation is too short
            if (codon_num == (nchar(translation) + 1)) { # may be stop codon
              if (ref.AA == "*") { # yep!
                trans.AA <- "*"
                if (debug) 
                  cat("  + SNP falls within the stop codon\n")
              } else {
                cat("SNP @", pos, ": annotated translation is too short\n")
              }
            }
          } else {
            trans.AA <- subseq(translation, codon_num, codon_num)
          }
          if (debug) 
            cat("  + trans.AA", trans.AA, "ref.AA", ref.AA, "snp.AA", snp.AA, "\n")
          oldAA_newAA <- paste(sep="/", ref.AA, snp.AA)
          oldcodon_newcodon <- paste(sep="/", ref.codon, snp.codon)
          CDS_size <- nchar(translation)
          # adjust effect
          if (ref.AA == snp.AA) 
            effect <- SNP_EFFECT["SYNONYMOUS_CODING"]
          else if (snp.AA == "*")
            effect <- SNP_EFFECT["STOP_GAINED"]
          else if (ref.AA == "*")
            effect <- SNP_EFFECT["STOP_LOST"]
          else 
            effect <- SNP_EFFECT["NON_SYNONYMOUS_CODING"]

          if (debug) 
            cat("  + effect", effect, "\n")

        } else if (CDS$strand == "-") {

          seq.CDS <- reverseComplement(subseq(genome, start=start(CDS), width=width(CDS)))
          CDS.pos <- end(CDS) - pos  # 0-based
          CDS.codon <- CDS.pos %/% 3  # 0-based
          codon_num <- CDS.codon + 1 # 1-based
          CDS.phase <- CDS.pos %% 3  # 0, 1, 2
          if (debug) 
            cat("  - CDS.pos", CDS.pos, "CDS.codon", CDS.codon, "CDS.phase", CDS.phase, "\n")
          seq.codon.pos <- end(CDS) - (CDS.codon * 3) - 2
          ref.codon.rc <- subseq(genome, start=seq.codon.pos, width=3)
          snp.codon.rc <- ref.codon.rc
          subseq(snp.codon.rc, start=(3 - CDS.phase), width=1) <- change
          ref.codon <- as.character(reverseComplement(ref.codon.rc))
          snp.codon <- as.character(reverseComplement(snp.codon.rc))
          if (debug) 
            cat("  - seq.codon.pos", seq.codon.pos, "ref.codon", ref.codon, "snp.codon", snp.codon, "\n")
          ref.AA <- GENETIC_CODE[ref.codon]
          snp.AA <- GENETIC_CODE[snp.codon]
          trans.AA <- "#"
          if (nchar(translation) < codon_num) { # translation is too short
            if (codon_num == (nchar(translation) + 1)) { # may be stop codon
              if (ref.AA == "*") { # yep!
                trans.AA <- "*"
                if (debug) 
                  cat("  - SNP falls within the stop codon\n")
              } else {
                cat("SNP @", pos, ": annotated translation is too short\n")
              }
            }
          } else {
            trans.AA <- subseq(translation, codon_num, codon_num)
          }
          if (debug) 
            cat("  - trans.AA", trans.AA, "ref.AA", ref.AA, "snp.AA", snp.AA, "\n")
          oldAA_newAA <- paste(sep="/", ref.AA, snp.AA)
          oldcodon_newcodon <- paste(sep="/", ref.codon, snp.codon)
          CDS_size <- nchar(translation)
          # adjust effect
          if (ref.AA == snp.AA) 
            effect <- SNP_EFFECT["SYNONYMOUS_CODING"]
          else if (snp.AA == "*")
            effect <- SNP_EFFECT["STOP_GAINED"]
          else if (ref.AA == "*")
            effect <- SNP_EFFECT["STOP_LOST"]
          else 
            effect <- SNP_EFFECT["NON_SYNONYMOUS_CODING"]

          if (debug) 
            cat("  - effect", effect, "\n")
        }
      }
      break
    }
  }

  # create sequence window around SNP
  w1start <- max(pos - 10, 1); w1width <- pos - w1start
  w2start <- pos + 1; w2width <- min(width(genome) - pos, 10)
  window1010bp <- paste(sep="|", 
                        as.character(subseq(genome, start=w1start, width=w1width)), 
                        change,
                        as.character(subseq(genome, start=w2start, width=w2width))) 
  note <- unique.no.NA(feat$product)
  if ((t1 <- unique.no.NA(feat$protein_id))[1] != "")
    note <- paste(sep="; ", note, paste(sep=":", collapse="; ", "protein_id", t1))
  if ((t1 <- unique.no.NA(feat$locus_tag))[1] != "")
    note <- paste(sep="; ", note, paste(sep=":", collapse="; ", "locus_tag", t1))

  if (verbose) {
    if (gene_name != "")
      cat(gene_name, " ")
    cat(effect, "\n")
  }

  ans <- data.frame(
    chromosome=I(names(genome)[1]),
    position=pos,
    ref=I(ref),
    change=I(change),
    changetype=I("SNP"),
    Homozygous=I("Hap"),
    changefreq=changefreq,
    gene_id=I(""),
    gene_name=I(gene_name),
    strand=paste(sep="/", strand),
    effect=I(effect),
    oldAA_newAA=I(oldAA_newAA),
    oldcodon_newcodon=I(oldcodon_newcodon),
    codon_num=codon_num,
    CDS_size=CDS_size,
    window1010bp=window1010bp,
    note=note)

  ####
  ans
}



indelEffect <-
function(snp, gff, genome, nHaplotypes=20, report.interval=10, ...)
{

  # snp is a (currently) VarScan-formatted data frame; gff is a BioC
  # import.gff-formatted file; genome is a BioC read.DNAStringSet-formatted
  # XStringSet.

  snp.ans <- NULL

  #          data.frame(chromosome=I(""),
  #                     position=0,
  #                     ref=I(""),
  #                     change=I(""),
  #                     changetype=I(""),
  #                     Homozygous=I(""),
  #                     changefreq=0,
  #                     gene_id=I(""),
  #                     gene_name=I(""),
  #                     strand=I(""),
  #                     effect=I(""),
  #                     oldAA_newAA=I(""),
  #                     oldcodon_newcodon=I(""),
  #                     codon_num=0,
  #                     CDS_size=0,
  #                     window1010bp=I(""),
  #                     note=I(""))[0, ]
  cat(nrow(snp), "indels\n")
  for (i in 1:nrow(snp)) {
    pos <- snp[i, "Position"]
    feat <- posFeaturesOverlap(pos, gff)

    # first, check to see if we have overlapping genes (+ and -, + and +, etc.)
    locus_tags <- unique.no.NA(feat$locus_tag[feat$codon_start == 1])
    if (length(locus_tags) > 0) {
      for (tag in locus_tags) {
        #thisfeat <- subset(feat, locus_tag == tag)
        rw <- feat$locus_tag == tag
        rw[is.na(rw)] <- FALSE
        thisfeat <- feat[rw, ]
        thisfeat <- rbind(thisfeat, subset(feat, type == "source"))
        thisIndelEffect <- singleIndelEffect(snp[i, ], thisfeat, genome, nHaplotypes, ...)
        snp.ans <- rbind(snp.ans, thisIndelEffect)
      }
    } else {
      thisIndelEffect <- singleIndelEffect(snp[i, ], feat, genome, nHaplotypes, ...)
      snp.ans <- rbind(snp.ans, thisIndelEffect)
    }
    if (! (i %% report.interval)) 
      cat(i, "indels processed, last at position", pos, "\n")
  }
  rownames(snp.ans) <- 1:nrow(snp.ans)
  ####
  snp.ans
}



singleIndelEffect <-
function(thisindel, feat, genome, nHaplotypes=20, debug=FALSE, verbose=FALSE)  
# thisindel is position etc. of snp, features is overlap of gff
{
  pos <- thisindel$Position
  ref <- as.character(thisindel$Ref)
  change <- thisindel$VarAllele
  refref <- as.character(subseq(genome, start=pos, width=1))
  changefreq <- thisindel$VarFreq
  effect.feat <- feat
  if ("CDS" %in% feat$type) {
    # at least one hit is coding
    effect <- SNP_EFFECT["CODING"]
    effect.feat <- subset(feat, type == "CDS")
  } else if ("tRNA" %in% feat$type) {
    effect <- SNP_EFFECT["TRNA"]
    effect.feat <- subset(feat, type == "tRNA")
  } else if ("rRNA" %in% feat$type) {
    effect <- SNP_EFFECT["RRNA"]
    effect.feat <- subset(feat, type == "rRNA")
  } else if ("gene" %in% feat$type) {
    effect <- SNP_EFFECT["INTRAGENIC_NON_CODING_RNA"]
    effect.feat <- subset(feat, type == "gene")
  } else if (nrow(feat) == 1 && feat$type == "source") {
    effect <- SNP_EFFECT["EXTRAGENIC"]
    effect.feat <- subset(feat, type == "source")
  }
  strand <- unique.no.NA(effect.feat$strand)
  # find the gene; this may be a group entry so we have to find the
  # first of the group
  gene_name <- ""
  if (effect != SNP_EFFECT["EXTRAGENIC"]) {
    if (! all(is.na(effect.feat$gene))) {
      # the matching feature has a gene name
      gene_name <- unique.no.NA(effect.feat$gene)
    } else {
      # we may be part of a group, find the gene via the group
      group.name <- unique.no.NA(effect.feat$group)
      if (group.name != "") {
        gfeat <- subset(feat, group == group.name)
        gene_name <- unique.no.NA(gfeat$gene)
      }
    }
  }

  oldAA_newAA <- ""
  oldcodon_newcodon <- ""
  codon_num <- 0
  CDS_size <- 0
  types <- unique.no.NA(feat$type)
  types <- types[! types %in% c("gene", "source")]
  if (verbose) 
    cat("indel", "@", pos, paste(sep="/", ref, change), paste(collapse=":", types), " ")

  if (effect == "CODING") {
    while (TRUE) {
      CDS <- subset(feat, type == "CDS")
      if (nrow(CDS) > 1) {
        cat("indel @", pos, ": can't handle SNPs in genes with >1 CDS\n")
        break
      }
      #if (CDS$strand != "+") {
      #  cat("  : can only handle + strand\n")
      #  break
      #}
      if (debug) 
        cat(" ", CDS$strand, "start", start(CDS), "end", end(CDS), "width", width(CDS), "\n")
      translation <- unique.no.NA(CDS$translation)
      if (translation == "") {
        cat("SNP @", pos, ": no annotated translation available\n")
      } else {  # use the existing translation

        if (CDS$strand == "+") {

          seq.CDS <- subseq(genome, start=start(CDS), width=width(CDS))
          CDS.pos <- pos - start(CDS) # 0-based
          CDS.codon <- CDS.pos %/% 3  # 0-based
          codon_num <- CDS.codon + 1 # 1-based
          CDS.phase <- CDS.pos %% 3  # 0, 1, 2
          if (debug) 
            cat("  + CDS.pos", CDS.pos, "CDS.codon", CDS.codon, "CDS.phase", CDS.phase, "\n")
          seq.codon.pos <- start(CDS) + (CDS.codon * 3)
          ref.codon <- as.character(subseq(genome, start=seq.codon.pos, width=3))
          snp.codon <- "###"
          substr(snp.codon, start=(CDS.phase + 1), stop=(CDS.phase + 1)) <- change
          if (debug)
            cat("  + seq.codon.pos", seq.codon.pos, "ref.codon", ref.codon, "snp.codon", snp.codon, "\n")
          ref.AA <- GENETIC_CODE[ref.codon]
          snp.AA <- "#"
          trans.AA <- "#"
          if (nchar(translation) < codon_num) { # translation is too short
            if (codon_num == (nchar(translation) + 1)) { # may be stop codon
              if (ref.AA == "*") { # yep!
                trans.AA <- "*"
                if (debug) 
                  cat("  + indel falls within the stop codon\n")
              } else {
                cat("indel @", pos, ": annotated translation is too short\n")
              }
            }
          } else {
            trans.AA <- subseq(translation, codon_num, codon_num)
          }
          if (debug) 
            cat("  + trans.AA", trans.AA, "ref.AA", ref.AA, "snp.AA", snp.AA, "\n")
          oldAA_newAA <- paste(sep="/", ref.AA, snp.AA)
          oldcodon_newcodon <- paste(sep="/", ref.codon, snp.codon)
          CDS_size <- nchar(translation)
          # adjust effect, not yet implemented for indels
          if (FALSE) {
            if (ref.AA == snp.AA)  {
              effect <- SNP_EFFECT["SYNONYMOUS_CODING"]
            } else if (snp.AA == "*") {
              effect <- SNP_EFFECT["STOP_GAINED"]
            } else if (ref.AA == "*") {
              effect <- SNP_EFFECT["STOP_LOST"]
            } else {
              effect <- SNP_EFFECT["NON_SYNONYMOUS_CODING"]
            }
          }

          if (debug) 
            cat("  + effect", effect, "\n")

        } else if (CDS$strand == "-") {

          seq.CDS <- reverseComplement(subseq(genome, start=start(CDS), width=width(CDS)))
          CDS.pos <- end(CDS) - pos  # 0-based
          CDS.codon <- CDS.pos %/% 3  # 0-based
          codon_num <- CDS.codon + 1 # 1-based
          CDS.phase <- CDS.pos %% 3  # 0, 1, 2
          if (debug) 
            cat("  - CDS.pos", CDS.pos, "CDS.codon", CDS.codon, "CDS.phase", CDS.phase, "\n")
          seq.codon.pos <- end(CDS) - (CDS.codon * 3) - 2
          ref.codon.rc <- subseq(genome, start=seq.codon.pos, width=3)
          snp.codon.rc <- ref.codon.rc
          subseq(snp.codon.rc, start=(3 - CDS.phase), width=1) <- change
          ref.codon <- as.character(reverseComplement(ref.codon.rc))
          snp.codon <- "###" # as.character(reverseComplement(snp.codon.rc))
          if (debug) 
            cat("  - seq.codon.pos", seq.codon.pos, "ref.codon", ref.codon, "snp.codon", snp.codon, "\n")
          ref.AA <- GENETIC_CODE[ref.codon]
          snp.AA <- "#" #GENETIC_CODE[snp.codon]
          trans.AA <- "#"
          if (nchar(translation) < codon_num) { # translation is too short
            if (codon_num == (nchar(translation) + 1)) { # may be stop codon
              if (ref.AA == "*") { # yep!
                trans.AA <- "*"
                if (debug) 
                  cat("  - indel falls within the stop codon\n")
              } else {
                cat("indel @", pos, ": annotated translation is too short\n")
              }
            }
          } else {
            trans.AA <- subseq(translation, codon_num, codon_num)
          }
          if (debug) 
            cat("  - trans.AA", trans.AA, "ref.AA", ref.AA, "snp.AA", snp.AA, "\n")
          oldAA_newAA <- paste(sep="/", ref.AA, snp.AA)
          oldcodon_newcodon <- paste(sep="/", ref.codon, snp.codon)
          CDS_size <- nchar(translation)
          # adjust effect
          if (FALSE) {
            if (ref.AA == snp.AA) {
              effect <- SNP_EFFECT["SYNONYMOUS_CODING"]
            } else if (snp.AA == "*") {
              effect <- SNP_EFFECT["STOP_GAINED"]
            } else if (ref.AA == "*") {
              effect <- SNP_EFFECT["STOP_LOST"]
            } else {
              effect <- SNP_EFFECT["NON_SYNONYMOUS_CODING"]
            }
          }

          if (debug) 
            cat("  - effect", effect, "\n")
        }
      }
      break
    }
  }

  # create sequence window around indel
  w1start <- max(pos - 10, 1); w1width <- pos - w1start
  w2start <- pos + 1; w2width <- min(width(genome) - pos, 10)
  window1010bp <- paste(sep="|", 
                        as.character(subseq(genome, start=w1start, width=w1width)), 
                        change,
                        as.character(subseq(genome, start=w2start, width=w2width))) 
  note <- unique.no.NA(feat$product)
  if ((t1 <- unique.no.NA(feat$protein_id))[1] != "")
    note <- paste(sep="; ", note, paste(sep=":", collapse="; ", "protein_id", t1))
  if ((t1 <- unique.no.NA(feat$locus_tag))[1] != "")
    note <- paste(sep="; ", note, paste(sep=":", collapse="; ", "locus_tag", t1))

  if (verbose) {
    if (gene_name != "")
      cat(gene_name, " ")
    cat(effect, "\n")
  }

  ans <- data.frame(
    chromosome=I(names(genome)[1]),
    position=pos,
    ref=I(ref),
    change=I(paste(sep="", "\'", change, "\'")),
    changetype=I("SNP"),
    Homozygous=I("Hap"),
    changefreq=changefreq,
    gene_id=I(""),
    gene_name=I(gene_name),
    strand=paste(sep="/", strand),
    effect=I(effect),
    oldAA_newAA=I(oldAA_newAA),
    oldcodon_newcodon=I(oldcodon_newcodon),
    codon_num=codon_num,
    CDS_size=CDS_size,
    window1010bp=window1010bp,
    note=note)

  ####
  ans
}
