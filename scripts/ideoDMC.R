# PiGx BSseq Pipeline.
#
# Copyright © 2017 Altuna Akalin <altuna.akalin@mdc-berlin.de>
# Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
# It has been adapted from
# https://gist.github.com/al2na/f002dafb1f3a5782f099#file-ideodmc-r
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' function for making ideogram for differential methylation values
#' requires methylKit, ggbio and GenomicRanges
#'
#' @example
#' library(BSgenome)
#' library("BSgenome.Hsapiens.UCSC.hg18")
#' chr.len = seqlengths(Hsapiens)  # get chromosome lengths
#' # remove X,Y,M and random chromosomes
#' chr.len = chr.len[grep("_|M|X|Y", names(chr.len), invert = T)] 
#' 
#' download.file("http://methylkit.googlecode.com/files/myDiff.rda", 
#'               destfile = "myDiff.rda")
#' load("myDiff.rda")
#' 
#' ideoDMC(myDiff, chrom.length = chr.len, difference = 25, qvalue = 0.01, 
#'        circos = TRUE, title = "test", hyper.col = "magenta", hypo.col = "green")
#'        
ideoDMC <- function(methylDiff.obj, chrom.length, difference = 25, 
                    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
                    hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  # TODO: rewrite the require statement as library::function calls to not pollute the reports environment

  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
                                                                     width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  
  hypo = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                        type = "hypo")
  hyper = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                         type = "hyper")
  
  g.per = as(hyper, "GRanges")
  seqlevels(g.per) = seqlevels(myIdeo)
  seqlengths(g.per)=(chrom.length)
  
  g.po = as(hypo, "GRanges")
  seqlevels(g.po) = seqlevels(myIdeo)
  seqlengths(g.po)=(chrom.length)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                                  radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                           size = 1, aes(x = midpoint, 
                                         y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
      scale_colour_manual(values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
                      vjust = 0, radius = 55, trackWidth = 7) + labs(title = title)
    
  } else {
    
    p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
    p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
                         aes(x = midpoint, 
                             y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
                                                                                          hypo.col)) + labs(title = title)
    # new alternative commented out
    #autoplot(c(g.po, g.per), layout = "karyogram", geom = "point", size = 0.65, 
    #aes(x = midpoint,y = meth.diff, color = id))  + scale_colour_manual(values = c(hyper.col, 
    #                                                                                        hypo.col)) + labs(title = title)
    
  }
}

ideoDMC_hyper_hypo <- function(hyper, hypo, chrom.length, 
                               circos = FALSE, 
                               title = "Differentially methylated cytosines", 
                               hyper.col = "magenta", 
                               hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
                                                                     width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  
  if(class(hypo) %in% c("methylDiffDB","methylDiff")) {
    hypo <- as(hypo, "GRanges")
  }
  if(class(hyper) %in% c("methylDiffDB","methylDiff")) {
    hyper <- as(hyper, "GRanges")
  }
  
  g.per = hyper
  seqlevels(g.per) = seqlevels(myIdeo)
  seqlengths(g.per)=(chrom.length)
  
  g.po = hypo
  seqlevels(g.po) = seqlevels(myIdeo)
  seqlengths(g.po)=(chrom.length)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                                  radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                           size = 1, aes(x = midpoint, 
                                         y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
      scale_colour_manual("Regions",values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
                      vjust = 0, radius = 55, trackWidth = 7) + labs(title = title)
    
  } else {
    
    # p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
    # p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
    #                      aes(x = midpoint, 
    #                          y = meth.diff, color = id)) + scale_colour_manual("Regions", values = c(hyper.col, 
    #                                                                                                  hypo.col)) + labs(title = title)
    # new alternative commented out
    d = c(g.per, g.po)
    p = autoplot(myIdeo, layout = "karyogram")
    p + layout_karyogram(d, geom = "point", size = 0.65,
                         aes(x = start,  y = meth.diff, color = id))+
      scale_colour_manual("", values = c(hyper.col, hypo.col))+
      labs(title = title) 
  }
}




