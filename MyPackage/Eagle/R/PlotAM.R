#' @title Visualisation of multiple locus association mapping results
#' @description    A plotting function that provides additional information on the significant 
#'     marker-trait associations found by \code{\link{AM}}
#' @param  AMobj  the (list) object obtained from running \code{\link{AM}}. 
#' @param  itnum  the iterate number of the model building process whose results are to be viewed.
#' @param  chr  either "All" for all chromosomes or the label of a specific chromosome.
#' @param  type either "Manhattan" or  "Score"  for a manhattan plot or a plot of the 
#'    score statistics across SNPs. 
#'
#' @details
#'  A function useful for viewing the strength of association across the whole genome and 
#'  how this association changes as the model is built. 
#'
#'  The  score statistics (\code{type="Score"}) or p-values 
#'  (\code{type="Manhattan"}) are plotted against the location of the SNPs.  The horizontal lines 
#'  denote the location of the SNPs found by \code{\link{AM}}. The number gives the order 
#'  in which they were found. 
#'
#'  A single chromosome or all (\code{chr="All"}) chromosomes can be viewed.  
#'
#'  By setting \code{itnum} to different values, how the score statistics or p-values 
#'  increase/decrease 
#'  over the model building process can be observed. 
#'
#' 
#' @examples
#'  \dontrun{
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #---------------
#'   # read the map 
#'   #---------------
#'   #
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   # read marker data
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   # read phenotype data
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   # Perform multiple-locus genome-wide association mapping 
#'   res <- AM(trait = 'y',
#'                            fformula=c("cov1 + cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj)
#'
#'  # Plotting the p-values from the first iteration of the module building process. 
#'  # You can see why Eagle has identified the SNP that is has. 
#'   PlotAM(AMobj=res, itnum=1)
#'
#'
#'  # Plotting the results from the final step of the model building process
#'  # By accounting for the effect of SNP in strong association with the trait, the 
#'  # strength of association changes across the genome. 
#'   PlotAM(AMobj=res, itnum=3)
#'
#'  }
#'
#' 
#' 
#' @seealso \code{\link{AM}}
#'
PlotAM <- function(AMobj=NULL, itnum=1, chr="All", type="Manhattan" )
{


 if(is.null(AMobj)){
    message(" PlotAM function requires AMobj object to be specified. This object is obtained by running AM().")
    return(NULL)
    }
 if(!is.list(AMobj)){
    message(" PlotAM function requires AMobj object to be a list object.")
    return(NULL)
   }

if (!is.integer(itnum)){

 if(itnum < 1 ||  (itnum > length(AMobj$outlierstat)  )){
    message(" The parameter itnum must have an integer value between ", 1, " and ",
         length(AMobj$outlierstat))
    return(NULL)
  }
}

 if (!(chr=="All")  &&  !any(chr %in%  AMobj$Chr) ){
    message(" The parameter chr does not match any of the chromosome labels. ")
    message(" The allowable chromosome labels are ", c("All", unique(AMobj$map[,2])))
    return(NULL)
 }



  # we do not have a map
     xindx <- 1:length( AMobj$outlierstat[[itnum]] )
     xvals <- xindx

     yvals <- AMobj$outlierstat[[itnum]]
     isit  <- IsItBigger(vals=AMobj$outlierstat, itnum=itnum )
     bigger <- isit$bigger


     percentagechange <- isit$percentagechange
     chrm <- rep(1, length(xindx))
     pos  <- xvals

   # map exisits
     if(!is.null(AMobj$map)){
        if(chr != "All"){
          # picking a single chrm to plot
          xindx <- which(as.character(AMobj$map[,2]) == chr)
          xvals <- AMobj$map[xindx, ncol(AMobj$map)]
          chrm <-  AMobj$map[xindx, 2]
          pos <- xvals
          yvals <-  AMobj$outlierstat[[as.numeric(itnum)]][xindx]
          isit <- IsItBigger(vals=AMobj$outlierstat, itnum=itnum, xindx=xindx )
          bigger <- isit$bigger
          percentagechange <- isit$percentagechange
       } else {
          # plotting all the chromosomes - more difficult
          # reordering based on chrm then map position
          oindx <- order(AMobj$map[,2], AMobj$map[, ncol(AMobj$map)])
          yvals <- AMobj$outlierstat[[as.numeric(itnum)]][oindx]  ## reordering yvals

          if( as.numeric(itnum)  > 1){
             bigger <- rep(""    , length(  AMobj$outlierstat[[as.numeric(itnum)]][oindx] ) )
             percentagechange <- rep(0, length(  AMobj$outlierstat[[as.numeric(itnum)]][oindx] ) )


            a <-  AMobj$outlierstat[[as.numeric(itnum)]][oindx]
            b <- AMobj$outlierstat[[as.numeric(itnum) - 1 ]][oindx]

             indx <- which(  a >  b )
             bigger[indx] <- "Increased value"
             percentagechange[indx] <-  (( b - a)) [indx]

             indx <- which(  a <=  b )
             bigger[indx] <- "Decreased value"
             percentagechange[indx] <-  (( a - b)) [indx]

          }


          mapordered <- AMobj$map[oindx,]
          # map position is within chrm, need cumulative postion. 
          chrms <- unique(mapordered[,2])
          xvals <- mapordered[, ncol(mapordered)]
          if (length(chrms) > 1){
            xvals <- rep(0, nrow(mapordered))
            indx <- which(mapordered[,2] == chrms[1])
            xvals[indx] <- mapordered[indx,ncol(mapordered)]
            genometot <- max(xvals)
             for(ii in chrms[-1]){
               indx <- which(mapordered[,2] == ii)
               xvals[indx] <- mapordered[indx, ncol(mapordered)] + genometot
               genometot <- max(xvals)
              }  ## end for
          } ## end if length(chrms)
           chrm <- mapordered[,2]
           pos <- mapordered[, ncol(mapordered)]
       }  # if else 


     }  ##  if(!is.null(map))

    xlabel <- "Map Position (bp)"
    if(is.null(AMobj$map ))
      xlabel <- "Column Position of SNP"

    ylabel <- "Score Statistic"
    if(type=="Manhattan")
       ylabel <- "-log10(p value)"


     # addition on SNP-trait positions on map
     if(length(AMobj$Chr) > 1){  ## first entry of list is always NA
         # found associations 
          found.chr <- AMobj$Chr[!is.na(AMobj$Chr)]
          found.pos <- AMobj$Pos[!is.na(AMobj$Pos)]
          found.label <- 1:length(found.chr)  ## used for annotation in plot
      }


     # place on -lgo10 scale if manhattan selected
     if(type=="Manhattan"){
       yvals[is.nan(yvals)] <- 0
       yvals[yvals < 0] <- 0  ## rounding error - very close to 0 when negative
       ts <- sqrt(yvals)
       pval <- 1 - pnorm(ts)
       logp <- -1*log10(pval)
       yvals <- logp
     }


     # create data frame for plotting 
     df <- data.frame(x=xvals, y=yvals, chrm=chrm, pos=pos, foundchr=FALSE, foundpos=FALSE, foundlabel=0 )

     # check for SNP findings from AM
     if (length(AMobj$Chr)>1){
          for(ii in 1:length(found.chr)){
             indx <- which(df$chrm == found.chr[ii])
             if(!is.null(indx))
                   df$foundchr[indx] <- TRUE
             indx <- which(df$pos == found.pos[ii])
             if(!is.null(indx))
                   df$foundpos[indx] <- TRUE
              indx <- which(df$foundchr & df$foundpos & df$foundlabel==0)
              if(!is.null(indx))
                   df$foundlabel[indx] <- ii
       }  ## end for ii
                                    
     }  ## end if length()       
     geomX <- with(df, x[foundchr&foundpos])
     geomLabels <- with(df, foundlabel[foundchr&foundpos])

     if(is.null(bigger)){
       p  <- ggplot(data=df, aes(x=df$x, y=df$y )) + geom_point()

     } else {

       df$Increase <- bigger ## used for color coding points that have increased/decreased from previous iteration 
       df$Percentagechange <- percentagechange

       p  <- ggplot(data=df, aes(x=df$x, y=df$y , color=df$Increase, size=df$Percentagechange) )  + geom_point() +  scale_color_manual(values=c("#3b5998","#cae1ff" ))
    }

  p <- p + theme_hc()
            p <- p + ylab(ylabel) + xlab(xlabel)
            p <- p +  theme(legend.title=element_blank())  ## no legend title
            p <- p + theme(legend.position="right")

            if(!is.null(geomX)){
              for(ii in geomX){
               yadj <- sample(seq(0.5,0.9,0.1), 1)
               p <- p + geom_vline(xintercept = ii, linetype="solid", color="#FFE4B5", size=1.5)
               p <- p + annotate("text", size=8, label=geomLabels[which(geomX==ii)] ,
                       x=(ii  - ( diff(range(df$x))*0.02) )   ,
                        y = max(df$y)*yadj )
              }  ## end for ii
              p <- p + scale_size(guide='none')
            } ## if !is.null
           p <- p + theme( axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
                           axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),
                           legend.text=element_text(size=14) )
           p <- p  + guides(colour = guide_legend(override.aes = list(size=8)))  ## changing point size in legend

           p <- p + theme(legend.position = 'bottom', legend.spacing.x = unit(0.5, 'cm'))





           if(chr=="All"){
              # entire genome
               txt2 <- "across all chromosomes"
               if(type=="Manhattan"){
                  txt1 <- " -log p value of the score statistic"
                  txt3 <- "-log p value"
               } else {
                 txt1 <- "score statistic"
                  txt3 <- txt1
               } ## inner else
            } else{
              # chrm selected
               txt2 <- paste("on chromosome", chr)
               if(type=="Manhattan"){
                  txt1 <- " -log p value of the score statistic"
                  txt3 <- "-log p value"
               } else {
                  txt1 <- "score statistic"
                  txt3 <- txt1

               } ## inner else
            }  ## end outer else






return(p)




}  ## end function PlotAM ... 


   IsItBigger <- function(vals, itnum, xindx=NULL){
   ## calculates percentage change in size (increase or decrease)
   bigger <- NULL
   percentagechange <- NULL
   if(as.numeric(itnum) > 1){
       # entire genome
       bigger <- rep("" , length(vals[[as.numeric(itnum)]]) ) ## initialize 
       percentagechange <- rep(0, length(vals[[as.numeric(itnum)]]) ) ## initialize
       a <-  vals[[as.numeric(itnum)]]
       b <- vals[[as.numeric(itnum) - 1 ]]
       # a > b       
       indx <- which( a  >  b )
       bigger[indx] <- "Increased value"
       percentagechange[indx] <- ( (a  -    b ) )[indx]
       # a < b
       indx <- which( a <=  b  )
       bigger[indx] <- "Decreased value"
       percentagechange[indx] <- ( (b  -    a ) )[indx]

       if(!is.null(xindx)){
          # reduce to chromosome
          bigger <- bigger[xindx]
          percentagechange <- percentagechange[xindx]
       }
     }
     res <- list(bigger=bigger, percentagechange=percentagechange)
     return(res)
  }




