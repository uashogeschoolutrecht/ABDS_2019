# Kamil Slowikowski
# Apr 11, 2010
#
# Display a Needleman-Wunsch alignment matrix with arrows
# showing the best possible alignments.

# Display an alignment matrix with arrows.
plotNeedlemanWunsch = function(seq1, seq2, match, mismatch, gap,
                               numbers=T, arrows=T, path=T) {
  q1 = needleman_wunsch(seq1, seq2, match, mismatch, gap)
  
  draw_alignment_matrix(q1$F, q1$Ptr, numbers=numbers, arrows=arrows)
  
  if (path) {
    bboxes(q1$Ptr)
    paths(q1$F, q1$Ptr)
  }
  
  d1 = dim(q1$Ptr)
  
  text(d1[2] / 2, d1[1] + 2, "Needleman-Wunsch")
  
  xs = seq(1, d1[2], d1[2] / 3)
  text(xs[1], d1[1] + 1, col='blue',  paste("match =", match))
  text(xs[2], d1[1] + 1, col='red',   paste("mismatch =", -mismatch))
  text(xs[3], d1[1] + 1, col='black', paste("gap =", -gap))
}

#############################################################################
# Author     : Dr. Liping Tong
# Website    : http://webpages.math.luc.edu/~ltong/
# 
# x and y are the two sequences to compare
# x is on the left and y is on the right
# m is the score for match
# s is the penalty for mismatch
# d is the penalty for alignment with a gap
needleman_wunsch = function(x,y,m,s,d)
{
  nx = nchar(x)
  xx = rep(0,nx)
  for(i in 1:nx)
    xx[i] = substr(x,start=i,stop=i)
  
  ny = nchar(y)
  yy = rep(0,ny)
  for(i in 1:ny) yy[i] = substr(y,start=i,stop=i)
  
  # initialize F and ptr
  F = ptr = matrix(0,nx+1,ny+1)
  dimnames(F) = dimnames(ptr) = list(c("",xx),c("",yy))
  for (i in 1:nx) F[i+1,1] = -i*d
  for (j in 1:ny) F[1,j+1] = -j*d
  
  # main iteration
  for (i in 1:nx)
  {
    for (j in 1:ny)
    {
      if (xx[i] == yy[j])
      {
        t1 = F[i,j] + m
      } else {
        t1 = F[i,j] - s
      }
      
      t2 = F[i,j+1] - d
      t3 = F[i+1,j] - d
      
      F[i+1,j+1] = tt = max(t1,t2,t3)
      
      if (t1 == tt) ptr[i+1,j+1] = ptr[i+1,j+1] + 2
      if (t2 == tt) ptr[i+1,j+1] = ptr[i+1,j+1] + 3
      if (t3 == tt) ptr[i+1,j+1] = ptr[i+1,j+1] + 4
    }
  }
  return(list(F = F, Ptr = ptr))
}

# Highlight cells in the grid with thick lines.
bbox = function(x, y) {
  rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, border='darkgrey', lwd=3)
}

# Create arrows between cells in the grid.
arrow1 = function(x0, y0, x1, y1, color) {
  arrows(x0, y0, x1, y1, col=color, length=0.07, angle=30, lwd=1.5)
}

# Up, Left, Up-Left
up = function(x, y, color) arrow1(x,       y + 0.3, x,       y + 0.6, color)
le = function(x, y, color) arrow1(x - 0.3, y,       x - 0.6, y,       color)
ul = function(x, y, color) arrow1(x - 0.3, y + 0.3, x - 0.6, y + 0.6, color)

bboxes = function(p) {
  rows = dim(p)[1]
  f1 = function(i, j) {
    y = rows - i
    bbox(j, y)
    if (p[i, j] %in% c(2, 5, 6, 9)) {
      f1(i - 1, j - 1)
    }
    if (p[i, j] %in% c(3, 5, 7, 9)) {
      f1(i - 1, j)
    }
    if (p[i, j] %in% c(4, 6, 7, 9)) {
      f1(i, j - 1)
    }
  }
  f1(dim(p)[1], dim(p)[2])
}

paths = function(v, p) {
  rows = dim(p)[1]
  f1 = function(i, j) {
    y = rows - i
    if (p[i, j] %in% c(2, 5, 6, 9)) {
      if (v[i, j] > v[i - 1, j - 1]) {
        ul(j, y, 'blue')
      } else {
        ul(j, y, 'red')
      }
      f1(i - 1, j - 1)
    }
    if (p[i, j] %in% c(3, 5, 7, 9)) {
      up(j, y, 'black')
      f1(i - 1, j)
    }
    if (p[i, j] %in% c(4, 6, 7, 9)) {
      le(j, y, 'black')
      f1(i, j - 1)
    }
  }
  f1(dim(p)[1], dim(p)[2])
}

# v is a matrix of values
# p is a matrix of pointers that represent arrow directions
draw_alignment_matrix = function(v, p, numbers=T, arrows=T)
{
  rows = dim(v)[1]
  cols = dim(v)[2]
  
  # Create a plot of appropriate size.
  par(mar=c(0, 0, 0, 0))
  plot('', axes=FALSE, xlab='', ylab='',
       xlim=c(-0.5, cols + 0.5),
       ylim=c(-0.5, rows + 2.5))
  
  # Draw a grid.
  segments(x0=seq(-0.5, cols + 0.5), y0=-0.5, y1=rows + 0.5, col='grey')
  segments(y0=seq(-0.5, rows + 0.5), x0=-0.5, x1=cols + 0.5, col='grey')
  
  # Print the letters of the sequences.
  for (i in 1:rows) text(0, rows - i, dimnames(v)[[1]][i], font=2)
  for (j in 1:cols) text(j, rows, dimnames(v)[[2]][j], font=2)
  
  #color = 'cornflowerblue'
  color = 'grey'
  
  for (i in 1:rows) {
    y = rows - i
    
    for (x in 1:cols) {
      if (numbers) {
        text(x, y, v[i, x])
      }
      if (arrows) {
        if (p[i, x] %in% c(2, 5, 6, 9)) { ul(x, y, color) }
        if (p[i, x] %in% c(3, 5, 7, 9)) { up(x, y, color) }
        if (p[i, x] %in% c(4, 6, 7, 9)) { le(x, y, color) }
      }
    }
  }
}
