library(RcppColors)
library(jacobi)
library(rgl)
library(Rvcg)

mesh <- MeshesOperations::torusMesh(1, 0.4, nu = 800, nv = 700)

Z <- apply(mesh$vb[-4L, ], 2L, function(xyz){
  a <- xyz[1]
  b <- xyz[2]
  c <- xyz[3]
  z <- wsigma(a + 1i* b, tau = (1i+c(crossprod(xyz)))/2)
})

color <- colorMap1(Z*Z, reverse = c(TRUE, FALSE, TRUE))
mesh$material <- list(color = color)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("lightgrey")
shade3d(mesh)



# # -- if you want an animation
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 120,
  duration = 1,
  dir = ".",
  movie = "zzpic",
  convert = FALSE,
  webshot = FALSE
)

command <- "gifski --fps=9 --frames=zzpic*.png -o SigmaCyclide.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
