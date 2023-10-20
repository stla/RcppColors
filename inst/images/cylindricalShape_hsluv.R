library(rgl)
library(cgalMeshes)
library(RcppColors)

f <- function(θ, z){
  ρ <- (1 - 0.25*z*z) * (1 + 0.5*sin(1.5*pi*z) + 0.3*cos(5*θ)) 
  x <- ρ*cos(θ)
  y <- ρ*sin(θ)
  rbind(x, y, z)  
}

fcolor <- Vectorize(function(θ, z) {
  d <- θ * 180/pi
  h <- if(d < 180) 2*d else 360-2*(d-180)
  s <- 90
  l <- 25 * (z + 2)
  hsluv(h, s, l)
})

mesh <- parametricMesh(
  f, c(0, 2*pi), c(-2, 2), periodic = c(TRUE, FALSE), 
  nu = 512L, nv = 512L, fcolor = fcolor
)
mesh <- addNormals(mesh)

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -90, zoom = 0.6)
shade3d(mesh)

# animation ###
movie3d(spin3d(axis = c(0, 0, 1), rpm = 60),
        duration = 1, fps = 60,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/60)

pngs <- Sys.glob("zzpic*.png")
gifski::gifski(
  pngs,
  "cylindricalShape_hsluv.gif",
  width = 512, height = 512,
  delay = 1/11
)
file.remove(pngs)
