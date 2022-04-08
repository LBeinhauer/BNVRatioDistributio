### Ratio Distribution (Normal) in R ###



# Probability Density Functions of:

# uncorrelated central normal ratio #

# pdf of bivariate uncorrelated gaussian with means 0
p_uncor_cent_xy <- function(x, y){
  (1/sqrt(2*pi)) * exp(-(1/2)*x^2) * (1/sqrt(2*pi)) * exp(-(1/2)*y^2)
} 


p_uncor_cent_xy(0, 0)

# pdf of ratio Z = X/Y
p_uncor_cent <- function(z){
  (1/(pi))/(1 + z^2)
}

p_uncor_cent(0)



# uncorrelated non-central normal ratio

# pdf of ratio Z = X/Y, with non-zero means
p_uncor <- function(z, sig_y, sig_x, mu_x, mu_y){
  
  a <- sqrt((1/(sig_x^2))*(z^2) + (1/(sig_y^2)))
  
  b <- sqrt((mu_x/(sig_x^2))*z + (mu_y/(sig_y^2)))
  
  c <- ((mu_x^2)/(sig_x^2)) + ((mu_y^2)/(sig_y^2))
  
  d <- exp(((b^2) - (c*(a^2)))/(2*(a^2)))
  
  p <- ((b*d)/(a^3)) * (1/(sqrt(2*pi)*sig_x*sig_y)) * (pnorm((b/a)) - pnorm(-(b/a))) + 
    (1/((a^2)*pi*sig_x*sig_y)) * exp(-(c/2))
  
  p
}

p_uncor(1, 1, 1, 1, 1)



# correlated central normal ratio

# pdf of ratio Z = X/Y, with zero means and correlated variables
p_cor_cent <- function(z, sig_y, sig_x, cor){
  
  alph <- cor*(sig_x/sig_y)
  
  bet <-  (sig_x/sig_y) * sqrt(1 - cor)
  
  p <- (1/pi) * (bet/(((z-alph)^2) + (bet^2)))
  
  p
    
}

# p_cor_cent(1, 1, 1, .9140783)




# correlated non-central normal ratio (Katz Approx.)

z = (x + mu_x)/(y + mu_y)

xadj <- (x - cor*y*(sig_x/sig_y))

sig_xadj <- (sig_x * sqrt(1-cor))

# this means:

z = (xadj + (cor*y*(sig_x/sig_y)) + mu_x)/(y + mu_y)




# comples gaussian ratio distribution (Baxley)

# if ...







