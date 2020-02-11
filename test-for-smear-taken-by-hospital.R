
# test for difference in smear taken
# between hospitals
#


library(MASS)

# data from main.R

dat$taken <- dat$Smear != "Not taken"

tbl <- 
  table(dat$taken,
      dat$SiteID)

tbl %>% 
  prop.table(margin = 2) %>%
  round(2)


d <- data.frame(A = c(4,34),
                B = c(0,17),
                C = c(0,3),
                F = c(1,5),
                G = c(0,7),
                H = c(1,11),
                L = c(5,16),
                N = c(1,56),
                # R = c(0,0), # this breaks chisq.test()
                S = c(1,9))

chisq.test(d)
# chisq.test(tbl)
