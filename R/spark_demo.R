## Apache Spark Demo

library(sparklyr)
sc <- spark_connect(master = "local")
library(dplyr)
iris_tbl <- copy_to(sc, iris)
install.packages("nycflights13")
install.packages("Lahman")
flights_tbl <- copy_to(sc, nycflights13::flights, "flights")
batting_tbl <- copy_to(sc, Lahman::Batting, "batting")
src_tbls(sc)

## using dplyr 



## qou can also use SQL direcly