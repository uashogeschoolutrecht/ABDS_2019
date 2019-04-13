## get data from 
## Bike Share 
# https://s3.amazonaws.com/capitalbikeshare-data/index.html


#Loading the rvest package
library('rvest')

#Specifying the url for desired website to be scraped


# amazon open data buckets
# https://registry.opendata.aws/

library("aws.s3")

y <- get_bucket(bucket = "nyc-tlc")
y

y[[2]]

x <- get_object(y[[3]])

key <- y[[2]]$Key

#aws.s3::s3read_using(read.csv, object = key, bucket = "nyc-tlc")

?s3read_using

library(sparklyr)



## air quality data
##https://github.com/ropensci/ropenaq
## aws dataset (complete bucket) https://registry.opendata.aws/openaq/
## https://openaq.org/
## part of the data:
# devtools::install_github("ropensci/ropenaq")


library(ropenaq)
## relational data

countries <- aq_countries()
cities <- aq_cities(country = "NL")
## get data
air_nl <- aq_measurements(city = "Amsterdam", country = "NL", limit = 100)


#Lets try to read using sparklyr packages 
y <- get_bucket(bucket = "openaq-fetches")

i = 2

bucket_item <- y[[i]]

bucket_item$Bucket
bucket_item$Key

library(sparklyr)
library(dplyr)
sc <- spark_connect(master = "local")

usercsv_tbl <- 
  sparklyr::spark_read_csv(sc,name = "usercsvtlb",
                              path = 
                   file.path("s3a",
                             bucket_item$Bucket,
                             bucket_item$Key),
                 memory = FALSE)
