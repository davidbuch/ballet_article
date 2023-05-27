# The raw data loaded in this script come from
# http://icg.port.ac.uk/~nicholb/edsgc/other.html
edsgc_ascii_colnames <- c("RA",
                          "DEC",
                          "XMIN",
                          "XMAX",
                          "YMIN",
                          "YMAX",
                          "AREA",
                          "IMAX",
                          "COSMAGCAL",
                          "ISKY",
                          "IXCEN",
                          "IYCEN",
                          "UMAJAX",
                          "UMINAX",
                          "UTHETA",
                          "IMAJAX",
                          "IMINAX",
                          "POSANGLE",
                          "CORMAGCAL",
                          "SIGMA",
                          "IDSEQ",
                          "LOGAREA",
                          "GEOM",
                          "GEOMLOG",
                          "SPARE1",
                          "SPARE2",
                          "SPARE3"
                          )

conn <- file("data/raw_data/full_edsgc_ascii/full_edsgc_ascii.dat",open="r")
galaxy_data <- as.data.frame(scan(conn, 
                            what = c(as.list(rep(0,24)),as.list(rep("",3)))), 
                             col.names = edsgc_ascii_colnames) 
close(conn)
conn <- file("data/raw_data/full_edsgc_ascii/full_edsgc_ascii2.dat",open="r")
galaxy_data_2 <- as.data.frame(scan(conn, 
                              what = c(as.list(rep(0,24)),as.list(rep("",3)))),
                               col.names = edsgc_ascii_colnames) 
close(conn)
galaxy_data <- rbind(galaxy_data, galaxy_data_2)
rm(galaxy_data_2)
galaxy_data <- galaxy_data %>% 
  filter((RA > 0) & (RA < 40/60) & (DEC < -28) & (DEC > -38))

write.csv(galaxy_data, 
          "data/intermediate_data/primary/edsgc.csv", 
          row.names = FALSE)


galaxy_clusters <- read.fwf("data/raw_data/edcc_table_from_paper.data", 
                            skip = 3, header = FALSE,
                            col.names = c("ID", "RA", "DEC", "m1", "m3", "m10", 
                                          "n_c", "n_b", "r_A", "field", "run", 
                                          "deb", "abell"),
                            widths = c(7,13,13,7,7,7,5,5,8,7,5,3,10))
galaxy_clusters <- lapply(galaxy_clusters, trimws) %>% as.data.frame()
galaxy_clusters <- galaxy_clusters %>% 
  mutate(ID_qmark = grepl('\\?',ID),
         ID = as.numeric(gsub("\\?", "", ID)),
         RA = gsub(" +", ":", RA),
         DEC = gsub(" +", ":", DEC),
         m1 = as.numeric(m1),
         m3 = as.numeric(m3),
         m10 = as.numeric(gsub('\\*', '', m10)),
         n_c_asterisk = grepl('\\*', n_c),
         n_c = as.numeric(gsub('\\*', '', n_c)),
         n_b = as.numeric(n_b),
         r_A = as.numeric(r_A)
  )

galaxy_clusters$RA <- sapply(str_split(galaxy_clusters$RA, ":"), 
                             \(x) as.numeric(x[1]) + 
                               as.numeric(x[2]) / 60 + 
                               as.numeric(x[3]) / (60 * 60))

galaxy_clusters$DEC <- sapply(str_split(galaxy_clusters$DEC, ":"), 
                              \(x) as.numeric(x[1]) - 
                                as.numeric(x[2]) / 60 - 
                                as.numeric(x[3]) / (60 * 60))

galaxy_clusters <- galaxy_clusters %>%
  filter((RA > 0) & (RA < 40/60) & (DEC < -28) & (DEC > -38))

write.csv(galaxy_clusters, 
          "data/intermediate_data/primary/edsgc_clusters.csv",
          row.names = FALSE)
