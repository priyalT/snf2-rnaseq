library(readr)

#step 1: try to make 1 count matrix by reading one file: genes x column 2 only
#step 2: figure out how to make multiple dataframes
#step 3: join all the dataframes to make a
#        genes X samples dataframe, samples = 12

counts  <- as.data.frame(
            read_delim("results/alignment/star/ERR458502_ReadsPerGene.out.tab",
                       delim = "\t", col_names = TRUE))
counts <- as.data.frame(counts)
typeof(counts)
