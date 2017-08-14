##Takes the file containing the names of the samples and an ibs matrix and labels which two samples are being compared at each point in the matrix 
wildintrogressions <- read.table(file="/Users/connordepies/WildIntrogression/data/plink.mibs", header = FALSE, sep = "\t")
wildintrogresssions.names <- read.table(file = "/Users/connordepies/WildIntrogression/data/names.txt")
colnames(x=wildintrogressions) <- wildintrogresssions.names$V1
rownames(x=wildintrogressions) <- wildintrogresssions.names$V1
write.table(x=wildintrogressions, file="/Users/connordepies/WildIntrogression/data/labledintrogressions.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

