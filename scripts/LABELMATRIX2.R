# Connor Depies August, 15, 2017
##Takes the file containing the names of the samples and an ibs matrix and labels which two samples are being compared at each point in the matrix 
wildintrogressions <- read.table(file="/Users/connordepies/WildIntrogression/data/onlywild.mibs", header = FALSE, sep = "\t")
wildintrogressions.names <- read.table(file = "/Users/connordepies/WildIntrogression/data/wildnames.txt")
colnames(x=wildintrogressions) <- wildintrogressions.names$V1
rownames(x=wildintrogressions) <- wildintrogressions.names$V1
write.table(x=wildintrogressions, file="/Users/connordepies/WildIntrogression/data/labledwildonly.txt", sep ="\t", row.names = TRUE, col.names = NA, quote = FALSE)