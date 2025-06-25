library(arrow)
library(stats)

# suppose report-lib.parquet is a predicted lib converted to .parquet
lib<-read_parquet('E:/NML/Rd05_Alignment_LibConversion/output_library_predicted.parquet',as_data_frame = F)
libr<-read_parquet('E:/NML/Rd05_Alignment_LibConversion/output_library_predicted.parquet')

# suppose C:/Raw/MultiColumn/regular.parquet is a report based on a regular run
df<-read_parquet('E:/NML/Rd05_AlignmentSearch/report.parquet')

prs <- intersect(libr$Precursor.Id, df$Precursor.Id)
x <- libr[match(prs, libr$Precursor.Id),]
y <- df[match(prs, df$Precursor.Id),]

fit <- loess(y$RT ~ x$RT)
pred <- predict(fit, libr$RT)

libr$RT = pred
libr <- libr[!is.na(libr$RT),] # predict() returns NA for points that would require extrapolation - just remove those from the library
alibr <- arrow_table(libr, schema = lib$schema)
write_parquet(alibr, 'E:/NML/Rd05_Alignment_LibConversion/aligned.parquet')
