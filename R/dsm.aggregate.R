dsm.aggregate <- function(predicted.values=predicted.values, prediction.layers=prediction.layers,
                          path) 
#
#   dsm.aggregate - purpose: sum predicted abundances over user-specified polygons
#
#   Arguments:
#         predicted.values        - array of predicted values produced by predict.gam()
#         prediction.layers       - known as 'prediction.layers.dat.r' that will be derived from
#                                   the prediction grid layer; one row for each record in the prediction
#                                   grid, entries are ID of polygon being predicted over; 0 means
#                                   don't predict [manufactured internally within Distance]
#         path                    - needed to send R looking in the temporary directory for files
#         prediction.layer.names  - for each layer being predicted over; a labelling dataframe, 
#                                   consisting of an ID and the label associated with that ID
#                                    [manufactured internally within Distance]   
#
#    Value:
#         labelled.result - dataframe with labelled layer-specific abundances
#     Functions used:  aggregate()
{
#     Cat a bookmark for the benefit of the Distance results file
cat('\tResponse Surface/Prediction: Aggregate sums\t\n')


#         proactively take action in the case of missing values coming back from predict
#           use 'merge' rather than merely data.frame; use all.y=TRUE so all cells defined in
#           prediction.layers will be represented
#     maddeningly, predict.gam and predict.glm return 
#                        their predicted values in slightly different forms
 if (is.null(dimnames(predicted.values))) {
    predicted.values.df <- data.frame(ID=as.integer(unlist(names(predicted.values))), values=predicted.values)
 } else {
    predicted.values.df <- data.frame(ID=as.integer(unlist(dimnames(predicted.values))), values=predicted.values)
 }
   predict.with.ids <- merge(predicted.values.df, prediction.layers, by="ID", all.y=TRUE)
# predict.with.ids <- data.frame(predicted.values, prediction.layers)
 file.prefix <- "pred.lyr.dat.r."
 #  Step through as many files as there are layers over to be predicted
 for (i in 2:length(names(prediction.layers)))
 {
   non.labelled.result <- aggregate(predict.with.ids, list(ID=prediction.layers[,i]), sum, na.rm=TRUE)
   table.name <- paste(path, file.prefix, names(prediction.layers)[i], sep="")
   temporary.label.object <- read.table(file=table.name, header=TRUE, sep="\t")
   labelled.result <- merge(temporary.label.object, non.labelled.result[,c(1,3)], by="ID")
   cat('Aggregated predicted values:')
   print(write.table(labelled.result[,2:3], row.names=FALSE, quote=FALSE))
 }
}
