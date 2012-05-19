\name{dsm.aggregate}
\alias{dsm.aggregate}
\title{Sum predicted abundances over user-specified polygons}
\usage{
  dsm.aggregate(predicted.values = predicted.values,
    prediction.layers = prediction.layers, path)
}
\arguments{
  \item{predicted.values}{array of predicted values
  produced by predict.gam()}

  \item{prediction.layers}{known as
  'prediction.layers.dat.r' that will be derived from the
  prediction grid layer; one row for each record in the
  prediction grid, entries are ID of polygon being
  predicted over; 0 means don't predict (manufactured
  internally within Distance).}

  \item{path}{needed to send R looking in the temporary
  directory for files}

  \item{prediction.layer.names}{for each layer being
  predicted over; a labelling dataframe, consisting of an
  ID and the label associated with that ID (manufactured
  internally within Distance).}
}
\value{
  dataframe with labelled layer-specific abundances
}
\description{
  Sum predicted abundances over user-specified polygons
}

