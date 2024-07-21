# Obermoser: quickly copy paste factor names in makeContrasts
# This is for training finger, ps() can be used for anything in Obermoser

ps <- function(x, grep = NA) {
  if (!is.na(grep)) {
    return(paste("(",paste(x[grep(grep,x)], collapse = " + "), ")/",length(x[grep(grep,x)]), sep = ""))
  } else {
    return(paste("(",paste(x, collapse = " + "), ")/", length(x), sep = ""))
  }
}

## For training_set_finger
one <- colnames(design)[grep("168",colnames(design))]
two <- colnames(design)[grep("0",colnames(design))]
three <- colnames(design)[grep("1.5",colnames(design))]
four <- colnames(design)[grep("3\\.",colnames(design))]
five <- colnames(design)[grep("\\.6",colnames(design))]
six <- colnames(design)[grep("9",colnames(design))]
seven <- colnames(design)[grep("12",colnames(design))]
eight <- colnames(design)[grep("15",colnames(design))]
nine <- colnames(design)[grep("24",colnames(design))]
ten <- colnames(design)[grep("36",colnames(design))]
eleven <- colnames(design)[grep("48",colnames(design))]

## For training_set_vein
one <- colnames(design)[grep("\\.\\.7",colnames(design))]
two <- colnames(design)[grep("\\.0\\.",colnames(design))]
three <- colnames(design)[grep("\\.1\\.",colnames(design))]
four <- colnames(design)[grep("3\\.",colnames(design))]
five <- colnames(design)[grep("[FM]\\.7",colnames(design))]
six <- colnames(design)[grep("10",colnames(design))]
seven <- colnames(design)[grep("14",colnames(design))]
eight <- colnames(design)[grep("21",colnames(design))]
nine <- colnames(design)[grep("28",colnames(design))]

## For test_set_finger
one <- colnames(design)[grep("\\.\\.7",colnames(design))]
two <- colnames(design)[grep("\\.0\\.",colnames(design))]
three <- colnames(design)[grep("[FM]\\.7",colnames(design))]

## For test_set_vein
one <- colnames(design)[grep("\\.\\.7",colnames(design))]
two <- colnames(design)[grep("\\.0\\.[FSP]",colnames(design))]
three <- colnames(design)[grep("\\.0\\.5",colnames(design))]
four <- colnames(design)[grep("\\.1\\.",colnames(design))]
five <- colnames(design)[grep("3\\.",colnames(design))]
six <- colnames(design)[grep("[FM]\\.7",colnames(design))]
seven <- colnames(design)[grep("10",colnames(design))]
eight <- colnames(design)[grep("28",colnames(design))]