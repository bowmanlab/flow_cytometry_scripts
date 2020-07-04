# flow_cytometry_scripts
Base scripts for flow cytometry analysis.

This repository has basic scripts for analyzing data from our Guava 11HT, generalizable to most flow cytometry instruments.  Some older scripts from our defunt CyFlow Space are also present in the old_cyflow directory.

fcm_model.r - Constructs a model to classify FCM data from the Guava 11HT.

fcm_predict.r - Classifies data according to a model constructed by fcm_model.r.

Current models and associated figures can be found as .Rdata and .pdf files, respectively.
