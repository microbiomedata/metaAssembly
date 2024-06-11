#!/bin/bash
#


# zip bundle.zip *wdl
WD=/global/cfs/cdirs/m3408/aim2/dev/kli_training/mg_annotation/

curl --netrc -X POST "https://nmdc-cromwell.freeddns.org:8443/api/workflows/v1" \
	-H "accept: application/json" \
	-H "Content-Type: multipart/form-data" \
	-F "workflowSource=@jgi_assembly.wdl" \
	-F "workflowInputs=@input.json;type=application/json"  \
    -F "workflowDependencies=@imports.zip" \
	-F "labels=@labels.json;type=application/json"

